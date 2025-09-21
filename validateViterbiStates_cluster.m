% validateViterbiStates_cluster.m
%
% This is the second step in the CLUSTER analysis pipeline.
% It validates HMM state assignments by calculating FRET for each Viterbi-decoded
% state segment, plotting overlaid histograms, and displaying transition counts.
% It should be run AFTER runViterbi_cluster.m.
%

function validateViterbiStates_cluster()
    clearvars -except break_debug; close all; clc;
    fprintf('--- Starting Viterbi State Validation (Cluster Workflow) ---\n');

    %% 0. Setup Paths
    project_base_path = fileparts(mfilename('fullpath'));
    paths.hmm_core = fullfile(project_base_path, 'hmm_core_functions');
    if exist(paths.hmm_core, 'dir'), addpath(genpath(paths.hmm_core)); end
    original_pwd = pwd;

    try
        %% 1. Select Processed Cluster Run Directory and Load Necessary Files
        fprintf('\n--- 1. Selecting Processed Cluster HMM Run Directory ---\n');
        selected_run_path = uigetdir(pwd, 'Select the Cluster HMM Run directory (containing Viterbi & Model results)');
        if isequal(selected_run_path, 0), disp('User cancelled. Exiting.'); return; end
        fprintf('Processing results from: %s\n', selected_run_path);
        cd(selected_run_path);

        % --- Load Processed HMM Parameters (sorted) ---
        if ~exist('HMM_parm.mat', 'file'), error('VALIDATE_CLUSTER: HMM_parm.mat not found. Please process the HMM run first.'); end
        load('HMM_parm.mat', 'Obs', 'Trans', 'K'); 
        Nstates = size(Trans, 1);
        [~, max_idx] = max(abs(K(:))); [row, col] = ind2sub(size(K), max_idx);
        if K(row, col) ~= 0 && Trans(row, col) ~= 0
             dt_model_sec = Trans(row,col) / K(row,col);
        else, error('Cannot infer dt from HMM_parm.mat.'); end

        % --- Load Est.mat to get original model index ---
        if ~exist('Est.mat', 'file'), error('VALIDATE_CLUSTER: Est.mat not found.'); end
        load('Est.mat', 'ass');
        chosen_model_idx_str = num2str(ass);

        % --- Load Viterbi Paths and Transition Analysis ---
        if ~exist('Viterbi_decoded_states_Q0_cluster.mat', 'file'), error('VALIDATE_CLUSTER: Viterbi results not found. Please run runViterbi_cluster.m first.'); end
        load('Viterbi_decoded_states_Q0_cluster.mat', 'Q0');
        if ~exist('Viterbi_transition_analysis_cluster.mat', 'file'), error('VALIDATE_CLUSTER: Transition analysis not found.'); end
        load('Viterbi_transition_analysis_cluster.mat', 'trans_num');

        % --- Load Raw Photon Data ---
        if ~exist('data.mat', 'file'), error('VALIDATE_CLUSTER: data.mat not found.'); end
        load('data.mat', 'data');
        
        %% 2. User Input (Minimum Photons per Segment)
        fprintf('\n--- 2. Getting Parameters for Segment FRET Calculation ---\n');
        prompt_min_photons = {'Enter minimum photons per segment to calculate FRET (e.g., 5 or 10):'};
        answer_min_photons = inputdlg(prompt_min_photons, 'Segment FRET Parameter', [1 60], {'5'});
        if isempty(answer_min_photons), disp('User cancelled. Exiting.'); cd(original_pwd); return; end
        min_photons_per_segment_for_fret = str2double(answer_min_photons{1});
        
        %% 3. Process Trajectories: Segment FRET and Count Transitions
        fprintf('\n--- 3. Segmenting Trajectories and Calculating FRET ---\n');
        fret_values_per_segment_by_state = cell(1, Nstates);
        
        num_trajectories = length(data);
        h_wait_val = waitbar(0, 'Processing trajectories for FRET segments...');
        
        for i_traj = 1:num_trajectories
            waitbar(i_traj/num_trajectories, h_wait_val);
            
            trajectory_photon_symbols = data{i_traj}(:,2);
            viterbi_path = Q0{i_traj};

            if isempty(viterbi_path) || isempty(trajectory_photon_symbols) || length(viterbi_path) ~= length(trajectory_photon_symbols), continue; end

            path_changes_indices = find([true; diff(viterbi_path(:)) ~= 0; true]);
            for k_segment = 1:(length(path_changes_indices)-1)
                segment_start_idx = path_changes_indices(k_segment);
                segment_end_idx = path_changes_indices(k_segment+1) - 1;
                segment_state = viterbi_path(segment_start_idx);
                segment_photons = trajectory_photon_symbols(segment_start_idx:segment_end_idx);

                if length(segment_photons) >= min_photons_per_segment_for_fret
                    n_acceptor = nnz(segment_photons == 2);
                    n_donor = nnz(segment_photons == 1);
                    if (n_acceptor + n_donor) > 0
                        fret_for_segment = n_acceptor / (n_acceptor + n_donor);
                        fret_values_per_segment_by_state{segment_state}(end+1) = fret_for_segment;
                    end
                end
            end
        end
        if ishandle(h_wait_val), close(h_wait_val); end
        
        % --- Display total transitions ---
        fprintf('\n--- Total Transition Counts ---\n');
        fprintf('Total number of transitions found across all trajectories: %d\n', sum(trans_num));

        %% 4. Plot Overlaid Histograms of Per-Segment FRET Values
        fprintf('\n--- 4. Plotting Per-Segment FRET Histograms ---\n');
        fig_hist = figure('Name', 'Per-Segment FRET by Viterbi State (Cluster)', 'Position', [200, 200, 800, 600]);
        ax_hist = axes(fig_hist);
        hold(ax_hist, 'on');

        hist_edges = -0.025:0.025:1.025;
        plot_colors = lines(Nstates);
        legend_handles = []; legend_texts = {};

        % Load original, unsorted obsmat for comparison
        load('HMM_output.mat', 'obsmat');
        model_obsmat_for_comparison = obsmat{1,Nstates}{ass};

        for i_state = 1:Nstates
            fret_values = fret_values_per_segment_by_state{i_state};
            
            if ~isempty(fret_values)
                h_bar = histogram(ax_hist, fret_values, hist_edges, 'Normalization', 'pdf', ...
                                  'FaceColor', plot_colors(i_state,:), 'EdgeColor', 'k', 'FaceAlpha', 0.7);
                legend_handles(end+1) = h_bar;
                avg_fret = mean(fret_values); std_fret = std(fret_values);
                legend_text = sprintf('State %d Segments (N=%d, E=%.2f \\pm %.2f)', ...
                                      i_state, length(fret_values), avg_fret, std_fret);
            else
                h_dummy = plot(ax_hist, NaN, NaN, 'Color', plot_colors(i_state,:), 'LineWidth', 2);
                legend_handles(end+1) = h_dummy;
                legend_text = sprintf('State %d (No valid segments)', i_state);
            end
            
            model_fret = model_obsmat_for_comparison(i_state, 2); % P(Acceptor|State)
            plot(ax_hist, [model_fret, model_fret], get(ax_hist,'YLim'), ...
                 'LineStyle', '--', 'Color', plot_colors(i_state,:)*0.6, 'LineWidth', 3, 'HandleVisibility', 'off');
            legend_texts{end+1} = [legend_text sprintf(', Model P(A|S)=%.2f', model_fret)];
        end

        hold(ax_hist, 'off');
        xlabel(ax_hist, 'Segment FRET Efficiency (Calculated from D/A Symbols)');
        ylabel(ax_hist, 'Probability Density'); xlim(ax_hist, [-0.05, 1.05]);
        
        title_str1 = sprintf('Per-Segment FRET by Viterbi State (Model #%s, Min %d photons/seg)', chosen_model_idx_str, min_photons_per_segment_for_fret);
        title_str2 = sprintf('Total Transitions = %d, HMM dt=%.1e s', sum(trans_num), dt_model_sec);
        title(ax_hist, {title_str1, title_str2}, 'Interpreter', 'none');
        
        if ~isempty(legend_handles), legend(ax_hist, legend_handles, legend_texts, 'Location', 'best', 'Interpreter', 'tex'); end
        grid(ax_hist, 'on'); box(ax_hist, 'on');

        %% 5. Saving Results
        fprintf('\n--- 5. Saving Validation Results ---\n');
        saveas(fig_hist, 'State_Segment_FRET_Histograms_Overlaid_Cluster.png');
        
        save_struct.fret_values_per_segment_by_state = fret_values_per_segment_by_state;
        save_struct.min_photons_per_segment = min_photons_per_segment_for_fret;
        save_struct.total_transitions = sum(trans_num);
        save('State_Segment_FRET_Data_Cluster.mat', '-struct', 'save_struct');
        fprintf('Validation histograms and data saved to .png and .mat files.\n');

        fprintf('\n--- Viterbi State Validation (Cluster) Finished ---\n');
        cd(original_pwd);

    catch ME_validate
        fprintf(2,'ERROR in validateViterbiStates_cluster: %s\n', ME_validate.message);
        disp(ME_validate.getReport());
        if exist('h_wait_val', 'var') && ishandle(h_wait_val), close(h_wait_val); end
        cd(original_pwd);
    end
end