% validateViterbiStates_local.m
% Validates HMM state assignments by calculating FRET for each Viterbi-decoded
% state segment, plotting overlaid histograms, and counting state transitions.

function validateViterbiStates_local()
    clearvars -except break_debug; close all; clc;
    fprintf('--- Starting Viterbi State Validation (FRET Histograms & Transitions) ---\n');

    %% 0. Setup Paths
    project_base_path = fileparts(mfilename('fullpath'));
    paths.results_base = fullfile(project_base_path, 'HMM_Local_Results');
    paths.hmm_core = fullfile(project_base_path, 'hmm_core_functions');
    if exist(paths.hmm_core, 'dir'), addpath(genpath(paths.hmm_core)); end
    original_pwd = pwd;

    try
        %% 1. Select HMM Run Directory and Load Necessary Files
        % ... (This section remains IDENTICAL to the previous version) ...
        fprintf('\n--- 1. Selecting HMM Run Directory and Loading Data ---\n');
        if ~exist(paths.results_base, 'dir'), error('VALIDATE_VITERBI: Base results directory not found: %s', paths.results_base); end
        selected_run_path = uigetdir(paths.results_base, 'Select the HMM_Run_* directory containing Viterbi & Model results');
        if isequal(selected_run_path, 0), disp('User cancelled. Exiting.'); return; end
        fprintf('Processing results from: %s\n', selected_run_path);
        cd(selected_run_path);

        model_params_loaded = []; Nstates = NaN; dt_model_sec = NaN; chosen_model_idx_str = 'N/A';
        if exist('Best_HMM_Model_Local.mat', 'file')
            fprintf('Loading Best_HMM_Model_Local.mat...\n');
            load('Best_HMM_Model_Local.mat', 'best_model_overall', 'config_params_to_save', 'best_idx');
            model_params_loaded = best_model_overall;
            if isfield(config_params_to_save, 'Nstates') && isfield(config_params_to_save, 'dt_analysis_sec')
                Nstates = config_params_to_save.Nstates;
                dt_model_sec = config_params_to_save.dt_analysis_sec;
            else, error('VALIDATE_VITERBI: Nstates or dt_analysis_sec not found in config_params_to_save.'); end
            chosen_model_idx_str = num2str(best_idx);
            fprintf('INFO: Using Model #%s. Nstates = %d, dt_model_sec = %.2e s.\n', chosen_model_idx_str, Nstates, dt_model_sec);
        else, error('VALIDATE_VITERBI: Best_HMM_Model_Local.mat not found. Ensure runHMM_local & processHMM_results_local have run.'); end

        if isempty(model_params_loaded) || ~isfield(model_params_loaded, 'obsmat') || isempty(Nstates) || isnan(Nstates)
            error('VALIDATE_VITERBI: Loaded model parameters or Nstates are invalid/empty.');
        end
        model_obsmat_for_comparison = model_params_loaded.obsmat;

        fprintf('Loading Viterbi_decoded_states_Q0.mat...\n');
        if ~exist('Viterbi_decoded_states_Q0.mat', 'file'), error('VALIDATE_VITERBI: Viterbi_decoded_states_Q0.mat not found. Ensure runViterbi_local has run.'); end
        load('Viterbi_decoded_states_Q0.mat', 'Q0');
        if isempty(Q0), error('VALIDATE_VITERBI: Decoded states (Q0) are empty.'); end

        fprintf('Loading data_re_for_hmm.mat...\n');
        if ~exist('data_re_for_hmm.mat', 'file'), error('VALIDATE_VITERBI: data_re_for_hmm.mat not found.'); end
        load('data_re_for_hmm.mat', 'data_re_for_hmm');
        if isempty(data_re_for_hmm), error('VALIDATE_VITERBI: Trajectory data (data_re_for_hmm) is empty.'); end

        if length(Q0) ~= length(data_re_for_hmm)
            error('VALIDATE_VITERBI: Mismatch between number of trajectories in Q0 (%d) and data_re_for_hmm (%d).', length(Q0), length(data_re_for_hmm));
        end


        %% 2. User Input (Minimum Photons per Segment)
        % ... (This section remains IDENTICAL to the previous version) ...
        fprintf('\n--- 2. Getting Parameters for Segment FRET Calculation ---\n');
        prompt_min_photons = {'Enter minimum photons per segment to calculate FRET (e.g., 5 or 10):'};
        dlgtitle_min_photons = 'Segment FRET Parameter'; definput_min_photons = {'5'};
        answer_min_photons = inputdlg(prompt_min_photons, dlgtitle_min_photons, [1 60], definput_min_photons);
        if isempty(answer_min_photons), disp('User cancelled. Exiting.'); cd(original_pwd); return; end
        try
            min_photons_per_segment_for_fret = str2double(answer_min_photons{1});
            if isnan(min_photons_per_segment_for_fret) || min_photons_per_segment_for_fret < 1, error('Invalid minimum photons (must be >= 1).'); end
        catch ME_minp_input, cd(original_pwd); error('Error parsing minimum photons: %s', ME_minp_input.message); end
        fprintf('INFO: Minimum photons per segment for FRET calculation: %d\n', min_photons_per_segment_for_fret);


        %% 3. Process Trajectories: Segment FRET and Count Transitions
        fprintf('\n--- 3. Segmenting Trajectories, Calculating FRET, and Counting Transitions ---\n');
        fret_values_per_segment_by_state = cell(1, Nstates);
        for k_st_init = 1:Nstates, fret_values_per_segment_by_state{k_st_init} = []; end

        % Initialize transition counters (assuming Nstates=2 for specific 1->2, 2->1)
        % For general Nstates, you would use a transition_counts_matrix(Nstates, Nstates)
        total_transitions_1_to_2 = 0;
        total_transitions_2_to_1 = 0;
        % You could also count total time spent in each state if desired
        % total_photons_in_state = zeros(1, Nstates);

        num_trajectories = length(data_re_for_hmm);
        h_wait_val = waitbar(0, 'Processing trajectories for FRET & Transitions...');
        cleanupObj_waitbar_val = onCleanup(@() closeValidationWaitbar(h_wait_val));

        for i_traj = 1:num_trajectories
            if mod(i_traj, round(num_trajectories/100)+1)==0 || i_traj==1 || i_traj==num_trajectories
                waitbar(i_traj/num_trajectories, h_wait_val, sprintf('Trajectory %d/%d', i_traj, num_trajectories));
            end
            
            trajectory_photon_symbols = data_re_for_hmm{i_traj}(:,2);
            viterbi_path = Q0{i_traj};

            if isempty(trajectory_photon_symbols) || isempty(viterbi_path) || ...
               length(trajectory_photon_symbols) ~= length(viterbi_path) || ...
               length(trajectory_photon_symbols) < 1
                continue;
            end

            % --- Calculate FRET for Segments ---
            path_changes_indices = find([true; diff(viterbi_path(:)) ~= 0; true]);
            for k_segment = 1:(length(path_changes_indices)-1)
                segment_start_photon_idx = path_changes_indices(k_segment);
                segment_end_photon_idx = path_changes_indices(k_segment+1) - 1;
                current_segment_state = viterbi_path(segment_start_photon_idx);
                segment_photon_symbols = trajectory_photon_symbols(segment_start_photon_idx:segment_end_photon_idx);
                num_photons_in_segment = length(segment_photon_symbols);

                if num_photons_in_segment >= min_photons_per_segment_for_fret
                    donor_photons_in_segment = sum(segment_photon_symbols == 1);
                    acceptor_photons_in_segment = sum(segment_photon_symbols == 2);
                    total_DA_in_segment = donor_photons_in_segment + acceptor_photons_in_segment;
                    if total_DA_in_segment > 0
                        fret_for_segment = acceptor_photons_in_segment / total_DA_in_segment;
                        if current_segment_state >= 1 && current_segment_state <= Nstates
                            fret_values_per_segment_by_state{current_segment_state}(end+1) = fret_for_segment; %#ok<AGROW>
                        else
                            warning('VALIDATE_VITERBI: Traj %d, invalid Viterbi state %d.', i_traj, current_segment_state);
                        end
                    end
                end
            end
            
            % --- Count Transitions (for Nstates=2, specifically 1->2 and 2->1) ---
            if Nstates == 2 && length(viterbi_path) > 1
                for t_photon = 1:(length(viterbi_path)-1)
                    current_state = viterbi_path(t_photon);
                    next_state = viterbi_path(t_photon+1);
                    if current_state == 1 && next_state == 2
                        total_transitions_1_to_2 = total_transitions_1_to_2 + 1;
                    elseif current_state == 2 && next_state == 1
                        total_transitions_2_to_1 = total_transitions_2_to_1 + 1;
                    end
                end
            end
            % For general Nstates, you'd update a matrix:
            % if Nstates > 1 && length(viterbi_path) > 1
            %     for t_photon = 1:(length(viterbi_path)-1)
            %         current_s = viterbi_path(t_photon);
            %         next_s = viterbi_path(t_photon+1);
            %         if current_s ~= next_s && current_s >=1 && current_s <=Nstates && next_s >=1 && next_s <=Nstates
            %             transition_counts_matrix(current_s, next_s) = transition_counts_matrix(current_s, next_s) + 1;
            %         end
            %     end
            % end
        end % End of trajectory loop
        
        fprintf('\n--- Transition Counts (for Nstates=2 model) ---\n');
        if Nstates == 2
            fprintf('Total Transitions State 1 -> State 2: %d\n', total_transitions_1_to_2);
            fprintf('Total Transitions State 2 -> State 1: %d\n', total_transitions_2_to_1);
        else
            fprintf('Transition counting for 1->2 and 2->1 is specific to Nstates=2.\n');
            % If you implement general transition_counts_matrix, display it here.
        end


        %% 4. Plot Overlaid Histograms of Per-Segment FRET Values
        % ... (This section for plotting histograms remains IDENTICAL to the previous version,
        %      including the legend fix with \\pm) ...
        fprintf('\n--- 4. Plotting Per-Segment FRET Histograms ---\n');
        fig_segment_fret_hist = figure('Name', 'Per-Segment FRET by Viterbi State', 'NumberTitle','off', 'Color', 'w', 'Position', [200, 200, 800, 600]);
        ax_hist = axes(fig_segment_fret_hist);
        hold(ax_hist, 'on');

        hist_edges_fret_segment = -0.025:0.025:1.025;
        plot_colors = lines(Nstates);
        legend_handles = [];
        legend_texts = {};

        for current_state_idx = 1:Nstates
            fret_values_for_this_state = fret_values_per_segment_by_state{current_state_idx};
            base_legend_text = ''; 

            if ~isempty(fret_values_for_this_state)
                h_bar = histogram(ax_hist, fret_values_for_this_state, hist_edges_fret_segment, 'Normalization', 'pdf', ...
                                  'FaceColor', plot_colors(current_state_idx,:), 'EdgeColor', 'k', 'FaceAlpha', 0.65);
                legend_handles(end+1) = h_bar; %#ok<AGROW>
                avg_fret = mean(fret_values_for_this_state);
                std_fret = std(fret_values_for_this_state);
                base_legend_text = sprintf('State %d Segments (N=%d, E=%.2f \\pm %.2f)', ...
                                            current_state_idx, length(fret_values_for_this_state), avg_fret, std_fret);
                fprintf('HMM State %d: Mean Segment FRET = %.3f, Std = %.3f, N_segments = %d\n', ...
                        current_state_idx, avg_fret, std_fret, length(fret_values_for_this_state));
            else
                fprintf('HMM State %d: No segments with >= %d photons found.\n', current_state_idx, min_photons_per_segment_for_fret);
                h_dummy = plot(ax_hist, NaN, NaN, 'Color', plot_colors(current_state_idx,:), 'LineWidth', 2);
                legend_handles(end+1) = h_dummy;
                base_legend_text = sprintf('State %d (No valid segments)', current_state_idx);
            end
            
            if current_state_idx <= size(model_obsmat_for_comparison, 1) && size(model_obsmat_for_comparison,2) >= 2
                model_char_value = model_obsmat_for_comparison(current_state_idx, 2);
                plot(ax_hist, [model_char_value, model_char_value], get(ax_hist,'YLim'), ...
                     'LineStyle', '--', 'Color', plot_colors(current_state_idx,:)*0.6, 'LineWidth', 2.5, ...
                     'HandleVisibility', 'off');
                base_legend_text = [base_legend_text sprintf(', Model P(A|S)=%.2f', model_char_value)];
            end
            legend_texts{end+1} = base_legend_text; %#ok<AGROW>
        end

        hold(ax_hist, 'off');
        xlabel(ax_hist, 'Segment FRET Efficiency (Calculated from D/A Symbols)');
        ylabel(ax_hist, 'Probability Density');
        xlim(ax_hist, [-0.05, 1.05]);
                % ... (after hold(ax_hist, 'off'); and xlabel, ylabel, xlim) ...

        % Prepare the base title string
        base_title_str = sprintf('Per-Segment FRET by Viterbi State (Model #%s, Min %d photons/seg, HMM dt=%.1e s)', ...
                                 chosen_model_idx_str, min_photons_per_segment_for_fret, dt_model_sec);
        
        % Add transition counts to the title if Nstates == 2
        if Nstates == 2
            transitions_title_str = sprintf('Transitions: S1->S2 = %d, S2->S1 = %d', ...
                                            total_transitions_1_to_2, total_transitions_2_to_1);
            % Combine with a newline: \n
            title_str = {base_title_str, transitions_title_str}; % Use a cell array for multi-line title
        else
            title_str = base_title_str; % Single line title if not 2 states
        end
        
        title(ax_hist, title_str, 'Interpreter', 'none'); % 'Interpreter', 'none' is fine for this content
        
        if ~isempty(legend_handles) && ~isempty(legend_texts)
            legend(ax_hist, legend_handles, legend_texts, 'Location', 'NorthEast', 'Interpreter', 'tex');
        end
        % ...
        grid(ax_hist, 'on'); box(ax_hist, 'on'); set(ax_hist, 'FontSize', 10);

        try
            saveas(fig_segment_fret_hist, 'State_Segment_FRET_Histograms_Overlaid.png');
            % Save transition counts if Nstates == 2
            if Nstates == 2
                save('State_Segment_FRET_Data.mat', 'fret_values_per_segment_by_state', ...
                     'min_photons_per_segment_for_fret', 'Nstates', 'dt_model_sec', ...
                     'model_params_loaded', 'chosen_model_idx_str', ...
                     'total_transitions_1_to_2', 'total_transitions_2_to_1', '-append', '-v7.3'); % Append transition counts
            else
                 save('State_Segment_FRET_Data.mat', 'fret_values_per_segment_by_state', ...
                     'min_photons_per_segment_for_fret', 'Nstates', 'dt_model_sec', ...
                     'model_params_loaded', 'chosen_model_idx_str', '-append', '-v7.3'); % Save without specific transition counts
            end
            fprintf('State-segment FRET histograms and data saved/updated.\n');
        catch ME_save_hist
            warning(ME_save_hist.identifier, 'Could not save FRET segment histogram/data: %s', ME_save_hist.message);
        end


        fprintf('\n--- Viterbi State Validation Finished ---\n');
        cd(original_pwd);
        fprintf('Returned to original directory: %s\n', original_pwd);

    catch ME_validate
        fprintf(2,'ERROR in validateViterbiStates_local: %s\n', ME_validate.message);
        disp(ME_validate.getReport());
        if exist('h_wait_val', 'var') && ishandle(h_wait_val), close(h_wait_val); end
        if exist('original_pwd', 'var') && isfolder(original_pwd) && ~strcmp(pwd, original_pwd)
            cd(original_pwd);
            fprintf('Returned to original directory %s after error.\n', original_pwd);
        end
    end

    function closeValidationWaitbar(h_input)
        if exist('h_input', 'var') && ishandle(h_input), close(h_input); end
    end
end