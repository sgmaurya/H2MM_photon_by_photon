% runTPT_and_VisualizeBinned_Advanced.m
%
% This is the main script to run. It performs three actions:
%   1. Asks the user to select a TPT method and its parameters.
%   2. Calculates TPT boundaries for all trajectories using the chosen method.
%   3. Launches an interactive viewer to display the binned data with the
%      calculated TPTs as shaded regions.
%
function runTPT_and_VisualizeBinned_Advanced()
    clearvars -except break_debug; close all; clc;
    fprintf('--- Advanced TPT Visualization Launcher ---\n');

    original_pwd = pwd;
    try
        %% 1. Select HMM Run Directory and Load All Data
        fprintf('\n--- 1. Selecting HMM Run Directory ---\n');
        
        selected_run_path = uigetdir(pwd, 'Select the HMM_Run_* directory with TPT results');
        if isequal(selected_run_path, 0), disp('User cancelled. Exiting.'); return; end
        fprintf('Processing results from: %s\n', selected_run_path);
        cd(selected_run_path);

        % --- Load ALL necessary files ---
        fprintf('Loading all necessary result files...\n');
        load('Best_HMM_Model_Local.mat', 'best_model_overall', 'config_params_to_save');
        load('Viterbi_decoded_states_Q0.mat', 'Q0');
        load('data_re_for_hmm.mat', 'data_re_for_hmm');
        load('AA.mat', 'AA');
        load('Viterbi_transition_analysis.mat', 'ind_trans', 'trans_num');
        
        dt_model_sec = config_params_to_save.dt_analysis_sec;
        Nstates = config_params_to_save.Nstates;
        
        %% 2. Choose TPT Method, Parameters, and Calculate Boundaries
        fprintf('\n--- 2. Choosing TPT Method and Calculating Boundaries ---\n');
        
        obsmat = best_model_overall.obsmat;
        [~, sorted_idx] = sort(obsmat(:,2));
        state_map(sorted_idx(1)) = 1; state_map(sorted_idx(2)) = 2;
        P_Red_L = obsmat(sorted_idx(1), 2); P_Red_H = obsmat(sorted_idx(2), 2);
        
        method_options = {'Standard LLR', 'Iterative LLR (Convergent)', 'Standard Bayesian', 'Iterative Bayesian (Convergent)'};
        [indx, tf] = listdlg('ListString', method_options, 'SelectionMode','single', 'Name','Visualization Method', 'PromptString','Choose TPT method to visualize:');
        if ~tf, disp('User cancelled.'); return; end
        analysis_method = method_options{indx};
        params = struct('method', analysis_method);
    
        switch analysis_method
            case 'Standard LLR'
                ans_lambda = inputdlg({'Enter LLR Threshold (Lambda):'}, 'Input', [1 50], {'4.6'});
                if isempty(ans_lambda), return; end
                params.lambda_thresh = str2double(ans_lambda{1});
            case 'Iterative LLR (Convergent)'
                ans_iter = inputdlg({'LLR Threshold (Lambda):', 'Max Iterations:', 'Tolerance (photons):'}, 'Input', [1 60], {'2.09', '10', '2'});
                if isempty(ans_iter), return; end
                params.lambda_thresh = str2double(ans_iter{1}); params.max_iter = str2double(ans_iter{2}); params.conv_tol_photons = str2double(ans_iter{3});
            case 'Standard Bayesian'
                ans_prob = inputdlg({'Posterior Probability Threshold (0-1):'}, 'Input', [1 50], {'0.89'});
                if isempty(ans_prob), return; end
                params.prob_thresh = str2double(ans_prob{1});
            case 'Iterative Bayesian (Convergent)'
                ans_iter = inputdlg({'Posterior Probability Threshold (0-1):', 'Max Iterations:', 'Tolerance (photons):'}, 'Input', [1 60], {'0.89', '10', '2'});
                if isempty(ans_iter), return; end
                params.prob_thresh = str2double(ans_iter{1}); params.max_iter = str2double(ans_iter{2}); params.conv_tol_photons = str2double(ans_iter{3});
        end
        
        num_trajectories = length(data_re_for_hmm);
        tpt_boundaries_per_traj = cell(1, num_trajectories);
        
        h_wait = waitbar(0, 'Calculating TPT boundaries for visualization...');
        for i_traj = 1:num_trajectories
            waitbar(i_traj / num_trajectories, h_wait);
            viterbi_path_orig = Q0{i_traj}; photon_stream = data_re_for_hmm{i_traj};
            if isempty(viterbi_path_orig) || length(viterbi_path_orig) < 2, continue; end
            viterbi_path_mapped = state_map(viterbi_path_orig);
            transition_indices = find(diff(viterbi_path_mapped) ~= 0) + 1;
            boundaries_this_traj = [];
            
            for k_trans = 1:length(transition_indices)
                idx_change_viterbi = transition_indices(k_trans);
                state_before=viterbi_path_mapped(idx_change_viterbi-1); state_after=viterbi_path_mapped(idx_change_viterbi);
                if state_before == 1 && state_after == 2, type = 'LH';
                elseif state_before == 2 && state_after == 1, type = 'HL';
                else, continue; end
                
                photon_colors = photon_stream(:,2); [idx_exit, idx_entry] = deal(NaN);

                switch params.method
                    case 'Standard LLR'
                        idx_exit = findExitPoint_Standard(photon_colors, idx_change_viterbi, P_Red_H, P_Red_L, params.lambda_thresh, type);
                        idx_entry = findEntryPoint_Standard(photon_colors, idx_change_viterbi, P_Red_H, P_Red_L, params.lambda_thresh, type);
                    case 'Iterative LLR (Convergent)'
                        current_center_idx = idx_change_viterbi;
                        for iter = 1:params.max_iter
                            exit_pass = findExitPoint_Standard(photon_colors, current_center_idx, P_Red_H, P_Red_L, params.lambda_thresh, type);
                            entry_pass = findEntryPoint_Standard(photon_colors, current_center_idx, P_Red_H, P_Red_L, params.lambda_thresh, type);
                            if isnan(exit_pass) || isnan(entry_pass), break; end
                            new_center_idx = round((exit_pass + entry_pass) / 2);
                            if abs(new_center_idx - current_center_idx) <= params.conv_tol_photons, idx_exit = exit_pass; idx_entry = entry_pass; break; end
                            current_center_idx = new_center_idx;
                        end
                        if isnan(idx_exit), idx_exit = exit_pass; idx_entry = entry_pass; end
                    case 'Standard Bayesian'
                        idx_exit = findExitPoint_Bayesian(photon_colors, idx_change_viterbi, P_Red_H, P_Red_L, params.prob_thresh, type);
                        idx_entry = findEntryPoint_Bayesian(photon_colors, idx_change_viterbi, P_Red_H, P_Red_L, params.prob_thresh, type);
                    case 'Iterative Bayesian (Convergent)'
                        current_center_idx = idx_change_viterbi;
                        for iter = 1:params.max_iter
                            exit_pass = findExitPoint_Bayesian(photon_colors, current_center_idx, P_Red_H, P_Red_L, params.prob_thresh, type);
                            entry_pass = findEntryPoint_Bayesian(photon_colors, current_center_idx, P_Red_H, P_Red_L, params.prob_thresh, type);
                            if isnan(exit_pass) || isnan(entry_pass), break; end
                            new_center_idx = round((exit_pass + entry_pass) / 2);
                            if abs(new_center_idx - current_center_idx) <= params.conv_tol_photons, idx_exit = exit_pass; idx_entry = entry_pass; break; end
                            current_center_idx = new_center_idx;
                        end
                        if isnan(idx_exit), idx_exit = exit_pass; idx_entry = entry_pass; end
                end

                if ~isnan(idx_exit) && ~isnan(idx_entry) && idx_entry > idx_exit
                    boundaries_this_traj(end+1, :) = [idx_exit, idx_entry];
                end
            end
            tpt_boundaries_per_traj{i_traj} = boundaries_this_traj;
        end
        if ishandle(h_wait), close(h_wait); end

        %% 3. Call the Visualization Function
        fprintf('\n--- 3. Launching Interactive Trajectory Viewer ---\n');
        
        data_for_plotting = data_re_for_hmm;
        for i = 1:length(data_re_for_hmm)
            if ~isempty(data_re_for_hmm{i})
                data_for_plotting{i}(:,1) = data_re_for_hmm{i}(:,1) * dt_model_sec;
            end
        end

        Traj_Vet_local_with_TPT(data_for_plotting, Q0, ind_trans, trans_num, AA, dt_model_sec, Nstates, tpt_boundaries_per_traj);

        fprintf('\n--- Visualization Finished ---\n');
        cd(original_pwd);

    catch ME_TPT_Vis
        fprintf(2,'\nERROR in runTPT_and_VisualizeBinned_Advanced: %s\n', ME_TPT_Vis.message);
        disp(ME_TPT_Vis.getReport());
        if exist('h_wait', 'var') && ishandle(h_wait), close(h_wait); end
        cd(original_pwd);
    end
end

% --- START OF MODIFIED VIEWER FUNCTION ---
function Traj_Vet_local_with_TPT(data_tv, Q0_tv, ind_trans_tv, trans_num_tv, AA_tv, dt_tv, Nstates_tv, tpt_boundaries_per_traj)
    
    prompt_tv = {'Enter time resolution for binning traces (us):', ...
                 'Enter maximum number of trajectories to display:', ...
                 'Enter sorting option (1=Transitions, 2=Flux, 3=None, 4=Length, 5=Specific State):'};
    dlgtitle_tv = 'Trajectory Visualization Input';
    definput_tv = {'5', '10', '1'}; 
    answer_tv = inputdlg(prompt_tv, dlgtitle_tv, [1 60], definput_tv);
    if isempty(answer_tv), disp('Trajectory visualization cancelled by user.'); return; end
    Time_Res_us = str2double(answer_tv{1});
    num_traj_to_display_max = str2double(answer_tv{2});
    sort_option = str2double(answer_tv{3});
    
    num_total_trajectories = length(data_tv);
    switch sort_option
        case 1, [~, sorted_indices] = sort(trans_num_tv, 'descend');
        case 2
            Av_rate_of_photons = cellfun(@(x) size(x,1)/((x(end,1)-x(1,1))+dt_tv), data_tv);
            [~, sorted_indices] = sort(Av_rate_of_photons, 'descend');
        case 3, sorted_indices = 1:num_total_trajectories;
        case 4
            burst_length_sec = cellfun(@(x) (x(end,1)-x(1,1))+dt_tv, data_tv);
            [~, sorted_indices] = sort(burst_length_sec, 'descend');
        case 5
            st_ans = inputdlg('Enter HMM state(s) to select (e.g., 1 or 1,2):', 'Select State(s)', [1 50], {'1'});
            if isempty(st_ans), return; end
            states_to_find = str2num(st_ans{1}); %#ok<ST2NM>
            has_state = cellfun(@(q) any(ismember(unique(q), states_to_find)), Q0_tv);
            sorted_indices = find(has_state);
        otherwise, sorted_indices = 1:num_total_trajectories;
    end
    Trans_list_indices = sorted_indices;

    % <<< FIX: Check if any trajectories were found BEFORE creating the figure >>>
    if isempty(Trans_list_indices)
        msgbox('No trajectories match the selected sorting criteria. Please try another option.', 'No Trajectories Found', 'warn');
        fprintf('\nWARNING: No trajectories found matching the sort criteria. Visualization cancelled.\n');
        return; % Exit the function gracefully
    end
    
    num_traj_to_plot_actual = min(num_traj_to_display_max, length(Trans_list_indices));
    fprintf('Found %d trajectories to display. Showing top %d.\n', length(Trans_list_indices), num_traj_to_plot_actual);
    
    fig_traj_vet = figure('Name', 'Trajectory Visualization with Viterbi & TPTs', 'Color', 'w', 'Position', [50, 50, 1200, 800]);

    for i_plot_loop = 1:num_traj_to_plot_actual
        clf(fig_traj_vet, 'reset');
        idx_in_data = Trans_list_indices(i_plot_loop);
        
        current_trajectory_photons = data_tv{idx_in_data};
        if isempty(current_trajectory_photons), continue; end
        Time_axis_photons_sec = current_trajectory_photons(:,1); 
        Time_axis_photons_usec = Time_axis_photons_sec * 1e6;
        Photon_Symbols = current_trajectory_photons(:,2);
        
        X_bins_usec = Time_axis_photons_usec(1) : Time_Res_us : Time_axis_photons_usec(end);
        if isempty(X_bins_usec) || length(X_bins_usec)<2, X_bins_usec = [Time_axis_photons_usec(1), Time_axis_photons_usec(1)+Time_Res_us]; end
        if X_bins_usec(end) < Time_axis_photons_usec(end), X_bins_usec(end+1)=X_bins_usec(end)+Time_Res_us; end
        X_bins_plot_centers = X_bins_usec(1:end-1) + Time_Res_us/2;
        Donor_binned = histcounts(Time_axis_photons_usec(Photon_Symbols == 1), X_bins_usec);
        Acceptor_binned = histcounts(Time_axis_photons_usec(Photon_Symbols == 2), X_bins_usec);
        Total_photons_binned = Donor_binned + Acceptor_binned;
        FRET_binned = Acceptor_binned ./ Total_photons_binned;
        FRET_binned(isnan(FRET_binned) | isinf(FRET_binned)) = 0;

        ax1 = subplot(2,1,1);
        hold(ax1, 'on');
        legend_handles = []; legend_labels = {};
        
        boundaries = tpt_boundaries_per_traj{idx_in_data};
        if ~isempty(boundaries)
            t_start_us = data_tv{idx_in_data}(boundaries(:, 1), 1) * 1e6;
            t_end_us   = data_tv{idx_in_data}(boundaries(:, 2), 1) * 1e6;
            for k=1:length(t_start_us)
                h_patch = patch(ax1, [t_start_us(k), t_end_us(k), t_end_us(k), t_start_us(k)], [-0.1, -0.1, 1.1, 1.1], ...
                      [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                if k == 1, h_tpt = h_patch; end
            end
            legend_handles(end+1) = h_tpt; legend_labels{end+1} = 'TPT Region';
        end
        
        h_fret = plot(ax1, X_bins_plot_centers, FRET_binned, '-b');
        h_viterbi = stairs(ax1, Time_axis_photons_usec, AA_tv(2, Q0_tv{idx_in_data}), '-k', 'LineWidth', 2);
        legend_handles = [legend_handles, h_fret, h_viterbi];
        legend_labels = [legend_labels, 'Binned FRET', 'Viterbi Path'];

        if ~isempty(ind_trans_tv{idx_in_data}), plot(ax1, [Time_axis_photons_usec(ind_trans_tv{idx_in_data})'; Time_axis_photons_usec(ind_trans_tv{idx_in_data})'], [-0.1; 1.1], '--m', 'HandleVisibility', 'off'); end
        hold(ax1, 'off'); grid on; ylim([-0.1, 1.1]);
        ylabel(ax1, 'FRET Efficiency / P(A|State)');
        title(ax1, sprintf('Trajectory #%d (Original Index: %d)', i_plot_loop, idx_in_data));
        legend(ax1, legend_handles, legend_labels, 'Location', 'best');

        ax2 = subplot(2,1,2);
        hold(ax2, 'on');
        if ~isempty(boundaries)
            t_start_us = data_tv{idx_in_data}(boundaries(:, 1), 1) * 1e6;
            t_end_us   = data_tv{idx_in_data}(boundaries(:, 2), 1) * 1e6;
            y_max = max(Total_photons_binned) * 1.1; if y_max < 5, y_max=5; end
            for k=1:length(t_start_us)
                patch(ax2, [t_start_us(k), t_end_us(k), t_end_us(k), t_start_us(k)], [0, 0, y_max, y_max], ...
                      [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
        end
        plot(ax2, X_bins_plot_centers, Donor_binned, '-g', 'DisplayName', 'Donor');
        plot(ax2, X_bins_plot_centers, Acceptor_binned, '-r', 'DisplayName', 'Acceptor');
        plot(ax2, X_bins_plot_centers, Total_photons_binned, '--c', 'DisplayName', 'Total');
        hold(ax2, 'off'); grid on;
        y_max_plot = max(Total_photons_binned) * 1.15; if y_max_plot < 10, y_max_plot = 10; end
        ylim(ax2, [0, y_max_plot]);
        
        ylabel(ax2, sprintf('Photons / %g us', Time_Res_us));
        xlabel(ax2, 'Time (\mus)');
        legend(ax2, 'Location','best');
        linkaxes([ax1, ax2], 'x'); xlim(ax1, [Time_axis_photons_usec(1), Time_axis_photons_usec(end)]);

        sgtitle(fig_traj_vet, sprintf('Displaying Trajectory %d of %d (from sorted list)', i_plot_loop, num_traj_to_plot_actual), 'FontSize', 10, 'FontWeight', 'bold');
        if i_plot_loop < num_traj_to_plot_actual, fprintf('Displaying traj %d. Press key to continue...\n', i_plot_loop); pause; end
    end
    fprintf('Displayed last selected trajectory.\n');
end

%% --- ALL TPT SUB-FUNCTIONS ---
% (Identical to the functions in calculateTPT_Advanced_Experimental.m)
function idx = findEntryPoint_Standard(colors, start_idx, pR_H, pR_L, lambda, type)
LLR=0; idx=NaN; if strcmp(type,'LH'), lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); else, lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); end; for i=start_idx:length(colors), if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end; if LLR>lambda, idx=i; return; end; end; end
function idx = findExitPoint_Standard(colors, start_idx, pR_H, pR_L, lambda, type)
LLR=0; idx=NaN; if strcmp(type,'LH'), lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); else, lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); end; for i=(start_idx-1):-1:1, if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end; if LLR>lambda, idx=i+1; return; end; end; end
function idx = findEntryPoint_Bayesian(colors, start_idx, pR_H, pR_L, p_thresh, type)
LLR=0; idx=NaN; if strcmp(type,'LH'), lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); else, lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); end; for i=start_idx:length(colors), if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end; prob = 1 / (1 + exp(-LLR)); if prob > p_thresh, idx=i; return; end; end; end
function idx = findExitPoint_Bayesian(colors, start_idx, pR_H, pR_L, p_thresh, type)
LLR=0; idx=NaN; if strcmp(type,'LH'), lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); else, lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); end; for i=(start_idx-1):-1:1, if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end; prob = 1 / (1 + exp(-LLR)); if prob > p_thresh, idx=i+1; return; end; end; end