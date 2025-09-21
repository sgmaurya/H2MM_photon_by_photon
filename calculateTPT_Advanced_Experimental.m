% calculateTPT_Advanced_Experimental.m
%
% This is the final, unified analysis script for experimental data.
% It combines the robust file handling of the original script with the four
% advanced TPT calculation methods developed during the validation project.
%
function calculateTPT_Advanced_Experimental()
    clearvars -except break_debug; close all; clc;
    fprintf('--- Starting Advanced TPT Calculation for Experimental Data ---\n');

    original_pwd = pwd;
    try
        %% 1. Select HMM Run Directory and Load All Necessary Files
        fprintf('\n--- 1. Selecting HMM Run Directory and Loading Data ---\n');
        
        project_base_path = fileparts(mfilename('fullpath'));
        results_base_path = fullfile(project_base_path, 'HMM_Local_Results');
        if ~exist(results_base_path, 'dir'), results_base_path = pwd; end
        
        selected_run_path = uigetdir(results_base_path, 'Select the HMM_Run_* directory from your local results');
        if isequal(selected_run_path, 0), disp('User cancelled. Exiting.'); return; end
        cd(selected_run_path);
        fprintf('Processing results from: %s\n', selected_run_path);

        fprintf('Loading HMM model, Viterbi path, and photon data...\n');
        load('Best_HMM_Model_Local.mat', 'best_model_overall', 'config_params_to_save', 'best_idx');
        load('Viterbi_decoded_states_Q0.mat', 'Q0');
        load('data_re_for_hmm.mat', 'data_re_for_hmm');
        data = data_re_for_hmm;
        model_params = best_model_overall; config_params = config_params_to_save; chosen_model_idx = best_idx;
        
        %% 2. Extract and Prepare Model Parameters
        fprintf('\n--- 2. Preparing HMM Model Parameters ---\n');
        Nstates = config_params.Nstates;
        if Nstates ~= 2, error('This script is designed for 2-state models.'); end
        dt_model_sec = config_params.dt_analysis_sec;
        obsmat = model_params.obsmat;
        [~, sorted_idx] = sort(obsmat(:,2));
        state_map = zeros(1, Nstates);
        state_map(sorted_idx(1)) = 1; state_map(sorted_idx(2)) = 2;
        P_Red_L = obsmat(sorted_idx(1), 2); P_Green_L = obsmat(sorted_idx(1), 1);
        P_Red_H = obsmat(sorted_idx(2), 2); P_Green_H = obsmat(sorted_idx(2), 1);
        fprintf('State Remapping Complete: Low FRET (State 1) P(A)=%.3f, High FRET (State 2) P(A)=%.3f\n', P_Red_L, P_Red_H);

        %% 3. Choose TPT Calculation Method and Parameters
        fprintf('\n--- 3. Choosing TPT Calculation Method ---\n');
        method_options = {'Standard LLR', 'Iterative LLR (Convergent)', 'Standard Bayesian', 'Iterative Bayesian (Convergent)'};
        [indx, tf] = listdlg('ListString', method_options, 'SelectionMode','single', 'Name','Method Selection', 'PromptString','Choose TPT Analysis Method:','ListSize',[250, 100]);
        if ~tf, disp('User cancelled.'); cd(original_pwd); return; end
        analysis_method = method_options{indx};
        params = struct('method', analysis_method);
    
        switch analysis_method
            case 'Standard LLR'
                ans_lambda = inputdlg({'Enter LLR Threshold (Lambda):'}, 'Input', [1 50], {'4.6'});
                if isempty(ans_lambda), cd(original_pwd); return; end
                params.lambda_thresh = str2double(ans_lambda{1});
                title_param_str = sprintf('Lambda=%.2f', params.lambda_thresh);
            case 'Iterative LLR (Convergent)'
                ans_iter = inputdlg({'LLR Threshold (Lambda):', 'Max Iterations:', 'Tolerance (photons):'}, 'Input', [1 60], {'4.6', '10', '2'});
                if isempty(ans_iter), cd(original_pwd); return; end
                params.lambda_thresh = str2double(ans_iter{1}); params.max_iter = str2double(ans_iter{2}); params.conv_tol_photons = str2double(ans_iter{3});
                title_param_str = sprintf('Iter. LLR, Lambda=%.2f', params.lambda_thresh);
            case 'Standard Bayesian'
                ans_prob = inputdlg({'Posterior Probability Threshold (0-1):'}, 'Input', [1 50], {'0.95'});
                if isempty(ans_prob), cd(original_pwd); return; end
                params.prob_thresh = str2double(ans_prob{1});
                title_param_str = sprintf('P_{thresh}=%.2f', params.prob_thresh);
            case 'Iterative Bayesian (Convergent)'
                ans_iter = inputdlg({'Posterior Probability Threshold (0-1):', 'Max Iterations:', 'Tolerance (photons):'}, 'Input', [1 60], {'0.95', '10', '2'});
                if isempty(ans_iter), cd(original_pwd); return; end
                params.prob_thresh = str2double(ans_iter{1}); params.max_iter = str2double(ans_iter{2}); params.conv_tol_photons = str2double(ans_iter{3});
                title_param_str = sprintf('Iter. Bayes, P_{thresh}=%.2f', params.prob_thresh);
        end

        %% 4. Main Calculation Loop
        fprintf('\n--- 4. Calculating TPTs using ''%s'' ---\n', params.method);
        all_TPTs_sec_HL = []; all_TPTs_sec_LH = [];
        num_trajectories = length(data); h_wait = waitbar(0, 'Calculating TPTs...');
        
        for i_traj = 1:num_trajectories
            waitbar(i_traj/num_trajectories, h_wait);
            viterbi_path_orig = Q0{i_traj}; photon_stream = data{i_traj};
            if isempty(viterbi_path_orig) || length(viterbi_path_orig) < 2, continue; end
            viterbi_path_mapped = state_map(viterbi_path_orig);
            transition_indices = find(diff(viterbi_path_mapped) ~= 0) + 1;
            
            for k_trans = 1:length(transition_indices)
                idx_change_viterbi = transition_indices(k_trans);
                state_before = viterbi_path_mapped(idx_change_viterbi - 1); state_after  = viterbi_path_mapped(idx_change_viterbi);
                if state_before == 1 && state_after == 2, type = 'LH'; elseif state_before == 2 && state_after == 1, type = 'HL'; else, continue; end
                
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
                    tpt_sec = (photon_stream(idx_entry, 1) - photon_stream(idx_exit, 1)) * dt_model_sec;
                    if tpt_sec >= 0, if strcmp(type, 'HL'), all_TPTs_sec_HL(end+1) = tpt_sec; else, all_TPTs_sec_LH(end+1) = tpt_sec; end; end
                end
            end
        end
        if ishandle(h_wait), close(h_wait); end

        %% 5. Plotting and Saving Results
        fprintf('\n--- 5. Plotting and Saving TPT Results ---\n');
        if isempty(all_TPTs_sec_HL) && isempty(all_TPTs_sec_LH), msgbox('No valid TPTs found.', 'Warning', 'warn'); cd(original_pwd); return; end
        
        plot_style = questdlg('Choose histogram style:', 'Plot Style', 'Overlaid (Transparent)', 'Outlines (Stairs)', 'Overlaid (Transparent)');
        if isempty(plot_style), cd(original_pwd); return; end
        bin_width_ans = inputdlg('Enter bin width in microseconds:', 'Bin Width', [1 40], {'5'});
        if isempty(bin_width_ans), cd(original_pwd); return; end
        bin_width_us = str2double(bin_width_ans{1});

        fig_tpt_hist = figure('Name', 'Advanced TPT Distribution', 'Color', 'w', 'Position', [300, 300, 900, 600]);
        ax = axes(fig_tpt_hist); hold(ax, 'on');
        mean_tpt_us_HL = mean(all_TPTs_sec_HL)*1e6; mean_tpt_us_LH = mean(all_TPTs_sec_LH)*1e6;
        num_tpts_HL = length(all_TPTs_sec_HL); num_tpts_LH = length(all_TPTs_sec_LH);
        color_HL = [0.2 0.5 0.8]; color_LH = [0.9 0.4 0.1];

        switch plot_style
            case 'Overlaid (Transparent)', if ~isempty(all_TPTs_sec_HL), histogram(ax, all_TPTs_sec_HL*1e6, 'BinWidth', bin_width_us, 'Normalization', 'pdf', 'FaceColor', color_HL, 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'DisplayName', 'High -> Low'); end; if ~isempty(all_TPTs_sec_LH), histogram(ax, all_TPTs_sec_LH*1e6, 'BinWidth', bin_width_us, 'Normalization', 'pdf', 'FaceColor', color_LH, 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'DisplayName', 'Low -> High'); end
            case 'Outlines (Stairs)', if ~isempty(all_TPTs_sec_HL), histogram(ax, all_TPTs_sec_HL*1e6, 'BinWidth', bin_width_us, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', color_HL, 'LineWidth', 2.5, 'DisplayName', 'High -> Low'); end; if ~isempty(all_TPTs_sec_LH), histogram(ax, all_TPTs_sec_LH*1e6, 'BinWidth', bin_width_us, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', color_LH, 'LineWidth', 2.5, 'DisplayName', 'Low -> High'); end
        end
        if ~isempty(all_TPTs_sec_HL), xline(ax, mean_tpt_us_HL, '--', 'Color', color_HL, 'LineWidth', 2, 'DisplayName', sprintf('Mean H->L = %.2f us', mean_tpt_us_HL), 'Label', 'Mean H->L'); end
        if ~isempty(all_TPTs_sec_LH), xline(ax, mean_tpt_us_LH, '--', 'Color', color_LH, 'LineWidth', 2, 'DisplayName', sprintf('Mean L->H = %.2f us', mean_tpt_us_LH), 'Label', 'Mean L->H'); end
        hold(ax, 'off');
        
        xlabel(ax, 'Transition Path Time (\mus)'); ylabel(ax, 'Probability Density');
        grid(ax, 'on'); box(ax, 'on'); set(ax, 'FontSize', 11); legend(ax, 'Location', 'NorthEast', 'FontSize', 10);
        title_str1 = sprintf('TPT from ''%s'' (Model #%d, %s)', params.method, chosen_model_idx, title_param_str);
        title_str2 = sprintf('N_{H->L} = %d | N_{L->H} = %d', num_tpts_HL, num_tpts_LH);
        title(ax, {title_str1, title_str2}, 'Interpreter', 'none');

        fprintf('Saving TPT histogram figure and data...\n');
        safe_method_name = strrep(params.method, ' ', '_'); safe_method_name = strrep(safe_method_name, '(', ''); safe_method_name = strrep(safe_method_name, ')', '');
        saveas(fig_tpt_hist, sprintf('TPT_Distribution_%s.png', safe_method_name));
        save_struct.TPTs_sec_HL = all_TPTs_sec_HL; save_struct.TPTs_sec_LH = all_TPTs_sec_LH;
        save_struct.params = params;
        save(sprintf('TPT_Results_%s.mat', safe_method_name), '-struct', 'save_struct', '-v7.3');
        fprintf('Results saved to .mat and .png files with method name.\n');
        
    catch ME_TPT_Adv
        fprintf(2,'\nERROR in calculateTPT_Advanced_Experimental: %s\n', ME_TPT_Adv.message);
        disp(ME_TPT_Adv.getReport());
        if exist('h_wait', 'var') && ishandle(h_wait), close(h_wait); end
    end
    cd(original_pwd);
end

%% --- SUB-FUNCTIONS for TPT calculation ---
function idx = findEntryPoint_Standard(colors, start_idx, pR_H, pR_L, lambda, type)
LLR=0; idx=NaN; if strcmp(type,'LH'), lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); else, lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); end; for i=start_idx:length(colors), if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end; if LLR>lambda, idx=i; return; end; end; end
function idx = findExitPoint_Standard(colors, start_idx, pR_H, pR_L, lambda, type)
LLR=0; idx=NaN; if strcmp(type,'LH'), lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); else, lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); end; for i=(start_idx-1):-1:1, if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end; if LLR>lambda, idx=i+1; return; end; end; end
function idx = findEntryPoint_Bayesian(colors, start_idx, pR_H, pR_L, p_thresh, type)
LLR=0; idx=NaN; if strcmp(type,'LH'), lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); else, lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); end; for i=start_idx:length(colors), if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end; prob = 1 / (1 + exp(-LLR)); if prob > p_thresh, idx=i; return; end; end; end
function idx = findExitPoint_Bayesian(colors, start_idx, pR_H, pR_L, p_thresh, type)
LLR=0; idx=NaN; if strcmp(type,'LH'), lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); else, lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); end; for i=(start_idx-1):-1:1, if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end; prob = 1 / (1 + exp(-LLR)); if prob > p_thresh, idx=i+1; return; end; end; end