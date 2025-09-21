% calculateTPT_Advanced_Cluster.m
%
% This script extends the TPT calculation framework to large, cluster-analyzed
% datasets. It combines the data loading structure from 'calculateTPT_cluster.m'
% with the four advanced analysis methods implemented in
% 'calculateTPT_Advanced_Experimental.m'.
%
% This script serves as the final analysis step for calculating TPTs from
% cluster-based HMM results, offering a choice between four statistically
% distinct methods to measure the period of ambiguity during molecular
% transitions.
%
% Methods available:
%   1. Standard LLR: Viterbi-centered search with a log-likelihood ratio threshold.
%   2. Iterative LLR: A data-driven method that removes Viterbi bias by
%      iteratively re-calculating the transition midpoint.
%   3. Standard Bayesian: Viterbi-centered search with a posterior probability threshold.
%   4. Iterative Bayesian: A robust, data-driven method combining the iterative
%      midpoint search with a Bayesian stopping condition.
%

function calculateTPT_Advanced_Cluster()
    clearvars -except break_debug; close all; clc;
    fprintf('--- Starting Advanced TPT Calculation for Cluster HMM Results ---\n');

    original_pwd = pwd;
    try
        %% 1. Select Processed Cluster Run Directory and Load All Data
        fprintf('\n--- 1. Selecting Processed Cluster HMM Run Directory ---\n');
        selected_run_path = uigetdir(pwd, 'Select the Cluster HMM Run directory with Viterbi results');
        if isequal(selected_run_path, 0), disp('User cancelled. Exiting.'); return; end
        cd(selected_run_path);
        fprintf('Processing results from: %s\n', selected_run_path);

        % --- Load Processed HMM Parameters from Cluster Output ---
        fprintf('Loading HMM parameters, Viterbi paths, and raw data...\n');
        if ~exist('HMM_parm.mat', 'file'), error('CALC_TPT_ADV_CLUSTER: HMM_parm.mat not found.'); end
        load('HMM_parm.mat', 'Obs', 'Trans', 'K'); 
        
        if ~exist('Est.mat', 'file'), error('CALC_TPT_ADV_CLUSTER: Est.mat not found.'); end
        load('Est.mat', 'ass', 'B');
        chosen_model_idx = ass;

        % --- Load Viterbi Paths and Raw Photon Data ---
        if ~exist('Viterbi_decoded_states_Q0_cluster.mat', 'file'), error('CALC_TPT_ADV_CLUSTER: Viterbi results file not found.'); end
        load('Viterbi_decoded_states_Q0_cluster.mat', 'Q0');
        
        if ~exist('data.mat', 'file'), error('CALC_TPT_ADV_CLUSTER: data.mat not found.'); end
        load('data.mat', 'data');
        
        %% 2. Extract and Prepare Model Parameters
        fprintf('\n--- 2. Preparing HMM Model Parameters ---\n');
        Nstates = size(Trans, 1);
        if Nstates ~= 2, error('TPT calculation is configured for 2-state models only.'); end
        
        % Infer dt from the model kinetics, as in the original cluster script
        [~, max_idx] = max(abs(K(:))); [row, col] = ind2sub(size(K), max_idx);
        if K(row, col)~=0 && Trans(row, col)~=0
            dt_model_sec = Trans(row,col)/K(row,col);
            fprintf('INFO: Inferred dt = %.2e seconds from model kinetics.\n', dt_model_sec);
        else
            error('Cannot infer dt from HMM_parm.mat. Check model files.'); 
        end
        
        % --- Map states to Low-FRET (1) and High-FRET (2) ---
        % The cluster pipeline saves Obs as [P(Red|L); P(Red|H)] and B as the
        % original state indices corresponding to [L, H].
        P_Red_L = Obs(1); P_Green_L = 1 - P_Red_L;
        P_Red_H = Obs(2); P_Green_H = 1 - P_Red_H;
        
        state_map = zeros(1, max(B)); % Ensure state_map is large enough
        state_map(B(1)) = 1; % Map original low-FRET state index to 1
        state_map(B(2)) = 2; % Map original high-FRET state index to 2
        
        fprintf('State Remapping Complete:\n');
        fprintf('  - Original HMM state %d -> Low FRET (State 1) with P(Acceptor)=%.3f\n', B(1), P_Red_L);
        fprintf('  - Original HMM state %d -> High FRET (State 2) with P(Acceptor)=%.3f\n', B(2), P_Red_H);
        
        %% 3. Choose TPT Calculation Method and Parameters
        fprintf('\n--- 3. Choosing TPT Calculation Method ---\n');
        method_options = {'Standard LLR', 'Iterative LLR (Hazy Viterbi)', 'Standard Bayesian', 'Iterative Bayesian (Hazy Viterbi)'};
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
            case 'Iterative LLR (Hazy Viterbi)'
                ans_iter = inputdlg({'LLR Threshold (Lambda):', 'Max Iterations:', 'Tolerance (photons):'}, 'Input', [1 60], {'4.6', '10', '2'});
                if isempty(ans_iter), cd(original_pwd); return; end
                params.lambda_thresh = str2double(ans_iter{1}); params.max_iter = str2double(ans_iter{2}); params.conv_tol_photons = str2double(ans_iter{3});
                title_param_str = sprintf('Iter. LLR, Lambda=%.2f', params.lambda_thresh);
            case 'Standard Bayesian'
                ans_prob = inputdlg({'Posterior Probability Threshold (0-1):'}, 'Input', [1 50], {'0.99'});
                if isempty(ans_prob), cd(original_pwd); return; end
                params.prob_thresh = str2double(ans_prob{1});
                title_param_str = sprintf('P_{thresh}=%.2f', params.prob_thresh);
            case 'Iterative Bayesian (Hazy Viterbi)'
                ans_iter = inputdlg({'Posterior Probability Threshold (0-1):', 'Max Iterations:', 'Tolerance (photons):'}, 'Input', [1 60], {'0.99', '10', '2'});
                if isempty(ans_iter), cd(original_pwd); return; end
                params.prob_thresh = str2double(ans_iter{1}); params.max_iter = str2double(ans_iter{2}); params.conv_tol_photons = str2double(ans_iter{3});
                title_param_str = sprintf('Iter. Bayes, P_{thresh}=%.2f', params.prob_thresh);
        end

        %% 4. Main Calculation Loop
        fprintf('\n--- 4. Calculating TPTs using ''%s'' ---\n', params.method);
        all_TPTs_sec_HL = []; all_TPTs_sec_LH = [];
        num_trajectories = length(data); 
        h_wait = waitbar(0, 'Calculating TPTs...');
        
        for i_traj = 1:num_trajectories
            waitbar(i_traj/num_trajectories, h_wait, sprintf('Processing Trajectory %d/%d', i_traj, num_trajectories));
            viterbi_path_orig = Q0{i_traj}; 
            photon_stream = data{i_traj};
            if isempty(viterbi_path_orig) || isempty(photon_stream) || length(viterbi_path_orig) < 2, continue; end
            
            % Map the original Viterbi path (e.g., states 4 and 7) to our analysis
            % path (states 1 and 2 for Low/High FRET)
            viterbi_path_mapped = state_map(viterbi_path_orig);
            transition_indices = find(diff(viterbi_path_mapped) ~= 0) + 1;
            
            for k_trans = 1:length(transition_indices)
                idx_change_viterbi = transition_indices(k_trans);
                
                state_before = viterbi_path_mapped(idx_change_viterbi - 1); 
                state_after  = viterbi_path_mapped(idx_change_viterbi);
                
                if state_before == 1 && state_after == 2
                    type = 'LH'; % Low -> High
                elseif state_before == 2 && state_after == 1
                    type = 'HL'; % High -> Low
                else
                    continue; 
                end
                
                photon_colors = photon_stream(:,2); 
                [idx_exit, idx_entry] = deal(NaN);

                switch params.method
                    case 'Standard LLR'
                        idx_exit = findExitPoint_Standard(photon_colors, idx_change_viterbi, P_Red_H, P_Red_L, params.lambda_thresh, type);
                        idx_entry = findEntryPoint_Standard(photon_colors, idx_change_viterbi, P_Red_H, P_Red_L, params.lambda_thresh, type);
                    
                    case 'Iterative LLR (Hazy Viterbi)'
                        current_center_idx = idx_change_viterbi;
                        exit_pass = NaN; entry_pass = NaN;
                        for iter = 1:params.max_iter
                            exit_pass = findExitPoint_Standard(photon_colors, current_center_idx, P_Red_H, P_Red_L, params.lambda_thresh, type);
                            entry_pass = findEntryPoint_Standard(photon_colors, current_center_idx, P_Red_H, P_Red_L, params.lambda_thresh, type);
                            if isnan(exit_pass) || isnan(entry_pass) || entry_pass <= exit_pass, break; end
                            new_center_idx = round((exit_pass + entry_pass) / 2);
                            if abs(new_center_idx - current_center_idx) <= params.conv_tol_photons
                                idx_exit = exit_pass; idx_entry = entry_pass; 
                                break; 
                            end
                            current_center_idx = new_center_idx;
                        end
                        if isnan(idx_exit), idx_exit = exit_pass; idx_entry = entry_pass; end % Store last valid result if not converged

                    case 'Standard Bayesian'
                        idx_exit = findExitPoint_Bayesian(photon_colors, idx_change_viterbi, P_Red_H, P_Red_L, params.prob_thresh, type);
                        idx_entry = findEntryPoint_Bayesian(photon_colors, idx_change_viterbi, P_Red_H, P_Red_L, params.prob_thresh, type);

                    case 'Iterative Bayesian (Hazy Viterbi)'
                        current_center_idx = idx_change_viterbi;
                        exit_pass = NaN; entry_pass = NaN;
                        for iter = 1:params.max_iter
                            exit_pass = findExitPoint_Bayesian(photon_colors, current_center_idx, P_Red_H, P_Red_L, params.prob_thresh, type);
                            entry_pass = findEntryPoint_Bayesian(photon_colors, current_center_idx, P_Red_H, P_Red_L, params.prob_thresh, type);
                            if isnan(exit_pass) || isnan(entry_pass) || entry_pass <= exit_pass, break; end
                            new_center_idx = round((exit_pass + entry_pass) / 2);
                            if abs(new_center_idx - current_center_idx) <= params.conv_tol_photons
                                idx_exit = exit_pass; idx_entry = entry_pass;
                                break;
                            end
                            current_center_idx = new_center_idx;
                        end
                        if isnan(idx_exit), idx_exit = exit_pass; idx_entry = entry_pass; end % Store last valid result
                end

                if ~isnan(idx_exit) && ~isnan(idx_entry) && idx_entry > idx_exit
                    time_exit_dt_units = photon_stream(idx_exit, 1);
                    time_entry_dt_units = photon_stream(idx_entry, 1);
                    tpt_sec = (time_entry_dt_units - time_exit_dt_units) * dt_model_sec;

                    if tpt_sec >= 0
                        if strcmp(type, 'HL'), all_TPTs_sec_HL(end+1) = tpt_sec; 
                        else, all_TPTs_sec_LH(end+1) = tpt_sec; end
                    end
                end
            end
        end
        if ishandle(h_wait), close(h_wait); end

        %% 5. Plotting and Saving Results
        fprintf('\n--- 5. Plotting and Saving TPT Results ---\n');
        if isempty(all_TPTs_sec_HL) && isempty(all_TPTs_sec_LH)
            msgbox('No valid TPTs found with the selected method and parameters.', 'Warning', 'warn'); 
            cd(original_pwd); 
            return; 
        end
        
        plot_style = questdlg('Choose a histogram style:', 'Plot Style', 'Overlaid (Transparent)', 'Outlines (Stairs)', 'Overlaid (Transparent)');
        if isempty(plot_style), disp('Plotting cancelled.'); cd(original_pwd); return; end
        bin_width_ans = inputdlg('Enter bin width in microseconds (e.g., 5):', 'Bin Width', [1 40], {'5'});
        if isempty(bin_width_ans), disp('Plotting cancelled.'); cd(original_pwd); return; end
        bin_width_us = str2double(bin_width_ans{1});

        fig_tpt_hist = figure('Name', 'Advanced TPT Distribution (Cluster)', 'Color', 'w', 'Position', [300, 300, 900, 600]);
        ax = axes(fig_tpt_hist); hold(ax, 'on');
        
        mean_tpt_us_HL = mean(all_TPTs_sec_HL)*1e6; 
        mean_tpt_us_LH = mean(all_TPTs_sec_LH)*1e6;
        num_tpts_HL = length(all_TPTs_sec_HL); 
        num_tpts_LH = length(all_TPTs_sec_LH);
        color_HL = [0.2 0.5 0.8]; % Blue
        color_LH = [0.9 0.4 0.1]; % Orange

        switch plot_style
            case 'Overlaid (Transparent)'
                if ~isempty(all_TPTs_sec_HL), histogram(ax, all_TPTs_sec_HL*1e6, 'BinWidth', bin_width_us, 'Normalization', 'pdf', 'FaceColor', color_HL, 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'DisplayName', 'High -> Low'); end
                if ~isempty(all_TPTs_sec_LH), histogram(ax, all_TPTs_sec_LH*1e6, 'BinWidth', bin_width_us, 'Normalization', 'pdf', 'FaceColor', color_LH, 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'DisplayName', 'Low -> High'); end
            case 'Outlines (Stairs)'
                if ~isempty(all_TPTs_sec_HL), histogram(ax, all_TPTs_sec_HL*1e6, 'BinWidth', bin_width_us, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', color_HL, 'LineWidth', 2.5, 'DisplayName', 'High -> Low'); end
                if ~isempty(all_TPTs_sec_LH), histogram(ax, all_TPTs_sec_LH*1e6, 'BinWidth', bin_width_us, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', color_LH, 'LineWidth', 2.5, 'DisplayName', 'Low -> High'); end
        end
        
        if ~isempty(all_TPTs_sec_HL) && ~isnan(mean_tpt_us_HL)
             xline(ax, mean_tpt_us_HL, '--', 'Color', color_HL, 'LineWidth', 2, 'DisplayName', sprintf('Mean H->L = %.2f us', mean_tpt_us_HL), 'Label', 'Mean H->L');
        end
        if ~isempty(all_TPTs_sec_LH) && ~isnan(mean_tpt_us_LH)
            xline(ax, mean_tpt_us_LH, '--', 'Color', color_LH, 'LineWidth', 2, 'DisplayName', sprintf('Mean L->H = %.2f us', mean_tpt_us_LH), 'Label', 'Mean L->H');
        end
        hold(ax, 'off');
        
        xlabel(ax, 'Transition Path Time (\mus)'); ylabel(ax, 'Probability Density');
        grid(ax, 'on'); box(ax, 'on'); set(ax, 'FontSize', 11); 
        legend(ax, 'Location', 'NorthEast', 'FontSize', 10);
        
        title_str1 = sprintf('TPT from ''%s'' (Model #%d, %s)', params.method, chosen_model_idx, title_param_str);
        title_str2 = sprintf('N_{H->L} = %d | N_{L->H} = %d', num_tpts_HL, num_tpts_LH);
        title(ax, {title_str1, title_str2}, 'Interpreter', 'none'); % 'none' interpreter is safer

        fprintf('Saving TPT histogram figure and data...\n');
        safe_method_name = strrep(params.method, ' ', '_'); 
        safe_method_name = strrep(safe_method_name, '(', ''); 
        safe_method_name = strrep(safe_method_name, ')', '');
        
        saveas(fig_tpt_hist, sprintf('TPT_Distribution_Cluster_%s.png', safe_method_name));
        
        save_struct.TPTs_sec_HL = all_TPTs_sec_HL; 
        save_struct.TPTs_sec_LH = all_TPTs_sec_LH;
        save_struct.params_used = params;
        save_struct.model_info.chosen_model_index = chosen_model_idx;
        save_struct.model_info.dt_model_sec = dt_model_sec;
        save_struct.model_info.P_Red_H = P_Red_H;
        save_struct.model_info.P_Red_L = P_Red_L;

        save(sprintf('TPT_Results_Cluster_%s.mat', safe_method_name), '-struct', 'save_struct', '-v7.3');
        fprintf('Results saved to .mat and .png files with method name suffix.\n');
        
        fprintf('\n--- Advanced TPT Calculation for Cluster Data Finished ---\n');

    catch ME_TPT_Adv_Cluster
        fprintf(2,'\nERROR in calculateTPT_Advanced_Cluster: %s\n', ME_TPT_Adv_Cluster.message);
        disp(ME_TPT_Adv_Cluster.getReport());
        if exist('h_wait', 'var') && ishandle(h_wait), close(h_wait); end
    end
    cd(original_pwd);
end

%% --- SUB-FUNCTIONS for TPT calculation ---
% These functions calculate the boundaries of the transition path.
% Note: P(Green) is calculated internally as 1-P(Red).

function idx = findEntryPoint_Standard(colors, start_idx, pR_H, pR_L, lambda, type)
    LLR=0; idx=NaN; 
    if strcmp(type,'LH'), lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); 
    else, lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); end
    for i=start_idx:length(colors)
        if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end
        if LLR>lambda, idx=i; return; end
    end
end

function idx = findExitPoint_Standard(colors, start_idx, pR_H, pR_L, lambda, type)
    LLR=0; idx=NaN; 
    if strcmp(type,'LH'), lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); 
    else, lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); end
    for i=(start_idx-1):-1:1
        if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end
        if LLR>lambda, idx=i+1; return; end
    end
end

function idx = findEntryPoint_Bayesian(colors, start_idx, pR_H, pR_L, p_thresh, type)
    LLR=0; idx=NaN; 
    if strcmp(type,'LH'), lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); 
    else, lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); end
    for i=start_idx:length(colors)
        if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end
        prob = 1 / (1 + exp(-LLR)); 
        if prob > p_thresh, idx=i; return; end
    end
end

function idx = findExitPoint_Bayesian(colors, start_idx, pR_H, pR_L, p_thresh, type)
    LLR=0; idx=NaN; 
    if strcmp(type,'LH'), lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); 
    else, lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); end
    for i=(start_idx-1):-1:1
        if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end
        prob = 1 / (1 + exp(-LLR));
        if prob > p_thresh, idx=i+1; return; end
    end
end