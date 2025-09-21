% calculateTPT_cluster.m (v2 - Enhanced Plotting)
%
% This is the third step in the CLUSTER analysis pipeline.
%
% v2 Corrections:
%   - Adds detailed statistics (N counts, mean TPTs) to the histogram title.
%   - Adds vertical lines for the mean TPTs to the plot.
%

function calculateTPT_cluster()
    clearvars -except break_debug; close all; clc;
    fprintf('--- Starting TPT Calculation from Cluster HMM Results ---\n');

    original_pwd = pwd;
    try
        %% 1. Select Processed Cluster Run Directory and Load All Data
        fprintf('\n--- 1. Selecting Processed Cluster HMM Run Directory ---\n');
        selected_run_path = uigetdir(pwd, 'Select the Cluster HMM Run directory with Viterbi results');
        if isequal(selected_run_path, 0), disp('User cancelled. Exiting.'); return; end
        fprintf('Processing results from: %s\n', selected_run_path);
        cd(selected_run_path);

        % --- Load Processed HMM Parameters ---
        if ~exist('HMM_parm.mat', 'file'), error('CALC_TPT_CLUSTER: HMM_parm.mat not found.'); end
        load('HMM_parm.mat', 'Obs', 'Trans', 'K'); 
        Nstates = size(Trans, 1);
        if Nstates ~= 2, error('TPT calculation is configured for 2-state models only.'); end
        
        [~, max_idx] = max(abs(K(:))); [row, col] = ind2sub(size(K), max_idx);
        if K(row, col)~=0 && Trans(row, col)~=0, dt_model_sec=Trans(row,col)/K(row,col);
        else, error('Cannot infer dt from HMM_parm.mat.'); end

        if ~exist('Est.mat', 'file'), error('CALC_TPT_CLUSTER: Est.mat not found.'); end
        load('Est.mat', 'ass', 'B');
        chosen_model_idx = ass;

        % --- Load Viterbi Paths and Raw Data ---
        if ~exist('Viterbi_decoded_states_Q0_cluster.mat', 'file'), error('CALC_TPT_CLUSTER: Viterbi results not found.'); end
        load('Viterbi_decoded_states_Q0_cluster.mat', 'Q0');
        if ~exist('data.mat', 'file'), error('CALC_TPT_CLUSTER: data.mat not found.'); end
        load('data.mat', 'data');
        
        %% 2. User Input for LLR Algorithm and Parameter Setup
        fprintf('\n--- 2. Setting LLR Algorithm Parameters ---\n');
        ans_lambda = inputdlg({'Enter LLR Threshold (Lambda_thresh, e.g., log(100)):'}, 'LLR Parameter', [1 60], {'4.6'});
        if isempty(ans_lambda), disp('Cancelled. Exiting.'); return; end
        lambda_thresh = str2double(ans_lambda{1});

        P_Red_L = Obs(1); P_Green_L = 1 - P_Red_L;
        P_Red_H = Obs(2); P_Green_H = 1 - P_Red_H;
        state_map(B(1)) = 1; state_map(B(2)) = 2;

        %% 3. Main Loop: Iterate Through Trajectories and Transitions
        fprintf('\n--- 3. Calculating TPTs for all Transitions ---\n');
        all_TPTs_sec_HL = []; all_TPTs_sec_LH = [];
        
        num_trajectories = length(data);
        h_wait = waitbar(0, 'Calculating TPTs...');
        
        for i_traj = 1:num_trajectories
            waitbar(i_traj/num_trajectories, h_wait, sprintf('Trajectory %d/%d', i_traj, num_trajectories));
            viterbi_path_orig = Q0{i_traj}; photon_stream = data{i_traj};
            if isempty(viterbi_path_orig) || length(viterbi_path_orig) < 2, continue; end
            viterbi_path_mapped = state_map(viterbi_path_orig);
            transition_indices = find(diff(viterbi_path_mapped) ~= 0) + 1;
            for k_trans = 1:length(transition_indices)
                idx_change=transition_indices(k_trans);
                if viterbi_path_mapped(idx_change-1)==1&&viterbi_path_mapped(idx_change)==2, type='LH';
                elseif viterbi_path_mapped(idx_change-1)==2&&viterbi_path_mapped(idx_change)==1, type='HL';
                else, continue; end
                idx_exit=findExitPoint(photon_stream(:,2),idx_change,P_Red_H,P_Green_H,P_Red_L,P_Green_L,lambda_thresh,type);
                idx_entry=findEntryPoint(photon_stream(:,2),idx_change,P_Red_H,P_Green_H,P_Red_L,P_Green_L,lambda_thresh,type);
                if ~isnan(idx_exit)&&~isnan(idx_entry)&&idx_entry>idx_exit
                    tpt_sec=(photon_stream(idx_entry,1)-photon_stream(idx_exit,1))*dt_model_sec;
                    if tpt_sec>=0
                        if strcmp(type,'HL'), all_TPTs_sec_HL(end+1)=tpt_sec;
                        else, all_TPTs_sec_LH(end+1)=tpt_sec; end
                    end
                end
            end
        end
        if ishandle(h_wait), close(h_wait); end

        %% 4. Plotting and Saving Results (Enhanced)
        fprintf('\n--- 4. Plotting and Saving TPT Results ---\n');
        if isempty(all_TPTs_sec_HL) && isempty(all_TPTs_sec_LH), warning('No valid TPTs found.'); cd(original_pwd); return; end
        
        fig_tpt = figure('Name', 'TPT Distribution (Cluster)', 'Color', 'w', 'Position', [300, 300, 900, 600]);
        ax = axes(fig_tpt);
        hold(ax, 'on');
        bin_width_us_ans = inputdlg('Enter bin width for TPT histogram (us):', 'Bin Width', [1 50], {'3'});
        bin_width_us = str2double(bin_width_us_ans{1});
        
        color_HL = [0.2 0.55 0.8]; % Blue
        color_LH = [0.9 0.45 0.1]; % Orange

        if ~isempty(all_TPTs_sec_HL), histogram(ax, all_TPTs_sec_HL * 1e6, 'BinWidth', bin_width_us, 'Normalization', 'pdf', 'FaceColor', color_HL, 'EdgeColor', 'k', 'FaceAlpha', 0.7, 'DisplayName', 'High -> Low'); end
        if ~isempty(all_TPTs_sec_LH), histogram(ax, all_TPTs_sec_LH * 1e6, 'BinWidth', bin_width_us, 'Normalization', 'pdf', 'FaceColor', color_LH, 'EdgeColor', 'k', 'FaceAlpha', 0.7, 'DisplayName', 'Low -> High'); end
        
        % <<< START OF CORRECTION >>>
        mean_tpt_us_HL = mean(all_TPTs_sec_HL) * 1e6;
        mean_tpt_us_LH = mean(all_TPTs_sec_LH) * 1e6;
        num_tpts_HL = length(all_TPTs_sec_HL);
        num_tpts_LH = length(all_TPTs_sec_LH);

        if ~isnan(mean_tpt_us_HL), xline(ax, mean_tpt_us_HL, '--', 'Color', color_HL*0.9, 'LineWidth', 2, 'DisplayName', sprintf('Mean H->L = %.2f us', mean_tpt_us_HL)); end
        if ~isnan(mean_tpt_us_LH), xline(ax, mean_tpt_us_LH, '--', 'Color', color_LH*0.9, 'LineWidth', 2, 'DisplayName', sprintf('Mean L->H = %.2f us', mean_tpt_us_LH)); end
        
        hold(ax, 'off'); xlabel(ax, 'Transition Path Time (\mus)'); ylabel(ax, 'Probability Density');
        grid(ax, 'on'); box(ax, 'on'); legend(ax);
        
        title_str1 = sprintf('TPT Distribution from LLR (Model #%d, Lambda=%.2f)', chosen_model_idx, lambda_thresh);
        title_str2 = sprintf('H->L: N=%d, Mean=%.2f \\mus | L->H: N=%d, Mean=%.2f \\mus', num_tpts_HL, mean_tpt_us_HL, num_tpts_LH, mean_tpt_us_LH);
        title(ax, {title_str1, title_str2}, 'Interpreter', 'tex');
        % <<< END OF CORRECTION >>>
        
        saveas(fig_tpt, 'TPT_Distribution_LLR_Cluster.png');
        
        save_struct.all_TPTs_sec_HL = all_TPTs_sec_HL;
        save_struct.all_TPTs_sec_LH = all_TPTs_sec_LH;
        save_struct.lambda_thresh = lambda_thresh;
        save('TPT_Analysis_Results_Cluster.mat', '-struct', 'save_struct');
        fprintf('Results saved to TPT_Analysis_Results_Cluster.mat and .png\n');

        fprintf('\n--- TPT Calculation Finished ---\n');
        cd(original_pwd);

    catch ME_TPT
        fprintf(2,'\nERROR in calculateTPT_cluster: %s\n', ME_TPT.message);
        disp(ME_TPT.getReport());
        if exist('h_wait', 'var') && ishandle(h_wait), close(h_wait); end
        cd(original_pwd);
    end
end
% ... (LLR Sub-functions remain unchanged) ...
function idx_exit = findExitPoint(photon_colors, idx_change, pR_H, pG_H, pR_L, pG_L, lambda_thresh, transition_type), LLR=0; idx_exit=NaN; if strcmp(transition_type,'LH'),lr=log(pR_L/pR_H);lg=log(pG_L/pG_H);else,lr=log(pR_H/pR_L);lg=log(pG_H/pG_L);end;for i=(idx_change-1):-1:1,if photon_colors(i)==2,LLR=LLR+lr;else,LLR=LLR+lg;end;if LLR>lambda_thresh,idx_exit=i+1;return;end;end;end
function idx_entry = findEntryPoint(photon_colors, idx_change, pR_H, pG_H, pR_L, pG_L, lambda_thresh, transition_type), LLR=0;idx_entry=NaN;if strcmp(transition_type,'LH'),lr=log(pR_H/pR_L);lg=log(pG_H/pG_L);else,lr=log(pR_L/pR_H);lg=log(pG_L/pG_H);end;for i=idx_change:length(photon_colors),if photon_colors(i)==2,LLR=LLR+lr;else,LLR=LLR+lg;end;if LLR>lambda_thresh,idx_entry=i;return;end;end;end