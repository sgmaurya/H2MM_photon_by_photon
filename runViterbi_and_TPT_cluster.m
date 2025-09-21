% runViterbi_and_TPT_cluster.m
%
% This is the main analysis script for cluster data. It should be run AFTER
% the HMM results have been processed (e.g., using your Recoloring_Param script).
%
% This script performs two main tasks:
%   1. Runs the Viterbi algorithm on all trajectories using the processed model.
%   2. Calculates the Transition Path Times (TPTs) for all transitions using
%      the LLR method.
%

function runViterbi_and_TPT_cluster()
    clearvars -except break_debug; close all; clc;
    fprintf('--- Starting Cluster Viterbi Decoding & TPT Calculation ---\n');

    original_pwd = pwd;
    try
        %% 1. Select Processed Cluster Run Directory and Load Data
        fprintf('\n--- 1. Selecting Processed Cluster HMM Run Directory ---\n');
        selected_run_path = uigetdir(pwd, 'Select the Cluster HMM Run directory (containing HMM_parm.mat, data.mat, etc.)');
        if isequal(selected_run_path, 0), disp('User cancelled. Exiting.'); return; end
        fprintf('Processing results from: %s\n', selected_run_path);
        cd(selected_run_path);

        % --- Load Processed HMM Parameters ---
        fprintf('Loading HMM_parm.mat...\n');
        if ~exist('HMM_parm.mat', 'file'), error('VITERBI_CLUSTER_ERROR: HMM_parm.mat not found. Please process the HMM run first.'); end
        load('HMM_parm.mat', 'Obs', 'Trans', 'Prior', 'K'); % Obs=Sorted FRET, Trans=Sorted transmat
        
        % --- Infer dt from the saved K matrix (most robust method) ---
        % K = (Trans - eye(N)) / dt  =>  dt = (Trans(i,j)) / K(i,j) for off-diagonal
        Nstates = size(Trans, 1);
        [max_val, max_idx] = max(abs(K(:))); % Find largest off-diagonal rate
        [row, col] = ind2sub(size(K), max_idx);
        if K(row, col) ~= 0 && Trans(row, col) ~=0
             dt_model_sec = Trans(row,col) / K(row,col);
             fprintf('INFO: Inferred HMM dt = %.2e s from HMM_parm.mat.\n', dt_model_sec);
        else
            error('Cannot infer dt from HMM_parm.mat as K matrix is zero. Please process HMM results again.');
        end

        % --- Construct the original, unsorted model for Viterbi ---
        fprintf('Loading Est.mat to reconstruct original model for Viterbi...\n');
        if ~exist('Est.mat', 'file'), error('VITERBI_CLUSTER_ERROR: Est.mat not found.'); end
        load('Est.mat', 'Est_ass', 'Est_B'); % Est_B is the sorting index
        
        % The 'ass' index refers to the model choice from the original HMM_output.mat
        % Load the original HMM_output to get the unsorted model
        load('HMM_output.mat', 'obsmat', 'transmat', 'prior');
        model_params.obsmat = obsmat{1,Nstates}{Est_ass};
        model_params.transmat = transmat{1,Nstates}{Est_ass};
        model_params.prior = prior{1,Nstates}{Est_ass};
        if size(model_params.prior, 2) > 1, model_params.prior = model_params.prior'; end

        % --- Load Raw Photon Data ---
        fprintf('Loading data.mat...\n');
        if ~exist('data.mat', 'file'), error('VITERBI_CLUSTER_ERROR: data.mat not found.'); end
        load('data.mat', 'data');
        if ~iscell(data) || isempty(data), error('VITERBI_CLUSTER_ERROR: Trajectory data is empty or invalid.'); end

        %% 2. Execute Viterbi Algorithm on All Trajectories
        fprintf('\n--- 2. Executing Viterbi Algorithm ---\n');
        
        num_trajectories = length(data);
        Q0 = cell(1, num_trajectories);
        h_wait_viterbi = waitbar(0, 'Running Viterbi...');
        
        for i = 1:num_trajectories
            waitbar(i/num_trajectories, h_wait_viterbi, sprintf('Viterbi: Traj %d/%d', i, num_trajectories));
            current_traj = data{i};
            if isempty(current_traj) || size(current_traj,1) < 1, Q0{i} = []; continue; end
            
            times_in_dt_units = current_traj(:,1);
            observations = current_traj(:,2);
            obslik = multinomial_prob(observations, model_params.obsmat);

            inter_photon_durations_sec = [];
            if size(times_in_dt_units,1) > 1
                inter_photon_durations_dt_units = diff(times_in_dt_units);
                inter_photon_durations_sec = inter_photon_durations_dt_units * dt_model_sec;
                inter_photon_durations_sec(inter_photon_durations_sec < (dt_model_sec/2)) = dt_model_sec;
            end
            
            % Using a standard Viterbi function, assuming it's in the path
            Q0{i} = viterbi_path(model_params.prior, model_params.transmat, obslik, inter_photon_durations_sec, dt_model_sec);
        end
        if ishandle(h_wait_viterbi), close(h_wait_viterbi); end
        
        save('Viterbi_decoded_states_Q0_cluster.mat','Q0');
        fprintf('Viterbi decoded states saved to Viterbi_decoded_states_Q0_cluster.mat\n');

        %% 3. Calculate Transition Path Times using LLR
        fprintf('\n--- 3. Calculating TPTs for all Transitions ---\n');
        if Nstates ~= 2, error('TPT calculation is currently configured only for 2-state models.'); end

        % The HMM_parm.mat 'Obs' variable is already sorted by FRET
        P_Red_L = Obs(1); P_Green_L = 1 - P_Red_L;
        P_Red_H = Obs(2); P_Green_H = 1 - P_Red_H;
        
        ans_lambda = inputdlg({'Enter LLR Threshold (Lambda_thresh, e.g., log(100)):'}, 'LLR Parameter', [1 60], {'4.6'});
        if isempty(ans_lambda), error('User cancelled Lambda input.'); end
        lambda_thresh = str2double(ans_lambda{1});
        fprintf('INFO: LLR evidence threshold set to %.2f\n', lambda_thresh);

        all_TPTs_sec_HL = []; all_TPTs_sec_LH = [];
        h_wait_tpt = waitbar(0, 'Calculating TPTs...');
        
        % The sorting index Est_B maps: Est_B(1)=original_idx_of_LowFRET, Est_B(2)=original_idx_of_HighFRET
        state_map = zeros(1, Nstates);
        state_map(Est_B(1)) = 1; % Low FRET
        state_map(Est_B(2)) = 2; % High FRET

        for i_traj = 1:num_trajectories
            waitbar(i_traj/num_trajectories, h_wait_tpt, sprintf('TPT: Traj %d/%d', i_traj, num_trajectories));
            
            viterbi_path_orig = Q0{i_traj};
            photon_stream = data{i_traj};
            if isempty(viterbi_path_orig) || length(viterbi_path_orig) < 2, continue; end
            
            viterbi_path_mapped = state_map(viterbi_path_orig);
            transition_indices = find(diff(viterbi_path_mapped) ~= 0) + 1;
            
            for k_trans = 1:length(transition_indices)
                idx_change = transition_indices(k_trans);
                if viterbi_path_mapped(idx_change-1)==1 && viterbi_path_mapped(idx_change)==2, transition_type='LH';
                elseif viterbi_path_mapped(idx_change-1)==2 && viterbi_path_mapped(idx_change)==1, transition_type='HL';
                else, continue; end
                
                idx_exit = findExitPoint(photon_stream(:,2), idx_change, P_Red_H, P_Green_H, P_Red_L, P_Green_L, lambda_thresh, transition_type);
                idx_entry = findEntryPoint(photon_stream(:,2), idx_change, P_Red_H, P_Green_H, P_Red_L, P_Green_L, lambda_thresh, transition_type);

                if ~isnan(idx_exit) && ~isnan(idx_entry) && idx_entry > idx_exit
                    tpt_sec = (photon_stream(idx_entry, 1) - photon_stream(idx_exit, 1)) * dt_model_sec;
                    if tpt_sec >= 0
                        if strcmp(transition_type, 'HL'), all_TPTs_sec_HL(end+1) = tpt_sec;
                        else, all_TPTs_sec_LH(end+1) = tpt_sec; end
                    end
                end
            end
        end
        if ishandle(h_wait_tpt), close(h_wait_tpt); end

        save_struct.all_TPTs_sec_HL = all_TPTs_sec_HL;
        save_struct.all_TPTs_sec_LH = all_TPTs_sec_LH;
        save('TPT_Analysis_Results_Cluster.mat', '-struct', 'save_struct');
        fprintf('TPT analysis results saved to TPT_Analysis_Results_Cluster.mat\n');

        %% 4. Generate TPT Histogram Plot
        fprintf('\n--- 4. Generating Final TPT Histogram Plot ---\n');
        figure('Name', 'TPT Distribution (Cluster)', 'Color', 'w', 'Position', [300, 300, 900, 600]);
        hold on;
        if ~isempty(all_TPTs_sec_HL), histogram(all_TPTs_sec_HL * 1e6, 'Normalization', 'pdf', 'DisplayName', 'High -> Low', 'FaceAlpha', 0.7); end
        if ~isempty(all_TPTs_sec_LH), histogram(all_TPTs_sec_LH * 1e6, 'Normalization', 'pdf', 'DisplayName', 'Low -> High', 'FaceAlpha', 0.7); end
        hold off; xlabel('Transition Path Time (\mus)'); ylabel('Probability Density');
        grid on; box on; legend;
        title_str = sprintf('TPT Distribution from LLR (Model #%d, Lambda=%.2f, N_{Traj}=%d)', ...
                            Est_ass, lambda_thresh, num_trajectories);
        title(title_str, 'Interpreter', 'none');
        saveas(gcf, 'TPT_Distribution_LLR_Cluster.png');

        fprintf('\n--- Cluster Viterbi & TPT Analysis Finished ---\n');
        cd(original_pwd);

    catch ME_Viterbi_TPT
        fprintf(2,'\nERROR in runViterbi_and_TPT_cluster: %s\n', ME_Viterbi_TPT.message);
        disp(ME_Viterbi_TPT.getReport());
        if exist('h_wait_viterbi', 'var') && ishandle(h_wait_viterbi), close(h_wait_viterbi); end
        if exist('h_wait_tpt', 'var') && ishandle(h_wait_tpt), close(h_wait_tpt); end
        cd(original_pwd);
    end
end

%% LLR Sub-functions
function idx_exit = findExitPoint(photon_colors, idx_change, pR_H, pG_H, pR_L, pG_L, lambda_thresh, transition_type)
    LLR=0; idx_exit=NaN;
    if strcmp(transition_type, 'LH'), log_update_red=log(pR_L/pR_H); log_update_green=log(pG_L/pG_H);
    else, log_update_red=log(pR_H/pR_L); log_update_green=log(pG_H/pG_L); end
    for i=(idx_change-1):-1:1
        if photon_colors(i)==2, LLR=LLR+log_update_red; else, LLR=LLR+log_update_green; end
        if LLR > lambda_thresh, idx_exit=i+1; return; end
    end
end

function idx_entry = findEntryPoint(photon_colors, idx_change, pR_H, pG_H, pR_L, pG_L, lambda_thresh, transition_type)
    LLR=0; idx_entry=NaN;
    if strcmp(transition_type, 'LH'), log_update_red=log(pR_H/pR_L); log_update_green=log(pG_H/pG_L);
    else, log_update_red=log(pR_L/pR_H); log_update_green=log(pG_L/pG_H); end
    for i=idx_change:length(photon_colors)
        if photon_colors(i)==2, LLR=LLR+log_update_red; else, LLR=LLR+log_update_green; end
        if LLR > lambda_thresh, idx_entry=i; return; end
    end
end