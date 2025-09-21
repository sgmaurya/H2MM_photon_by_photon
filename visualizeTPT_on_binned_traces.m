% visualizeTPT_on_binned_traces.m
%
% LAUNCHER SCRIPT
% This script loads all necessary data from an HMM run folder, calculates
% the TPT boundaries for every trajectory, and then calls the modified
% viewer 'Traj_Vet_local_with_TPT.m' to display the results.
%

function visualizeTPT_on_binned_traces()
    clearvars -except break_debug; close all; clc;
    fprintf('--- Launcher for TPT Visualization on Binned Traces ---\n');

    original_pwd = pwd;
    try
        %% 1. Select HMM Run Directory and Load All Data
        fprintf('\n--- 1. Selecting HMM Run Directory ---\n');
        
        project_base_path = fileparts(mfilename('fullpath'));
        results_base_path = fullfile(project_base_path, 'HMM_Local_Results');
        if ~exist(results_base_path, 'dir'), results_base_path = pwd; end
        
        selected_run_path = uigetdir(results_base_path, 'Select the HMM_Run_* directory with TPT results');
        if isequal(selected_run_path, 0), disp('User cancelled. Exiting.'); return; end
        fprintf('Processing results from: %s\n', selected_run_path);
        cd(selected_run_path);

        % --- Load ALL necessary files ---
        fprintf('Loading results...\n');
        tpt_results = load('TPT_Analysis_Results_Local.mat');
        load('Viterbi_decoded_states_Q0.mat', 'Q0');
        load('data_re_for_hmm.mat', 'data_re_for_hmm');
        load('AA.mat', 'AA');
        load('Viterbi_transition_analysis.mat', 'ind_trans', 'trans_num');

        data = data_re_for_hmm;
        dt_model_sec = tpt_results.dt_model_sec;
        Nstates = size(tpt_results.model_params_used.obsmat_orig, 1);
        
        %% 2. Re-calculate TPT Boundaries (Photon Indices)
        fprintf('\n--- 2. Locating TPT Boundaries for Visualization ---\n');
        
        P_Red_H = tpt_results.model_params_used.P_Red_H;
        P_Green_H = 1 - P_Red_H;
        P_Red_L = tpt_results.model_params_used.P_Red_L;
        P_Green_L = 1 - P_Red_L;
        lambda_thresh = tpt_results.lambda_thresh;
        obsmat_orig = tpt_results.model_params_used.obsmat_orig;
        
        [~, sorted_idx] = sort(obsmat_orig(:,2));
        state_map(sorted_idx(1)) = 1; state_map(sorted_idx(2)) = 2;

        num_trajectories = length(data);
        tpt_boundaries_per_traj = cell(1, num_trajectories);
        
        h_wait = waitbar(0, 'Re-calculating TPT boundaries for plotting...');
        for i_traj = 1:num_trajectories
            waitbar(i_traj / num_trajectories, h_wait);
            
            viterbi_path_orig = Q0{i_traj};
            photon_stream = data{i_traj};
            if isempty(viterbi_path_orig) || length(viterbi_path_orig) < 2, continue; end
            
            viterbi_path_mapped = state_map(viterbi_path_orig);
            transition_indices = find(diff(viterbi_path_mapped) ~= 0) + 1;
            
            boundaries_this_traj = [];
            for k_trans = 1:length(transition_indices)
                idx_change = transition_indices(k_trans);
                state_before=viterbi_path_mapped(idx_change-1); state_after=viterbi_path_mapped(idx_change);
                
                transition_type = '';
                if state_before == 1 && state_after == 2, transition_type = 'LH';
                elseif state_before == 2 && state_after == 1, transition_type = 'HL';
                else, continue; end
                
                idx_exit = findExitPoint(photon_stream(:,2), idx_change, P_Red_H, P_Green_H, P_Red_L, P_Green_L, lambda_thresh, transition_type);
                idx_entry = findEntryPoint(photon_stream(:,2), idx_change, P_Red_H, P_Green_H, P_Red_L, P_Green_L, lambda_thresh, transition_type);

                if ~isnan(idx_exit) && ~isnan(idx_entry) && idx_entry > idx_exit
                    boundaries_this_traj(end+1, :) = [idx_exit, idx_entry];
                end
            end
            tpt_boundaries_per_traj{i_traj} = boundaries_this_traj;
        end
        if ishandle(h_wait), close(h_wait); end

        %% 3. Call the Visualization Function
        fprintf('\n--- 3. Launching Interactive Trajectory Viewer ---\n');
        
        if ~exist('Traj_Vet_local_with_TPT.m', 'file')
            error('The viewer script "Traj_Vet_local_with_TPT.m" was not found in the path.');
        end
        
        % Call the modified viewer and pass the TPT boundaries to it
        Traj_Vet_local_with_TPT(data, Q0, ind_trans, trans_num, AA, dt_model_sec, Nstates, tpt_boundaries_per_traj);

        fprintf('\n--- Visualization Finished ---\n');
        cd(original_pwd);

    catch ME_TPT_Vis
        fprintf(2,'\nERROR in visualizeTPT_on_binned_traces: %s\n', ME_TPT_Vis.message);
        disp(ME_TPT_Vis.getReport());
        if exist('h_wait', 'var') && ishandle(h_wait), close(h_wait); end
        cd(original_pwd);
    end
end


%% LLR Sub-functions (needed for boundary calculation)
function idx_exit = findExitPoint(photon_colors, idx_change, pR_H, pG_H, pR_L, pG_L, lambda_thresh, transition_type)
    LLR = 0;
    idx_exit = NaN;
    if strcmp(transition_type, 'LH'), log_update_red=log(pR_L/pR_H); log_update_green=log(pG_L/pG_H);
    else, log_update_red=log(pR_H/pR_L); log_update_green=log(pG_H/pG_L); end
    for i=(idx_change-1):-1:1
        if photon_colors(i)==2, LLR=LLR+log_update_red; else, LLR=LLR+log_update_green; end
        if LLR > lambda_thresh, idx_exit=i+1; return; end
    end
end

function idx_entry = findEntryPoint(photon_colors, idx_change, pR_H, pG_H, pR_L, pG_L, lambda_thresh, transition_type)
    LLR = 0;
    idx_entry = NaN;
    if strcmp(transition_type, 'LH'), log_update_red=log(pR_H/pR_L); log_update_green=log(pG_H/pG_L);
    else, log_update_red=log(pR_L/pR_H); log_update_green=log(pG_L/pG_H); end
    for i=idx_change:length(photon_colors)
        if photon_colors(i)==2, LLR=LLR+log_update_red; else, LLR=LLR+log_update_green; end
        if LLR > lambda_thresh, idx_entry=i; return; end
    end
end