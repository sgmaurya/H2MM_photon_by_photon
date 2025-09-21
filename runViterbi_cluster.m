% runViterbi_cluster.m (v4 - Fixed AA Matrix for Visualization)
%
% This is the first step in the CLUSTER analysis pipeline.
%
% v4 Corrections:
%   - Prepares a correctly formatted AA matrix for the Traj_Vet_local viewer
%     to prevent the "AA_tv(1,:) is not [1:Nstates_tv]" warning and ensure
%     Viterbi levels are plotted with the correct absolute FRET values.
%

function runViterbi_cluster()
    clearvars -except break_debug; close all; clc;
    fprintf('--- Starting Cluster Viterbi Path Decoding ---\n');

    %% 0. Setup Paths
    project_base_path = fileparts(mfilename('fullpath'));
    paths.hmm_core = fullfile(project_base_path, 'hmm_core_functions');
    if exist(paths.hmm_core, 'dir')
        addpath(genpath(paths.hmm_core));
        fprintf('INFO: Added hmm_core_functions path: %s\n', paths.hmm_core);
    else
        warning('VITERBI_CLUSTER_WARN: hmm_core_functions directory not found at %s. Dependent functions may not work.', paths.hmm_core);
    end
    original_pwd = pwd;

    try
        %% 1. Select Processed Cluster Run Directory and Load Data
        fprintf('\n--- 1. Selecting Processed Cluster HMM Run Directory ---\n');
        selected_run_path = uigetdir(pwd, 'Select the Cluster HMM Run directory (containing HMM_parm.mat, data.mat, etc.)');
        if isequal(selected_run_path, 0), disp('User cancelled. Exiting.'); return; end
        fprintf('Processing results from: %s\n', selected_run_path);
        cd(selected_run_path);

        if ~exist('HMM_parm.mat', 'file'), error('VITERBI_CLUSTER_ERROR: HMM_parm.mat not found.'); end
        load('HMM_parm.mat', 'Trans', 'K'); 
        Nstates = size(Trans, 1);
        
        [~, max_idx] = max(abs(K(:))); [row, col] = ind2sub(size(K), max_idx);
        if K(row, col) ~= 0 && Trans(row, col) ~=0
             dt_model_sec = Trans(row,col) / K(row,col);
             fprintf('INFO: Inferred HMM dt = %.2e s from HMM_parm.mat.\n', dt_model_sec);
        else, error('Cannot infer dt from HMM_parm.mat.'); end

        if ~exist('Est.mat', 'file'), error('VITERBI_CLUSTER_ERROR: Est.mat not found.'); end
        loaded_est = load('Est.mat', 'ass'); 
        ass = loaded_est.ass; 
        
        load('HMM_output.mat', 'obsmat', 'transmat', 'prior');
        model_params.obsmat = obsmat{1,Nstates}{ass};
        model_params.transmat = transmat{1,Nstates}{ass};
        model_params.prior = prior{1,Nstates}{ass};
        if size(model_params.prior, 2) > 1, model_params.prior = model_params.prior(:); end
        
        fprintf('Successfully reconstructed original HMM model from guess #%d.\n', ass);

        if ~exist('data.mat', 'file'), error('VITERBI_CLUSTER_ERROR: data.mat not found.'); end
        load('data.mat', 'data');

        %% 2. Execute Viterbi Algorithm on All Trajectories
        fprintf('\n--- 2. Executing Viterbi Algorithm ---\n');
        
        viterbi_function_name = 'viterbi_path_PhotnoByPhoton2_local';
        if ~exist([viterbi_function_name, '.m'], 'file'), error('Viterbi function ''%s.m'' not found.', viterbi_function_name); end
        if ~exist('multinomial_prob.m', 'file'), error('Dependency ''multinomial_prob.m'' not found.'); end

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
            
            Q0{i} = feval(viterbi_function_name, model_params.prior, model_params.transmat, obslik, inter_photon_durations_sec, dt_model_sec);
        end
        if ishandle(h_wait_viterbi), close(h_wait_viterbi); end
        
        save('Viterbi_decoded_states_Q0_cluster.mat','Q0');
        fprintf('Viterbi decoded states saved to Viterbi_decoded_states_Q0_cluster.mat\n');

        %% 3. Finding Number of Transitions
        fprintf('\n--- 3. Analyzing Transitions ---\n');
        ind_trans = cell(1, num_trajectories); trans_num = zeros(1, num_trajectories);
        for i = 1:num_trajectories
            q_temp = Q0{i}; if isempty(q_temp), continue; end; T_q = length(q_temp); ind_temp_transitions = [];
            if T_q > 1, for t_q = 1:(T_q-1), if q_temp(t_q) ~= q_temp(t_q+1), ind_temp_transitions(end+1) = t_q; end; end; end
            ind_trans{i} = ind_temp_transitions; trans_num(i) = length(ind_temp_transitions);
        end
        save('Viterbi_transition_analysis_cluster.mat','ind_trans','trans_num');
        fprintf('Transition analysis saved to Viterbi_transition_analysis_cluster.mat\n');

        %% 4. Call Trajectory Visualization
        fprintf('\n--- 4. Visualizing Trajectories (Optional) ---\n');
        
        ans_trajvet = questdlg('Display Viterbi paths with Traj_Vet_local?', 'Trajectory Visualization', 'Yes', 'No', 'Yes');
        if strcmp(ans_trajvet, 'Yes')
            if exist('Traj_Vet_local.m', 'file')
                fprintf('Preparing data and running Traj_Vet_local...\n');
                
                % --- START OF CORRECTION ---
                
                % Load the AA matrix saved by Recoloring_Param.m
                load('AA.mat', 'AA');
                
                % Create a new AA matrix specifically for the viewer
                % that conforms to the viewer's expected format.
                AA_for_plotting = zeros(2, Nstates);
                AA_for_plotting(1,:) = 1:Nstates; % First row is [1 2] to satisfy the warning
                
                % "Unsort" the FRET values to match the original state indices
                unsorted_fret_values(AA(1,:)) = AA(2,:);
                AA_for_plotting(2,:) = unsorted_fret_values;
                
                % Prepare time data for plotting (normalized to start at t=0)
                data_for_plotting = data;
                for i = 1:length(data)
                    if ~isempty(data{i}) && size(data{i}, 1) > 0
                        time_in_sec = data{i}(:,1) * dt_model_sec;
                        normalized_time_sec = time_in_sec - time_in_sec(1);
                        data_for_plotting{i}(:,1) = normalized_time_sec;
                    end
                end
                
                % Call the viewer with the correctly formatted data
                Traj_Vet_local(data_for_plotting, Q0, ind_trans, trans_num, AA_for_plotting, dt_model_sec, Nstates);
                
                % --- END OF CORRECTION ---
                
            else
                warning('Traj_Vet_local.m not found on path. Skipping visualization.');
            end
        else, fprintf('Skipping trajectory visualization.\n'); end

        fprintf('\n--- Cluster Viterbi Path Decoding Finished ---\n');
        cd(original_pwd);

    catch ME_Viterbi
        fprintf(2,'\nERROR in runViterbi_cluster: %s\n', ME_Viterbi.message);
        disp(ME_Viterbi.getReport());
        if exist('h_wait_viterbi', 'var') && ishandle(h_wait_viterbi), close(h_wait_viterbi); end
        cd(original_pwd);
    end
end