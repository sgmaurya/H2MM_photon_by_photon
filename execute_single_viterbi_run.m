% execute_single_viterbi_run.m
% A fully automated, non-interactive worker function to run the Viterbi
% algorithm on the results of a single HMM run. The visualization step
% is controlled by an input flag.

function execute_single_viterbi_run(results_path, run_visualization_flag)
    
    % --- DIARY SETUP ---
    log_base_dir = fullfile(results_path, 'Logs');
    if ~exist(log_base_dir, 'dir'), mkdir(log_base_dir); end
    timestamp_for_diary = datestr(now, 'yyyy-mm-dd_HHMMSS');
    diary_filename = fullfile(log_base_dir, sprintf('viterbi_run_log_%s.txt', timestamp_for_diary));
    if exist(diary_filename, 'file'), delete(diary_filename); end
    diary(diary_filename);
    cleanupObj_diary = onCleanup(@() diary('off'));
    fprintf('--- Starting Automated Single Viterbi Path Decoding ---\n');
    fprintf('Processing results in: %s\n', results_path);
    
    original_pwd = pwd;
    cleanupObj_cd = onCleanup(@() cd(original_pwd)); % Ensure we always return
    
    try
        %% 1. Load Best Model & Data
        fprintf('\n--- 1. Loading Model and Data ---\n');
        cd(results_path); % Change to the results directory
        
        if exist('Best_HMM_Model_Local.mat', 'file')
            load('Best_HMM_Model_Local.mat', 'best_model_overall', 'config_params_to_save');
            model_params = best_model_overall;
        else
            error('VITERBI_LOAD_ERROR: Best_HMM_Model_Local.mat not found.');
        end

        load('data_re_for_hmm.mat', 'data_re_for_hmm');
        
        Nstates = config_params_to_save.Nstates;
        dt_model_sec = config_params_to_save.dt_analysis_sec;
        
        %% 2. Extract Model Parameters
        fprintf('\n--- 2. Extracting Model Parameters ---\n');
        model_prior = model_params.prior(:)'; % Ensure row vector
        model_obsmat = model_params.obsmat;
        model_transmat = model_params.transmat;
        
        AA(1,:) = 1:Nstates;
        AA(2,:) = model_obsmat(:,2)'; % P(Acceptor|State)
        save('AA.mat', 'AA');
        
        %% 3. Execute Viterbi Algorithm
        fprintf('\n--- 3. Executing Viterbi Algorithm ---\n');
        
        num_trajectories = length(data_re_for_hmm);
        Q0 = cell(1, num_trajectories);
        
        parfor i = 1:num_trajectories
            current_trajectory_data = data_re_for_hmm{i};
            if isempty(current_trajectory_data), continue; end
            
            times_in_dt_units = current_trajectory_data(:,1);
            observations = current_trajectory_data(:,2);
            obslik = multinomial_prob(observations, model_obsmat);

            actual_inter_photon_durations_sec = [];
            if size(times_in_dt_units,1) > 1
                inter_photon_durations_in_dt_units = diff(times_in_dt_units);
                actual_inter_photon_durations_sec = inter_photon_durations_in_dt_units * dt_model_sec;
                actual_inter_photon_durations_sec(actual_inter_photon_durations_sec < (dt_model_sec/2)) = dt_model_sec;
            end
            Q0{i} = viterbi_path_PhotnoByPhoton2_local(model_prior, model_transmat, obslik, actual_inter_photon_durations_sec, dt_model_sec);
        end
        save('Viterbi_decoded_states_Q0.mat','Q0');
        fprintf('Viterbi decoded states saved to Viterbi_decoded_states_Q0.mat\n');

        %% 4. Analyze Transitions
        fprintf('\n--- 4. Analyzing Transitions ---\n');
        ind_trans = cell(1, num_trajectories);
        trans_num = zeros(1, num_trajectories);
        for i = 1:num_trajectories
            q_temp = Q0{i};
            if ~isempty(q_temp) && length(q_temp) > 1
                ind_trans{i} = find(diff(q_temp) ~= 0);
                trans_num(i) = length(ind_trans{i});
            end
        end
        save('Viterbi_transition_analysis.mat','ind_trans','trans_num');
        fprintf('Transition analysis saved.\n');

        %% 5. Call Trajectory Visualization (if flagged)
        if run_visualization_flag
            fprintf('\n--- 5. Visualizing Trajectories ---\n');
            data_for_plotting = cell(size(data_re_for_hmm));
            for i = 1:length(data_re_for_hmm)
                if ~isempty(data_re_for_hmm{i})
                    temp_data = data_re_for_hmm{i};
                    temp_data(:,1) = temp_data(:,1) * dt_model_sec; % Convert time to seconds
                    data_for_plotting{i} = temp_data;
                end
            end
            
            if exist('Traj_Vet_local.m', 'file')
                Traj_Vet_local(data_for_plotting, Q0, ind_trans, trans_num, AA, dt_model_sec, Nstates);
            else
                warning('TRAJVET_NOT_FOUND: Traj_Vet_local.m not found. Skipping visualization.');
            end
        else
            fprintf('\n--- 5. Skipping Trajectory Visualization as requested ---\n');
        end
        
        fprintf('\n--- Viterbi Decoding Finished Successfully ---\n');
        
    catch ME_viterbi
        fprintf(2, '\n--- ERROR in Automated Viterbi Run ---\n');
        fprintf(2, 'Error occurred in folder: %s\n', results_path);
        rethrow(ME_viterbi);
    end
end