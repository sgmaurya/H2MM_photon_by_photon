function runViterbi_local()
clearvars -except break_debug; close all; clc;
fprintf('--- Starting Local Viterbi Path Decoding ---\n');

      
%% 0. Setup Paths
project_base_path = fileparts(mfilename('fullpath'));
fprintf('INFO: Viterbi project base path: %s\n', project_base_path);

paths.hmm_core = fullfile(project_base_path, 'hmm_core_functions');
paths.results_base = fullfile(project_base_path, 'HMM_Local_Results');
addpath(project_base_path);
if exist(paths.hmm_core, 'dir'), addpath(genpath(paths.hmm_core)); fprintf('INFO: Added hmm_core_functions path: %s\n', paths.hmm_core);
else, warning('VITERBI_PATH_WARN: hmm_core_functions directory not found at %s', paths.hmm_core); end
original_pwd = pwd;

try
    %% 1. Select HMM Run Directory and Load Best Model & Data
    fprintf('\n--- 1. Selecting HMM Run Directory and Loading Model/Data ---\n');
    if ~exist(paths.results_base, 'dir'), error('VITERBI_SETUP_ERROR: Base results directory not found.'); end
    selected_run_path = uigetdir(paths.results_base, 'Select the HMM_Run_* directory');
    if isequal(selected_run_path, 0), disp('User cancelled. Exiting.'); return; end
    fprintf('Using results from: %s\n', selected_run_path);
    cd(selected_run_path);
    
    % Load model
    if exist('Best_HMM_Model_Local.mat', 'file')
        fprintf('Loading Best_HMM_Model_Local.mat...\n');
        load('Best_HMM_Model_Local.mat', 'best_model_overall', 'config_params_to_save', 'best_idx');
        model_params = best_model_overall; chosen_model_index = best_idx;
        fprintf('Using best model from guess #%d.\n', chosen_model_index);
    else
        fprintf('Best_HMM_Model_Local.mat not found. Loading All_HMM_Results_Collected_Local.mat...\n');
        if ~exist('All_HMM_Results_Collected_Local.mat', 'file'), error('VITERBI_LOAD_ERROR: All_HMM_Results_Collected_Local.mat not found.'); end
        load('All_HMM_Results_Collected_Local.mat', 'all_hmm_results_collected', 'all_final_log_likelihoods', 'config_params_to_save');
        valid_ll_indices = find(~isnan(all_final_log_likelihoods));
        if isempty(valid_ll_indices), error('VITERBI_LOAD_ERROR: No valid HMM models.'); end
        figure('Name', 'Select Model'); plot(1:length(all_final_log_likelihoods), all_final_log_likelihoods, 'o-'); xlabel('Model Index'); ylabel('LogL'); title('Select Model'); grid on;
        [max_ll_val, auto_idx] = max(all_final_log_likelihoods);
        ans_model = inputdlg(sprintf('Enter model index (e.g., %d for max LogL %.2f):', auto_idx, max_ll_val),'Choose Model', [1 50], {num2str(auto_idx)});
        if isempty(ans_model), error('VITERBI_USER_CANCEL: Model selection cancelled.'); end
        chosen_model_index = str2double(ans_model{1});
        if isnan(chosen_model_index) || chosen_model_index < 1 || chosen_model_index > length(all_hmm_results_collected) || isnan(all_final_log_likelihoods(chosen_model_index))
            error('VITERBI_INPUT_ERROR: Invalid model index.');
        end
        model_params = all_hmm_results_collected{chosen_model_index};
        fprintf('Using selected model from guess #%d.\n', chosen_model_index);
        if ishandle(gcf), close(gcf); end
    end

    % Load data
    fprintf('Loading data_re_for_hmm.mat...\n');
    if ~exist('data_re_for_hmm.mat', 'file'), error('VITERBI_LOAD_ERROR: data_re_for_hmm.mat not found.'); end
    load('data_re_for_hmm.mat', 'data_re_for_hmm'); % data_re_for_hmm{i}(:,1) = time in dt_analysis_sec units, {i}(:,2) = symbols
    
    % Load config parameters including dt_analysis_sec
    if isempty(config_params_to_save) || ~isfield(config_params_to_save, 'Nstates') || ~isfield(config_params_to_save, 'dt_analysis_sec')
        error('VITERBI_LOAD_ERROR: config_params_to_save incomplete.');
    end
    Nstates_from_config = config_params_to_save.Nstates;
    dt_model_sec = config_params_to_save.dt_analysis_sec; % This is the HMM model's dt (e.g., 1e-12 s)
    fprintf('INFO: HMM model dt (dt_model_sec) = %.2e seconds.\n', dt_model_sec);
    fprintf('INFO: Assuming photon arrival times in data_re_for_hmm{i}(:,1) are in units of this dt_model_sec.\n');

    if isempty(model_params) || ~isfield(model_params, 'transmat'), error('VITERBI_LOAD_ERROR: Loaded model_params invalid.'); end
    if isempty(data_re_for_hmm), error('VITERBI_LOAD_ERROR: data_re_for_hmm is empty.'); end
    save('LastViterbiModelUsed.mat','chosen_model_index', 'model_params', 'config_params_to_save');

    %% 2. Extract Model Parameters and Prepare AA matrix for Visualization
    fprintf('\n--- 2. Extracting Model Parameters & Preparing AA for Visualization ---\n');
    model_prior = model_params.prior(:);
    if size(model_prior,2) > 1, model_prior = model_prior'; end
    model_obsmat = model_params.obsmat;
    model_transmat = model_params.transmat;

    AA_to_save = zeros(2, Nstates_from_config);
    AA_to_save(1,:) = 1:Nstates_from_config;
    if Nstates_from_config == size(model_obsmat,1) && size(model_obsmat,2) >= 2
        plot_values_for_states = model_obsmat(:,2)';
        fprintf('INFO: Using P(AcceptorSymbol|State) from model_obsmat for Viterbi path plotting levels.\n');
        if Nstates_from_config > 0, fprintf('       HMM State 1 will be plotted at %.4f\n', plot_values_for_states(1)); end
        if Nstates_from_config >= 2, fprintf('       HMM State 2 will be plotted at %.4f\n', plot_values_for_states(2)); end
    else, error('AA_CREATION_ERROR: model_obsmat dimensions mismatch.'); end
    AA_to_save(2,:) = plot_values_for_states;
    AA = AA_to_save;
    save('AA.mat', 'AA');
    fprintf('AA.mat saved (AA(1,:)=HMM states; AA(2,:)=P(A|State) for plotting).\n');

    %% 3. Execute Viterbi on each trajectory
    fprintf('\n--- 3. Executing Viterbi Algorithm ---\n');
    if ~exist('viterbi_path_PhotnoByPhoton2_local.m', 'file'), error('VITERBI_FUNC_MISSING: viterbi_path_PhotnoByPhoton2_local.m missing.'); end
    if ~exist('multinomial_prob.m', 'file'), error('VITERBI_FUNC_MISSING: multinomial_prob.m missing.'); end

    num_trajectories = length(data_re_for_hmm);
    Q0 = cell(1, num_trajectories);
    h_wait = waitbar(0, 'Running Viterbi...');
    cleanupObj_waitbar = onCleanup(@() closeViterbiWaitbar(h_wait));

    for i = 1:num_trajectories
        if mod(i, round(num_trajectories/100)+1) == 0 || i == num_trajectories || i == 1
             waitbar(i/num_trajectories, h_wait, sprintf('Viterbi: Traj %d/%d', i, num_trajectories));
        end

        current_trajectory_data_units = data_re_for_hmm{i};
        if isempty(current_trajectory_data_units) || size(current_trajectory_data_units,1) < 1
            Q0{i} = []; continue;
        end

        times_in_dt_units = current_trajectory_data_units(:,1); % These are counts of dt_model_sec
        observations = current_trajectory_data_units(:,2);
        obslik = multinomial_prob(observations, model_obsmat);

        actual_inter_photon_durations_sec = [];
        if size(times_in_dt_units,1) > 1
            inter_photon_durations_in_dt_units = diff(times_in_dt_units);
            % Convert durations from dt_model_sec units to actual seconds for Viterbi function
            actual_inter_photon_durations_sec = inter_photon_durations_in_dt_units * dt_model_sec;
            % Ensure minimum duration is at least dt_model_sec if very small (e.g., if diff was 0 or 1 dt unit)
            actual_inter_photon_durations_sec(actual_inter_photon_durations_sec < (dt_model_sec/2)) = dt_model_sec;
        end
        Q0{i} = viterbi_path_PhotnoByPhoton2_local(model_prior, model_transmat, obslik, actual_inter_photon_durations_sec, dt_model_sec);
    end
    save('Viterbi_decoded_states_Q0.mat','Q0');
    fprintf('Viterbi decoded states saved to Viterbi_decoded_states_Q0.mat\n');

    %% 4. Finding number of transitions
    fprintf('\n--- 4. Analyzing Transitions ---\n');
    ind_trans = cell(1, num_trajectories); trans_num = zeros(1, num_trajectories);
    for i = 1:num_trajectories
        q_temp = Q0{i}; if isempty(q_temp), continue; end; T_q = length(q_temp); ind_temp_transitions = [];
        if T_q > 1, for t_q = 1:(T_q-1), if q_temp(t_q) ~= q_temp(t_q+1), ind_temp_transitions(end+1) = t_q; end; end; end %#ok<AGROW>
        ind_trans{i} = ind_temp_transitions; trans_num(i) = length(ind_temp_transitions);
    end
    save('Viterbi_transition_analysis.mat','ind_trans','trans_num');
    fprintf('Transition analysis saved.\n');

    %% 5. Call Trajectory Visualization
    fprintf('\n--- 5. Visualizing Trajectories (Optional) ---\n');
    
    % Prepare data for Traj_Vet_local: times must be in SECONDS
    data_for_traj_vet_plotting = cell(size(data_re_for_hmm));
    for i_traj = 1:length(data_re_for_hmm)
        if ~isempty(data_re_for_hmm{i_traj})
            temp_data = data_re_for_hmm{i_traj};
            temp_data(:,1) = temp_data(:,1) * dt_model_sec; % Convert time column from dt_units to seconds
            data_for_traj_vet_plotting{i_traj} = temp_data;
        else
            data_for_traj_vet_plotting{i_traj} = [];
        end
    end
    
    ans_trajvet_type = questdlg('Display Viterbi paths with Traj_Vet_local?', 'Trajectory Visualization', 'Yes', 'No', 'Yes');
    if strcmp(ans_trajvet_type, 'Yes')
        if exist('Traj_Vet_local.m', 'file')
            fprintf('Running Traj_Vet_local...\n');
            % Pass dt_model_sec as the 'dt_tv' argument to Traj_Vet_local
            Traj_Vet_local(data_for_traj_vet_plotting, Q0, ind_trans, trans_num, AA, dt_model_sec, Nstates_from_config);
        else, warning('TRAJVET_NOT_FOUND: Traj_Vet_local.m not found.'); end
    else, fprintf('Skipping trajectory visualization.\n'); end

    fprintf('\n--- Local Viterbi Path Decoding Finished ---\n');
    cd(original_pwd); fprintf('Returned to original directory: %s\n', original_pwd);
catch ME_viterbi_local
    fprintf(2, 'ERROR in runViterbi_local: %s\n', ME_viterbi_local.message); disp(ME_viterbi_local.getReport());
    if exist('h_wait', 'var') && ishandle(h_wait), close(h_wait); end
    if exist('original_pwd', 'var') && isfolder(original_pwd) && ~strcmp(pwd, original_pwd), cd(original_pwd); fprintf('Returned to original directory after error.\n'); end
end
function closeViterbiWaitbar(h)
    if exist('h', 'var') && ishandle(h), close(h); end
end

end