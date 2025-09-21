% execute_single_hmm_run.m
% Worker function for running a single H2MM analysis. NON-INTERACTIVE.
% VERSION 10: FINAL & COMPLETE. Corrects the save command for 'All_HMM_Results_Collected_Local.mat'
%             to include all variables expected by the original processing script.

function execute_single_hmm_run(full_data_file_path, main_run_output_path, config)
    
    % --- DIARY SETUP ---
    log_base_dir = fullfile(main_run_output_path, 'Logs');
    if ~exist(log_base_dir, 'dir'), mkdir(log_base_dir); end
    timestamp_for_diary = datestr(now, 'yyyy-mm-dd_HHMMSS');
    diary_filename = fullfile(log_base_dir, sprintf('hmm_run_log_%s.txt', timestamp_for_diary));
    if exist(diary_filename, 'file'), delete(diary_filename); end
    diary(diary_filename);
    cleanupObj_diary = onCleanup(@() diary('off'));
    fprintf('--- Starting Automated Single HMM Run ---\n');
    
    try
        %% 0. Setup Paths
        project_base_path = fileparts(mfilename('fullpath'));
        paths.hmm_core = fullfile(project_base_path, 'hmm_core_functions'); 
        paths.helper_fcns = fullfile(project_base_path, 'helper_functions');
        addpath(genpath(paths.hmm_core)); 
        addpath(genpath(paths.helper_fcns));
        path_to_adapted_core_script = fullfile(paths.hmm_core, 'h2mm_wexac_ver1_local.m');

        %% 1. Unpack Config, Load and Prepare Data
        fprintf('\n--- 1. Loading and Preparing Data ---\n');
        loaded_data = load(full_data_file_path);
        data_re_for_hmm = loaded_data.data_for_hmm; % Assumes data is already sanitized
        dt_analysis_sec = loaded_data.dt;
        if length(data_re_for_hmm) > config.MAX_TRAJECTORIES
            data_re_for_hmm = data_re_for_hmm(1:config.MAX_TRAJECTORIES);
        end
        O_channels = 2;

        %% 2. Initial Guess Generation
        InitialGuess_all = InitialGuess_PhotonByPhoton(config.num_initial_guesses, config.Nstates, O_channels, dt_analysis_sec);
        
        %% 3. Save Config and Data
        config_params_to_save = config;
        config_params_to_save.dt_analysis_sec = dt_analysis_sec;
        Guessmodel_options = {'No restrictions (PhotonByPhoton)', 'Chain', 'three colors chain', 'other (e.g., special with fixS)'};
        config_params_to_save.guessName = Guessmodel_options{config.indx_guessmodel};
        save(fullfile(main_run_output_path, 'data_re_for_hmm.mat'), 'data_re_for_hmm', '-v7.3');
        save(fullfile(main_run_output_path, 'HMM_Run_Config.mat'), 'config_params_to_save');
        
        %% 4. Execute Main HMM Analysis
        fprintf('\n--- 4. Executing Main HMM Analysis ---\n');
        
        % Use a short, robust temporary path for parallel workers
        base_temp_dir = 'C:\HMM_Temp';
        if ~exist(base_temp_dir, 'dir'), mkdir(base_temp_dir); end
        timestamp_str = datestr(now, 'yyyymmdd_HHMMSS_FFF');
        temp_runs_base_path_main = fullfile(base_temp_dir, ['run_' timestamp_str]);
        mkdir(temp_runs_base_path_main);
        cleanupObj = onCleanup(@() rmdir(temp_runs_base_path_main, 's'));
        
        main_block_config.max_iter = config.max_iter_hmm;
        main_block_config.run_guesses_parallel_flag = config.run_guesses_parallel_flag;
        main_block_config.par_hmm_internal_flag_Estep = (config.run_guesses_parallel_flag == 0);
        main_block_config.obs_fix_flag = 1;
        main_block_config.analysis_mode_An = config.analysis_mode_An;
        main_block_config.MaxRunTime_cluster_hours = config.max_runtime_hmm_hours;
        
        [best_model_overall, all_hmm_results_collected, overall_hmm_execution_duration_sec, all_guess_durations_sec, all_final_log_likelihoods] = ...
            execute_hmm_run_block(data_re_for_hmm, InitialGuess_all, main_block_config, ...
                                  path_to_adapted_core_script, temp_runs_base_path_main);
        
        %% 5. Save Results
        fprintf('\n--- 5. Saving HMM Results ---\n');
        
        % <<< --- THIS IS THE FIX --- >>>
        % The previous version was missing several variables in this save command.
        % This new version saves the complete list, matching the original script's output.
        final_results_filename = fullfile(main_run_output_path, 'All_HMM_Results_Collected_Local.mat');
        save(final_results_filename, 'all_hmm_results_collected', 'all_final_log_likelihoods', 'InitialGuess_all', ...
             'config_params_to_save', 'all_guess_durations_sec', 'overall_hmm_execution_duration_sec', '-v7.3');
        fprintf('All results saved to: %s\n', final_results_filename);
        % <<< --- END FIX --- >>>

        if ~isempty(best_model_overall)
            [best_logL_main, best_idx_main] = max(all_final_log_likelihoods);
            best_initial_guess_struct_main = InitialGuess_all{best_idx_main};
            best_model_duration_sec_main = all_guess_durations_sec(best_idx_main);
            
            save(fullfile(main_run_output_path, 'Best_HMM_Model_Local.mat'), ...
                 'best_model_overall', 'best_idx_main', 'best_logL_main', ...
                 'best_initial_guess_struct_main', 'config_params_to_save', ...
                 'best_model_duration_sec_main', '-v7.3');
            fprintf('Best model saved to: Best_HMM_Model_Local.mat\n');
        else
            fprintf('WARNING: No valid best model found.\n');
        end
        fprintf('\n--- Single HMM Run Finished Successfully ---\n');
        
    catch ME_runhmm 
        fprintf(2, '\n--- ERROR in Single HMM Run ---\n');
        rethrow(ME_runhmm);
    end
end