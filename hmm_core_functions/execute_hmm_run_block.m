% --- START OF FILE execute_hmm_run_block.m ---
function [best_model_struct, all_results_collected_block, block_duration_sec, all_guess_durations_block, all_logL_block] = ...
    execute_hmm_run_block(data_for_block, initial_guesses_block, block_config, path_to_adapted_core_script, temp_base_for_runs_block)
    % Executes a block of HMM runs (either pilot or main).
    %
    % INPUTS:
    %   data_for_block:             Cell array of trajectory data for this block.
    %   initial_guesses_block:      Cell array of initial guess structs for this block.
    %   block_config:               Struct with parameters for this block, e.g.:
    %                               .max_iter
    %                               .run_guesses_parallel_flag (0 or 1, for parallelizing over guesses)
    %                               .par_hmm_internal_flag_Estep (0 or 1, for parallelizing E-step inside HMM)
    %                               .obs_fix_flag
    %                               .analysis_mode_An
    %                               .MaxRunTime_cluster_hours
    %                               .fixS_val (optional, if An=2)
    %   path_to_adapted_core_script: Full path to h2mm_wexac_ver1_local.m
    %   temp_base_for_runs_block:   Base path for temporary individual run directories for this block.
    %
    % OUTPUTS:
    %   best_model_struct:          The model struct (from all_results_collected_block) with the highest log-likelihood.
    %   all_results_collected_block: Cell array of all HMM result structs from this block.
    %   block_duration_sec:         Total duration for executing all guesses in this block.
    %   all_guess_durations_block:  Vector of durations for each individual guess.
    %   all_logL_block:             Vector of final log-likelihoods for each guess.

    num_guesses_block = length(initial_guesses_block);
    all_results_collected_block = cell(1, num_guesses_block);
    all_logL_block = NaN(1, num_guesses_block);
    all_guess_durations_block = NaN(1, num_guesses_block);

    fprintf('INFO (Block): Preparing %d iteration run configurations.\n', num_guesses_block);
    iteration_run_configs_block = cell(1, num_guesses_block);
    for i_cfg = 1:num_guesses_block
        run_cfg = struct();
        run_cfg.temp_dir_base = temp_base_for_runs_block;
        run_cfg.script_to_run_path = path_to_adapted_core_script;
        run_cfg.current_guess_index = i_cfg;
        run_cfg.all_initial_guesses = initial_guesses_block; % This is the set of guesses for THIS block
        run_cfg.trajectory_data = data_for_block;
        run_cfg.max_iter = block_config.max_iter;
        
        % Determine parallelization for the HMM E-step (par variable in h2mm_main)
        if block_config.run_guesses_parallel_flag == 1
            run_cfg.par_hmm_internal_flag = 0; % If guesses run in parallel, HMM E-step is serial
        else % Guesses run serially
            if isfield(block_config, 'par_hmm_internal_flag_Estep') && block_config.par_hmm_internal_flag_Estep == 1
                 run_cfg.par_hmm_internal_flag = 1; % HMM E-step runs in parallel
                 % fprintf('INFO (Block Guess %d config): Guesses serial, Inner HMM E-step PARALLEL.\n', i_cfg);
            else
                 run_cfg.par_hmm_internal_flag = 0; % HMM E-step runs serially
                 % fprintf('INFO (Block Guess %d config): Guesses serial, Inner HMM E-step SERIAL.\n', i_cfg);
            end
        end
                                               
        run_cfg.obs_fix_flag = block_config.obs_fix_flag;
        run_cfg.analysis_mode_An = block_config.analysis_mode_An;
        run_cfg.MaxRunTime_cluster_hours = block_config.MaxRunTime_cluster_hours;
        
        if block_config.analysis_mode_An == 2
            if ~isfield(block_config, 'fixS_val') || isempty(block_config.fixS_val)
                error('EXECUTE_HMM_BLOCK: An=2 but fixS_val is not provided in block_config.');
            end
            run_cfg.fixS_val = block_config.fixS_val;
            % Optional: if isfield(block_config, 'range_FRET_val'), run_cfg.range_FRET_val = block_config.range_FRET_val; end
        end
        iteration_run_configs_block{i_cfg} = run_cfg;
    end

    fprintf('INFO (Block): Starting HMM execution loop for %d guesses...\n', num_guesses_block);
    block_overall_start_time = tic;

    pool_shutdown_needed_block = false;
    % Decide if this block of guesses runs in parallel
    run_guesses_in_parallel_this_block = block_config.run_guesses_parallel_flag == 1;

    if run_guesses_in_parallel_this_block && ~isempty(ver('parallel'))
        current_pool_block = gcp('nocreate');
        if isempty(current_pool_block)
            fprintf('INFO (Block): Starting parallel pool for guesses...\n');
            try
                parpool; % Use default profile
                pool_shutdown_needed_block = true; 
            catch ME_parpool_start_block
                warning('EXECUTE_HMM_BLOCK:PARPOOL_START_FAILED', ...
                        'Failed to start parpool for block. Running guesses serially. Error: %s', ME_parpool_start_block.message);
                run_guesses_in_parallel_this_block = 0; % Fallback to serial
            end
        else
            fprintf('INFO (Block): Using existing parallel pool for guesses.\n');
        end
        
        if run_guesses_in_parallel_this_block % If parpool started or existed
            for i_cfg_par = 1:num_guesses_block % Mark that parpool is handled externally for local_hmm_runner
                iteration_run_configs_block{i_cfg_par}.control_parpool_externally = true;
            end
        end
    end

    if run_guesses_in_parallel_this_block
        parfor i_guess_par = 1:num_guesses_block
            worker_id_str = ''; 
            try currentTask = getCurrentTask(); if ~isempty(currentTask), worker_id_str = sprintf('(Worker %d)', currentTask.ID); end; catch; end
            fprintf('INFO (Block): Processing guess %d of %d %s (Parallel)...\n', i_guess_par, num_guesses_block, worker_id_str);
            current_run_cfg_par_block = iteration_run_configs_block{i_guess_par};
            try
                [model_out_par, ~, duration_par] = local_hmm_runner(current_run_cfg_par_block); 
                all_results_collected_block{i_guess_par} = model_out_par; 
                all_guess_durations_block(i_guess_par) = duration_par;
            catch ME_par_run_block
                fprintf('ERROR (Block): Parallel guess %d %s failed: %s\n', i_guess_par, worker_id_str, ME_par_run_block.message);
                % Consider logging ME_par_run_block.getReport()
                all_results_collected_block{i_guess_par} = struct(); % Ensure it's an empty struct on error
                all_guess_durations_block(i_guess_par) = NaN;
            end
        end
    else % Guesses run serially
        for i_guess_ser = 1:num_guesses_block
            fprintf('INFO (Block): Processing guess %d of %d (Serial)...\n', i_guess_ser, num_guesses_block);
            current_run_cfg_ser_block = iteration_run_configs_block{i_guess_ser};
            current_run_cfg_ser_block.control_parpool_externally = false; % local_hmm_runner might manage its own if needed AND if inner E-step is parallel
            try
                [model_out_ser, ~, duration_ser] = local_hmm_runner(current_run_cfg_ser_block);
                all_results_collected_block{i_guess_ser} = model_out_ser; 
                all_guess_durations_block(i_guess_ser) = duration_ser;
            catch ME_ser_run_block
                fprintf('ERROR (Block): Serial guess %d failed: %s\n', i_guess_ser, ME_ser_run_block.message);
                % Consider logging ME_ser_run_block.getReport()
                all_results_collected_block{i_guess_ser} = struct(); % Ensure it's an empty struct on error
                all_guess_durations_block(i_guess_ser) = NaN;
            end
        end
    end
    
    block_duration_sec = toc(block_overall_start_time);
    fprintf('INFO (Block): Total HMM execution time for this block of %d guesses: %.2f seconds.\n', num_guesses_block, block_duration_sec);
    
    for i_res = 1:num_guesses_block
        if ~isempty(all_results_collected_block{i_res}) && isstruct(all_results_collected_block{i_res}) && ...
           isfield(all_results_collected_block{i_res}, 'LL') && ~isempty(all_results_collected_block{i_res}.LL) && ...
           isnumeric(all_results_collected_block{i_res}.LL)
           all_logL_block(i_res) = all_results_collected_block{i_res}.LL(end);
        else
           fprintf('WARNING (Block): Results for guess %d were empty or lacked a valid LogLikelihood.\n', i_res);
           % all_logL_block(i_res) is already NaN
        end
    end

    best_model_struct = []; % Initialize to empty
    if any(~isnan(all_logL_block))
        [~, best_idx_block] = max(all_logL_block); % find index of max logL
        best_model_struct = all_results_collected_block{best_idx_block};
        fprintf('INFO (Block): Best model from this block is from guess %d with LogLikelihood = %f.\n', best_idx_block, all_logL_block(best_idx_block));
    else
        fprintf('WARNING (Block): No valid models found in this block.\n');
    end

    if pool_shutdown_needed_block
        fprintf('INFO (Block): Shutting down parallel pool for this block...\n');
        delete(gcp('nocreate'));
    end
end
% --- END OF FILE execute_hmm_run_block.m ---