% local_hmm_runner.m
% Manages the execution of a single HMM optimization run for one initial guess.
% Called by runHMM_local.m.

function [model_results, saved_output_filename_full_path, duration_sec] = local_hmm_runner(run_config)
    % run_config: A struct containing all necessary parameters for a single HMM run.
    % Expected fields in run_config:
    %   .temp_dir_base             % Base path for temporary directories (e.g., .../HMM_Local_Results/HMM_Run_.../temp_individual_runs/)
    %   .script_to_run_path        % Full path to h2mm_wexac_ver1_local.m
    %   .current_guess_index       % Index for the current initial guess
    %   .all_initial_guesses       % Cell array of all InitialGuess structs
    %   .trajectory_data           % Cell array of processed trajectory data (data_re)
    %   .max_iter                  % Max iterations for HMM
    %   .par_hmm_internal_flag     % 0 or 1, for the 'par' variable inside h2mm_main
    %   .obs_fix_flag              % For 'OBS' in varargin
    %   .analysis_mode_An          % Analysis mode (1 for normal, 2 for special)
    %   .MaxRunTime_cluster_hours  % Max runtime for HMM convergence check
    %   .fixS_val                  % (Optional) Value for 'fixed_state' if An=2
    %   .H_val                     % (Optional) Value for 'H' parameter
    %   .range_FRET_val            % (Optional) Value for 'range_FRET' parameter
    %   .control_parpool_externally% (Optional, boolean) If true, h2mm_wexac_ver1_local.m won't manage its own parpool

    % Declare globals that h2mm_wexac_ver1_local.m will use
    global LOCAL_JOB_INDEX_FOR_HMM_CORE;
    global LOCAL_PAR_CONTROL_FOR_HMM_CORE;

    % Construct the specific temporary directory for this guess's run
    temp_run_dir = fullfile(run_config.temp_dir_base, sprintf('guess_%04d_run', run_config.current_guess_index));
    
    % Initialize outputs
    saved_output_filename_full_path = '';
    model_results = struct(); 
    duration_sec = NaN; % Initialize duration

    original_dir = pwd; % Save current directory before changing
    fprintf('RUNNER_INFO (Guess %d): Original directory is: %s\n', run_config.current_guess_index, original_dir);
    
    runner_start_time = tic; % Start timer for this entire guess run

    try
        % Create the temporary directory if it doesn't exist
        if ~exist(temp_run_dir, 'dir')
            mkdir(temp_run_dir);
            fprintf('RUNNER_INFO (Guess %d): Created temp directory: %s\n', run_config.current_guess_index, temp_run_dir);
        end
        
        % Change to the temporary directory for this specific HMM run
        cd(temp_run_dir);
        fprintf('RUNNER_INFO (Guess %d): Changed to temp run directory: %s\n', run_config.current_guess_index, pwd);

        % 1. Prepare 'varargin' cell array for h2mm_wexac_ver1_local.m
        varargin_list = {};
        varargin_list{end+1} = 'data';
        varargin_list{end+1} = run_config.trajectory_data;
        varargin_list{end+1} = 'InitialGues'; 
        varargin_list{end+1} = run_config.all_initial_guesses; 
        varargin_list{end+1} = 'number_of_initialguess';
        varargin_list{end+1} = length(run_config.all_initial_guesses);
        varargin_list{end+1} = 'max_iter';
        varargin_list{end+1} = run_config.max_iter;
        varargin_list{end+1} = 'MaxRunTime';
        varargin_list{end+1} = run_config.MaxRunTime_cluster_hours;
        varargin_list{end+1} = 'OBS';
        varargin_list{end+1} = run_config.obs_fix_flag;
        varargin_list{end+1} = 'par';
        varargin_list{end+1} = run_config.par_hmm_internal_flag;

        if run_config.analysis_mode_An == 2 && isfield(run_config, 'fixS_val') && ~isempty(run_config.fixS_val)
            varargin_list{end+1} = 'fixed_state';
            varargin_list{end+1} = run_config.fixS_val;
        end
        if isfield(run_config, 'H_val') && ~isempty(run_config.H_val)
            varargin_list{end+1} = 'H';
            varargin_list{end+1} = run_config.H_val;
        end
        if isfield(run_config, 'range_FRET_val') && ~isempty(run_config.range_FRET_val)
            varargin_list{end+1} = 'range_FRET';
            varargin_list{end+1} = run_config.range_FRET_val;
        end
        
        varargin_to_save = varargin_list; 
        save('varargin.mat', 'varargin_to_save', '-v7.3'); 
        
        fprintf('RUNNER_DEBUG (Guess %d): Attempted to save varargin.mat in directory: %s\n', run_config.current_guess_index, pwd);
        if exist(fullfile(pwd, 'varargin.mat'), 'file')
            fprintf('RUNNER_DEBUG (Guess %d): varargin.mat CONFIRMED to exist in current directory.\n', run_config.current_guess_index);
        else
            fprintf('RUNNER_CRITICAL_ERROR (Guess %d): varargin.mat was NOT found in current directory (%s) after save!\n', run_config.current_guess_index, pwd);
            error('RUNNER_ERROR: varargin.mat creation failed in temp directory.'); 
        end

        % 2. Set up global variables for h2mm_wexac_ver1_local.m
        LOCAL_JOB_INDEX_FOR_HMM_CORE = run_config.current_guess_index;
        if isfield(run_config,'control_parpool_externally') && run_config.control_parpool_externally
            LOCAL_PAR_CONTROL_FOR_HMM_CORE = false; 
        else
            LOCAL_PAR_CONTROL_FOR_HMM_CORE = true;  
        end
        
        % 3. Run the modified core HMM script
        fprintf('RUNNER_INFO (Guess %d): Executing HMM core script: %s\n', run_config.current_guess_index, run_config.script_to_run_path);
        fprintf('RUNNER_INFO (Guess %d): Current directory before running core script is: %s\n', run_config.current_guess_index, pwd);
        
        if ~exist(run_config.script_to_run_path, 'file')
            error('RUNNER_ERROR (Guess %d): Core HMM script not found at specified path: %s', run_config.current_guess_index, run_config.script_to_run_path);
        end
        
        % Run the script by name, assuming its directory is on the path (added by runHMM_local.m)
        % This ensures it runs in the current 'temp_run_dir'.
        run('h2mm_wexac_ver1_local.m'); 

        fprintf('RUNNER_INFO (Guess %d): HMM core script execution finished.\n', run_config.current_guess_index);

        % 4. Load results from the file saved by h2mm_wexac_ver1_local.m
        raw_output_filename = sprintf('hmm_PhotonByPhotonJobArray%d.mat', run_config.current_guess_index);
        saved_output_filename_full_path = fullfile(pwd, raw_output_filename); 

        if exist(raw_output_filename, 'file')
            fprintf('RUNNER_INFO (Guess %d): Loading results from %s.\n', run_config.current_guess_index, raw_output_filename);
            loaded_results = load(raw_output_filename, 'LL_out_save', 'prior_out_save', 'transmat_out_save', 'obsmat_out_save', 'lastIteration_out_save', 'all_gamma_out_save', 'LastData_struct', 'InitialGuess_for_hmm');
            
            model_results.LL = loaded_results.LL_out_save;
            model_results.prior = loaded_results.prior_out_save;
            model_results.transmat = loaded_results.transmat_out_save;
            model_results.obsmat = loaded_results.obsmat_out_save;
            model_results.lastIteration = loaded_results.lastIteration_out_save;
            model_results.all_gamma = loaded_results.all_gamma_out_save; 
            model_results.LastDataFromRun = loaded_results.LastData_struct; 
            model_results.InitialGuessUsed = loaded_results.InitialGuess_for_hmm; 
        else
            warning('RUNNER_WARNING (Guess %d): Output file %s not found in %s. HMM run might have failed or saved elsewhere.', run_config.current_guess_index, raw_output_filename, pwd);
        end

        cd(original_dir);
        fprintf('RUNNER_INFO (Guess %d): Returned to original directory: %s.\n', run_config.current_guess_index, original_dir);
        
        cleanup_temp_files = false; 
        if cleanup_temp_files && exist(temp_run_dir, 'dir')
            fprintf('RUNNER_INFO (Guess %d): Cleaning up temp directory: %s\n', run_config.current_guess_index, temp_run_dir);
            try, rmdir(temp_run_dir, 's'); 
            catch ME_rmdir, warning('RUNNER_WARNING (Guess %d): Could not remove temp directory %s. Error: %s', run_config.current_guess_index, temp_run_dir, ME_rmdir.message); end
        end

    catch ME_runner
        fprintf('RUNNER_ERROR (Guess %d): An error occurred within local_hmm_runner: %s\n', run_config.current_guess_index, ME_runner.message);
        disp('Stack trace for local_hmm_runner error:');
        disp(ME_runner.getReport()); 
        
        if ~strcmp(pwd, original_dir)
            fprintf('RUNNER_INFO (Guess %d): Attempting to return to original directory %s after error.\n', run_config.current_guess_index, original_dir);
            cd(original_dir);
        end
    end

    duration_sec = toc(runner_start_time); % Calculate duration for this guess run
    fprintf('RUNNER_INFO (Guess %d): Execution time for this guess: %.2f seconds (%.2f minutes).\n', run_config.current_guess_index, duration_sec, duration_sec/60);

    clear global LOCAL_JOB_INDEX_FOR_HMM_CORE;
    clear global LOCAL_PAR_CONTROL_FOR_HMM_CORE;
end