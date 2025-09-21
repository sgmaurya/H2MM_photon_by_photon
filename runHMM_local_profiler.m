% --- START OF FILE runHMM_local_profiler.m ---

% runHMM_local_profiler.m
% Main script for running H2MM analysis locally, including optional pilot run.
% Modified to remove interactive dialogs for profiling purposes.

function runHMM_local_profiler() % <-- NAME CHANGED HERE
    clearvars -except break_debug; 
    close all; 
    clc;       
    
    % --- DIARY SETUP ---
    project_base_path_for_log = fileparts(mfilename('fullpath')); 
    log_base_dir = fullfile(project_base_path_for_log, 'HMM_Run_Logs');
    if ~exist(log_base_dir, 'dir'), mkdir(log_base_dir); end
    timestamp_for_diary = datestr(now, 'yyyy-mm-dd_HHMMSS'); % Used for diary and output folder names
    diary_filename = fullfile(log_base_dir, sprintf('hmm_run_log_profiler_%s.txt', timestamp_for_diary));
    if exist(diary_filename, 'file'), try delete(diary_filename); catch; end; end
    diary(diary_filename);
    cleanupObj_diary = onCleanup(@() diary('off'));
    fprintf('--- Starting Local HMM Analysis Pipeline (Profiler Mode) ---\n');
    fprintf('INFO: Command window output is being logged to: %s\n', diary_filename);
    % --- END DIARY SETUP ---

    try % Main try block
        %% --- PROFILER CONFIGURATION --- %%
        % Set all parameters here to avoid manual input during profiling.
        
        % 1. Set the full path to your SSSdata.mat file
        data_file_HMM = 'C:\Users\Satya\Documents\MATLAB\h2mm_wexac\MyLocalHMMSuite\DATA\SSSdata_spikesFiltered2_13-02-2023_15_13.mat'; % <-- IMPORTANT: SET YOUR FILE PATH HERE
        
        % 2. Main HMM Job Parameters
        main_run_params.Nstates = 2;
        main_run_params.num_initial_guesses = 8;
        main_run_params.dt_instrument_ns = 100; % Will be overwritten by time_unit from file if available
        main_run_params.analysis_mode_An = 1;   % 1=Normal, 2=Special
        main_run_params.run_guesses_parallel_flag = 1; % 1=Yes, 0=No
        main_run_params.max_iter_hmm = 50000;
        main_run_params.max_runtime_hmm_hours = 24;
        main_run_params.MAX_TRAJECTORIES = 1000; % Max trajectories to process

        % 3. HMM Model Topology Selection
        % 1='No restrictions (PhotonByPhoton)', 2='Chain', 3='three colors chain', 4='other (e.g., special with fixS)'
        selected_topology_index = 1; 

        % 4. Pilot Run Configuration
        run_pilot_hmm = false; % Set to 'true' to run the pilot, 'false' to skip it.
        pilot_run_params.subset_percent = 10;
        pilot_run_params.num_guesses = 5;
        pilot_run_params.max_iter = 1000;
        pilot_run_params.run_guesses_parallel = 0; % 1=Yes, 0=No

        % 5. 'fixS' parameter for 'other' topology (analysis_mode_An = 2)
        % This is only used if selected_topology_index is 4 or analysis_mode_An is 2
        fixS_value_for_special_mode = 2; 

        fprintf('--- PROFILER CONFIG LOADED ---\n');
        fprintf('Data file: %s\n', data_file_HMM);
        fprintf('Main run: %d states, %d guesses.\n', main_run_params.Nstates, main_run_params.num_initial_guesses);
        fprintf('Pilot run enabled: %s\n', mat2str(run_pilot_hmm));
        fprintf('--------------------------------\n');
        %% --- END PROFILER CONFIGURATION --- %%


        %% 0. Setup Paths and Configuration
        project_base_path = fileparts(mfilename('fullpath')); 
        fprintf('INFO: Project base path determined as: %s\n', project_base_path);
        paths.hmm_core = fullfile(project_base_path, 'hmm_core_functions'); 
        paths.initial_guess = fullfile(project_base_path, 'initial_guess_functions'); 
        paths.initial_guess_from_pilot = fullfile(project_base_path, 'initial_guess_from_pilot_functions'); % For InitialGuess_FromPilot.m
        paths.helper_fcns = fullfile(project_base_path, 'helper_functions'); % For execute_hmm_run_block.m
        paths.output_base = fullfile(project_base_path, 'HMM_Local_Results'); 
        
        % Create helper directories if they don't exist and add to path
        if ~isfolder(paths.hmm_core), error('Required directory not found: %s', paths.hmm_core); end
        if ~isfolder(paths.initial_guess), error('Required directory not found: %s', paths.initial_guess); end
        if ~isfolder(paths.initial_guess_from_pilot), mkdir(paths.initial_guess_from_pilot); fprintf('INFO: Created directory: %s\n', paths.initial_guess_from_pilot); end
        if ~isfolder(paths.helper_fcns), mkdir(paths.helper_fcns); fprintf('INFO: Created directory: %s\n', paths.helper_fcns); end

        original_normalise_path = fullfile(paths.hmm_core, 'normalise.m');
        local_normalise_path = fullfile(paths.hmm_core, 'normalise_local.m'); 
        if ~exist(local_normalise_path, 'file') && exist(original_normalise_path, 'file')
            fprintf('INFO: Attempting to create local version of normalise.m (normalise_local.m)...\n');
            copyfile(original_normalise_path, local_normalise_path);
            try
                content = fileread(local_normalise_path);
                content_modified = regexprep(content, 'repmatC\(', 'repmat('); 
                if ~strcmp(content, content_modified)
                    fid = fopen(local_normalise_path, 'w'); fprintf(fid, '%s', content_modified); fclose(fid);
                    fprintf('INFO: Successfully created and modified normalise_local.m.\n');
                else
                    fprintf('INFO: normalise_local.m created, but "repmatC(" was not found. Using as is.\n');
                end
            catch ME_norm, warning('NORMALISE_MODIFY_FAIL', 'Could not auto-modify %s. Error: %s', local_normalise_path, ME_norm.message); end
        elseif ~exist(local_normalise_path, 'file') && ~exist(original_normalise_path, 'file'), warning('NORMALISE_NOT_FOUND', 'normalise.m not found in %s.', paths.hmm_core);
        elseif exist(local_normalise_path, 'file'), fprintf('INFO: Using existing normalise_local.m.\n');
        else, fprintf('INFO: Using original normalise.m (as normalise_local.m was not created/found).\n'); end

        addpath(genpath(paths.hmm_core)); 
        addpath(genpath(paths.initial_guess)); 
        addpath(genpath(paths.initial_guess_from_pilot));
        addpath(genpath(paths.helper_fcns));
        fprintf('INFO: Added HMM function paths to MATLAB search path.\n');
        
        if ~exist(paths.output_base, 'dir'), mkdir(paths.output_base); end
        path_to_adapted_core_script = fullfile(paths.hmm_core, 'h2mm_wexac_ver1_local.m');
        if ~exist(path_to_adapted_core_script, 'file'), error('Adapted core HMM script h2mm_wexac_ver1_local.m not found in %s.', paths.hmm_core); end
        if ~exist('execute_hmm_run_block.m', 'file')
            warning('EXECUTE_HMM_BLOCK_NOT_FOUND', 'Helper function execute_hmm_run_block.m not found. Please ensure it is in a directory on the MATLAB path (e.g., %s).', paths.helper_fcns);
        end


        %% 1. Parameter Gathering & Data Loading
        fprintf('\n--- 1. Gathering Parameters and Loading Data ---\n');
        
        full_data_file_path = data_file_HMM;
        if ~exist(full_data_file_path, 'file'), error('Data file not found at the specified path: %s', full_data_file_path); end
        fprintf('Loading data from: %s\n', full_data_file_path);
        try data_loaded = load(full_data_file_path); 
        catch ME_load, diary off; error('Failed to load data file: %s. Error: %s', full_data_file_path, ME_load.message); end
        
        Sdata_for_analysis = {}; Sdata3_for_analysis = {}; burst_data_for_analysis = {}; time_unit_sec_from_file = [];
        if isfield(data_loaded, 'Sdata'), Sdata_for_analysis = extract_trajectory_cell('Sdata', data_loaded.Sdata); else, warning('Variable "Sdata" not found in loaded file.'); end
        if isfield(data_loaded, 'burst_data'), burst_data_for_analysis = extract_trajectory_cell('burst_data', data_loaded.burst_data); else, fprintf('INFO: Variable "burst_data" not found in loaded file.\n'); end
        if isfield(data_loaded, 'Sdata3'), Sdata3_for_analysis = extract_trajectory_cell('Sdata3', data_loaded.Sdata3); else, warning('Variable "Sdata3" not found in loaded file.'); end
        if isfield(data_loaded, 'time_unit'), time_unit_sec_from_file = data_loaded.time_unit; end
        
        default_deltaT_ns = main_run_params.dt_instrument_ns; % Use value from config
        if ~isempty(time_unit_sec_from_file) && isnumeric(time_unit_sec_from_file) && isscalar(time_unit_sec_from_file)
            default_deltaT_ns = round(time_unit_sec_from_file * 1e9); 
            fprintf('INFO: Instrument resolution updated to %d ns based on loaded file.\n', default_deltaT_ns);
        end
        
        try
            Nstates = main_run_params.Nstates;
            num_initial_guesses = main_run_params.num_initial_guesses;
            dt_instrument_ns = default_deltaT_ns; % Use the potentially updated value
            analysis_mode_An = main_run_params.analysis_mode_An;
            run_guesses_parallel_flag = main_run_params.run_guesses_parallel_flag;
            max_iter_hmm = main_run_params.max_iter_hmm;
            max_runtime_hmm_hours = main_run_params.max_runtime_hmm_hours;
            MAX_TRAJECTORIES = main_run_params.MAX_TRAJECTORIES;

            if any(isnan([Nstates, num_initial_guesses, dt_instrument_ns, analysis_mode_An, run_guesses_parallel_flag, max_iter_hmm, max_runtime_hmm_hours, MAX_TRAJECTORIES]))
                error('Invalid numeric input in HMM parameters.'); end
            if Nstates <=0 || num_initial_guesses <=0 || dt_instrument_ns <=0 || max_iter_hmm <=0 || max_runtime_hmm_hours <=0 || MAX_TRAJECTORIES <=0
                error('State count, guess count, resolution, iterations, runtime, and max trajectories must be positive.'); end
            if analysis_mode_An ~= 1 && analysis_mode_An ~= 2, error('Analysis mode must be 1 (Normal) or 2 (Special).'); end
        catch ME_parse, diary off; error('Error parsing HMM parameters: %s Ensure all inputs are valid numbers.', ME_parse.message); end
        
        dt_analysis_sec = dt_instrument_ns * 1e-9;
        if analysis_mode_An == 2, fprintf('INFO: Analysis Mode 2 (Special) selected for main run.\n'); end

        Guessmodel_options = [{'No restrictions (PhotonByPhoton)'}, {'Chain'}, {'three colors chain'}, {'other (e.g., special with fixS)'}];
        indx_guessmodel = selected_topology_index;
        guessName = Guessmodel_options{indx_guessmodel};
        fprintf('Selected HMM Topology for main run: %s\n', guessName);

        %% 1.A. Select Primary Data Source & Truncate (for main analysis)
        fprintf('\n--- 1.A. Selecting and Preparing Data Source for Main Analysis ---\n');
        data_source_for_hmm = {}; 
        if indx_guessmodel == 3 
            if ~isempty(Sdata3_for_analysis), data_source_for_hmm = Sdata3_for_analysis; fprintf('Using 3-color data source (Sdata3_for_analysis with %d trajectories).\n', length(data_source_for_hmm));
            else, diary off; error('3-color model selected, but Sdata3_for_analysis is empty or not found/processed correctly. Check SSSdata.mat.'); end
        else 
            if ~isempty(burst_data_for_analysis), data_source_for_hmm = burst_data_for_analysis; fprintf('Using "burst_data_for_analysis" as primary data source (%d trajectories).\n', length(data_source_for_hmm));
            elseif ~isempty(Sdata_for_analysis), data_source_for_hmm = Sdata_for_analysis; fprintf('Using "Sdata_for_analysis" as primary data source (burst_data was empty/invalid) (%d trajectories).\n', length(data_source_for_hmm));
            else, diary off; error('No suitable default data source (burst_data or Sdata) found or processed correctly. Check SSSdata.mat and ensure variables are correctly named and structured.'); end
        end
        
        if length(data_source_for_hmm) > MAX_TRAJECTORIES
            fprintf('Warning: Data source has %d trajectories, truncating to first %d as per user input for main analysis.\n', length(data_source_for_hmm), MAX_TRAJECTORIES);
            data_source_for_hmm = data_source_for_hmm(1:MAX_TRAJECTORIES);
        end
        if isempty(data_source_for_hmm), diary off; error('Data source for HMM is empty after selection/truncation. Please check your input SSSdata.mat file and processing logic.'); end

        %% 1.B. Preprocess Data (to data_re_for_hmm for main analysis) & Determine O_channels
        fprintf('\n--- 1.B. Preprocessing Data for Main HMM Analysis & Determining Observation Channels ---\n');
        data_re_for_hmm = cell(size(data_source_for_hmm));
        AbsBurstInitialTime = NaN(length(data_source_for_hmm), 1); 
        for i_traj = 1:length(data_source_for_hmm)
            current_traj_data = data_source_for_hmm{i_traj}; 
            if isempty(current_traj_data), data_re_for_hmm{i_traj} = zeros(0,2); continue; end
            if ~isnumeric(current_traj_data) || size(current_traj_data,2) ~= 2 || ndims(current_traj_data) > 2
                warning('RUNHMM_WARNING: Trajectory %d not Nx2 numeric matrix after extraction. Skipping. Size: %s, Class: %s', i_traj, mat2str(size(current_traj_data)), class(current_traj_data));
                data_re_for_hmm{i_traj} = zeros(0,2); continue;
            end
            photon_times_abs = current_traj_data(:,1); photon_channels = current_traj_data(:,2);
            [sorted_times, sort_idx] = sort(photon_times_abs); sorted_channels = photon_channels(sort_idx);
            if ~isempty(sorted_times)
                AbsBurstInitialTime(i_traj) = sorted_times(1);
                relative_times = sorted_times - sorted_times(1); 
                data_re_for_hmm{i_traj} = [relative_times, sorted_channels];
            else, data_re_for_hmm{i_traj} = zeros(0,2); end
        end
        valid_traj_idx = cellfun(@(x) ~isempty(x) && size(x,1)>0, data_re_for_hmm);
        data_re_for_hmm = data_re_for_hmm(valid_traj_idx);
        AbsBurstInitialTime = AbsBurstInitialTime(valid_traj_idx);
        if isempty(data_re_for_hmm), diary off; error('No valid trajectories remaining after preprocessing for main analysis. Check input data format.'); end
        fprintf('INFO: Main data preprocessing complete. %d valid trajectories for HMM analysis.\n', length(data_re_for_hmm));

        max_obs_val = 0;
        for k_traj_chk = 1:length(data_re_for_hmm)
            if ~isempty(data_re_for_hmm{k_traj_chk}) && size(data_re_for_hmm{k_traj_chk},2) >=2
                trajectory_channels = data_re_for_hmm{k_traj_chk}(:,2);
                if ~isempty(trajectory_channels), max_obs_val = max(max_obs_val, max(trajectory_channels)); end
            end
        end
        O_channels = max(1, round(max_obs_val));
        fprintf('INFO: Determined O_channels = %d from processed main data.\n', O_channels);
        
        fixS_for_special = []; % Initialize for main run config saving

        %% 1.5. Pilot Run (Optional)
        best_pilot_model_params = []; 
        pilot_run_performed_successfully = false;
        if run_pilot_hmm
            fprintf('\n--- 1.5. Configuring and Executing Pilot HMM Run ---\n');
                try
                    pilot_subset_percent = pilot_run_params.subset_percent;
                    pilot_num_guesses = pilot_run_params.num_guesses;
                    pilot_max_iter = pilot_run_params.max_iter;
                    pilot_run_guesses_parallel = pilot_run_params.run_guesses_parallel;

                    if isnan(pilot_subset_percent) || pilot_subset_percent <= 0 || pilot_subset_percent > 100 || ...
                       isnan(pilot_num_guesses) || pilot_num_guesses <= 0 || ...
                       isnan(pilot_max_iter) || pilot_max_iter <= 0 || ...
                       isnan(pilot_run_guesses_parallel) || (pilot_run_guesses_parallel ~= 0 && pilot_run_guesses_parallel ~= 1)
                        error('Invalid pilot parameters.');
                    end
                catch ME_pilot_parse, warning('PILOT_PARAM_ERROR', 'Pilot parameter parsing failed or invalid. Skipping pilot run. Error: %s', ME_pilot_parse.message); pilot_subset_percent = 0; end

                if pilot_subset_percent > 0
                    num_pilot_trajs = min(length(data_re_for_hmm), ceil(length(data_re_for_hmm) * pilot_subset_percent / 100));
                    rng('shuffle'); pilot_data_indices = randperm(length(data_re_for_hmm), num_pilot_trajs);
                    pilot_data_re_subset = data_re_for_hmm(pilot_data_indices);
                    
                    fprintf('INFO: Generating %d initial guesses for pilot run using PhotonByPhoton...\n', pilot_num_guesses);
                    pilot_InitialGuesses = InitialGuess_PhotonByPhoton(pilot_num_guesses, Nstates, O_channels, dt_analysis_sec);

                    pilot_block_config.max_iter = pilot_max_iter;
                    pilot_block_config.run_guesses_parallel_flag = pilot_run_guesses_parallel;
                    pilot_block_config.par_hmm_internal_flag_Estep = (pilot_run_guesses_parallel == 0); % Enable E-step parallel if pilot guesses are serial
                    pilot_block_config.obs_fix_flag = 1; 
                    pilot_block_config.analysis_mode_An = analysis_mode_An; 
                    pilot_block_config.MaxRunTime_cluster_hours = max(1, max_runtime_hmm_hours / 4); 
                    if analysis_mode_An == 2
                        if indx_guessmodel == 4 && isempty(fixS_for_special)
                             fixS_for_special_pilot = fixS_value_for_special_mode;
                             fprintf('INFO: Using pre-configured fixS value for pilot: %d\n', fixS_for_special_pilot);
                             
                             if ~isempty(fixS_for_special_pilot) && (fixS_for_special_pilot < 1 || fixS_for_special_pilot > Nstates), error('Invalid fixS for pilot.'); end
                             pilot_block_config.fixS_val = fixS_for_special_pilot;
                             if isempty(fixS_for_special), fixS_for_special = fixS_for_special_pilot; end
                        elseif ~isempty(fixS_for_special)
                            pilot_block_config.fixS_val = fixS_for_special;
                        else
                             if isempty(fixS_for_special)
                                 fixS_for_special_pilot_generic = fixS_value_for_special_mode;
                                 fprintf('INFO: Using pre-configured fixS value for pilot (An=2): %d\n', fixS_for_special_pilot_generic);
                                 
                                 if fixS_for_special_pilot_generic >= 1 && fixS_for_special_pilot_generic <= Nstates
                                     pilot_block_config.fixS_val = fixS_for_special_pilot_generic;
                                     if isempty(fixS_for_special), fixS_for_special = fixS_for_special_pilot_generic; end
                                 end
                             elseif ~isempty(fixS_for_special)
                                 pilot_block_config.fixS_val = fixS_for_special;
                             end
                        end
                    end
                    
                    pilot_temp_base = fullfile(paths.output_base, sprintf('HMM_PilotRun_Temp_%s', timestamp_for_diary));
                    if ~exist(pilot_temp_base, 'dir'), mkdir(pilot_temp_base); end

                    fprintf('INFO: Starting pilot HMM run: %d trajectories (%.1f%%), %d guesses, up to %d iter.\n', ...
                        length(pilot_data_re_subset), pilot_subset_percent, pilot_num_guesses, pilot_max_iter);
                    
                    [pilot_best_model_struct, pilot_all_results, pilot_overall_duration, ~, ~] = ...
                        execute_hmm_run_block(pilot_data_re_subset, pilot_InitialGuesses, pilot_block_config, ...
                                              path_to_adapted_core_script, pilot_temp_base);
                    
                    if ~isempty(pilot_best_model_struct) && isfield(pilot_best_model_struct, 'prior') && ~isempty(pilot_best_model_struct.prior)
                        best_pilot_model_params.prior = pilot_best_model_struct.prior;
                        best_pilot_model_params.transmat = pilot_best_model_struct.transmat;
                        best_pilot_model_params.obsmat = pilot_best_model_struct.obsmat;
                        pilot_run_performed_successfully = true;
                        fprintf('INFO: Pilot run completed successfully. Best model parameters obtained.\n');
                        pilot_run_output_dir = fullfile(paths.output_base, sprintf('HMM_Run_Pilot_Results_%s', timestamp_for_diary));
                        if ~exist(pilot_run_output_dir, 'dir'), mkdir(pilot_run_output_dir); end
                        save(fullfile(pilot_run_output_dir, 'Pilot_Run_Details.mat'), ...
                             'best_pilot_model_params', 'pilot_block_config', 'pilot_all_results', 'pilot_InitialGuesses', 'pilot_overall_duration', 'pilot_data_indices');
                        fprintf('INFO: Pilot run results saved to: %s\n', pilot_run_output_dir);
                    else, fprintf('WARNING: Pilot run did not yield a valid best model.\n'); end
                else, fprintf('INFO: Invalid pilot parameters or subset percentage. Skipping pilot run.\n'); end
        else, fprintf('INFO: Pilot run not selected by user.\n'); end

        %% 2. HMM Model Topology & Initial Guess Generation (Main Run)
        fprintf('\n--- 2. Generating Initial Guesses for Main HMM Run ---\n');
        InitialGuess_all = cell(1, num_initial_guesses); 
        
        if indx_guessmodel == 4 && isempty(fixS_for_special)
            if analysis_mode_An ~= 2, fprintf('Warning: "other" topology selected. Forcing Analysis Mode to 2 for main run.\n'); analysis_mode_An = 2; end
            fixS_for_special = fixS_value_for_special_mode;
            fprintf('INFO: Using pre-configured fixS for main run: %d\n', fixS_for_special);
            
            if fixS_for_special < 1 || fixS_for_special > Nstates, diary off; error('Invalid fixS value for main run.'); end
        elseif analysis_mode_An == 2 && isempty(fixS_for_special)
             fixS_temp = fixS_value_for_special_mode;
             fprintf('INFO: Using pre-configured fixS for main run (An=2): %d\n', fixS_temp);
             if fixS_temp >=1 && fixS_temp <= Nstates, fixS_for_special = fixS_temp; end
             
             if isempty(fixS_for_special)
                 fprintf('INFO: An=2 selected but no valid fixS provided for main run. HMM core might behave unexpectedly if it requires fixS.\n');
             end
        end


        if pilot_run_performed_successfully && ~isempty(best_pilot_model_params)
            fprintf('INFO: Using pilot results to generate %d initial guesses for main run.\n', num_initial_guesses);
            if ~exist('InitialGuess_FromPilot.m', 'file')
                warning('InitialGuess_FromPilot.m not found. Falling back to PhotonByPhoton. Please create InitialGuess_FromPilot.m in %s.', paths.initial_guess_from_pilot);
                InitialGuess_all = InitialGuess_PhotonByPhoton(num_initial_guesses, Nstates, O_channels, dt_analysis_sec);
            else
                perturb_config.rate_noise_factor = 0.25; perturb_config.obs_noise_level = 0.1; perturb_config.prior_noise_level = 0.1;
                InitialGuess_all = InitialGuess_FromPilot(best_pilot_model_params, num_initial_guesses, Nstates, O_channels, dt_analysis_sec, perturb_config);
            end
        else
            if run_pilot_hmm, fprintf('INFO: Pilot run attempted but failed. Generating guesses using base method.\n');
            else, fprintf('INFO: Pilot run not used. Generating %d initial guesses using base method.\n', num_initial_guesses); end
            
            switch indx_guessmodel
                case 1, InitialGuess_all = InitialGuess_PhotonByPhoton(num_initial_guesses, Nstates, O_channels, dt_analysis_sec);
                case 2, warning('RUNHMM:InitialGuess_chain needs update for O_channels, dt_analysis_sec.'); InitialGuess_all = InitialGuess_chain(num_initial_guesses, Nstates); % TODO
                case 3 
                    warning('RUNHMM:InitialGuess_3colors_chain needs update for dt_analysis_sec.'); 
                    if O_channels ~= 3, fprintf('INFO: Forcing O_channels to 3 for "three colors chain" model, determined O_channels=%d.\n', O_channels); end
                    InitialGuess_all = InitialGuess_3colors_chain(num_initial_guesses, Nstates); % TODO (pass 3 as O_channels, and dt_analysis_sec)
                case 4 
                    fprintf('INFO: Using "other" topology for main run. An=%d. fixS=%d.\n', analysis_mode_An, fixS_for_special);
                    warning('RUNHMM:InitialGuess_special needs update for O_channels.'); 
                    InitialGuess_all = InitialGuess_special(num_initial_guesses, Nstates, dt_analysis_sec, fixS_for_special); % TODO (pass O_channels)
                otherwise, diary off; error('Invalid guess model selection index.');
            end
        end
        fprintf('%d initial guesses generated for main run.\n', length(InitialGuess_all));

        %% 3. Create Output Folder & Save Data/Config (Main Run)
        fprintf('\n--- 3. Setting up Output Folders and Saving Data/Config for Main Run ---\n');
        output_run_folder_name = sprintf('HMM_Run_%dstates_%s_An%d_%s', Nstates, strrep(guessName, ' ', '_'), analysis_mode_An, timestamp_for_diary);
        main_run_output_path = fullfile(paths.output_base, output_run_folder_name);
        temp_runs_base_path_main = fullfile(main_run_output_path, 'temp_individual_runs_main'); 
        try mkdir(main_run_output_path); mkdir(temp_runs_base_path_main);
        catch ME_mkdir, diary off; error('Failed to create output directories for main run: %s', ME_mkdir.message); end
        fprintf('Main run output will be saved in: %s\n', main_run_output_path);
        
        config_params_to_save = struct('Nstates', Nstates, 'num_initial_guesses', num_initial_guesses, ...
            'dt_analysis_sec', dt_analysis_sec, 'O_channels_determined', O_channels, ...
            'analysis_mode_An', analysis_mode_An, 'run_guesses_parallel_flag', run_guesses_parallel_flag, ...
            'max_iter_hmm', max_iter_hmm, 'max_runtime_hmm_hours', max_runtime_hmm_hours, ...
            'MAX_TRAJECTORIES_used', MAX_TRAJECTORIES, 'indx_guessmodel', indx_guessmodel, 'guessName', guessName, ...
            'fixS_for_special', fixS_for_special, 'full_data_file_path', full_data_file_path, ... 
            'project_base_path', project_base_path, 'timestamp', timestamp_for_diary, ...
            'pilot_run_attempted', run_pilot_hmm, 'pilot_run_succeeded', pilot_run_performed_successfully);
        save(fullfile(main_run_output_path, 'HMM_Run_Config.mat'), 'config_params_to_save');
        fprintf('Main run configuration saved.\n');
        save(fullfile(main_run_output_path, 'data_re_for_hmm.mat'), 'data_re_for_hmm', 'AbsBurstInitialTime', '-v7.3');
        fprintf('Processed data (data_re_for_hmm) saved for main run.\n');

        %% 4. Execute Main HMM Analysis
        fprintf('\n--- 4. Executing Main HMM Analysis for %d Initial Guesses ---\n', length(InitialGuess_all));
        main_block_config.max_iter = max_iter_hmm;
        main_block_config.run_guesses_parallel_flag = run_guesses_parallel_flag;
        main_block_config.par_hmm_internal_flag_Estep = (run_guesses_parallel_flag == 0 && license('test','distrib_computing_toolbox')); % Enable E-step parallel only if main guesses are serial and toolbox available
        main_block_config.obs_fix_flag = 1;
        main_block_config.analysis_mode_An = analysis_mode_An;
        main_block_config.MaxRunTime_cluster_hours = max_runtime_hmm_hours;
        if analysis_mode_An == 2
            if isempty(fixS_for_special) && indx_guessmodel == 4 % "other" implies fixS is critical
                 error('Main Run: An=2 with "other" topology, but fixS_for_special is undefined.'); 
            elseif ~isempty(fixS_for_special)
                 main_block_config.fixS_val = fixS_for_special;
            end
        end
        
        [best_model_overall, all_hmm_results_collected, overall_hmm_execution_duration_sec, all_guess_durations_sec, all_final_log_likelihoods] = ...
            execute_hmm_run_block(data_re_for_hmm, InitialGuess_all, main_block_config, ...
                                  path_to_adapted_core_script, temp_runs_base_path_main);

        final_results_filename = fullfile(main_run_output_path, 'All_HMM_Results_Collected_Local.mat');
        save(final_results_filename, 'all_hmm_results_collected', 'all_final_log_likelihoods', 'InitialGuess_all', ...
             'config_params_to_save', 'all_guess_durations_sec', 'overall_hmm_execution_duration_sec', '-v7.3');
        fprintf('All main HMM results saved to: %s\n', final_results_filename);

        if ~isempty(best_model_overall)
            best_idx_main = find(all_final_log_likelihoods == max(all_final_log_likelihoods), 1, 'first'); % Robust way to get index
            best_logL_main = all_final_log_likelihoods(best_idx_main);
            best_initial_guess_struct_main = InitialGuess_all{best_idx_main};
            best_model_duration_sec_main = all_guess_durations_sec(best_idx_main);
            save(fullfile(main_run_output_path, 'Best_HMM_Model_Local.mat'), ...
                 'best_model_overall', 'best_idx_main', 'best_logL_main', ...
                 'best_initial_guess_struct_main', 'config_params_to_save', ...
                 'best_model_duration_sec_main', '-v7.3');
            fprintf('Best HMM model from main run saved.\n');
            if isfield(best_model_overall, 'prior') && ~isempty(best_model_overall.prior)
                fprintf('Best Main Model (Guess %d) Params:\n Prior: %s\n TransMat:\n', ...
                        best_idx_main, mat2str(best_model_overall.prior,3)); disp(best_model_overall.transmat);
                fprintf(' ObsMat:\n'); disp(best_model_overall.obsmat); 
                fprintf(' Iterations: %d, Duration: %.2f s\n', ...
                        best_model_overall.lastIteration, best_model_duration_sec_main);
            end
        else, fprintf('\nNo valid models found in the main run.\n'); end
        
        fprintf('\n--- Local HMM Analysis Pipeline Finished Successfully ---\n');

    catch ME_runhmm 
        fprintf('\n--- ERROR in HMM Analysis Pipeline ---\n');
        fprintf('Error ID: %s\nMessage: %s\n', ME_runhmm.identifier, ME_runhmm.message);
        fprintf('Stack Trace:\n'); disp(ME_runhmm.getReport('extended','hyperlinks','off'));
        fprintf('\n--- HMM Analysis Pipeline Terminated Due to Error ---\n');
    end
    
%% ——— Nested Helper Function for Data Extraction ———
    function final_traj_cell = extract_trajectory_cell(input_var_name, loaded_var)
        % ... (Content of extract_trajectory_cell as provided in the original script) ...
        final_traj_cell = {}; 
        if ~iscell(loaded_var) && ~isempty(loaded_var) 
            warning('HELPER_EXTRACT: Variable "%s" is not a cell array as expected. Returning empty.', input_var_name);
            return;
        elseif isempty(loaded_var) 
             fprintf('HELPER_EXTRACT: Variable "%s" is an empty cell. No trajectories to extract.\n', input_var_name);
            return;
        end

        temp_data = loaded_var;
        nesting_level = 0;
        max_nesting = 3; 

        while iscell(temp_data) && numel(temp_data) == 1 && iscell(temp_data{1}) && nesting_level < max_nesting
            temp_data = temp_data{1};
            nesting_level = nesting_level + 1;
        end
        if nesting_level > 0, fprintf('INFO: Un-nested "%s" %d level(s).\n', input_var_name, nesting_level); end

        if iscell(temp_data) && ~isempty(temp_data)
            first_element_peek = temp_data{1}; 
            
            if isnumeric(first_element_peek) && (isempty(first_element_peek) || size(first_element_peek,2) >= 2) 
                raw_trajectories = temp_data;
                processed_trajectories = cell(size(raw_trajectories)); 
                for k_traj = 1:numel(raw_trajectories)
                    current_raw_traj = raw_trajectories{k_traj};
                    if isnumeric(current_raw_traj) && (isempty(current_raw_traj) || size(current_raw_traj,2) >= 2)
                        if ~isempty(current_raw_traj)
                            processed_trajectories{k_traj} = current_raw_traj(:, 1:2); 
                        else
                            processed_trajectories{k_traj} = zeros(0,2); 
                        end
                    elseif ~isempty(current_raw_traj) 
                         warning('HELPER_EXTRACT: Trajectory element within "%s" at index %d is not numeric or has < 2 columns. Skipping. Class: %s, Size: %s', input_var_name, k_traj, class(current_raw_traj), mat2str(size(current_raw_traj)));
                         processed_trajectories{k_traj} = zeros(0,2); 
                    else 
                         processed_trajectories{k_traj} = zeros(0,2); 
                    end
                end
                valid_indices_filter = cellfun(@(x) ~isempty(x) || (isempty(x) && size(x,2)==2 && size(x,1)==0), processed_trajectories);
                final_traj_cell = processed_trajectories(valid_indices_filter);
                fprintf('INFO: "%s" (processed as direct cell array of matrices) -> %d valid trajectories.\n', input_var_name, length(final_traj_cell));

            elseif iscell(first_element_peek)
                fprintf('INFO: "%s" (after any un-nesting) appears to be a cell array of cell arrays. Concatenating inner cells...\n', input_var_name);
                try
                    concatenated_trajectories = {};
                    for k_outer = 1:numel(temp_data) 
                        inner_cell_array = temp_data{k_outer};
                        if iscell(inner_cell_array)
                            for k_inner = 1:numel(inner_cell_array) 
                                current_raw_traj_inner = inner_cell_array{k_inner};
                                if isnumeric(current_raw_traj_inner) && (isempty(current_raw_traj_inner) || size(current_raw_traj_inner,2) >= 2) 
                                    if ~isempty(current_raw_traj_inner)
                                        concatenated_trajectories{end+1} = current_raw_traj_inner(:, 1:2); %#ok<AGROW>
                                    else
                                        concatenated_trajectories{end+1} = zeros(0,2); %#ok<AGROW>
                                    end
                                elseif ~isempty(current_raw_traj_inner)
                                    warning('HELPER_EXTRACT: Inner cell element of "%s" at outer index %d, inner index %d is not numeric or has < 2 columns. Skipping. Class: %s, Size: %s',input_var_name, k_outer, k_inner, class(current_raw_traj_inner), mat2str(size(current_raw_traj_inner)));
                                end
                            end
                        elseif ~isempty(inner_cell_array) 
                             warning('HELPER_EXTRACT: Outer cell element of "%s" at index %d is not a cell array of trajectories. Skipping this batch. Class: %s',input_var_name, k_outer, class(inner_cell_array));
                        end
                    end
                    final_traj_cell = concatenated_trajectories;
                    fprintf('INFO: Concatenated inner cells from "%s" into %d total trajectories (taking first 2 columns).\n', input_var_name, length(final_traj_cell));
                catch ME_concat
                    warning('HELPER_EXTRACT: Failed to concatenate inner cells for "%s". Error: %s. Using empty.', input_var_name, ME_concat.message);
                    final_traj_cell = {};
                end
            else 
                if ~isempty(temp_data) && ~isempty(first_element_peek) 
                    warning('HELPER_EXTRACT: Format of "%s" (after any un-nesting) not recognized. First element class: %s. Expected numeric matrix (>=2 cols) or another cell.', input_var_name, class(first_element_peek));
                elseif isempty(temp_data)
                     fprintf('INFO: "%s" became an empty cell after un-nesting attempts.\n', input_var_name);
                end
                final_traj_cell = {};
            end
        elseif iscell(temp_data) && isempty(temp_data) 
            fprintf('INFO: "%s" (after any un-nesting) is an empty cell array.\n', input_var_name);
            final_traj_cell = {}; 
        else 
             warning('HELPER_EXTRACT: Variable "%s" is not a cell array after initial un-nesting attempts, or it is an unsupported structure. Final class: %s', input_var_name, class(temp_data));
             final_traj_cell = {};
        end
    end % End of nested helper function extract_trajectory_cell
end % --- END OF MAIN FUNCTION runHMM_local_profiler ---```