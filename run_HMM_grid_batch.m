% run_HMM_grid_batch.m
% A controller script to automate H2MM analysis for the entire data grid.
% VERSION 2: Corrected all Python syntax errors to native MATLAB code.

function run_HMM_grid_batch()
    clear; close all; clc;
    
    fprintf('--- Starting Automated HMM Batch Processing for Data Grid ---\n');
    
    %% --- USER CONFIGURATION ---
    % ---------------------------------------------------------------------
    % 1. DEFINE BASE PATH for the Photon Data Grid.
    base_data_path = 'Y:\Photon_Data_Grid';

    % 2. DEFINE HMM PARAMETERS
    config.Nstates = 2;
    config.num_initial_guesses = 8;
    config.max_iter_hmm = 50000;
    config.max_runtime_hmm_hours = 24;
    config.MAX_TRAJECTORIES = 20000; % Set high to process all bursts
    config.run_guesses_parallel_flag = 1; 
    config.indx_guessmodel = 1; 
    config.analysis_mode_An = 1; 
    config.fixS_for_special = []; 
    config.run_pilot = false;
    % ---------------------------------------------------------------------
    
    %% --- Script Core Logic ---
    worker_script_name = 'execute_single_hmm_run.m';
    if ~exist(worker_script_name, 'file'), error('Worker script "%s" not found.', worker_script_name); end
    
    % --- A. Find all the top-level 'TPT_Explore_*' folders ---
    tpt_system_folders_struct = dir(fullfile(base_data_path, 'TPT_Explore_*'));
    tpt_system_folders_struct = tpt_system_folders_struct([tpt_system_folders_struct.isdir]); % Filter for directories only
    
    if isempty(tpt_system_folders_struct)
        error('No "TPT_Explore_*" folders found in %s. Nothing to process.', base_data_path);
    end

    num_tpt_systems = length(tpt_system_folders_struct);
    fprintf('Found %d TPT system folders to process.\n', num_tpt_systems);

    % --- B. Outer loop to process each TPT system folder ---
    for i_tpt = 1:num_tpt_systems
        tpt_folder_info = tpt_system_folders_struct(i_tpt);
        tpt_folder_name = tpt_folder_info.name;
        tpt_folder_path = fullfile(tpt_folder_info.folder, tpt_folder_name);
        
        fprintf('\n%s\nProcessing TPT System %d/%d: %s\n%s\n', ...
                repmat('=', 1, 80), i_tpt, num_tpt_systems, tpt_folder_name, repmat('=', 1, 80));

        % Find all 'data_for_hmm.mat' files within this TPT system's subdirectories
        data_file_pattern = fullfile(tpt_folder_path, 'flux_*_hz', 'HMM_Ready_Data', 'data_for_hmm.mat');
        data_files = dir(data_file_pattern);
        
        if isempty(data_files)
            fprintf('Warning: No chunked data files found for TPT system %s. Skipping.\n', tpt_folder_name);
            continue;
        end
        
        num_flux_conditions = length(data_files);
        fprintf('Found %d flux conditions to analyze for this TPT system.\n', num_flux_conditions);
        
        % --- C. Inner loop to process each flux condition ---
        for i_file = 1:num_flux_conditions
            file_info = data_files(i_file);
            full_data_path = fullfile(file_info.folder, file_info.name);
            
            fprintf('\n--- Analyzing Flux Condition %d of %d: %s ---\n', i_file, num_flux_conditions, file_info.folder);
            
            % --- Create the descriptive output folder name ---
            Guessmodel_options = {'PBP', 'Chain', '3CChain', 'Special'};
            guessName_str = Guessmodel_options{config.indx_guessmodel};
            timestamp_str = datestr(now, 'yymmdd_HHMMSS'); % Use a longer timestamp to avoid collision

            output_folder_name = sprintf('HMM_Res_%ds_%s_An%d_%s', ...
                                         config.Nstates, guessName_str, ...
                                         config.analysis_mode_An, timestamp_str);

            output_dir = fullfile(file_info.folder, output_folder_name);
            if ~exist(output_dir, 'dir'), mkdir(output_dir); end
            
            try
                % Call the same worker script as before
                execute_single_hmm_run(full_data_path, output_dir, config);
                fprintf('--- Successfully completed HMM analysis. Results are in: %s ---\n', output_dir);
            catch ME
                fprintf(2, '--- ERROR processing file: %s ---\n', full_data_path);
                rethrow(ME);
            end
        end
    end
    
    fprintf('\n%s\n--- HMM Grid Batch Processing Complete ---\n', repmat('*', 1, 80));
end