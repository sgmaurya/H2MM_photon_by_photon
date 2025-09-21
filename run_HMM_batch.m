% run_HMM_batch.m
% Controller script to automate H2MM analysis for all flux-sweep data folders.
% VERSION 3: ADDS a timestamp to the output folder name to prevent overwriting results.

function run_HMM_batch()
    clear; close all; clc;
    
    fprintf('--- Starting Automated HMM Batch Processing ---\n');
    
    %% --- USER CONFIGURATION ---
    base_data_path = 'Y:\Flux_Sweep_Validation_Memory_Efficient\New folder';

    % Define HMM Parameters
    config.Nstates = 2;
    config.num_initial_guesses = 8;
    config.max_iter_hmm = 50000;
    config.max_runtime_hmm_hours = 24;
    config.MAX_TRAJECTORIES = 1000000; 
    config.run_guesses_parallel_flag = 1; 
    config.indx_guessmodel = 1; 
    config.analysis_mode_An = 1; 
    config.fixS_for_special = []; 
    config.run_pilot = false;
    % ... (pilot config can be added here if needed)
    % ---------------------------------------------------------------------
    
    %% --- Script Core Logic ---
    worker_script_name = 'execute_single_hmm_run.m';
    if ~exist(worker_script_name, 'file'), error('Worker script "%s" not found.', worker_script_name); end
    
    data_file_pattern = fullfile(base_data_path, 'flux_*_hz', 'HMM_Ready_Data', 'data_for_hmm.mat');
    data_files = dir(data_file_pattern);
    
    if isempty(data_files), error('No "data_for_hmm.mat" files found using pattern:\n%s', data_file_pattern); end
    
    num_files = length(data_files);
    fprintf('Found %d data files to process.\n', num_files);
    
    for i = 1:num_files
        file_info = data_files(i);
        full_data_path = fullfile(file_info.folder, file_info.name);
        
        fprintf('\n%s\n--- Processing File %d of %d ---\n%s\n', repmat('=', 1, 80), i, num_files, repmat('=', 1, 80));
        
        % --- MODIFICATION: ADD TIMESTAMP TO FOLDER NAME ---
        Guessmodel_options_abbrev = {'PBP', 'Chain', '3CChain', 'Special'};
        guessName_str_abbrev = Guessmodel_options_abbrev{config.indx_guessmodel};
        timestamp_str_abbrev = datestr(now, 'yymmdd_HHMM'); % Shorter timestamp

        output_folder_name = sprintf('HMM_Res_%ds_%s_An%d_%s', ...
                                     config.Nstates, ...
                                     guessName_str_abbrev, ...
                                     config.analysis_mode_An, ...
                                     timestamp_str_abbrev);
        % --- END MODIFICATION ---

        output_dir = fullfile(file_info.folder, output_folder_name);
        if ~exist(output_dir, 'dir'), mkdir(output_dir); end
        
        try
            execute_single_hmm_run(full_data_path, output_dir, config);
            fprintf('--- Successfully completed analysis. Results are in: %s ---\n', output_dir);
        catch ME
            fprintf(2, '--- ERROR processing file: %s ---\n', full_data_path);
            rethrow(ME);
        end
    end
    
    fprintf('\n--- Batch Processing Complete ---\n');

end