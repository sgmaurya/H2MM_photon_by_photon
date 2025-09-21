% process_HMM_batch_final.m
% Final controller script to fully automate post-processing.
% It finds all HMM results folders and calls the non-interactive worker
% 'execute_single_hmm_processing_auto.m' for each one.

function process_HMM_batch_final()
    clear; close all; clc;
    
    fprintf('--- Starting FULLY AUTOMATED HMM Results Post-Processing ---\n');
    
    %% --- USER CONFIGURATION ---
    base_results_path = 'Y:\Flux_Sweep_Validation_Memory_Efficient'; % Using the virtual drive
    % ---------------------------------------------------------------------

    %% --- Script Core Logic ---
    worker_script_name = 'execute_single_hmm_processing_auto.m';
    if ~exist(worker_script_name, 'file')
        error('The automated worker script "%s" was not found.', worker_script_name);
    end
    
    % Use the abbreviated name pattern to find the folders
    results_folder_pattern = fullfile(base_results_path, 'flux_*_hz', 'HMM_Ready_Data', 'HMM_Res_*');
    results_folders = dir(results_folder_pattern);
    results_folders = results_folders([results_folders.isdir]);
    
    if isempty(results_folders), error('No HMM results folders found using pattern:\n%s', results_folder_pattern); end
    
    num_folders = length(results_folders);
    fprintf('Found %d HMM results folders to process.\n', num_folders);
    
    for i = 1:num_folders
        folder_info = results_folders(i);
        full_results_path = fullfile(folder_info.folder, folder_info.name);
        
        fprintf('\n%s\n--- Processing Folder %d of %d: %s ---\n%s\n', ...
                repmat('=', 1, 80), i, num_folders, folder_info.name, repmat('=', 1, 80));
        
        try
            % Call the fully automated worker function
            execute_single_hmm_processing_auto(full_results_path);
            fprintf('--- Successfully completed post-processing for: %s ---\n', folder_info.name);
        catch ME
            fprintf(2, '\n--- ERROR processing folder: %s ---\n', full_results_path);
            rethrow(ME); % Stop the process if one folder fails
        end
    end
    
    fprintf('\n%s\n--- Batch Post-Processing Complete ---\n', repmat('*', 1, 80));
end