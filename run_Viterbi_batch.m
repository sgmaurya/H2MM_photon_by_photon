% run_Viterbi_batch.m
% Final controller script to fully automate the Viterbi decoding step.
% It finds all processed HMM results folders and calls the non-interactive
% worker 'execute_single_viterbi_run.m' for each one.

function run_Viterbi_batch()
    clear; close all; clc;
    
    fprintf('--- Starting FULLY AUTOMATED Viterbi Batch Decoding ---\n');
    
    %% --- USER CONFIGURATION ---
    % ---------------------------------------------------------------------
    % 1. DEFINE BASE PATH for the Flux Sweep Validation data.
    base_results_path = 'Y:\Flux_Sweep_Validation_Memory_Efficient'; % Using the virtual drive

    % 2. DECIDE WHETHER TO VISUALIZE TRAJECTORIES
    %    true  = Show the Traj_Vet_local plots for each run (pauses script)
    %    false = Run completely in the background, no plots (fully automated)
    config.run_visualization = false;
    % ---------------------------------------------------------------------

    %% --- Script Core Logic ---
    worker_script_name = 'execute_single_viterbi_run.m';
    if ~exist(worker_script_name, 'file')
        error('The automated worker script "%s" was not found.', worker_script_name);
    end
    
    % Use the abbreviated name pattern to find the folders
    results_folder_pattern = fullfile(base_results_path, 'flux_*_hz', 'HMM_Ready_Data', 'HMM_Res_*');
    results_folders = dir(results_folder_pattern);
    results_folders = results_folders([results_folders.isdir]);
    
    if isempty(results_folders), error('No HMM results folders found using pattern:\n%s', results_folder_pattern); end
    
    num_folders = length(results_folders);
    fprintf('Found %d processed HMM folders to run Viterbi on.\n', num_folders);
    
    for i = 1:num_folders
        folder_info = results_folders(i);
        full_results_path = fullfile(folder_info.folder, folder_info.name);
        
        fprintf('\n%s\n--- Running Viterbi on Folder %d of %d: %s ---\n%s\n', ...
                repmat('=', 1, 80), i, num_folders, folder_info.name, repmat('=', 1, 80));
        
        try
            % Call the fully automated worker function
            execute_single_viterbi_run(full_results_path, config.run_visualization);
            fprintf('--- Successfully completed Viterbi decoding for: %s ---\n', folder_info.name);

            % If visualizing, pause to allow user to see plots
            if config.run_visualization && i < num_folders
                input('Press Enter to close figures and continue to the next folder...', 's');
                close all;
            end
        catch ME
            fprintf(2, '\n--- ERROR running Viterbi on folder: %s ---\n', full_results_path);
            rethrow(ME);
        end
    end
    
    fprintf('\n%s\n--- Viterbi Batch Decoding Complete ---\n', repmat('*', 1, 80));
end