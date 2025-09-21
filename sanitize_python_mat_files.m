% sanitize_python_mat_files.m
% A dedicated pre-processing script to fix the Python-to-MATLAB "corrupt file"
% save issue. It loads each Python-generated data.mat, rebuilds the cell
% array natively in MATLAB, and saves a sanitized version.
% RUN THIS SCRIPT ONCE before starting your HMM analysis pipeline.

function sanitize_python_mat_files()
    clear; close all; clc;
    
    fprintf('--- Starting MATLAB Data Sanitization Pre-Processing ---\n');
    
    %% --- USER CONFIGURATION ---
    % ---------------------------------------------------------------------
    % DEFINE BASE PATH for the Flux Sweep Validation data.
    % This should point to the main folder containing all the 'flux_*_hz' subfolders.
    base_data_path = 'C:\Users\Satya\Documents\Randon_walk_energy barrier\Asymmetric_barrier_Final_LowBarrier_k1_15_x1_5_x2_6pt5_k2_5pt26_k3_8\Flux_Sweep_Validation_Memory_Efficient';
    % ---------------------------------------------------------------------
    
    %% --- Script Core Logic ---
    
    % Find all the original data files created by the Python chunking script
    data_file_pattern = fullfile(base_data_path, 'flux_*_hz', 'HMM_Ready_Data', 'data_for_hmm.mat');
    data_files = dir(data_file_pattern);
    
    if isempty(data_files)
        error('No "data_for_hmm.mat" files found using the pattern:\n%s\nPlease check your base_data_path.', data_file_pattern);
    end
    
    num_files = length(data_files);
    fprintf('Found %d data files to sanitize.\n', num_files);
    
    success_count = 0;
    
    % Loop through each file, load it, clean it, and save it.
    for i = 1:num_files
        file_info = data_files(i);
        full_data_path = fullfile(file_info.folder, file_info.name);
        
        fprintf('\n--- Processing File %d of %d: %s ---\n', i, num_files, full_data_path);
        
        try
            % Load the potentially problematic data
            fprintf('Loading original file...\n');
            loaded_data = load(full_data_path);
            data_re_for_hmm_loaded = loaded_data.data_for_hmm;
            dt = loaded_data.dt; % Preserve the dt variable
            
            % --- The "Deep Clean" Sanitization ---
            fprintf('Sanitizing data structure...\n');
            num_trajectories = numel(data_re_for_hmm_loaded);
            data_for_hmm = cell(size(data_re_for_hmm_loaded)); % Create new, clean variable
            for k = 1:num_trajectories
                if ~isempty(data_re_for_hmm_loaded{k})
                    data_for_hmm{k} = double(data_re_for_hmm_loaded{k});
                else
                    data_for_hmm{k} = [];
                end
            end
            clear data_re_for_hmm_loaded; % Free up memory
            
            % --- Save the Clean Data ---
            % We will overwrite the original file with the sanitized version.
            % This is safe because this is a one-time pre-processing step.
            fprintf('Saving sanitized file (overwriting original)...\n');
            save(full_data_path, 'data_for_hmm', 'dt', '-v7.3');
            
            fprintf('--- Successfully sanitized: %s ---\n', file_info.name);
            success_count = success_count + 1;
        catch ME
            fprintf(2, '--- ERROR sanitizing file: %s ---\n', full_data_path);
            fprintf(2, 'Error Message: %s\n', ME.message);
            fprintf(2, 'This file may have a deeper issue. The script will stop.\n');
            rethrow(ME);
        end
    end
    
    fprintf('\n%s\n--- Sanitization Complete ---\n', repmat('*', 1, 80));
    fprintf('Successfully processed: %d of %d files\n', success_count, num_files);
    fprintf('%s\n', repmat('*', 1, 80));

end