% convert_mat_to_txt_batch.m
% A utility script to solve the Python data loading issue.
% VERSION 3: FINAL & DEFINITIVE. Replaces 'writematrix' with the more
%            fundamental and robust 'dlmwrite' function to handle the
%            unusual 3D cell array data structure and bypass the
%            "Invalid parameter" error on all MATLAB versions.

function convert_mat_to_txt_batch()
    clear; close all; clc;
    
    fprintf('--- Starting Conversion of .mat files to .txt ---\_n');
    
    %% --- USER CONFIGURATION ---
    base_data_path = 'Y:\Photon_Data_Grid';
    % ---------------------------------------------------------------------

    %% --- Script Core Logic ---
    
    mat_file_pattern = fullfile(base_data_path, 'TPT_Explore_*', 'flux_*_hz', 'data.mat');
    mat_files = dir(mat_file_pattern);
    
    if isempty(mat_files)
        error('No "data.mat" files found using the pattern:\n%s', mat_file_pattern);
    end
    
    num_files = length(mat_files);
    fprintf('Found %d .mat files to convert to .txt.\n', num_files);
    
    h_wait = waitbar(0, 'Converting files...');
    successful_conversions = 0;
    
    for i = 1:num_files
        file_info = mat_files(i);
        full_mat_path = fullfile(file_info.folder, file_info.name);
        
        waitbar(i/num_files, h_wait, sprintf('Converting file %d of %d', i, num_files));
        
        try
            loaded_data = load(full_mat_path);
            
            % --- Robust Data Extraction ---
            % This correctly handles the strange 3D cell array format by
            % converting it to a standard numeric matrix.
            if iscell(loaded_data.data)
                photon_data_matrix = cell2mat(loaded_data.data);
            else
                photon_data_matrix = loaded_data.data;
            end
            
            % If the matrix is 3D (1xMx3), reshape it to 2D (Mx3)
            if ndims(photon_data_matrix) == 3
                photon_data_matrix = reshape(photon_data_matrix, [], 3);
            end
            
            output_txt_path = fullfile(file_info.folder, 'photon_stream.txt');
            
            % --- THIS IS THE FIX ---
            % Use the highly robust 'dlmwrite' function.
            % It is the workhorse for writing numeric matrices to text files.
            dlmwrite(output_txt_path, photon_data_matrix, 'delimiter', '\t', 'precision', 9);
            % --- END FIX ---
            
            successful_conversions = successful_conversions + 1;
        catch ME
            fprintf(2, '\n--- ERROR converting file: %s ---\n', full_mat_path);
            fprintf(2, 'Error Message: %s\n', ME.message);
        end
    end
    
    if ishandle(h_wait), close(h_wait); end
    
    fprintf('\n%s\n--- Conversion Complete! ---\n', repmat('*', 1, 80));
    fprintf('Successfully converted %d of %d files.\n', successful_conversions, num_files);
    fprintf('You can now run the Python chunking script.\n');
    fprintf('%s\n', repmat('*', 1, 80));

end