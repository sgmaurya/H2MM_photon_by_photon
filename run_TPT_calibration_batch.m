% run_TPT_calibration_batch.m
% Final controller script to fully automate the TPT calibration sweep.
% VERSION 2: Calls the new worker 'find_and_save_optimal_tpts.m' to not only
%            find the best thresholds but also save the corresponding TPT
%            distributions for later visualization.

function run_TPT_calibration_batch()
    clear; close all; clc;
    
    fprintf('--- Starting FULLY AUTOMATED TPT Calibration and Data Saving ---\n');
    
    %% --- USER CONFIGURATION ---
    base_results_path = 'Y:\Flux_Sweep_Validation_Memory_Efficient';
    GROUND_TRUTH_MEAN_TPT_US = 11.06;
    output_filename = 'TPT_Calibration_Summary.mat'; % Renamed for clarity
    % ---------------------------------------------------------------------

    %% --- Script Core Logic ---
    worker_script_name = 'find_and_save_optimal_tpts.m'; % <<< Use the new worker
    if ~exist(worker_script_name, 'file'), error('Worker script "%s" not found.', worker_script_name); end
    
    results_folder_pattern = fullfile(base_results_path, 'flux_*_hz', 'HMM_Ready_Data', 'HMM_Res_*');
    results_folders = dir(results_folder_pattern);
    results_folders = results_folders([results_folders.isdir]);
    
    if isempty(results_folders), error('No HMM results folders found.'); end
    
    num_folders = length(results_folders);
    fprintf('Found %d processed HMM folders to calibrate.\n', num_folders);
    
    calibration_results = table('Size', [num_folders, 3], ...
                                'VariableTypes', {'double', 'double', 'double'}, ...
                                'VariableNames', {'PhotonFlux_Hz', 'Optimal_Lambda', 'Optimal_P_Thresh'});
    
    for i = 1:num_folders
        folder_info = results_folders(i);
        full_results_path = fullfile(folder_info.folder, folder_info.name);
        
        fprintf('\n%s\n--- Calibrating Folder %d of %d ---\n', repmat('=', 1, 80), i, num_folders);
        
        try
            flux_str = regexp(folder_info.folder, 'flux_(\d+\.?\d*e\d+)_hz', 'tokens');
            photon_flux = str2double(flux_str{1}{1});
            
            % Call the new, more powerful worker function
            [opt_lambda, opt_p_thresh] = find_and_save_optimal_tpts(full_results_path, GROUND_TRUTH_MEAN_TPT_US);
            
            calibration_results.PhotonFlux_Hz(i) = photon_flux;
            calibration_results.Optimal_Lambda(i) = opt_lambda;
            calibration_results.Optimal_P_Thresh(i) = opt_p_thresh;

            fprintf('--- Successfully calibrated and saved TPT data for this flux ---\n');

        catch ME, rethrow(ME); end
    end
    
    calibration_results = sortrows(calibration_results, 'PhotonFlux_Hz');
    output_save_path = fullfile(base_results_path, output_filename);
    save(output_save_path, 'calibration_results');
    
    fprintf('\n%s\n--- TPT Calibration Sweep Complete ---\n', repmat('*', 1, 80));
    fprintf('Final summary table saved to: %s\n', output_save_path);
    disp(calibration_results);
end