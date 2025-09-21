% run_robustness_test_batch.m
% A controller script to fully automate the TPT robustness and specificity validations.
% FINAL VERSION: Defaults to 'SPECIFICITY' mode as requested.

function run_robustness_test_batch()
    clear; close all; clc;
    
    %% --- USER CONFIGURATION ---
    % ---------------------------------------------------------------------
    % 1. CHOOSE THE TEST MODE:
    %    'SPECIFICITY' - (Recommended) Starts from a COMPLETELY RANDOM point.
    %    'ROBUSTNESS'  - Starts from a point NEAR the true Viterbi changepoint.
    test_params.mode = 'SPECIFICITY'; 
    
    % 2. DEFINE BASE PATH for the Flux Sweep Validation data.
    base_results_path = 'Y:\Flux_Sweep_Validation_Memory_Efficient';

    % 3. DEFINE RANDOMIZATION WINDOW (only used for ROBUSTNESS mode)
    test_params.window = 50;
    % ---------------------------------------------------------------------

    %% --- Script Core Logic ---
    fprintf('--- Starting TPT Validation Sweep (Mode: %s) ---\n', test_params.mode);
    worker_script_name = 'execute_robustness_test.m';
    if ~exist(worker_script_name, 'file'), error('Worker script "%s" not found.', worker_script_name); end
    
    [gt_file, gt_path] = uigetfile('*.mat', 'Select the Ground Truth ..._for_MATLAB.mat file');
    if isequal(gt_file,0), disp('User cancelled.'); return; end
    load(fullfile(gt_path, gt_file), 'Inflection_Point');
    ground_truth_tpts_us = Inflection_Point;
    
    calib_file = fullfile(base_results_path, 'TPT_Calibration_Summary.mat');
    if ~exist(calib_file, 'file'), error('Calibration summary file not found: %s', calib_file); end
    load(calib_file, 'calibration_results');
    
    results_folder_pattern = fullfile(base_results_path, 'flux_*_hz', 'HMM_Ready_Data', 'HM*_Res_*');
    results_folders = dir(results_folder_pattern);
    results_folders = results_folders([results_folders.isdir]);
    
    if isempty(results_folders), error('No HMM results folders found.'); end
    
    fprintf('Found %d folders to process.\n', length(results_folders));
    
    for i = 1:length(results_folders)
        folder_info = results_folders(i);
        full_results_path = fullfile(folder_info.folder, folder_info.name);
        
        fprintf('\n%s\n--- Running Test on Folder %d of %d ---\n', repmat('=', 1, 80), i, length(results_folders));
        
        try
            flux_str = regexp(folder_info.folder, 'flux_(\d+\.?\d*e\d+)_hz', 'tokens');
            photon_flux = str2double(flux_str{1}{1});
            
            calib_row = abs(calibration_results.PhotonFlux_Hz - photon_flux) < 1e-9;
            if ~any(calib_row), warning('Could not find threshold for flux %.2e. Skipping.', photon_flux); continue; end
            
            test_params.p_thresh = calibration_results.Optimal_P_Thresh(calib_row);
            test_params.lambda = calibration_results.Optimal_Lambda(calib_row);
            
            fprintf('INFO: Using calibrated thresholds for flux %.2e Hz (P=%.3f, Lambda=%.3f)\n', ...
                    photon_flux, test_params.p_thresh, test_params.lambda);

            % Call the corrected worker script
            execute_robustness_test(full_results_path, ground_truth_tpts_us, test_params);
            
        catch ME, rethrow(ME); end
    end
    
    fprintf('\n%s\n--- Validation Batch Run Complete ---\n', repmat('*', 1, 80));
end