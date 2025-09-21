% validate_TPT_fixed_window.m
% The definitive validation script to test the effect of information density.
% It analyzes fixed-size time windows of continuous photon data across all
% fluxes to prove that convergence depends on flux, not burst length.

function validate_TPT_fixed_window()
    clearvars; close all; clc;
    
    %% --- USER CONFIGURATION ---
    % ---------------------------------------------------------------------
    % 1. DEFINE THE FIXED TIME WINDOW for analysis
    TIME_WINDOW_US = 500; % 500 microseconds

    % 2. DEFINE THE TEST MODE
    %    'SPECIFICITY' - (Recommended) Starts from a COMPLETELY RANDOM point in the window.
    TEST_MODE = 'SPECIFICITY';
    % ---------------------------------------------------------------------
    
    fprintf('--- Starting TPT Fixed Time Window Validation (Mode: %s) ---\n', TEST_MODE);
    
    %% --- File Loading ---
    base_results_path = uigetdir(pwd, 'Select the main validation folder (e.g., Flux_Sweep_Validation_Memory_Efficient)');
    if isequal(base_results_path, 0), disp('User cancelled.'); return; end
    
    [gt_file, gt_path] = uigetfile('*.mat', 'Select the Ground Truth ..._for_MATLAB.mat file');
    if isequal(gt_file,0), disp('User cancelled.'); return; end
    load(fullfile(gt_path, gt_file), 'Inflection_Point');
    ground_truth_tpts_us = Inflection_Point;
    
    calib_file = fullfile(base_results_path, 'TPT_Calibration_Summary.mat');
    if ~exist(calib_file, 'file'), error('Calibration summary file not found: %s', calib_file); end
    load(calib_file, 'calibration_results');
    
    %% --- Main Loop over Flux Conditions ---
    flux_folders = dir(fullfile(base_results_path, 'flux_*_hz'));
    flux_folders = flux_folders([flux_folders.isdir]);
    
    for i_flux = 1:length(flux_folders)
        flux_folder_info = flux_folders(i_flux);
        flux_folder_path = fullfile(flux_folder_info.folder, flux_folder_info.name);
        
        fprintf('\n%s\n--- Processing Flux Condition: %s ---\n%s\n', repmat('=', 1, 80), flux_folder_info.name, repmat('=', 1, 80));
        
        % --- Load continuous data for this flux ---
        data_file = fullfile(flux_folder_path, 'data.mat');
        if ~exist(data_file, 'file'), warning('Continuous data file "data.mat" not found. Skipping flux.'); continue; end
        load(data_file, 'data'); % Loads 'data' cell array of 20 long streams
        
        % --- Stitch all 20 trajectories into one massive stream ---
        fprintf('Stitching trajectories...');
        stitched_stream = stitch_trajectories(data);
        fprintf(' Done. Total duration: %.2f seconds.\n', stitched_stream.times_sec(end));
        
        % --- Load HMM model for this flux to get emission probabilities ---
        results_folder = dir(fullfile(flux_folder_path, 'HMM_Ready_Data', 'HMM_Res_*'));
        results_folder = results_folder([results_folder.isdir]);
        if isempty(results_folder), warning('No HMM results folder found. Skipping flux.'); continue; end
        load(fullfile(results_folder(1).folder, results_folder(1).name, 'Best_HMM_Model_Local.mat'), 'best_model_overall');
        obsmat = best_model_overall.obsmat;
        [~, sorted_idx] = sort(obsmat(:,2));
        P_Red_L = obsmat(sorted_idx(1), 2); P_Red_H = obsmat(sorted_idx(2), 2);
        
        % --- Find the calibrated thresholds for this flux ---
        flux_str = regexp(flux_folder_info.name, 'flux_(\d+\.?\d*e\d+)_hz', 'tokens');
        photon_flux = str2double(flux_str{1}{1});
        calib_row = abs(calibration_results.PhotonFlux_Hz - photon_flux) < 1e-9;
        p_thresh = calibration_results.Optimal_P_Thresh(calib_row);
        lambda_thresh = calibration_results.Optimal_Lambda(calib_row);
        
        %% --- Run the test on the fixed time windows ---
        [TPTs_sb, TPTs_ib, TPTs_sl, TPTs_il] = deal([]); % Standard/Iterative Bayes, Standard/Iterative LLR
        
        time_window_sec = TIME_WINDOW_US * 1e-6;
        num_windows = floor(stitched_stream.times_sec(end) / time_window_sec);
        
        for i_win = 1:num_windows
            start_time = (i_win - 1) * time_window_sec;
            end_time = i_win * time_window_sec;
            
            % Extract photons within this time window
            idx_in_window = find(stitched_stream.times_sec >= start_time & stitched_stream.times_sec < end_time);
            if length(idx_in_window) < 20, continue; end % Skip windows with too few photons
            
            window_photons.times = stitched_stream.times_sec(idx_in_window);
            window_photons.colors = stitched_stream.colors(idx_in_window);
            window_photons.gt_times = stitched_stream.gt_times_ns(idx_in_window);
            
            % Check if a true transition happened in this window
            if ~any(diff(floor(window_photons.gt_times / 1000)) ~= 0), continue; end % Simple check for change in ground truth
            
            % Generate ONE random start point
            start_idx = randi([10, length(window_photons.colors) - 10]);
            type = 'LH'; if rand > 0.5, type = 'HL'; end 
            
            % Run all four methods
            [e_sb, en_sb] = get_bayes_points(window_photons.colors, start_idx, P_Red_H, P_Red_L, p_thresh, type);
            if ~isnan(e_sb) && ~isnan(en_sb) && en_sb > e_sb, TPTs_sb(end+1) = (window_photons.times(en_sb) - window_photons.times(e_sb)); end
            [e_ib, en_ib] = find_tpt_iterative_bayes_forced(window_photons.colors, start_idx, P_Red_H, P_Red_L, p_thresh, type);
            if ~isnan(e_ib) && ~isnan(en_ib) && en_ib > e_ib, TPTs_ib(end+1) = (window_photons.times(en_ib) - window_photons.times(e_ib)); end
            [e_sl, en_sl] = get_llr_points(window_photons.colors, start_idx, P_Red_H, P_Red_L, lambda_thresh, type);
            if ~isnan(e_sl) && ~isnan(en_sl) && en_sl > e_sl, TPTs_sl(end+1) = (window_photons.times(en_sl) - window_photons.times(e_sl)); end
            [e_il, en_il] = find_tpt_iterative_llr_forced(window_photons.colors, start_idx, P_Red_H, P_Red_L, lambda_thresh, type);
            if ~isnan(e_il) && ~isnan(en_il) && en_il > e_il, TPTs_il(end+1) = (window_photons.times(en_il) - window_photons.times(e_il)); end
        end
        
        %% --- Plot and Save for this flux ---
        fig = figure('Name', ['Fixed Window Test: ' flux_folder_info.name], 'Visible', 'off');
        hold on; bin_width = 2.0;
        histogram(ground_truth_tpts_us, 'BinWidth', bin_width, 'Normalization', 'pdf', 'DisplayName', sprintf('Ground Truth (N=%d, Mean=%.2f us)', length(ground_truth_tpts_us), mean(ground_truth_tpts_us)), 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.6);
        histogram(TPTs_sb*1e6, 'BinWidth', bin_width, 'Normalization', 'pdf', 'DisplayName', sprintf('Std. Bayes (N=%d, Mean=%.2f us)', length(TPTs_sb), mean(TPTs_sb*1e6)), 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 2.0, 'LineStyle', ':');
        histogram(TPTs_sl*1e6, 'BinWidth', bin_width, 'Normalization', 'pdf', 'DisplayName', sprintf('Std. LLR (N=%d, Mean=%.2f us)', length(TPTs_sl), mean(TPTs_sl*1e6)), 'DisplayStyle', 'stairs', 'EdgeColor', [0.6 0 0], 'LineWidth', 2.0, 'LineStyle', ':');
        histogram(TPTs_ib*1e6, 'BinWidth', bin_width, 'Normalization', 'pdf', 'DisplayName', sprintf('Iter. Bayes (N=%d, Mean=%.2f us)', length(TPTs_ib), mean(TPTs_ib*1e6)), 'DisplayStyle', 'stairs', 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 2.5);
        histogram(TPTs_il*1e6, 'BinWidth', bin_width, 'Normalization', 'pdf', 'DisplayName', sprintf('Iter. LLR (N=%d, Mean=%.2f us)', length(TPTs_il), mean(TPTs_il*1e6)), 'DisplayStyle', 'stairs', 'EdgeColor', [0.4660 0.6740 0.1880], 'LineWidth', 2.5);
        hold off; grid on; box on;
        title({sprintf('Fixed Window (%d us) Test for Flux = %.1e Hz', TIME_WINDOW_US, photon_flux), 'Start Points: SPECIFICITY (Forced Result)'}, 'FontSize', 14);
        xlabel('Transition Path Time (\mus)'); ylabel('Probability Density');
        legend('show', 'Location', 'NorthEast'); xlim([0, 60]);
        
        output_dir = fullfile(flux_folder_path, 'HMM_Ready_Data', results_folder(1).name);
        saveas(fig, fullfile(output_dir, 'TPT_Fixed_Window_Test_Plot.png'));
        savefig(fig, fullfile(output_dir, 'TPT_Fixed_Window_Test_Plot.fig'));
        close(fig);
    end
    fprintf('\n--- Fixed Window Validation Batch Run Complete ---\n');
end

%% --- SUB-FUNCTIONS ---
function stitched = stitch_trajectories(data_cell)
    total_photons = sum(cellfun(@(c) size(c,1), data_cell));
    stitched.times_sec = zeros(total_photons, 1);
    stitched.colors = zeros(total_photons, 1);
    stitched.gt_times_ns = zeros(total_photons, 1);
    
    current_offset_sec = 0;
    current_idx = 1;
    for i = 1:length(data_cell)
        traj = data_cell{i};
        if isempty(traj), continue; end
        num_photons = size(traj, 1);
        
        end_idx = current_idx + num_photons - 1;
        stitched.times_sec(current_idx:end_idx) = traj(:,1) + current_offset_sec;
        stitched.colors(current_idx:end_idx) = traj(:,2);
        stitched.gt_times_ns(current_idx:end_idx) = traj(:,3); % Assuming 3rd column is ground truth time
        
        current_idx = end_idx + 1;
        current_offset_sec = stitched.times_sec(end_idx);
    end
end
% ... (Include all the find_tpt... and get...points helper functions from the previous script here) ...
function [e,en]=find_tpt_iterative_bayes_forced(c,s,ph,pl,p,t),max_iter=500;tol=2;cci=s;[e,en,ep,enp]=deal(NaN);for iter=1:max_iter,[ep,enp]=get_bayes_points(c,cci,ph,pl,p,t);if isnan(ep)||isnan(enp),e=NaN;en=NaN;return;end;nci=round((ep+enp)/2);if abs(nci-cci)<=tol,e=ep;en=enp;return;end;cci=nci;end;e=ep;en=enp;end
function [e,en]=find_tpt_iterative_llr_forced(c,s,ph,pl,l,t),max_iter=500;tol=2;cci=s;[e,en,ep,enp]=deal(NaN);for iter=1:max_iter,[ep,enp]=get_llr_points(c,cci,ph,pl,l,t);if isnan(ep)||isnan(enp),e=NaN;en=NaN;return;end;nci=round((ep+enp)/2);if abs(nci-cci)<=tol,e=ep;en=enp;return;end;cci=nci;end;e=ep;en=enp;end
function [idx_e,idx_en]=get_bayes_points(c,s,ph,pl,p,t),idx_e=findExitPoint_Bayesian(c,s,ph,pl,p,t);idx_en=findEntryPoint_Bayesian(c,s,ph,pl,p,t);end
function [idx_e,idx_en]=get_llr_points(c,s,ph,pl,l,t),idx_e=findExitPoint_Standard(c,s,ph,pl,l,t);idx_en=findEntryPoint_Standard(c,s,ph,pl,l,t);end
function idx=findEntryPoint_Standard(c,s,ph,pl,l,t),llr=0;idx=NaN;if strcmp(t,'LH'),lr=log(ph/pl);lg=log((1-ph)/(1-pl));else,lr=log(pl/ph);lg=log((1-pl)/(1-ph));end;for i=s:length(c),if c(i)==2,llr=llr+lr;else,llr=llr+lg;end;if llr>l,idx=i;return;end;end;end
function idx=findExitPoint_Standard(c,s,ph,pl,l,t),llr=0;idx=NaN;if strcmp(t,'LH'),lr=log(pl/ph);lg=log((1-pl)/(1-ph));else,lr=log(ph/pl);lg=log((1-ph)/(1-pl));end;for i=(s-1):-1:1,if c(i)==2,llr=llr+lr;else,llr=llr+lg;end;if llr>l,idx=i+1;return;end;end;end
function idx=findEntryPoint_Bayesian(c,s,ph,pl,p,t),llr=0;idx=NaN;if strcmp(t,'LH'),lr=log(ph/pl);lg=log((1-ph)/(1-pl));else,lr=log(pl/ph);lg=log((1-pl)/(1-ph));end;for i=s:length(c),if c(i)==2,llr=llr+lr;else,llr=llr+lg;end;if 1/(1+exp(-llr))>p,idx=i;return;end;end;end
function idx=findExitPoint_Bayesian(c,s,ph,pl,p,t),llr=0;idx=NaN;if strcmp(t,'LH'),lr=log(pl/ph);lg=log((1-pl)/(1-ph));else,lr=log(ph/pl);lg=log((1-ph)/(1-pl));end;for i=(s-1):-1:1,if c(i)==2,llr=llr+lr;else,llr=llr+lg;end;if 1/(1+exp(-llr))>p,idx=i+1;return;end;end;end