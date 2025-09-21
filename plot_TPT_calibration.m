% plot_TPT_calibration_final.m
% This script loads the results from the TPT calibration sweep and generates
% two separate, publication-quality plots:
%   1. Calibrated Bayesian Threshold (P_thresh) vs. Photon Flux
%   2. Calibrated LLR Threshold (Lambda) vs. Photon Flux

function plot_TPT_calibration()
    clear; close all; clc;
    
    %% 1. Select the Base Folder and Load the Results File
    fprintf('Please select the main validation folder containing the results file.\n');
    
    % Use uigetdir to ask the user to select the folder
    base_results_path = uigetdir(pwd, 'Select the main validation folder (e.g., Flux_Sweep_Validation_Memory_Efficient)');
    if isequal(base_results_path, 0)
        disp('User cancelled selection. Exiting.');
        return;
    end
    
    % Construct the full path to the results file
    results_filename = 'TPT_Calibration_Results.mat';
    full_results_path = fullfile(base_results_path, results_filename);
    
    if ~exist(full_results_path, 'file')
        error('Results file not found: %s\nPlease ensure the file exists in the selected directory.', full_results_path);
    end
    
    fprintf('Loading calibration data from: %s\n', full_results_path);
    load(full_results_path, 'calibration_results');
    
    % Extract data from the table for easy plotting
    photon_flux = calibration_results.PhotonFlux_Hz;
    optimal_lambda = calibration_results.Optimal_Lambda;
    optimal_p_thresh = calibration_results.Optimal_P_Thresh;
    
    %% 2. Create Plot 1: Iterative Bayesian Threshold
    fprintf('Generating plot for Bayesian Threshold (P_thresh)...\n');
    
    fig1 = figure('Name', 'Bayesian Threshold Calibration', 'Color', 'w', 'Position', [100, 200, 800, 550]);
    ax1 = axes(fig1);
    
    plot(ax1, photon_flux, optimal_p_thresh, 'o-', ...
        'Color', [0, 0.4470, 0.7410], ... % Blue
        'LineWidth', 2, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', [0, 0.4470, 0.7410]);
    hold(ax1, 'on');
    
    % Add a horizontal line to highlight the plateau
    plateau_p_thresh = optimal_p_thresh(end);
    yline(ax1, plateau_p_thresh, '--', sprintf('Plateau = %.3f', plateau_p_thresh), ...
        'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
    
    set(ax1, 'XScale', 'log');
    grid(ax1, 'on');
    box(ax1, 'on');
    xlabel(ax1, 'Photon Flux (Hz)');
    ylabel(ax1, 'Calibrated P_{thresh}');
    title(ax1, 'TPT Calibration: Optimal Bayesian Threshold vs. Photon Flux', 'FontSize', 14);
    ylim(ax1, [0.7, 1.0]); % Focus on the relevant range
    set(ax1, 'FontSize', 12);

    % Save the figure
    output_fig_name1 = 'Bayesian_Threshold_vs_Flux.png';
    saveas(fig1, fullfile(base_results_path, output_fig_name1));
    fprintf('Plot 1 successfully saved as: %s\n', output_fig_name1);
    
    %% 3. Create Plot 2: Iterative LLR Threshold
    fprintf('Generating plot for LLR Threshold (Lambda)...\n');
    
    fig2 = figure('Name', 'LLR Threshold Calibration', 'Color', 'w', 'Position', [950, 200, 800, 550]);
    ax2 = axes(fig2);
    
    plot(ax2, photon_flux, optimal_lambda, 's--', ...
        'Color', [0.8500, 0.3250, 0.0980], ... % Orange
        'LineWidth', 2, ...
        'MarkerSize', 8);
    hold(ax2, 'on');
        
    % Add a horizontal line to highlight the plateau
    plateau_lambda = optimal_lambda(end);
    yline(ax2, plateau_lambda, '--', sprintf('Plateau = %.2f', plateau_lambda), ...
        'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
        
    set(ax2, 'XScale', 'log');
    grid(ax2, 'on');
    box(ax2, 'on');
    xlabel(ax2, 'Photon Flux (Hz)');
    ylabel(ax2, 'Calibrated LLR Threshold (\Lambda)');
    title(ax2, 'TPT Calibration: Optimal LLR Threshold vs. Photon Flux', 'FontSize', 14);
    set(ax2, 'FontSize', 12);
    
    % Save the figure
    output_fig_name2 = 'LLR_Threshold_vs_Flux.png';
    saveas(fig2, fullfile(base_results_path, output_fig_name2));
    fprintf('Plot 2 successfully saved as: %s\n', output_fig_name2);
    
end