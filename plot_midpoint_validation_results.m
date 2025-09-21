% plot_midpoint_validation_results.m
%
% This script loads the data file created by the Python validation script
% ('midpoint_validation_results.mat') and generates three publication-quality,
% fully editable figures to visualize the results.
%
function plot_midpoint_validation_results()
    clearvars; close all; clc;

    %% 1. Load the Validation Data
    fprintf('Please select the ''midpoint_validation_results.mat'' file...\n');
    [file, path] = uigetfile('*.mat', 'Select Midpoint Validation Results');
    if isequal(file, 0), disp('User cancelled.'); return; end
    
    load(fullfile(path, file));

    %% 2. Calculate Statistics
    mean_viterbi_HL = mean(viterbi_errors_HL_us);
    std_viterbi_HL  = std(viterbi_errors_HL_us);
    mean_viterbi_LH = mean(viterbi_errors_LH_us);
    std_viterbi_LH  = std(viterbi_errors_LH_us);
    
    mean_bayes_HL = mean(bayesian_errors_HL_us);
    std_bayes_HL  = std(bayesian_errors_HL_us);
    mean_bayes_LH = mean(bayesian_errors_LH_us);
    std_bayes_LH  = std(bayesian_errors_LH_us);

    %% 3. Figure 1: Viterbi Midpoint Error Histogram
    figure('Name', 'Viterbi Midpoint Error', 'Color', 'w', 'Position', [100 400 700 500]);
    hold on;
    histogram(viterbi_errors_HL_us, 'BinWidth', 2, 'Normalization', 'pdf', 'FaceColor', [0.3 0.6 1.0], 'FaceAlpha', 0.7, 'DisplayName', sprintf('High -> Low (%.2f \\pm %.2f \\mus)', mean_viterbi_HL, std_viterbi_HL));
    histogram(viterbi_errors_LH_us, 'BinWidth', 2, 'Normalization', 'pdf', 'FaceColor', [1.0 0.5 0.3], 'FaceAlpha', 0.7, 'DisplayName', sprintf('Low -> High (%.2f \\pm %.2f \\mus)', mean_viterbi_LH, std_viterbi_LH));
    xline(0, '--k', 'LineWidth', 2.5, 'Label', 'Perfect Agreement');
    grid on; box on;
    xlabel('Viterbi Midpoint Error (\mus)');
    ylabel('Probability Density');
    title('Distribution of Viterbi Midpoint Errors');
    legend('show', 'Location', 'NorthWest');
    xlim([-50 50]); % Set fixed axis limits
    
    %% 4. Figure 2: Iterative Bayesian Midpoint Error Histogram
    figure('Name', 'Iterative Bayesian Midpoint Error', 'Color', 'w', 'Position', [850 400 700 500]);
    hold on;
    histogram(bayesian_errors_HL_us, 'BinWidth', 2, 'Normalization', 'pdf', 'FaceColor', [0.3 0.6 1.0], 'FaceAlpha', 0.7, 'DisplayName', sprintf('High -> Low (%.2f \\pm %.2f \\mus)', mean_bayes_HL, std_bayes_HL));
    histogram(bayesian_errors_LH_us, 'BinWidth', 2, 'Normalization', 'pdf', 'FaceColor', [1.0 0.5 0.3], 'FaceAlpha', 0.7, 'DisplayName', sprintf('Low -> High (%.2f \\pm %.2f \\mus)', mean_bayes_LH, std_bayes_LH));
    xline(0, '--k', 'LineWidth', 2.5, 'Label', 'Perfect Agreement');
    grid on; box on;
    xlabel('Iterative Bayesian Midpoint Error (\mus)');
    ylabel('Probability Density');
    title('Distribution of Corrected Bayesian Midpoint Errors');
    legend('show', 'Location', 'NorthWest');
    xlim([-50 50]); % Set fixed axis limits

    %% 5. Figure 3: The Combined Scatter Plot
    figure('Name', 'Correction of Viterbi Bias', 'Color', 'w', 'Position', [400 50 800 800]);
    hold on;
    scatter(viterbi_errors_HL_us, bayesian_errors_HL_us, 15, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [0.3 0.6 1.0], 'DisplayName', sprintf('High -> Low (N=%d)', numel(viterbi_errors_HL_us)));
    scatter(viterbi_errors_LH_us, bayesian_errors_LH_us, 15, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [1.0 0.5 0.3], 'DisplayName', sprintf('Low -> High (N=%d)', numel(viterbi_errors_LH_us)));
    
    lims = [min([xlim ylim]), max([xlim ylim])];
    plot(lims, lims, 'r--', 'LineWidth', 2, 'DisplayName', 'No Improvement (y=x)');
    
    axhline = refline(0, 0); axhline.Color = [0.2 0.2 0.2 0.8]; axhline.LineStyle = '--'; axhline.DisplayName = 'Perfect Bayesian Midpoint';
    axvline = line([0 0], lims, 'Color', [0.2 0.2 0.2 0.8], 'LineStyle', ':', 'DisplayName', 'Perfect Viterbi Midpoint');
    
    xlim(lims); ylim(lims);
    axis square; grid on; box on;
    xlabel('Viterbi Midpoint Error (\mus)');
    ylabel('Iterative Bayesian Midpoint Error (\mus)');
    title('Correction of Viterbi Bias by Iterative Bayesian Method');
    legend('show', 'Location', 'NorthWest');
    set(gca, 'FontSize', 11);
    
end