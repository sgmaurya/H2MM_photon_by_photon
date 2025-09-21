% plot_TPT_Comparison_Final_v2.m
%
% This is the definitive plotting script. It includes the fix for the 'LineStyle'
% error by setting the property on the histogram object after creation.
%
function plot_TPT_Comparison()
    clearvars; close all; clc;

    %% 1. Load the Foundational Ground Truth Data
    fprintf('--- Step 1: Select the Ground Truth Data ---\n');
    [gt_file, gt_path] = uigetfile('*.mat', 'Select the Ground Truth ..._for_MATLAB.mat file');
    if isequal(gt_file,0), disp('User cancelled.'); return; end
    full_gt_path = fullfile(gt_path, gt_file);
    fprintf('Loading Ground Truth data from: %s\n', full_gt_path);
    
    try
        load(full_gt_path);
        if exist('Inflection_Point', 'var')
            gt_tpts_us = Inflection_Point;
            clear Inflection_Point; 
        else
            errordlg('The required variable "Inflection_Point" was not found in the Ground Truth .mat file.', 'Variable Not Found');
            return;
        end
    catch ME
        errordlg(sprintf('Failed to load Ground Truth data: %s', ME.message), 'File Error');
        return;
    end
    
    %% 2. Ask How Many Measured Datasets to Compare
    prompt = {'How many MEASURED TPT result sets to compare? (1-4)'};
    dlgtitle = 'Number of Datasets for Comparison';
    definput = {'4'}; % Default to comparing all four methods
    answer = inputdlg(prompt, dlgtitle, [1 60], definput);
    if isempty(answer), disp('User cancelled.'); return; end
    
    try
        num_datasets = str2double(answer{1});
        if isnan(num_datasets) || ~ismember(num_datasets, 1:4), error('Please enter an integer between 1 and 4.'); end
    catch ME
        errordlg(sprintf('Invalid input: %s', ME.message), 'Input Error');
        return;
    end

    %% 3. Initialize Storage and Styles (with LineStyle fix)
    all_measured_data = cell(1, num_datasets);
    all_legends = cell(1, num_datasets);
    
    % --- FIX: Ensure all 'Stairs' styles have a 'LineStyle' field defined ---
    plot_styles = { ...
        struct('Type', 'Stairs', 'EdgeColor', [0 0.45 0.74], 'LineWidth', 2.5, 'LineStyle', '-'), ... % Blue Solid Outline
        struct('Type', 'Stairs', 'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 3.0, 'LineStyle', '--'), ... % Black Dashed Outline
        struct('Type', 'Stairs', 'EdgeColor', [0.47 0.67 0.19], 'LineWidth', 2.5, 'LineStyle', '-'), ... % Green Solid Outline
        struct('Type', 'Stairs', 'EdgeColor', [0.49 0.18 0.56], 'LineWidth', 2.5, 'LineStyle', ':') ... % Purple Dotted Outline
    };

    %% 4. Loop to Load Each Measured Dataset
    for i = 1:num_datasets
        [file_name, path_name] = uigetfile('*.mat', sprintf('Select Measured Dataset %d of %d', i, num_datasets));
        if isequal(file_name,0), disp('User cancelled.'); return; end
        full_path = fullfile(path_name, file_name);
        fprintf('Loading Measured Dataset %d from: %s\n', i, full_path);
        
        loaded_data = load(full_path);
        
        if isfield(loaded_data, 'TPTs_sec_HL')
            all_measured_data{i} = [loaded_data.TPTs_sec_HL, loaded_data.TPTs_sec_LH] * 1e6;
        else
            warning('Could not find TPT data in %s. Skipping this file.', file_name);
            continue;
        end
        
        % Automatically Generate Legend
        mean_val = mean(all_measured_data{i});
        params = loaded_data.params;
        switch params.method
            case 'Standard LLR', legend_str = sprintf('Std. LLR (\\lambda=%.2f, Mean=%.2f us)', params.lambda_thresh, mean_val);
            case 'Iterative LLR (Convergent)', legend_str = sprintf('Iter. LLR (\\lambda=%.2f, Mean=%.2f us)', params.lambda_thresh, mean_val);
            case 'Standard Bayesian', legend_str = sprintf('Std. Bayesian (P_{thresh}=%.2f, Mean=%.2f us)', params.prob_thresh, mean_val);
            case 'Iterative Bayesian (Convergent)', legend_str = sprintf('Iter. Bayesian (P_{thresh}=%.2f, Mean=%.2f us)', params.prob_thresh, mean_val);
            otherwise, legend_str = sprintf('Unknown Method (Mean=%.2f us)', mean_val);
        end
        all_legends{i} = legend_str;
    end

    %% 5. Create the Final Overlay Plot (with LineStyle fix)
    fprintf('\n--- Generating final comparison plot ---\n');
    figure('Position', [100, 100, 1000, 700]);
    hold on;

    bin_width_ans = inputdlg('Enter bin width in microseconds:', 'Bin Width', [1 40], {'2.0'});
    if isempty(bin_width_ans), return; end
    bin_width = str2double(bin_width_ans{1});

    % Plot 1: The Ground Truth (as a solid, filled histogram)
    histogram(gt_tpts_us, 'BinWidth', bin_width, 'Normalization', 'pdf', ...
              'DisplayName', sprintf('Ground Truth (Mean = %.2f us)', mean(gt_tpts_us)), ...
              'FaceColor', [0.85 0.33 0.1], 'FaceAlpha', 0.6, 'EdgeColor', 'none');

    % Plot 2...N: The Measured Datasets
    for i = 1:num_datasets
        if isempty(all_measured_data{i}), continue; end
        current_style = plot_styles{i};
        
        % --- FIX: Capture the histogram handle (h) ---
        h = histogram(all_measured_data{i}, 'BinWidth', bin_width, 'Normalization', 'pdf', ...
                      'DisplayName', all_legends{i}, ...
                      'DisplayStyle', 'stairs', 'EdgeColor', current_style.EdgeColor, ...
                      'LineWidth', current_style.LineWidth);
                      
        % --- FIX: Set the LineStyle property on the handle ---
        h.LineStyle = current_style.LineStyle;
    end

    hold off;
    grid on; box on;
    title('Final Validation: Ground Truth vs. Advanced TPT Measurements', 'FontSize', 16);
    xlabel('Transition Path Time (\mus)', 'FontSize', 14);
    ylabel('Probability Density', 'FontSize', 14);
    legend('show', 'Location', 'NorthEast', 'FontSize', 11, 'Interpreter', 'tex');
    set(gca, 'FontSize', 12);

    fprintf('\nPlot generation complete. Use "File > Save As" to save the figure.\n');
end