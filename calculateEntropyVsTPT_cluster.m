% calculateEntropyVsTPT_cluster.m
%
% This is a new analysis step that runs AFTER calculateTPT_cluster.m.
% It loads the TPT distributions and calculates the coarse-grained apparent
% entropy production (deltaS) as a function of the transition path time (t),
% based on the theory presented in Degünther et al., PNAS (2024).
%
% Formula used: deltaS(t) = k_B * ln [ P_forward(t) / P_reverse(t) ]
%   - 'forward' is assumed to be the Low -> High transition.
%   - 'reverse' is assumed to be the High -> Low transition.
%   - k_B is set to 1, so the result is in units of k_B.
%

function calculateEntropyVsTPT_cluster()
    clearvars -except break_debug; close all; clc;
    fprintf('--- Starting Entropy vs. TPT Calculation from Cluster HMM Results ---\n');

    original_pwd = pwd;
    try
        %% 1. Select Processed Cluster Run Directory and Load TPT Data
        fprintf('\n--- 1. Selecting Processed Cluster HMM Run Directory ---\n');
        selected_run_path = uigetdir(pwd, 'Select the Cluster HMM Run directory with TPT results');
        if isequal(selected_run_path, 0), disp('User cancelled. Exiting.'); return; end
        fprintf('Processing results from: %s\n', selected_run_path);
        cd(selected_run_path);

        % --- Load the TPT analysis results ---
        if ~exist('TPT_Analysis_Results_Cluster.mat', 'file')
            error('CALC_ENTROPY: TPT_Analysis_Results_Cluster.mat not found. Please run calculateTPT_cluster.m first.');
        end
        fprintf('Loading TPT data...\n');
        load('TPT_Analysis_Results_Cluster.mat', 'all_TPTs_sec_HL', 'all_TPTs_sec_LH', 'lambda_thresh');

        %% 2. Calculate Probability Densities from TPT Distributions
        fprintf('\n--- 2. Calculating Probability Densities for Forward and Reverse Transitions ---\n');
        
        % Convert to microseconds for easier handling
        tpts_us_forward_LH = all_TPTs_sec_LH * 1e6; % Low -> High is the forward process
        tpts_us_reverse_HL = all_TPTs_sec_HL * 1e6; % High -> Low is the reverse process

        % Let user define binning, must be consistent for both histograms
        bin_width_us_ans = inputdlg('Enter bin width for TPT histograms (us):', 'Bin Width', [1 50], {'3'});
        if isempty(bin_width_us_ans), disp('Cancelled. Exiting.'); return; end
        bin_width_us = str2double(bin_width_us_ans{1});

        % Define common bin edges to ensure PDFs are comparable
        max_tpt = max([tpts_us_forward_LH, tpts_us_reverse_HL]);
        bin_edges_us = 0 : bin_width_us : (ceil(max_tpt / bin_width_us) * bin_width_us);
        bin_centers_us = bin_edges_us(1:end-1) + bin_width_us/2;

        % Get histogram counts
        counts_forward_LH = histcounts(tpts_us_forward_LH, bin_edges_us);
        counts_reverse_HL = histcounts(tpts_us_reverse_HL, bin_edges_us);
        
        % Normalize counts to get Probability Density Functions (PDFs)
        % PDF = Counts / (Total Number of Events * Bin Width)
        pdf_forward_LH = counts_forward_LH / (sum(counts_forward_LH) * bin_width_us);
        pdf_reverse_HL = counts_reverse_HL / (sum(counts_reverse_HL) * bin_width_us);

        %% 3. Calculate Apparent Entropy Production (deltaS)
        fprintf('\n--- 3. Calculating Apparent Entropy Production (deltaS) vs. TPT ---\n');
        
        % To avoid division by zero or log(0), only calculate deltaS for bins
        % where both forward and reverse transitions were observed.
        valid_bins_mask = (pdf_forward_LH > 0) & (pdf_reverse_HL > 0);
        
        tpt_for_plot_us = bin_centers_us(valid_bins_mask);
        pdf_forward_for_calc = pdf_forward_LH(valid_bins_mask);
        pdf_reverse_for_calc = pdf_reverse_HL(valid_bins_mask);

        % Calculate deltaS, assuming k_B = 1
        deltaS_values = log(pdf_forward_for_calc ./ pdf_reverse_for_calc);

        fprintf('Calculation complete. Found %d time bins with valid data for both transitions.\n', length(tpt_for_plot_us));

        %% 4. Plotting and Saving Results
        fprintf('\n--- 4. Plotting and Saving Entropy Results ---\n');
        
        fig_deltaS = figure('Name', 'Entropy vs. TPT (Cluster)', 'Color', 'w', 'Position', [400, 300, 900, 600]);
        ax = axes(fig_deltaS);
        
        plot(ax, tpt_for_plot_us, deltaS_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        
        xlabel(ax, 'Transition Path Time (\mus)');
        ylabel(ax, 'Apparent Entropy Production \DeltaS (k_B units)');
        title_str1 = sprintf('Apparent Entropy Production vs. TPT (Lambda=%.2f)', lambda_thresh);
        title_str2 = 'L->H transitions are treated as Forward, H->L as Reverse';
        title(ax, {title_str1, title_str2});
        grid(ax, 'on');
        box(ax, 'on');
        
        % Add a horizontal line at deltaS = 0 for reference (equilibrium)
        hold(ax, 'on');
        yline(ax, 0, '--k', 'LineWidth', 1.5, 'DisplayName', 'Equilibrium (\DeltaS=0)');
        hold(ax, 'off');
        
        saveas(fig_deltaS, 'Entropy_vs_TPT_LLR_Cluster.png');
        
        save('Entropy_vs_TPT_Results.mat', 'tpt_for_plot_us', 'deltaS_values');
        fprintf('Results saved to Entropy_vs_TPT_Results.mat and .png\n');
        
        % ... inside Section 4, after the initial plot command ...

hold(ax, 'on');

% --- Add Moving Average ---
window_size = 5; % Number of points to average over. Can be adjusted.
moving_avg = movmean(deltaS_values, window_size);
plot(ax, tpt_for_plot_us, moving_avg, '-r', 'LineWidth', 3, 'DisplayName', sprintf('Moving Avg. (N=%d)', window_size));

% --- Add Shaded "Driven" Region ---
% Let's define the reliable region as where the TPT histogram is, say, >10% of its peak
% This is a more robust way to define it.
[max_counts_fwd, ~] = max(counts_forward_LH);
[max_counts_rev, ~] = max(counts_reverse_HL);
reliable_bins_mask = (counts_forward_LH > 0.1*max_counts_fwd) & (counts_reverse_HL > 0.1*max_counts_rev);
reliable_t_start = bin_centers_us(find(reliable_bins_mask, 1, 'first'));
reliable_t_end = bin_centers_us(find(reliable_bins_mask, 1, 'last'));

ylim_vals = get(ax, 'YLim');
patch(ax, [reliable_t_start, reliable_t_end, reliable_t_end, reliable_t_start], ...
     [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
     'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Statistically Reliable Region');

% --- Calculate and Display Mean Entropy in Reliable Region ---
mean_deltaS_reliable = mean(deltaS_values(reliable_bins_mask(valid_bins_mask)));
text_str = sprintf('<\\DeltaS> = %.2f k_B', mean_deltaS_reliable);
text(ax, reliable_t_end, ylim_vals(2)*0.9, text_str, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
     'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w');


legend(ax, 'show');
hold(ax, 'off');

% Now saveas and the rest of the script...

        fprintf('\n--- Entropy Calculation Finished ---\n');
        cd(original_pwd);

    catch ME_Entropy
        fprintf(2,'\nERROR in calculateEntropyVsTPT_cluster: %s\n', ME_Entropy.message);
        disp(ME_Entropy.getReport());
        cd(original_pwd);
    end
end