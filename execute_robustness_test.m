% execute_robustness_test_v2.m
% A new, powerful worker script for TPT validation.
% FINAL VERSION (Hypothesis 2): This version is modified to test the "failure mode".
% The iterative methods will now ALWAYS return a result, even if they don't
% converge, allowing us to see the distribution of non-converged TPTs.

function execute_robustness_test(results_path, ground_truth_tpts_us, test_params)
    
    % --- Load Data ---
    load(fullfile(results_path, 'Best_HMM_Model_Local.mat'), 'best_model_overall', 'config_params_to_save');
    load(fullfile(results_path, 'Viterbi_decoded_states_Q0.mat'), 'Q0');
    load(fullfile(results_path, 'data_re_for_hmm.mat'), 'data_re_for_hmm');
    
    % --- Prepare Model ---
    dt_model_sec = config_params_to_save.dt_analysis_sec;
    obsmat = best_model_overall.obsmat;
    [~, sorted_idx] = sort(obsmat(:,2));
    state_map(sorted_idx(1)) = 1; state_map(sorted_idx(2)) = 2;
    P_Red_L = obsmat(sorted_idx(1), 2); P_Red_H = obsmat(sorted_idx(2), 2);

    % --- Unpack Test Parameters ---
    test_mode = test_params.mode;
    rand_window = test_params.window;
    p_thresh = test_params.p_thresh;
    lambda_thresh = test_params.lambda;
    
    %% --- Main Calculation Loop ---
    fprintf('INFO: Running test in ''%s'' mode (Forced Result).\n', test_mode);
    TPTs_std_bayes = []; TPTs_iter_bayes = [];
    TPTs_std_llr = []; TPTs_iter_llr = [];
    
    for i_traj = 1:length(data_re_for_hmm)
        viterbi_path_orig = Q0{i_traj};
        photon_stream = data_re_for_hmm{i_traj};
        if isempty(viterbi_path_orig) || length(viterbi_path_orig) < 20, continue; end
        
        viterbi_path_mapped = state_map(viterbi_path_orig);
        true_transition_indices = find(diff(viterbi_path_mapped) ~= 0) + 1;
        
        start_indices = [];
        if strcmp(test_mode, 'ROBUSTNESS')
            for k=1:length(true_transition_indices), idx_viterbi = true_transition_indices(k); low_b = max(2, idx_viterbi - rand_window); high_b = min(length(viterbi_path_mapped) - 1, idx_viterbi + rand_window); if low_b < high_b, start_indices(end+1) = randi([low_b, high_b]); end; end
        elseif strcmp(test_mode, 'SPECIFICITY')
            for k=1:length(true_transition_indices), start_indices(end+1) = randi([10, length(viterbi_path_mapped) - 10]); end
        end
        
        for k_start = 1:length(start_indices)
            idx_random_start = start_indices(k_start);
            type = 'LH'; if rand > 0.5, type = 'HL'; end 
            colors = photon_stream(:,2);
            
            [exit_sb, entry_sb] = get_bayes_points(colors, idx_random_start, P_Red_H, P_Red_L, p_thresh, type);
            if ~isnan(exit_sb) && ~isnan(entry_sb) && entry_sb > exit_sb, TPTs_std_bayes(end+1) = (photon_stream(entry_sb, 1) - photon_stream(exit_sb, 1)) * dt_model_sec; end
            [exit_ib, entry_ib] = find_tpt_iterative_bayes_forced(colors, idx_random_start, P_Red_H, P_Red_L, p_thresh, type);
            if ~isnan(exit_ib) && ~isnan(entry_ib) && entry_ib > exit_ib, TPTs_iter_bayes(end+1) = (photon_stream(entry_ib, 1) - photon_stream(exit_ib, 1)) * dt_model_sec; end
            [exit_sl, entry_sl] = get_llr_points(colors, idx_random_start, P_Red_H, P_Red_L, lambda_thresh, type);
            if ~isnan(exit_sl) && ~isnan(entry_sl) && entry_sl > exit_sl, TPTs_std_llr(end+1) = (photon_stream(entry_sl, 1) - photon_stream(exit_sl, 1)) * dt_model_sec; end
            [exit_il, entry_il] = find_tpt_iterative_llr_forced(colors, idx_random_start, P_Red_H, P_Red_L, lambda_thresh, type);
            if ~isnan(exit_il) && ~isnan(entry_il) && entry_il > exit_il, TPTs_iter_llr(end+1) = (photon_stream(entry_il, 1) - photon_stream(exit_il, 1)) * dt_model_sec; end
        end
    end

    %% --- Plotting and Saving ---
    % ... (This section is unchanged) ...
    fig = figure('Name', ['TPT Validation: ' test_mode], 'Visible', 'off');
    hold on;
    bin_width = 2.0; 
    histogram(ground_truth_tpts_us, 'BinWidth', bin_width, 'Normalization', 'pdf', 'DisplayName', sprintf('Ground Truth (N=%d, Mean=%.2f us)', length(ground_truth_tpts_us), mean(ground_truth_tpts_us)), 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.6);
    histogram(TPTs_std_bayes*1e6, 'BinWidth', bin_width, 'Normalization', 'pdf', 'DisplayName', sprintf('Std. Bayes [Random] (N=%d, Mean=%.2f us)', length(TPTs_std_bayes), mean(TPTs_std_bayes*1e6)), 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 2.0, 'LineStyle', ':');
    histogram(TPTs_std_llr*1e6, 'BinWidth', bin_width, 'Normalization', 'pdf', 'DisplayName', sprintf('Std. LLR [Random] (N=%d, Mean=%.2f us)', length(TPTs_std_llr), mean(TPTs_std_llr*1e6)), 'DisplayStyle', 'stairs', 'EdgeColor', [0.6 0 0], 'LineWidth', 2.0, 'LineStyle', ':');
    histogram(TPTs_iter_bayes*1e6, 'BinWidth', bin_width, 'Normalization', 'pdf', 'DisplayName', sprintf('Iter. Bayes [Random] (N=%d, Mean=%.2f us)', length(TPTs_iter_bayes), mean(TPTs_iter_bayes*1e6)), 'DisplayStyle', 'stairs', 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 2.5);
    histogram(TPTs_iter_llr*1e6, 'BinWidth', bin_width, 'Normalization', 'pdf', 'DisplayName', sprintf('Iter. LLR [Random] (N=%d, Mean=%.2f us)', length(TPTs_iter_llr), mean(TPTs_iter_llr*1e6)), 'DisplayStyle', 'stairs', 'EdgeColor', [0.4660 0.6740 0.1880], 'LineWidth', 2.5);
    hold off; grid on; box on;
    title({['TPT Validation Test: ' test_mode ' (Forced Result)'], sprintf('Start Points: %s', strrep(test_mode, '_', ' '))}, 'FontSize', 14);
    xlabel('Transition Path Time (\mus)'); ylabel('Probability Density');
    legend('show', 'Location', 'NorthEast', 'FontSize', 10); set(gca, 'FontSize', 11);
    xlim([0, max(quantile(ground_truth_tpts_us, 0.99), 60)]);
    file_suffix = ['_random_viterbi_' lower(test_mode) '_forced.mat'];
    save(fullfile(results_path, ['TPT_Validation_Data' file_suffix]), 'TPTs_*', 'ground_truth_tpts_us', 'test_params');
    saveas(fig, fullfile(results_path, ['TPT_Validation_Plot' strrep(file_suffix, '.mat', '.png')]));
    savefig(fig, fullfile(results_path, ['TPT_Validation_Plot' strrep(file_suffix, '.mat', '.fig')]));
    close(fig);
end

%% --- SUB-FUNCTIONS WITH "FORCED RESULT" LOGIC ---
% <<< --- THIS IS THE KEY MODIFICATION --- >>>
function [e,en] = find_tpt_iterative_bayes_forced(c,s,ph,pl,p,t)
    max_iter=10; tol=2; cci=s; [e,en,ep,enp]=deal(NaN);
    for iter=1:max_iter
        [ep,enp]=get_bayes_points(c,cci,ph,pl,p,t);
        if isnan(ep)||isnan(enp), e=NaN; en=NaN; return; end
        nci=round((ep+enp)/2);
        if abs(nci-cci)<=tol, e=ep; en=enp; return; end % Converged
        cci=nci;
    end
    % If loop finishes, IT'S A NON-CONVERGENCE.
    % We now FORCE a result by returning the LAST calculated points.
    e=ep; en=enp;
end
function [e,en] = find_tpt_iterative_llr_forced(c,s,ph,pl,l,t)
    max_iter=500; tol=2; cci=s; [e,en,ep,enp]=deal(NaN);
    for iter=1:max_iter
        [ep,enp]=get_llr_points(c,cci,ph,pl,l,t);
        if isnan(ep)||isnan(enp), e=NaN; en=NaN; return; end
        nci=round((ep+enp)/2);
        if abs(nci-cci)<=tol, e=ep; en=enp; return; end % Converged
        cci=nci;
    end
    % If loop finishes, IT'S A NON-CONVERGENCE.
    % We now FORCE a result by returning the LAST calculated points.
    e=ep; en=enp;
end
% <<< --- END MODIFICATION --- >>>

function [idx_e,idx_en]=get_bayes_points(c,s,ph,pl,p,t),idx_e=findExitPoint_Bayesian(c,s,ph,pl,p,t);idx_en=findEntryPoint_Bayesian(c,s,ph,pl,p,t);end
function [idx_e,idx_en]=get_llr_points(c,s,ph,pl,l,t),idx_e=findExitPoint_Standard(c,s,ph,pl,l,t);idx_en=findEntryPoint_Standard(c,s,ph,pl,l,t);end
function idx=findEntryPoint_Standard(c,s,ph,pl,l,t),llr=0;idx=NaN;if strcmp(t,'LH'),lr=log(ph/pl);lg=log((1-ph)/(1-pl));else,lr=log(pl/ph);lg=log((1-pl)/(1-ph));end;for i=s:length(c),if c(i)==2,llr=llr+lr;else,llr=llr+lg;end;if llr>l,idx=i;return;end;end;end
function idx=findExitPoint_Standard(c,s,ph,pl,l,t),llr=0;idx=NaN;if strcmp(t,'LH'),lr=log(pl/ph);lg=log((1-pl)/(1-ph));else,lr=log(ph/pl);lg=log((1-ph)/(1-pl));end;for i=(s-1):-1:1,if c(i)==2,llr=llr+lr;else,llr=llr+lg;end;if llr>l,idx=i+1;return;end;end;end
function idx=findEntryPoint_Bayesian(c,s,ph,pl,p,t),llr=0;idx=NaN;if strcmp(t,'LH'),lr=log(ph/pl);lg=log((1-ph)/(1-pl));else,lr=log(pl/ph);lg=log((1-pl)/(1-ph));end;for i=s:length(c),if c(i)==2,llr=llr+lr;else,llr=llr+lg;end;if 1/(1+exp(-llr))>p,idx=i;return;end;end;end
function idx=findExitPoint_Bayesian(c,s,ph,pl,p,t),llr=0;idx=NaN;if strcmp(t,'LH'),lr=log(pl/ph);lg=log((1-pl)/(1-ph));else,lr=log(ph/pl);lg=log((1-ph)/(1-pl));end;for i=(s-1):-1:1,if c(i)==2,llr=llr+lr;else,llr=llr+lg;end;if 1/(1+exp(-llr))>p,idx=i+1;return;end;end;end