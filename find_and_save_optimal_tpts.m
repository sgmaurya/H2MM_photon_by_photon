% find_and_save_optimal_tpts.m
% A powerful worker function that finds the optimal Lambda and P_thresh
% to match a ground truth, then saves the results.
% VERSION 5: FINAL. Corrects the method name strings used in the calibration
%            sweep to match the helper function, fixing the "flat line" bug.

function [optimal_lambda, optimal_p_thresh] = find_and_save_optimal_tpts(results_path, ground_truth_tpt_us)
    
    %% 1. Load Data and Prepare Model
    % ... (This section is correct) ...
    load(fullfile(results_path, 'Best_HMM_Model_Local.mat'), 'best_model_overall', 'config_params_to_save');
    load(fullfile(results_path, 'Viterbi_decoded_states_Q0.mat'), 'Q0');
    load(fullfile(results_path, 'data_re_for_hmm.mat'), 'data_re_for_hmm');
    dt_model_sec = config_params_to_save.dt_analysis_sec;
    obsmat = best_model_overall.obsmat;
    [~, sorted_idx] = sort(obsmat(:,2));
    state_map(sorted_idx(1)) = 1; state_map(sorted_idx(2)) = 2;
    P_Red_L = obsmat(sorted_idx(1), 2); P_Red_H = obsmat(sorted_idx(2), 2);

    %% 2. Pre-calculate All Transitions
    % ... (This section is correct) ...
    transitions = struct('type', {}, 'center_idx', {}, 'photon_colors', {}, 'photon_times_dt_units', {});
    for i_traj = 1:length(data_re_for_hmm)
        viterbi_path_orig = Q0{i_traj};
        if isempty(viterbi_path_orig) || length(viterbi_path_orig) < 2, continue; end
        viterbi_path_mapped = state_map(viterbi_path_orig);
        transition_indices = find(diff(viterbi_path_mapped) ~= 0) + 1;
        for k_trans = 1:length(transition_indices)
            idx_change = transition_indices(k_trans);
            state_before = viterbi_path_mapped(idx_change - 1); state_after  = viterbi_path_mapped(idx_change);
            if state_before == 1 && state_after == 2, trans_type = 'LH';
            elseif state_before == 2 && state_after == 1, trans_type = 'HL';
            else, continue; end
            new_trans.type = trans_type; new_trans.center_idx = idx_change;
            new_trans.photon_colors = data_re_for_hmm{i_traj}(:,2);
            new_trans.photon_times_dt_units = data_re_for_hmm{i_traj}(:,1);
            transitions(end+1) = new_trans;
        end
    end
    if isempty(transitions), optimal_lambda = NaN; optimal_p_thresh = NaN; return; end
    
    %% 3. Calibration Sweep (with corrected method names)
    fprintf('INFO: Finding optimal thresholds by matching mean TPT to %.2f us...\n', ground_truth_tpt_us);
    prob_thresholds = linspace(0.70, 0.995, 40);
    
    % <<< --- FIX #1: Use the full, correct method name string --- >>>
    mean_tpts_bayes = arrayfun(@(p) mean([calculate_all_tpts('Iterative Bayesian (Convergent)', transitions, p, P_Red_H, P_Red_L, dt_model_sec).tpt])*1e6, prob_thresholds);
    
    [~, best_idx_bayes] = min(abs(mean_tpts_bayes - ground_truth_tpt_us));
    optimal_p_thresh = prob_thresholds(best_idx_bayes);

    lambda_thresholds = log(prob_thresholds ./ (1 - prob_thresholds));
    
    % <<< --- FIX #2: Use the full, correct method name string --- >>>
    mean_tpts_llr = arrayfun(@(l) mean([calculate_all_tpts('Iterative LLR (Convergent)', transitions, l, P_Red_H, P_Red_L, dt_model_sec).tpt])*1e6, lambda_thresholds);
    
    [~, best_idx_llr] = min(abs(mean_tpts_llr - ground_truth_tpt_us));
    optimal_lambda = lambda_thresholds(best_idx_llr);
    
    %% 4. Recalculate and Save TPT Distributions
    % ... (This section is correct and was already using the right names) ...
    fprintf('INFO: Recalculating TPTs with optimal thresholds and saving results...\n');
    params_ib.method = 'Iterative Bayesian (Convergent)'; params_ib.prob_thresh = optimal_p_thresh;
    [TPTs_sec_HL_ib, TPTs_sec_LH_ib] = calculate_all_tpts_by_type(params_ib.method, transitions, optimal_p_thresh, P_Red_H, P_Red_L, dt_model_sec);
    save_results_file(results_path, params_ib, TPTs_sec_HL_ib, TPTs_sec_LH_ib);
    params_sb.method = 'Standard Bayesian'; params_sb.prob_thresh = optimal_p_thresh;
    [TPTs_sec_HL_sb, TPTs_sec_LH_sb] = calculate_all_tpts_by_type(params_sb.method, transitions, optimal_p_thresh, P_Red_H, P_Red_L, dt_model_sec);
    save_results_file(results_path, params_sb, TPTs_sec_HL_sb, TPTs_sec_LH_sb);
    params_il.method = 'Iterative LLR (Convergent)'; params_il.lambda_thresh = optimal_lambda;
    [TPTs_sec_HL_il, TPTs_sec_LH_il] = calculate_all_tpts_by_type(params_il.method, transitions, optimal_lambda, P_Red_H, P_Red_L, dt_model_sec);
    save_results_file(results_path, params_il, TPTs_sec_HL_il, TPTs_sec_LH_il);
    params_sl.method = 'Standard LLR'; params_sl.lambda_thresh = optimal_lambda;
    [TPTs_sec_HL_sl, TPTs_sec_LH_sl] = calculate_all_tpts_by_type(params_sl.method, transitions, optimal_lambda, P_Red_H, P_Red_L, dt_model_sec);
    save_results_file(results_path, params_sl, TPTs_sec_HL_sl, TPTs_sec_LH_sl);
end

% --- All helper functions below are correct and unchanged ---
function save_results_file(results_path, params, TPTs_sec_HL, TPTs_sec_LH)
    safe_name = strrep(params.method, ' ', '_'); safe_name = strrep(safe_name, '(', ''); safe_name = strrep(safe_name, ')', '');
    filename = sprintf('TPT_Results_%s.mat', safe_name);
    save_struct.params = params; save_struct.TPTs_sec_HL = TPTs_sec_HL; save_struct.TPTs_sec_LH = TPTs_sec_LH;
    save(fullfile(results_path, filename), '-struct', 'save_struct', '-v7.3');
    fprintf('Saved: %s\n', filename);
end
function [tpts_hl, tpts_lh] = calculate_all_tpts_by_type(method, transitions, threshold, pR_H, pR_L, dt)
    tpts_hl = []; tpts_lh = []; all_tpts = calculate_all_tpts(method, transitions, threshold, pR_H, pR_L, dt);
    for i=1:length(all_tpts), if strcmp(all_tpts(i).type, 'HL'), tpts_hl(end+1) = all_tpts(i).tpt; else, tpts_lh(end+1) = all_tpts(i).tpt; end; end
end
function tpt_results = calculate_all_tpts(method, transitions, threshold, pR_H, pR_L, dt)
    tpt_results = struct('tpt', {}, 'type', {}); iter_params.max_iter = 10; iter_params.conv_tol_photons = 2;
    for k = 1:length(transitions)
        trans = transitions(k); current_center_idx = trans.center_idx; [idx_exit, idx_entry] = deal(NaN);
        switch method
            case {'Iterative Bayesian (Convergent)', 'Iterative LLR (Convergent)'}
                for iter = 1:iter_params.max_iter
                    if contains(method, 'Bayesian'), [exit_pass, entry_pass] = get_bayes_points(trans, current_center_idx, pR_H, pR_L, threshold);
                    else, [exit_pass, entry_pass] = get_llr_points(trans, current_center_idx, pR_H, pR_L, threshold); end
                    if isnan(exit_pass) || isnan(entry_pass), break; end
                    new_center_idx = round((exit_pass + entry_pass) / 2);
                    if abs(new_center_idx - current_center_idx) <= iter_params.conv_tol_photons, idx_exit = exit_pass; idx_entry = entry_pass; break; end
                    current_center_idx = new_center_idx;
                end
                if isnan(idx_exit), idx_exit = exit_pass; idx_entry = entry_pass; end
            case 'Standard Bayesian', [idx_exit, idx_entry] = get_bayes_points(trans, current_center_idx, pR_H, pR_L, threshold);
            case 'Standard LLR', [idx_exit, idx_entry] = get_llr_points(trans, current_center_idx, pR_H, pR_L, threshold);
        end
        if ~isnan(idx_exit) && ~isnan(idx_entry) && idx_entry > idx_exit
            tpt_sec = (trans.photon_times_dt_units(idx_entry) - trans.photon_times_dt_units(idx_exit)) * dt;
            if tpt_sec >= 0, tpt_results(end+1) = struct('tpt', tpt_sec, 'type', trans.type); end
        end
    end
    if isempty(tpt_results), tpt_results = struct('tpt',{NaN},'type',{''}); end
end
function [idx_e, idx_en] = get_bayes_points(t,s,ph,pl,p), idx_e=findExitPoint_Bayesian(t.photon_colors,s,ph,pl,p,t.type); idx_en=findEntryPoint_Bayesian(t.photon_colors,s,ph,pl,p,t.type); end
function [idx_e, idx_en] = get_llr_points(t,s,ph,pl,l), idx_e=findExitPoint_Standard(t.photon_colors,s,ph,pl,l,t.type); idx_en=findEntryPoint_Standard(t.photon_colors,s,ph,pl,l,t.type); end
function idx = findEntryPoint_Standard(colors, start_idx, pR_H, pR_L, lambda, type), LLR=0; idx=NaN; if strcmp(type,'LH'), lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); else, lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); end; for i=start_idx:length(colors), if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end; if LLR>lambda, idx=i; return; end; end; end
function idx = findExitPoint_Standard(colors, start_idx, pR_H, pR_L, lambda, type), LLR=0; idx=NaN; if strcmp(type,'LH'), lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); else, lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); end; for i=(start_idx-1):-1:1, if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end; if LLR>lambda, idx=i+1; return; end; end; end
function idx = findEntryPoint_Bayesian(colors, start_idx, pR_H, pR_L, p_thresh, type), LLR=0; idx=NaN; if strcmp(type,'LH'), lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); else, lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); end; for i=start_idx:length(colors), if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end; prob = 1 / (1 + exp(-LLR)); if prob > p_thresh, idx=i; return; end; end; end
function idx = findExitPoint_Bayesian(colors, start_idx, pR_H, pR_L, p_thresh, type), LLR=0; idx=NaN; if strcmp(type,'LH'), lr=log(pR_L/pR_H); lg=log((1-pR_L)/(1-pR_H)); else, lr=log(pR_H/pR_L); lg=log((1-pR_H)/(1-pR_L)); end; for i=(start_idx-1):-1:1, if colors(i)==2, LLR=LLR+lr; else, LLR=LLR+lg; end; prob = 1 / (1 + exp(-LLR)); if prob > p_thresh, idx=i+1; return; end; end; end