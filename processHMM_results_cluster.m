% processHMM_results_cluster.m (v2 - Corrected dt Handling)
%
% This script processes the output from a cluster-based HMM run.
% It reads the complex HMM_output.mat file, helps the user select the best
% model from multiple initial guesses, and crucially, PROMPTS the user for
% the correct HMM time resolution (dt) to ensure accurate rate calculations
% and recoloring.
%

function processHMM_results_cluster()
    clearvars -except break_debug; close all; clc;
    fprintf('--- Starting HMM Results Processing (Cluster Workflow) ---\n');

    original_pwd = pwd; 
    
    try 
        %% 1. Select Cluster HMM Run Directory and Load Results
        fprintf('\n--- 1. Selecting Cluster HMM Run Directory and Loading Results ---\n');
        
        selected_run_path = uigetdir(pwd, 'Select the Cluster HMM Run directory (containing HMM_output.mat, data.mat)');
        if isequal(selected_run_path, 0), disp('User cancelled directory selection. Exiting.'); return; end 
        
        fprintf('Processing results from: %s\n', selected_run_path);
        cd(selected_run_path);

        fprintf('Loading HMM_output.mat...\n');
        if ~exist('HMM_output.mat', 'file'), error('PROCESS_CLUSTER: HMM_output.mat not found.'); end
        loaded_vars_hmm = load('HMM_output.mat');

        fprintf('Loading data.mat...\n');
        if ~exist('data.mat', 'file'), error('PROCESS_CLUSTER: data.mat not found.'); end
        loaded_data = load('data.mat');
        data = loaded_data.data;

        %% 2. Unpack Cluster HMM Output and Find Best Model
        fprintf('\n--- 2. Unpacking Cluster Output and Finding Best Model ---\n');

        Nstates = loaded_vars_hmm.Nstates;
        obsmat_iterations = loaded_vars_hmm.obsmat{1,Nstates};
        transmat_iterations = loaded_vars_hmm.transmat{1,Nstates};
        prior_iterations = loaded_vars_hmm.prior{1,Nstates};
        num_guesses = length(obsmat_iterations);
        
        final_LL_for_selection = -inf(1, num_guesses);
        if isfield(loaded_vars_hmm, 'LL_end')
            ll_cell = loaded_vars_hmm.LL_end{1, Nstates};
            if iscell(ll_cell) && all(cellfun(@isnumeric, ll_cell))
                ll_vec = cell2mat(ll_cell);
                final_LL_for_selection(1:length(ll_vec)) = ll_vec;
            end
        else, error('PROCESS_CLUSTER: No usable LogLikelihood vector found.'); end

        figure('Name', 'Log-Likelihoods vs Initial Guess (Cluster)');
        plot(1:num_guesses, final_LL_for_selection, 'o-', 'MarkerSize', 8, 'LineWidth', 1.5);
        xlabel('Initial Guess #'); ylabel('Final Log-Likelihood');
        title('Log-Likelihoods of Cluster HMM Runs'); grid on;
        
        [max_ll_value, auto_selected_ass] = max(final_LL_for_selection);
        
        prompt_ass = {sprintf('Automatically selected Guess #%d (LogL: %.4f).\nEnter chosen guess # (or press Enter for auto):', auto_selected_ass, max_ll_value)};
        answer_ass = inputdlg(prompt_ass, 'Select Best Guess', [1 60], {num2str(auto_selected_ass)});
        if isempty(answer_ass) || isempty(answer_ass{1}), ass = auto_selected_ass;
        else, ass = str2double(answer_ass{1}); end
        
        best_model_params_struct.obsmat = obsmat_iterations{ass};
        best_model_params_struct.transmat = transmat_iterations{ass};
        best_model_params_struct.prior = prior_iterations{ass};
        if isempty(best_model_params_struct.obsmat), error('Selected model #%d has empty parameters.', ass); end

        % <<< START OF CORRECTION >>>
        % CRITICAL STEP: The 'dt' used for HMM training MUST be known.
        % It is not saved in HMM_output.mat, so we must ask the user.
        dt = NaN;
        found_dt = false;
        if isfield(loaded_vars_hmm, 'dt_analysis_sec'), dt = loaded_vars_hmm.dt_analysis_sec; found_dt = true;
        elseif isfield(loaded_vars_hmm, 'dt'), dt = loaded_vars_hmm.dt; found_dt = true; end

        if ~found_dt || isnan(dt)
            ans_dt = inputdlg({'CRITICAL: dt not found in HMM_output.mat. \nEnter the HMM time resolution (dt in seconds) used for TRAINING (e.g., 1e-6 for microseconds):'}, ...
                              'Input HMM Training dt', [1 100], {'1e-6'});
            if isempty(ans_dt), error('User cancelled dt input. Cannot proceed.'); end
            dt = str2double(ans_dt{1});
            if isnan(dt) || dt <= 0, error('Invalid dt entered.'); end
        end
        fprintf('INFO: Using dt = %.2e s for rate calculations and recoloring.\n', dt);
        % <<< END OF CORRECTION >>>

        %% 3. Perform Recoloring and Save Sorted Parameters
        fprintf('\n--- 3. Performing Recoloring and Saving HMM Parameters ---\n');

        model.Nstates = Nstates; 
        model.obsmat = best_model_params_struct.obsmat;
        model.prior = best_model_params_struct.prior;   
        model.transmat = best_model_params_struct.transmat;

        Realization_number = 5; 
        [~, ~, Recolored_FRET] = GenSimRecolored_hist(data, model, Realization_number); 

        calcFRET_inline = @(d_traj) nnz(d_traj(:,2)==2)/size(d_traj,1); 
        fret_data_all_traj = cellfun(calcFRET_inline, data, 'UniformOutput', false);
        fret_data = [fret_data_all_traj{:}]; 
        fret_data(isnan(fret_data)) = []; 
        Recolored_FRET_all = cat(2, Recolored_FRET{:});

        fig_recolor_hist = figure('Name', 'FRET Histograms (Data vs Recolored)', 'Position', [150, 150, 700, 500]);
        bin_width = 0.02; edges = -0.01:bin_width:1.01;
        histogram(fret_data, 'BinEdges', edges, 'Normalization', 'pdf', 'FaceColor', [0.2 0.5 0.8], 'DisplayName', 'Original Data');
        hold on;
        h_recoloring = histcounts(Recolored_FRET_all, 'BinEdges', edges, 'Normalization', 'pdf');
        stairs(edges, [h_recoloring, h_recoloring(end)], 'LineWidth', 2, 'Color', 'k', 'DisplayName', sprintf('Recolored Model (%d realizations)', Realization_number));
        
        obsmat_chosen = model.obsmat; prior_chosen = model.prior; transmat_chosen = model.transmat;
        [Obs_sort, B_sort_idx] = sort(obsmat_chosen(:,2)'); 
        Prio_sort = prior_chosen(B_sort_idx)'; 
        K_matrix_orig_style = (transmat_chosen - eye(Nstates)) / dt; 
        K_sorted_orig_style = K_matrix_orig_style(B_sort_idx, :); K_sorted_orig_style = K_sorted_orig_style(:, B_sort_idx); 
        K = K_sorted_orig_style; 
        Trans_sorted = transmat_chosen(B_sort_idx, :); Trans_sorted = Trans_sorted(:, B_sort_idx);
        Trans = Trans_sorted; Obs = Obs_sort; 

        [eigenvect, ~] = eig(Trans'); [~, unit_eigenval_idx] = min(abs(diag(eig(Trans')) - 1));
        Prior = real(eigenvect(:, unit_eigenval_idx) / sum(eigenvect(:, unit_eigenval_idx)));
        
        AA = [B_sort_idx; Obs]; save('AA.mat','AA'); 
        save('HMM_parm.mat','Obs','Trans','Prior','K','Prio_sort'); 
        Est_K = K; Est_B = B_sort_idx; Est_Prio_sort = Prio_sort; Est_Obs_sort = Obs_sort; Est_aVNorm = Prior; Est_AAA = Trans; Est_ass = ass;      
        save('Est.mat', 'Est_K', 'Est_B', 'Est_Prio_sort', 'Est_Obs_sort', 'Est_aVNorm', 'Est_AAA', 'Est_ass');
        fprintf('Sorted HMM parameters saved.\n');
        
        best_model_cluster = model;
        best_idx_cluster = ass;
        config_params_cluster.Nstates = Nstates;
        config_params_cluster.dt_analysis_sec = dt;
        save('Best_HMM_Model_Cluster.mat', 'best_model_cluster', 'config_params_cluster', 'best_idx_cluster');
        fprintf('Best model saved to Best_HMM_Model_Cluster.mat.\n');
        
        y_limits = get(gca, 'YLim');
        for i_state = 1:length(Obs), plot([Obs(i_state), Obs(i_state)], y_limits, 'r--', 'LineWidth', 1); end
        hold off; xlabel('FRET Efficiency'); ylabel('Probability Density');
        title(sprintf('FRET Distribution: Data vs. Recolored Model (Guess #%d)', ass));
        legend('Location', 'NorthEast'); grid on;
        saveas(fig_recolor_hist, 'FRET_Histogram_Recolored_Cluster.png');
        
        %% 4. Call HMMgraph
        fprintf('\n--- 4. Generating HMM Graph ---\n');
        try, HMMgraph;
        catch ME_graph, warning('HMMGraphFail', 'Failed to generate HMM graph: %s', ME_graph.message); end 

        fprintf('\n--- HMM Results Processing (Cluster) Finished ---\n');
        cd(original_pwd); 

    catch ME_main_process 
        fprintf('ERROR in processHMM_results_cluster: %s\n', ME_main_process.message);
        disp(ME_main_process.getReport());
        if exist('original_pwd', 'var') && isfolder(original_pwd), cd(original_pwd); end
    end 
end