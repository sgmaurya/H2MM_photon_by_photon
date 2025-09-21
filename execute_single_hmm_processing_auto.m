% execute_single_hmm_processing_auto.m
% A fully automated, non-interactive worker function to process HMM results.
% VERSION 3: FINAL & COMPLETE. This version's saving and plotting logic has been
%            restored to be IDENTICAL to the original processHMM_results_local.m,
%            ensuring all Viterbi-ready files are created and plots look the same.

function execute_single_hmm_processing_auto(results_path)
    
    % --- DIARY SETUP ---
    log_base_dir = fullfile(results_path, 'Logs');
    if ~exist(log_base_dir, 'dir'), mkdir(log_base_dir); end
    timestamp_for_diary = datestr(now, 'yyyy-mm-dd_HHMMSS');
    diary_filename = fullfile(log_base_dir, sprintf('hmm_processing_log_%s.txt', timestamp_for_diary));
    if exist(diary_filename, 'file'), delete(diary_filename); end
    diary(diary_filename);
    cleanupObj_diary = onCleanup(@() diary('off'));
    fprintf('--- Starting Automated Single HMM Post-Processing ---\n');
    fprintf('Processing results in: %s\n', results_path);
    
    try
        %% 0. Setup Paths
        project_base_path = fileparts(mfilename('fullpath'));
        paths.hmm_core = fullfile(project_base_path, 'hmm_core_functions'); 
        addpath(genpath(paths.hmm_core));
        
        %% 1. Load All Necessary Files
        fprintf('\n--- 1. Loading HMM Results and Data ---\n');
        
        loaded_results = load(fullfile(results_path, 'All_HMM_Results_Collected_Local.mat'));
        loaded_config = load(fullfile(results_path, 'HMM_Run_Config.mat'));
        loaded_data = load(fullfile(results_path, 'data_re_for_hmm.mat'));

        all_final_log_likelihoods = loaded_results.all_final_log_likelihoods;
        config_params_to_save = loaded_config.config_params_to_save;
        data = loaded_data.data_re_for_hmm;
        dt = config_params_to_save.dt_analysis_sec;
        Nstates = config_params_to_save.Nstates;

        %% 2. Automatically Select Best Model and Process
        fprintf('\n--- 2. Finding Best Model, Processing, and Saving All Files ---\n');
        
        [max_ll_value, best_idx] = max(all_final_log_likelihoods);
        fprintf('AUTOMATIC SELECTION: Best model is Guess #%d (LogL: %.4f).\n', best_idx, max_ll_value);

        best_model_params_struct = loaded_results.all_hmm_results_collected{best_idx};
        model.obsmat = best_model_params_struct.obsmat;
        model.prior = best_model_params_struct.prior;
        model.transmat = best_model_params_struct.transmat;
        model.Nstates = Nstates;
        
        % --- VERBATIM LOGIC FROM YOUR ORIGINAL SCRIPT ---
        % This block creates all the sorted parameters needed for saving and plotting.
        [Obs, B_sort_idx] = sort(model.obsmat(:,2)');
        Prio_sort = model.prior(B_sort_idx)';
        K_matrix_orig_style = (model.transmat - eye(Nstates)) / dt;
        K_sorted_orig_style = K_matrix_orig_style(B_sort_idx, :); 
        K_sorted_orig_style = K_sorted_orig_style(:, B_sort_idx);
        K = K_sorted_orig_style;
        Trans_sorted = model.transmat(B_sort_idx, :); 
        Trans_sorted = Trans_sorted(:, B_sort_idx);
        Trans = Trans_sorted;
        
        [eigenvect, eigenval_diag] = eig(Trans');
        [~, unit_eigenval_idx] = min(abs(diag(eigenval_diag) - 1)); 
        Prior = real(eigenvect(:, unit_eigenval_idx) / sum(eigenvect(:, unit_eigenval_idx)));
        
        AA = [B_sort_idx; Obs];
        Est_K = K; Est_B = B_sort_idx; Est_Prio_sort = Prio_sort; Est_Obs_sort = Obs; % Est_Obs_sort was missing
        Est_aVNorm = Prior; Est_AAA = Trans; Est_ass = best_idx;
        % --- END VERBATIM LOGIC ---
        
        % --- VERBATIM SAVE COMMANDS FROM YOUR ORIGINAL SCRIPT ---
        % This ensures all necessary files are created for the Viterbi step.
        save(fullfile(results_path, 'AA.mat'), 'AA');
        save(fullfile(results_path, 'HMM_parm.mat'), 'Obs', 'Trans', 'Prior', 'K', 'Prio_sort');
        save(fullfile(results_path, 'Est.mat'), 'Est_K', 'Est_B', 'Est_Prio_sort', 'Est_Obs_sort', 'Est_aVNorm', 'Est_AAA', 'Est_ass');
        
        best_model_overall = model; % 'model' now contains the Nstates field
        save(fullfile(results_path, 'Best_HMM_Model_Local.mat'), 'best_model_overall', 'config_params_to_save', 'best_idx');
        fprintf('All parameter and model files for Viterbi analysis have been saved.\n');

        %% 3. Generate and Save Plots (Identical to Original)
        fprintf('\n--- 3. Generating and Saving Plots ---\n');
        
        Realization_number = 5; 
        [~, ~, Recolored_FRET] = GenSimRecolored_hist(data, model, Realization_number);
        
        % --- VERBATIM PLOTTING CODE FROM YOUR ORIGINAL SCRIPT ---
        calcFRET_inline = @(d_traj) nnz(d_traj(:,2)==2)/size(d_traj,1); 
        fret_data_all_traj = cellfun(calcFRET_inline, data, 'UniformOutput', false);
        fret_data = [fret_data_all_traj{:}]; 
        fret_data(isnan(fret_data)) = []; 
        Recolored_FRET_all = cat(2, Recolored_FRET{:});

        fig_recolor_hist = figure('Name', 'FRET Histograms (Data vs Recolored)', 'Visible', 'off');
        bin_width = 0.02; edges = -0.01:bin_width:1.01;
        histogram(fret_data, 'BinEdges', edges, 'Normalization', 'pdf', 'FaceColor', [0.2 0.5 0.8], 'DisplayName', 'Original Data');
        hold on;
        h_recoloring = histcounts(Recolored_FRET_all, 'BinEdges', edges, 'Normalization', 'pdf');
        stairs(edges, [h_recoloring, h_recoloring(end)], 'LineWidth', 2, 'Color', 'k', 'DisplayName', sprintf('Recolored Model (%d realizations)', Realization_number));
        
        y_limits = get(gca, 'YLim');
        for i_state = 1:length(Obs), plot([Obs(i_state), Obs(i_state)], y_limits, 'r--', 'LineWidth', 1); end
        hold off;
        xlabel('FRET Efficiency'); ylabel('Probability Density');
        title(sprintf('FRET Distribution: Data vs. Recolored Model (Guess #%d)', best_idx));
        legend('Location', 'NorthEast'); grid on;
        % --- END VERBATIM PLOTTING CODE ---

        saveas(fig_recolor_hist, fullfile(results_path, 'FRET_Histogram_Recolored.png'));
        close(fig_recolor_hist);
        fprintf('Recoloring plot saved as PNG file.\n');
        
        fprintf('\n--- HMM Results Processing Finished Successfully ---\n');
        
    catch ME_process
        fprintf(2, '\n--- ERROR in Automated Post-Processing ---\n');
        fprintf(2, 'Error occurred in folder: %s\n', results_path);
        rethrow(ME_process);
    end
end