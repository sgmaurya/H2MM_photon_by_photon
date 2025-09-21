% --- START OF CORRECTED FILE processHMM_results_local.m ---

% processHMM_results_local.m
% This script processes the output from runHMM_local.m,
% finds the best model, performs recoloring, and generates graphs.
% MODIFIED: Can now analyze both main runs and pilot runs.
% CORRECTED: Now saves 'Best_HMM_Model_Local.mat' for direct use by runViterbi_local.m

function processHMM_results_local()
    clearvars -except break_debug;
    close all;
    clc;

    fprintf('--- Starting HMM Results Processing ---\n');

    %% 0. Setup Paths 
    project_base_path = fileparts(mfilename('fullpath')); 
    paths.hmm_core = fullfile(project_base_path, 'hmm_core_functions'); 
    paths.initial_guess = fullfile(project_base_path, 'initial_guess_functions'); 
    
    if exist(paths.hmm_core, 'dir'), addpath(genpath(paths.hmm_core)); end
    if exist(paths.initial_guess, 'dir'), addpath(genpath(paths.initial_guess)); end

    original_pwd = pwd; 
    
    try 
        %% 1. Select HMM Run Directory and Load Results
        fprintf('\n--- 1. Selecting HMM Run Directory and Loading Results ---\n');
        
        results_base_path = fullfile(project_base_path, 'HMM_Local_Results');
        if ~exist(results_base_path, 'dir')
            error('Base results directory not found: %s\nPlease run runHMM_local.m first.', results_base_path);
        end
        
        % --- MODIFIED PROMPT ---
        dialog_prompt = 'Select the HMM_Run_* or HMM_Run_Pilot_Results_* directory to process';
        selected_run_path = uigetdir(results_base_path, dialog_prompt);
        if isequal(selected_run_path, 0), disp('User cancelled directory selection. Exiting.'); return; end 
        
        fprintf('Processing results from: %s\n', selected_run_path);
        
        cd(selected_run_path); 
        fprintf('Changed current directory to: %s\n', selected_run_path);

        % The rest of the script uses standardized filenames, so it works for both run types.
        fprintf('Loading All_HMM_Results_Collected_Local.mat...\n');
        % Check for pilot-run specific variable names and adapt if necessary for compatibility
        try
            loaded_results = load('All_HMM_Results_Collected_Local.mat');
            if isfield(loaded_results, 'all_hmm_results_collected_pilot') % From Pilot Run
                all_hmm_results_collected = loaded_results.all_hmm_results_collected_pilot;
                all_final_log_likelihoods = loaded_results.all_final_log_likelihoods_pilot;
                InitialGuess_all = loaded_results.InitialGuess_all_pilot;
                config_params_to_save = loaded_results.config_params_to_save_pilot;
                 fprintf('INFO: Loaded pilot run result variables.\n');
            else % From Main Run
                all_hmm_results_collected = loaded_results.all_hmm_results_collected;
                all_final_log_likelihoods = loaded_results.all_final_log_likelihoods;
                InitialGuess_all = loaded_results.InitialGuess_all;
                config_params_to_save = loaded_results.config_params_to_save;
                fprintf('INFO: Loaded main run result variables.\n');
            end
        catch ME_load_results
            error('Failed to load or parse All_HMM_Results_Collected_Local.mat. Error: %s', ME_load_results.message);
        end

        fprintf('Loading data_re_for_hmm.mat...\n');
        loaded_data = load('data_re_for_hmm.mat');
        if isfield(loaded_data, 'data_re_for_hmm_pilot') % From Pilot Run
            data = loaded_data.data_re_for_hmm_pilot;
        else % From Main Run
            data = loaded_data.data_re_for_hmm;
        end

        Nstates = config_params_to_save.Nstates;
        dt = config_params_to_save.dt_analysis_sec; 

        if isempty(all_hmm_results_collected)
            error('No HMM results found. Aborting.');
        end

        %% 2. Adapt and Run Recoloring_Param Logic
        fprintf('\n--- 2. Finding Best Model and Performing Recoloring ---\n');
        
        num_guesses = length(all_final_log_likelihoods);
        valid_ll_indices = find(~isnan(all_final_log_likelihoods));

        if isempty(valid_ll_indices)
            error('All HMM runs resulted in NaN LogLikelihood. Cannot proceed.');
        end
        
        fig_ll_plot = figure('Name', 'Log-Likelihoods vs Initial Guess');
        plot(1:num_guesses, all_final_log_likelihoods, 'o-', 'MarkerSize', 8, 'LineWidth', 1.5);
        xlabel('Initial Guess #'); ylabel('Final Log-Likelihood');
        title('Log-Likelihoods of HMM Runs'); grid on;
        
        [max_ll_value, auto_selected_ass] = max(all_final_log_likelihoods);
        
        prompt_ass = {sprintf('Automatically selected Guess #%d (LogL: %.4f).\nEnter chosen guess # (or press Enter for auto):', auto_selected_ass, max_ll_value)};
        answer_ass = inputdlg(prompt_ass, 'Select Best Guess', [1 60], {num2str(auto_selected_ass)});

        if isempty(answer_ass) || isempty(answer_ass{1}), ass = auto_selected_ass;
        else, ass = str2double(answer_ass{1}); end
        
        best_model_params_struct = all_hmm_results_collected{ass};
        if isempty(best_model_params_struct) || ~isfield(best_model_params_struct, 'obsmat')
            error('Selected best guess #%d does not contain valid model parameters.', ass);
        end

        model.Nstates = Nstates; 
        model.obsmat = best_model_params_struct.obsmat;
        model.prior = best_model_params_struct.prior;   
        model.transmat = best_model_params_struct.transmat;

        Realization_number = 5; 
        [~, photon_color, Recolored_FRET] = GenSimRecolored_hist(data, model, Realization_number); 

        calcFRET_inline = @(d_traj) nnz(d_traj(:,2)==2)/size(d_traj,1); 
        fret_data_all_traj = cellfun(calcFRET_inline, data, 'UniformOutput', false);
        fret_data = [fret_data_all_traj{:}]; 
        fret_data(isnan(fret_data)) = []; 
        Recolored_FRET_all = cat(2, Recolored_FRET{:});

        fig_recolor_hist = figure('Name', 'FRET Histograms (Data vs Recolored)', 'Color', 'w', 'Position', [150, 150, 700, 500]);
        bin_width = 0.02; edges = -0.01:bin_width:1.01;
        histogram(fret_data, 'BinEdges', edges, 'Normalization', 'pdf', 'FaceColor', [0.2 0.5 0.8], 'DisplayName', 'Original Data');
        hold on;
        h_recoloring = histcounts(Recolored_FRET_all, 'BinEdges', edges, 'Normalization', 'pdf');
        stairs(edges, [h_recoloring, h_recoloring(end)], 'LineWidth', 2, 'Color', 'k', 'DisplayName', sprintf('Recolored Model (%d realizations)', Realization_number));
        
        obsmat_chosen = model.obsmat; prior_chosen_from_hmm = model.prior; transmat_chosen = model.transmat;
        [Obs_sort, B_sort_idx] = sort(obsmat_chosen(:,2)'); 
        Prio_sort = prior_chosen_from_hmm(B_sort_idx)'; 
        K_matrix_orig_style = (transmat_chosen - eye(Nstates)) / dt; 
        K_sorted_orig_style = K_matrix_orig_style(B_sort_idx, :); K_sorted_orig_style = K_sorted_orig_style(:, B_sort_idx); 
        K = K_sorted_orig_style; 
        Trans_sorted = transmat_chosen(B_sort_idx, :); Trans_sorted = Trans_sorted(:, B_sort_idx);
        Trans = Trans_sorted; Obs = Obs_sort; 

        try
            [eigenvect, eigenval_diag] = eig(Trans');
            [~, unit_eigenval_idx] = min(abs(diag(eigenval_diag) - 1)); 
            Prior = real(eigenvect(:, unit_eigenval_idx) / sum(eigenvect(:, unit_eigenval_idx)));
        catch ME_eig, warning(ME_eig.identifier, 'Could not calculate steady-state: %s', ME_eig.message); Prior = Prio_sort'; end
        
        AA = [B_sort_idx; Obs];
        save('AA.mat','AA'); 
        save('HMM_parm.mat','Obs','Trans','Prior','K','Prio_sort'); 
        Est_K = K; Est_B = B_sort_idx; Est_Prio_sort = Prio_sort; Est_Obs_sort = Obs_sort; Est_aVNorm = Prior; Est_AAA = Trans; Est_ass = ass;      
        save('Est.mat', 'Est_K', 'Est_B', 'Est_Prio_sort', 'Est_Obs_sort', 'Est_aVNorm', 'Est_AAA', 'Est_ass');
        fprintf('Sorted HMM parameters saved.\n');
        
        % --- START: ADDED CODE TO SAVE THE BEST MODEL FOR VITERBI ---
        % Create variables with the names expected by runViterbi_local.m
        best_model_overall = model; % 'model' contains the selected prior, obsmat, and transmat
        best_idx = ass;             % 'ass' is the index of the best guess you selected
        
        % Save the file runViterbi_local expects
        fprintf('Saving best model for Viterbi processing to Best_HMM_Model_Local.mat...\n');
        save('Best_HMM_Model_Local.mat', 'best_model_overall', 'config_params_to_save', 'best_idx');
        % --- END: ADDED CODE ---
        
        y_limits = get(gca, 'YLim');
        for i_state = 1:length(Obs), plot([Obs(i_state), Obs(i_state)], y_limits, 'r--', 'LineWidth', 1); end
        hold off;
        xlabel('FRET Efficiency'); ylabel('Probability Density');
        title(sprintf('FRET Distribution: Data vs. Recolored Model (Guess #%d)', ass));
        legend('Location', 'NorthEast'); grid on;
        saveas(fig_recolor_hist, 'FRET_Histogram_Recolored.png');
        
        %% 3. Call HMMgraph
        fprintf('\n--- 3. Generating HMM Graph ---\n');
        try 
            % --- THIS IS THE FIX ---
            % HMMgraph.m is a SCRIPT, not a function. It expects variables to
            % exist in the workspace with specific names. We will create them
            % here, just before calling it, to ensure compatibility.
            
            % It expects 'indx', but our config variable is 'indx_guessmodel'
            indx = config_params_to_save.indx_guessmodel; 
            
            % It expects 'guessName', but our config variable is 'guessName'
            % (This one is already correct, but we'll include it for clarity)
            guessName = config_params_to_save.guessName;

            % Save the temporary flag file that HMMgraph.m will load.
            % We save it in the current directory.
            save('HMMtypeFlag.mat', 'indx', 'guessName');
            % --- END FIX ---
            
            HMMgraph; % Call the original, unmodified graphing script
            
            % Clean up the temporary file
            if exist('HMMtypeFlag.mat', 'file'), delete('HMMtypeFlag.mat'); end
            
        catch ME_graph
            warning('HMMGraphFail', 'Failed to generate HMM graph: %s', ME_graph.message); 
            disp(ME_graph.getReport()); % Display full error for better debugging
        end 

        fprintf('\n--- HMM Results Processing Finished ---\n');
        cd(original_pwd); 

    catch ME_main_process 
        fprintf('ERROR in processHMM_results_local: %s\n', ME_main_process.message);
        disp(ME_main_process.getReport());
        if exist('original_pwd', 'var') && isfolder(original_pwd), cd(original_pwd); end
    end 
end

% --- END OF CORRECTED FILE processHMM_results_local.m ---