function runHMM_Nstates_Scan()
    % runHMM_Nstates_Scan.m
    % Wraps HMM analysis to scan over a range of Nstates, calculate AIC/BIC,
    % and help determine the optimal number of states.

    clearvars -except break_debug; 
    close all;                     
    clc;                           
    
    % --- DIARY SETUP ---
    diary_filename_final = ''; 
    cleanupObj_diary = [];     

    try 
        project_base_path_for_log = fileparts(mfilename('fullpath')); 
        log_base_dir = fullfile(project_base_path_for_log, 'HMM_Run_Logs');
        if ~exist(log_base_dir, 'dir'), mkdir(log_base_dir); end
        timestamp_for_diary = datestr(now, 'yyyy-mm-dd_HHMMSS');
        diary_filename_attempt = fullfile(log_base_dir, sprintf('hmm_nstates_scan_log_%s.txt', timestamp_for_diary));
        
        if exist(diary_filename_attempt, 'file')
            fprintf('INFO: Attempting to delete existing diary file: %s\n', diary_filename_attempt);
            try, delete(diary_filename_attempt); fprintf('INFO: Successfully deleted old diary file.\n');
            catch ME_del_diary_init, warning('DIARY_SETUP: Could not delete existing diary file: %s. Error: %s.', diary_filename_attempt, getReport(ME_del_diary_init, 'basic')); end
        end
        
        diary(diary_filename_attempt);
        diary_filename_final = diary_filename_attempt; 
        cleanupObj_diary = onCleanup(@() diary('off'));
        
        fprintf('--- Starting Nstates Scan for HMM Analysis ---\n'); 
        fprintf('INFO: Command window output is being logged to: %s\n', diary_filename_final);

    catch ME_diary_setup
        warning('CRITICAL_DIARY_ERROR: Failed during diary setup. Error: %s. Logging to diary might not be active.', getReport(ME_diary_setup, 'basic'));
        if strcmp(get(0,'Diary'),'on'), diary('off'); end
    end
    % --- End of Diary Setup ---

    try % MAIN TRY BLOCK FOR THE REST OF THE SCRIPT
        %% 0. Setup Paths
        project_base_path = fileparts(mfilename('fullpath'));
        fprintf('INFO: Project base path: %s\n', project_base_path);

        paths.hmm_core = fullfile(project_base_path, 'hmm_core_functions');
        paths.initial_guess = fullfile(project_base_path, 'initial_guess_functions');
        paths.output_base = fullfile(project_base_path, 'HMM_Nstates_Scan_Results');

        if ~isfolder(paths.hmm_core), error('Dir "hmm_core_functions" not found: %s', paths.hmm_core); end
        if ~isfolder(paths.initial_guess), error('Dir "initial_guess_functions" not found: %s', paths.initial_guess); end
        addpath(genpath(paths.hmm_core)); addpath(genpath(paths.initial_guess));
        fprintf('INFO: Added HMM function paths.\n');
        if ~exist(paths.output_base, 'dir'), mkdir(paths.output_base); end
        path_to_adapted_core_script = fullfile(paths.hmm_core, 'h2mm_wexac_ver1_local.m');
        if ~exist(path_to_adapted_core_script, 'file'), error('Core script h2mm_wexac_ver1_local.m not found in %s.', paths.hmm_core); end

        %% 1. User Inputs for Nstates Scan
        fprintf('\n--- 1. Gathering Scan Parameters ---\n');
        [file_name_data, path_name_data] = uigetfile('*.mat', 'Select SSSdata.mat', pwd);
        if isequal(file_name_data,0), error('File selection cancelled.'); end
        full_data_file_path = fullfile(path_name_data, file_name_data);
        fprintf('Loading data from: %s\n', full_data_file_path);
        try data_loaded = load(full_data_file_path); catch ME_load, error('Failed to load data: %s. Err: %s', full_data_file_path, ME_load.message); end
        
        Sdata_for_analysis={}; Sdata3_for_analysis={}; burst_data_for_analysis={};
        if isfield(data_loaded,'Sdata'),Sdata_for_analysis=extract_trajectory_cell_scan('Sdata',data_loaded.Sdata);else fprintf('INFO:"Sdata" not found.\n');end
        if isfield(data_loaded,'burst_data'),burst_data_for_analysis=extract_trajectory_cell_scan('burst_data',data_loaded.burst_data);else fprintf('INFO:"burst_data" not found.\n');end
        if isfield(data_loaded,'Sdata3'),Sdata3_for_analysis=extract_trajectory_cell_scan('Sdata3',data_loaded.Sdata3);else fprintf('INFO:"Sdata3" not found.\n');end
        time_unit_sec_from_file=[]; if isfield(data_loaded,'time_unit'),time_unit_sec_from_file=data_loaded.time_unit;end
        default_deltaT_ns=100; if ~isempty(time_unit_sec_from_file)&&isnumeric(time_unit_sec_from_file)&&isscalar(time_unit_sec_from_file)&&time_unit_sec_from_file>0,default_deltaT_ns=round(time_unit_sec_from_file*1e9);end
        
        dlg_title_scan_params = 'Core HMM Parameters for Nstates Scan';
        prompt_scan_params = {'#InitialGuesses/Nstate:','Resol.(ns):','AnalysisMode(1/2):','GuessesParallel(1/0):','MaxIterHMM:','MaxRuntimeHMM(hr):','MaxTraj.:','NstatesScan(e.g.2:6):','OBSmatUpdate(1/0):'};
        dims_scan_params = [1 60];
        definput_scan_params = {'10',num2str(default_deltaT_ns),'1','0','25000','24','10000','2:6','1'}; % MaxIter 25k
        answer_scan_params = inputdlg(prompt_scan_params,dlg_title_scan_params,dims_scan_params,definput_scan_params);
        if isempty(answer_scan_params), error('Scan parameters input cancelled.'); end

        try
            num_initial_guesses_per_N=str2double(answer_scan_params{1}); dt_instrument_ns=str2double(answer_scan_params{2}); analysis_mode_An=str2double(answer_scan_params{3}); run_guesses_parallel_flag=str2double(answer_scan_params{4}); max_iter_hmm=str2double(answer_scan_params{5}); max_runtime_hmm_hours=str2double(answer_scan_params{6}); MAX_TRAJECTORIES=str2double(answer_scan_params{7}); Nstates_str=answer_scan_params{8}; obs_learn_flag_for_k=str2double(answer_scan_params{9});
            if any(isnan([num_initial_guesses_per_N,dt_instrument_ns,analysis_mode_An,run_guesses_parallel_flag,max_iter_hmm,max_runtime_hmm_hours,MAX_TRAJECTORIES,obs_learn_flag_for_k])),error('Numeric input error.');end
            if num_initial_guesses_per_N<=0||dt_instrument_ns<=0||max_iter_hmm<=0||max_runtime_hmm_hours<=0||MAX_TRAJECTORIES<=0,error('Counts,resol,iter,runtime,maxtraj must be >0.');end
            if analysis_mode_An~=1&&analysis_mode_An~=2,error('AnalysisMode must be 1/2.');end; if run_guesses_parallel_flag~=0&&run_guesses_parallel_flag~=1,error('ParallelFlag must be 0/1.');end; if obs_learn_flag_for_k~=0&&obs_learn_flag_for_k~=1,error('OBSUpdateFlag must be 0/1.');end
            if contains(Nstates_str,':'),prts=strsplit(Nstates_str,':');if length(prts)~=2,error('Nstates range err.');end;Nstates_values=str2double(prts{1}):str2double(prts{2});else Nstates_values=str2num(Nstates_str);end %#ok<ST2NM>
            if isempty(Nstates_values)||any(isnan(Nstates_values))||any(Nstates_values<1)||any(mod(Nstates_values,1)~=0),error('Invalid Nstates_values string.');end
        catch ME_prs, error('Error parsing scan params: %s',ME_prs.message); end
        fprintf('Will scan Nstates: %s\n', mat2str(Nstates_values));
        Guessmodel_options=[{'PhotonByPhoton(NoRestr.)'},{'Chain'},{'3colorsChain'},{'other(special)'}]; [indx_guessmodel,tf_gm]=listdlg('ListString',Guessmodel_options,'SelectionMode','single','Name','HMM Topology','ListSize',[200 100]);
        if ~tf_gm,error('Topology selection cancelled.');end; guessName_scan=Guessmodel_options{indx_guessmodel}; fprintf('Selected HMM Topology: %s\n',guessName_scan);
        dt_analysis_sec=dt_instrument_ns*1e-9;
        
        ts_scan_folder = timestamp_for_diary; 
        main_scan_output_folder_name = sprintf('NScan_%s_An%d_%s', strrep(guessName_scan, ' ', '_'), analysis_mode_An, ts_scan_folder);
        main_scan_output_path = fullfile(paths.output_base, main_scan_output_folder_name);
        if ~exist(main_scan_output_path,'dir'),mkdir(main_scan_output_path);else fprintf('INFO: Output folder %s already exists.\n',main_scan_output_path);end
        fprintf('Full Nstates scan results will be saved in: %s\n', main_scan_output_path);

        %% 2. Data Preprocessing
        fprintf('\n--- 2. Preprocessing Data ---\n');
        data_src_scan={}; if indx_guessmodel==3,if ~isempty(Sdata3_for_analysis)&&iscell(Sdata3_for_analysis),data_src_scan=Sdata3_for_analysis;fprintf('Using Sdata3 (%d traj)\n',length(data_src_scan));else error('3-color Sdata3 invalid.');end; else if ~isempty(burst_data_for_analysis)&&iscell(burst_data_for_analysis),data_src_scan=burst_data_for_analysis;fprintf('Using burst_data (%d traj)\n',length(data_src_scan));elseif ~isempty(Sdata_for_analysis)&&iscell(Sdata_for_analysis),data_src_scan=Sdata_for_analysis;fprintf('Using Sdata (%d traj)\n',length(data_src_scan));else error('No suitable data source.');end;end
        if length(data_src_scan)>MAX_TRAJECTORIES,fprintf('Truncating data to %d\n',MAX_TRAJECTORIES);data_src_scan=data_src_scan(1:MAX_TRAJECTORIES);end
        if isempty(data_src_scan),error('Data source empty post-trunc.');end
        data_re_for_hmm=cell(size(data_src_scan));AbsBurstInitialTime=NaN(length(data_src_scan),1);N_photons_total=0;
        for i_trj=1:length(data_src_scan),curr_trj_d=data_src_scan{i_trj};if isempty(curr_trj_d),data_re_for_hmm{i_trj}=zeros(0,2);continue;end;if ~isnumeric(curr_trj_d)||ndims(curr_trj_d)>2||size(curr_trj_d,2)<2,warning('PREP: Bad traj %d',i_trj);data_re_for_hmm{i_trj}=zeros(0,2);continue;end;pt=curr_trj_d(:,1);pc=curr_trj_d(:,2);[st,si]=sort(pt);sc=pc(si);if~isempty(st),AbsBurstInitialTime(i_trj)=st(1);rt=st-st(1);data_re_for_hmm{i_trj}=[rt,sc];N_photons_total=N_photons_total+size(rt,1);else data_re_for_hmm{i_trj}=zeros(0,2);end;end
        vidx=cellfun(@(x)ismatrix(x)&&((isempty(x)&&size(x,2)==2)||(~isempty(x)&&size(x,1)>0)),data_re_for_hmm);data_re_for_hmm=data_re_for_hmm(vidx);AbsBurstInitialTime=AbsBurstInitialTime(vidx);
        if isempty(data_re_for_hmm),error('No valid traj after preproc.');end;fprintf('INFO:Preproc done. %d valid traj,%d photons.\n',length(data_re_for_hmm),N_photons_total);
        
        data_re_save_filepath = fullfile(main_scan_output_path, 'data_re_for_hmm.mat');
        try
            save(data_re_save_filepath, 'data_re_for_hmm', 'AbsBurstInitialTime', 'N_photons_total', '-v7.3');
            fprintf('INFO: Saved processed trajectory data to: %s\n', data_re_save_filepath);
        catch ME_save_data_re, warning('Failed to save data_re_for_hmm.mat. Error: %s', getReport(ME_save_data_re, 'basic')); end
        
        scan_config_params=struct('Nstates_values_scanned',Nstates_values,'num_initial_guesses_per_N',num_initial_guesses_per_N,'dt_analysis_sec',dt_analysis_sec,'analysis_mode_An',analysis_mode_An,'run_guesses_parallel_flag',run_guesses_parallel_flag,'max_iter_hmm',max_iter_hmm,'max_runtime_hmm_hours',max_runtime_hmm_hours,'MAX_TRAJECTORIES_used',MAX_TRAJECTORIES,'indx_guessmodel',indx_guessmodel,'guessName_scan',guessName_scan,'obs_learn_flag_for_k',obs_learn_flag_for_k,'full_data_file_path',full_data_file_path,'project_base_path',project_base_path,'N_photons_total',N_photons_total,'timestamp_scan_start',ts_scan_folder);
        overall_cfg_fp=fullfile(main_scan_output_path,'NScan_Overall_Config.mat');
        try save(overall_cfg_fp,'scan_config_params'); fprintf('INFO:Saved NScan_Overall_Config.mat to %s\n',main_scan_output_path);
        catch ME_scfg2, warning('Failed overall_cfg save: %s',getReport(ME_scfg2,'basic')); end

        %% 3. Loop Over Nstates_values
        fprintf('\n--- 3. Starting Nstates Scan Loop ---\n');
        results_summary = table('Size',[length(Nstates_values),6],'VariableTypes',repmat({'double'},1,6),'VariableNames',{'Nstates','Best_LL','Num_Params_k','AIC','BIC','Mean_Duration_Per_Guess_sec'});
        all_best_models_per_Nstates = cell(1,length(Nstates_values));
        
        fprintf('DEBUG: --- About to start the i_scan loop. Number of Nstates to scan: %d ---\n', length(Nstates_values));
        if ~exist('Nstates_values','var')||isempty(Nstates_values), error('CRITICAL: Nstates_values not defined or empty before main loop!'); end
        fprintf('DEBUG: Nstates_values confirmed right before i_scan loop: %s, length = %d\n', mat2str(Nstates_values), length(Nstates_values));

        for i_scan = 1:length(Nstates_values)
            fprintf('DEBUG: --- INSIDE i_scan LOOP: Scan Iteration %d of %d ---\n', i_scan, length(Nstates_values));
            current_Nstates = Nstates_values(i_scan);
            fprintf('\n========================================================\n');
            fprintf('INFO: Starting HMM analysis for Nstates = %d\n', current_Nstates);
            fprintf('========================================================\n');

            nstate_run_folder_name=sprintf('Nstates_%02d_Run',current_Nstates); 
            cur_n_out_p=fullfile(main_scan_output_path,nstate_run_folder_name); 
            cur_n_tmp_p=fullfile(cur_n_out_p,'temp_individual_runs');
            if ~exist(cur_n_out_p,'dir'),mkdir(cur_n_out_p);end; 
            if ~exist(cur_n_tmp_p,'dir'),mkdir(cur_n_tmp_p);end
            fprintf('Output for Nstates=%d in: %s\n',current_Nstates,cur_n_out_p);

            cur_fixS=[]; 
            if analysis_mode_An==2&&indx_guessmodel==4 
                if current_Nstates>1
                    ans_fs=inputdlg(sprintf('N%d:fixS(1-%d):',current_Nstates,current_Nstates),'FixS for Nstate',[1 40],{num2str(min(2,current_Nstates))});
                    if isempty(ans_fs)||isempty(str2double(ans_fs{1}))||isnan(str2double(ans_fs{1}))
                        disp('fixS input cancelled. Skipping this Nstates run.'); 
                        results_summary.Nstates(i_scan)=current_Nstates; results_summary.Best_LL(i_scan)=NaN; all_best_models_per_Nstates{i_scan}=[]; continue; 
                    end
                    cur_fixS=str2double(ans_fs{1});
                    if cur_fixS<1||cur_fixS>current_Nstates
                        warning('Bad fixS val. Skipping Nstates=%d.',current_Nstates); results_summary.Nstates(i_scan)=current_Nstates; results_summary.Best_LL(i_scan)=NaN; all_best_models_per_Nstates{i_scan}=[]; continue; 
                    end
                else 
                    cur_fixS=[]; fprintf('INFO: Nstates=1, no fixS for "other" topology.\n'); 
                end
            end
            
            IG_cur_N_orig={}; 
            try 
                switch indx_guessmodel
                    case 1, IG_cur_N_orig=InitialGuess_PhotonByPhoton(num_initial_guesses_per_N,current_Nstates);
                    case 2, IG_cur_N_orig=InitialGuess_chain(num_initial_guesses_per_N,current_Nstates);
                    case 3, IG_cur_N_orig=InitialGuess_3colors_chain(num_initial_guesses_per_N,current_Nstates);
                    case 4, IG_cur_N_orig=InitialGuess_special(num_initial_guesses_per_N,current_Nstates,dt_analysis_sec,cur_fixS);
                    otherwise, error('Internal error: Invalid indx_guessmodel during loop.');
                end
            catch ME_ig_gen
                 warning('Failed to generate initial guesses for Nstates=%d, Topology=%s. Error: %s. Skipping this Nstates.', current_Nstates, guessName_scan, ME_ig_gen.message);
                 results_summary.Nstates(i_scan)=current_Nstates; results_summary.Best_LL(i_scan)=NaN; all_best_models_per_Nstates{i_scan}=[]; continue;
            end

            if isempty(IG_cur_N_orig) || ~iscell(IG_cur_N_orig) || all(cellfun(@isempty, IG_cur_N_orig))
                warning('Nstates=%d: No valid initial guesses were generated. Skipping this Nstates.', current_Nstates);
                results_summary.Nstates(i_scan)=current_Nstates; results_summary.Best_LL(i_scan)=NaN; all_best_models_per_Nstates{i_scan}=[];
                continue;
            end
            IG_cur_N_to_use = IG_cur_N_orig(~cellfun(@isempty, IG_cur_N_orig)); 
            effective_num_guesses_this_N = length(IG_cur_N_to_use);

            if effective_num_guesses_this_N == 0
                 warning('Nstates=%d: No valid (non-empty) initial guesses after filtering. Skipping this Nstates.', current_Nstates);
                 results_summary.Nstates(i_scan)=current_Nstates; results_summary.Best_LL(i_scan)=NaN; all_best_models_per_Nstates{i_scan}=[];
                 continue;
            end
            if effective_num_guesses_this_N < num_initial_guesses_per_N
                fprintf('WARNING: Nstates=%d: Using %d valid initial guesses (expected %d).\n', current_Nstates, effective_num_guesses_this_N, num_initial_guesses_per_N);
            end

            if current_Nstates>0 
                debug_log_filename = fullfile(cur_n_out_p, sprintf('InitialGuesses_N%02d_DebugLog.txt', current_Nstates));
                fid_debug = fopen(debug_log_filename, 'w');
                if fid_debug == -1, warning('Could not open debug log for N%d: %s.',current_Nstates, debug_log_filename); fid_debug = 1; end
                fprintf(fid_debug, '--- DEBUG: Initial Guesses (Original Fields from Generator) for Nstates = %d ---\n', current_Nstates);
                fprintf(fid_debug, 'Timestamp: %s\nTopology: %s\nNumber of effective guesses logged (up to 3): %d\n\n', datestr(now), guessName_scan, min(effective_num_guesses_this_N,3) );
                num_guesses_to_print_debug = min(effective_num_guesses_this_N, 3);
                if isempty(IG_cur_N_to_use), fprintf(fid_debug, 'WARNING: No valid original guesses to log!\n');
                else
                    for i_ig_debug = 1:num_guesses_to_print_debug
                        fprintf(fid_debug, 'Original Guess Struct %d (for Nstates=%d):\n', i_ig_debug, current_Nstates);
                        co_g_d = IG_cur_N_to_use{i_ig_debug}; 
                        if isfield(co_g_d,'prior0'), fprintf(fid_debug,'prior0:\n%s\n',mat2str(co_g_d.prior0',4)); if abs(sum(co_g_d.prior0)-1)>1e-6, fprintf(fid_debug,' WARN prior0 sum=%f\n',sum(co_g_d.prior0));end; else fprintf(fid_debug,' WARN: Field prior0 missing from struct!\n');end
                        if isfield(co_g_d,'transmat0'), fprintf(fid_debug,'transmat0:\n'); for r_tr_dbg=1:size(co_g_d.transmat0,1), fprintf(fid_debug,'  %s\n',mat2str(co_g_d.transmat0(r_tr_dbg,:),4)); end; rs_tr_dbg=sum(co_g_d.transmat0,2); if any(abs(rs_tr_dbg-1)>1e-6),fprintf(fid_debug,' WARN TR0 sums: %s\n',mat2str(rs_tr_dbg',4));end; else fprintf(fid_debug,' WARN: Field transmat0 missing from struct!\n');end
                        if isfield(co_g_d,'obsmat0'), fprintf(fid_debug,'obsmat0:\n'); for r_e_dbg=1:size(co_g_d.obsmat0,1), fprintf(fid_debug,'  %s\n',mat2str(co_g_d.obsmat0(r_e_dbg,:),4)); end; rs_e_dbg=sum(co_g_d.obsmat0,2); if any(abs(rs_e_dbg-1)>1e-6),fprintf(fid_debug,' WARN E0 sums: %s\n',mat2str(rs_e_dbg',4));end; else fprintf(fid_debug,' WARN: Field obsmat0 missing from struct!\n');end
                        fprintf(fid_debug, '--------------------------------------\n\n');
                    end
                end
                if fid_debug ~= 1, fclose(fid_debug); fprintf('Initial guess debug log (original fields) saved to: %s\n', debug_log_filename);
                else, fprintf(1,'--- END OF DEBUG OUTPUT TO COMMAND WINDOW FOR Nstates = %d ---\n', current_Nstates); end
            end
            
            iter_rcfg_cN=cell(1,effective_num_guesses_this_N);
            for ic=1:effective_num_guesses_this_N
                rcfg=struct(); rcfg.temp_dir_base=cur_n_tmp_p; rcfg.script_to_run_path=path_to_adapted_core_script; rcfg.current_guess_index=ic;
                rcfg.all_initial_guesses=IG_cur_N_to_use; 
                rcfg.trajectory_data=data_re_for_hmm; rcfg.max_iter=max_iter_hmm;
                rcfg.obs_fix_flag=obs_learn_flag_for_k; rcfg.analysis_mode_An=analysis_mode_An; rcfg.MaxRunTime_cluster_hours=max_runtime_hmm_hours;
                if analysis_mode_An==2 && ~isempty(cur_fixS), rcfg.fixS_val=cur_fixS; end
                iter_rcfg_cN{ic}=rcfg;
            end

            % Variables for HMM execution results for current Nstate
            all_hmm_results_this_N = cell(1,effective_num_guesses_this_N); 
            all_guess_durations_this_N = NaN(1,effective_num_guesses_this_N);
            overall_hmm_exec_start_time_this_N = tic; 
            
            pool_shutdown_needed_this_N = false; 
            effective_parallel_flag_this_N = run_guesses_parallel_flag;

            if run_guesses_parallel_flag == 1 && ~isempty(ver('parallel'))
                cp_current = gcp('nocreate');
                if isempty(cp_current)
                    fprintf('Nstates=%d: Attempting to start parallel pool for guesses...\n', current_Nstates);
                    try parpool; pool_shutdown_needed_this_N=true; fprintf('Nstates=%d: Parallel pool started.\n', current_Nstates);
                    catch ME_parpool_currentN
                        warning('NStates=%d: Parpool startup failed: %s. Running guesses serially.', current_Nstates, getReport(ME_parpool_currentN,'basic'));
                        effective_parallel_flag_this_N=0;
                    end
                else, fprintf('Nstates=%d: Using existing parallel pool for guesses.\n', current_Nstates); end
            elseif run_guesses_parallel_flag == 1 && isempty(ver('parallel'))
                fprintf('Nstates=%d: Parallel exec requested, but Parallel Toolbox not found. Running serially.\n', current_Nstates);
                effective_parallel_flag_this_N=0;
            end

            for i_cfg_update_currentN = 1:effective_num_guesses_this_N
                if effective_parallel_flag_this_N == 1
                    iter_rcfg_cN{i_cfg_update_currentN}.control_parpool_externally=true;
                    iter_rcfg_cN{i_cfg_update_currentN}.par_hmm_internal_flag=0; 
                else 
                    iter_rcfg_cN{i_cfg_update_currentN}.control_parpool_externally=false;
                    iter_rcfg_cN{i_cfg_update_currentN}.par_hmm_internal_flag=0; 
                end
            end

            if effective_parallel_flag_this_N == 1
                parfor igp_currentN =1:effective_num_guesses_this_N
                    worker_id_currentN=''; try tsk_curN=getCurrentTask();if~isempty(tsk_curN),worker_id_currentN=sprintf('(W%d)',tsk_curN.ID);end;catch;end
                    fprintf('N%d,G%d/%d %s\n',current_Nstates,igp_currentN,effective_num_guesses_this_N,worker_id_currentN);
                    cfg_p_loop_currentN = iter_rcfg_cN{igp_currentN};
                    try [rp_loop_cN,~,dp_loop_cN]=local_hmm_runner(cfg_p_loop_currentN); 
                        all_hmm_results_this_N{igp_currentN}=rp_loop_cN; 
                        all_guess_durations_this_N(igp_currentN)=dp_loop_cN;
                    catch MEplr_cN
                        fprintf('ERR P N%d G%d %s:%s\n',current_Nstates,igp_currentN,worker_id_currentN,getReport(MEplr_cN,'basic'));
                        all_hmm_results_this_N{igp_currentN}=struct();
                    end
                end
            else 
                for igs_currentN=1:effective_num_guesses_this_N
                    fprintf('N%d,G%d/%d(S)\n',current_Nstates,igs_currentN,effective_num_guesses_this_N);
                    cfg_s_loop_currentN = iter_rcfg_cN{igs_currentN};
                    try [rs_loop_cN,~,ds_loop_cN]=local_hmm_runner(cfg_s_loop_currentN);
                        all_hmm_results_this_N{igs_currentN}=rs_loop_cN; 
                        all_guess_durations_this_N(igs_currentN)=ds_loop_cN;
                    catch MEslr_cN
                        fprintf('ERR S N%d G%d:%s\n',current_Nstates,igs_currentN,getReport(MEslr_cN,'basic'));
                        all_hmm_results_this_N{igs_currentN}=struct();
                    end
                end
            end
            overall_duration_this_N_val=toc(overall_hmm_exec_start_time_this_N);
            fprintf('N%d:TotalHMMtime %dG:%.2fs\n',current_Nstates,effective_num_guesses_this_N,overall_duration_this_N_val);
            if pool_shutdown_needed_this_N && ~isempty(gcp('nocreate')), try delete(gcp('nocreate'));fprintf('N%d:ParpoolShut\n',current_Nstates);catch ME_dps_cN,warning('N%d:ErrShutParpool:%s',current_Nstates,getReport(ME_dps_cN,'basic'));end;end

            % --- Process results for current_Nstates ---
            final_LLs_this_N = NaN(1,effective_num_guesses_this_N); 
            for irll_cN=1:effective_num_guesses_this_N
                res_s_cN = all_hmm_results_this_N{irll_cN}; 
                if~isempty(res_s_cN)&&isstruct(res_s_cN)&&isfield(res_s_cN,'LL')&&~isempty(res_s_cN.LL)&&isnumeric(res_s_cN.LL)&&isfinite(res_s_cN.LL(end))
                    final_LLs_this_N(irll_cN)=res_s_cN.LL(end);
                end
            end
            best_LL_this_N=-Inf; best_model_this_N=struct(); best_idx_this_N=NaN; 
            if any(isfinite(final_LLs_this_N))
                [best_LL_this_N,best_idx_this_N]=max(final_LLs_this_N);
                if best_idx_this_N>0 && best_idx_this_N<=length(all_hmm_results_this_N) && ~isempty(all_hmm_results_this_N{best_idx_this_N}) && isstruct(all_hmm_results_this_N{best_idx_this_N})
                    best_model_this_N=all_hmm_results_this_N{best_idx_this_N};
                else 
                    best_model_this_N=struct(); best_LL_this_N=NaN; fprintf('W:N%d:BestLL found,model invalid\n',current_Nstates);
                end
                all_best_models_per_Nstates{i_scan}=best_model_this_N;
                
                % Define variables to be saved for this Nstate's results file
                % These names MUST match what processHMM_results_local.m expects to load
                cur_N_all_hmm_res = all_hmm_results_this_N;       % Corrected variable name for saving
                cur_N_all_LLs     = final_LLs_this_N;             % Corrected variable name for saving
                % IG_cur_N_orig is already correctly named from InitialGuess call
                iter_rcfg_cur_N   = iter_rcfg_cN;                 % Use the iterated config name
                cur_N_all_g_dur_s = all_guess_durations_this_N;   % Corrected variable name
                cur_N_ovr_hmm_dur = overall_duration_this_N_val;  % Corrected variable name
                cfg_n_sv_val      = struct('N',current_Nstates,'fixS',cur_fixS,'effective_num_guesses',effective_num_guesses_this_N); % Use consistent name

                n_res_fn_currentN=fullfile(cur_n_out_p,sprintf('All_HMM_Results_N%02d.mat',current_Nstates));
                try 
                    save(n_res_fn_currentN, ...
                        'cur_N_all_hmm_res', 'cur_N_all_LLs', ...
                        'IG_cur_N_orig', 'iter_rcfg_cur_N', ...
                        'cur_N_all_g_dur_s', 'cur_N_ovr_hmm_dur', ...
                        'cfg_n_sv_val', '-v7.3');
                    fprintf('N%d:%dG res saved:%s\n',current_Nstates,effective_num_guesses_this_N,n_res_fn_currentN);
                catch ME_snr_cN
                    warning('N%d:FailSaveAllResN%02d.Err:%s',current_Nstates,current_Nstates,getReport(ME_snr_cN,'basic'));
                end
            else 
                fprintf('W:N%d:No valid HMM res (all LLs NaN/Inf)\n',current_Nstates);
                all_best_models_per_Nstates{i_scan}=[]; best_LL_this_N=NaN;
            end

            n_obs_sym_final=2;
            try if~isempty(best_model_this_N)&&isstruct(best_model_this_N)&&isfield(best_model_this_N,'obsmat0')&&~isempty(best_model_this_N.obsmat0),n_obs_sym_final=size(best_model_this_N.obsmat0,2);elseif~isempty(IG_cur_N_to_use)&&~isempty(IG_cur_N_to_use{1})&&isfield(IG_cur_N_to_use{1},'obsmat0')&&~isempty(IG_cur_N_to_use{1}.obsmat0),n_obs_sym_final=size(IG_cur_N_to_use{1}.obsmat0,2);end;catch ME_nos_final,warning('NoGet n_obs_sym N%d.Def=2.Err:%s',current_Nstates,getReport(ME_nos_final,'basic'));n_obs_sym_final=2;end
            k_T_final=0;if current_Nstates>=1,k_P_f=current_Nstates-1;if current_Nstates==1,k_P_f=0;end;k_TR_f=current_Nstates*(current_Nstates-1);if current_Nstates==1,k_TR_f=0;end;k_E_f=0;if obs_learn_flag_for_k==1,k_E_f=current_Nstates*(n_obs_sym_final-1);end;k_T_final=k_P_f+k_TR_f+k_E_f;end
            AIC_final=NaN;BIC_final=NaN;if isfinite(best_LL_this_N)&&N_photons_total>0&&k_T_final>=0,AIC_final=-2*best_LL_this_N+2*k_T_final;BIC_final=-2*best_LL_this_N+k_T_final*log(N_photons_total);end
            results_summary.Nstates(i_scan)=current_Nstates;results_summary.Best_LL(i_scan)=best_LL_this_N;results_summary.Num_Params_k(i_scan)=k_T_final;results_summary.AIC(i_scan)=AIC_final;results_summary.BIC(i_scan)=BIC_final;results_summary.Mean_Duration_Per_Guess_sec(i_scan)=nanmean(all_guess_durations_this_N);
            fprintf('N%d:BestLL=%.2f,k=%d,AIC=%.2f,BIC=%.2f,MeanDur=%.1fs\n',current_Nstates,best_LL_this_N,k_T_final,AIC_final,BIC_final,nanmean(all_guess_durations_this_N));
            fprintf('DEBUG: --- END OF i_scan LOOP BODY for Nstates = %d ---\n', current_Nstates);
        end % End loop over Nstates_values
        fprintf('DEBUG: --- FINISHED i_scan LOOP ---\n');
        
        %% 4. Finalize and Plot Nstates Scan Results
        fprintf('\n--- 4. Finalizing Nstates Scan ---\n');
        disp(results_summary);
        summary_fp=fullfile(main_scan_output_path,'NScan_Results_Summary.mat');try save(summary_fp,'results_summary','all_best_models_per_Nstates','scan_config_params');writetable(results_summary,fullfile(main_scan_output_path,'NScan_Results_Summary.csv'));fprintf('INFO:Scan summary saved to NScan_Results_Summary files in %s\n',main_scan_output_path);catch ME_ssum,warning('Fail save summary .mat/csv.Err:%s',getReport(ME_ssum,'basic'));end
        try fig_ab=figure('Name','NScan:AIC&BIC','Pos',[100,100,900,600],'Vis','off');yyaxis left;plot(results_summary.Nstates,results_summary.AIC,'bo-','LineWidth',1.5,'MarkerSize',8,'Disp','AIC');ylabel('AIC');hold on;[min_a,idx_a]=min(results_summary.AIC);min_aN=NaN;if~isempty(min_a)&&isfinite(min_a)&&idx_a>0&&idx_a<=length(results_summary.Nstates),min_aN=results_summary.Nstates(idx_a);plot(min_aN,min_a,'bp','MS',12,'MarkerFaceColor','b','Disp',sprintf('MinAIC(N=%d)',min_aN));text(min_aN,min_a,sprintf(' MinAIC(N=%d)',min_aN),'Vert','bot','Horiz','right');end;yyaxis right;plot(results_summary.Nstates,results_summary.BIC,'ro-','LineWidth',1.5,'MarkerSize',8,'Disp','BIC');ylabel('BIC');[min_b,idx_b]=min(results_summary.BIC);min_bN=NaN;if~isempty(min_b)&&isfinite(min_b)&&idx_b>0&&idx_b<=length(results_summary.Nstates),min_bN=results_summary.Nstates(idx_b);plot(min_bN,min_b,'rp','MS',12,'MarkerFaceColor','r','Disp',sprintf('MinBIC(N=%d)',min_bN));text(min_bN,min_b,sprintf(' MinBIC(N=%d)',min_bN),'Vert','top','Horiz','right');end;hold off;xlabel('Nstates');title(sprintf('AIC&BIC (Data:%s)',strrep(file_name_data,'_','\_')));legh=findobj(gca,'-regexp','DisplayName','[^'']');if~isempty(legh),legend(legh,'Loc','best','Interp','none');end;grid on;pfp_f=fullfile(main_scan_output_path,'NScan_AIC_BIC_plot.fig');pfp_p=fullfile(main_scan_output_path,'NScan_AIC_BIC_plot.png');saveas(fig_ab,pfp_f);saveas(fig_ab,pfp_p);if exist('fig_ab','var')&&isvalid(fig_ab)&&isa(fig_ab,'matlab.ui.Figure'),close(fig_ab);end;fprintf('INFO:AIC/BIC plot saved to %s\n',main_scan_output_path);catch ME_plt,warning('Fail gen/save AIC/BIC plot.Err:%s',getReport(ME_plt,'basic'));end
        if exist('min_aN','var')&&isfinite(min_aN),fprintf('NScan done.OptN_AIC:%d.\n',min_aN);end;if exist('min_bN','var')&&isfinite(min_bN),fprintf('NScan done.OptN_BIC:%d.\n',min_bN);end;fprintf('All scan res/plots saved in: %s\n',main_scan_output_path);
        fprintf('\n--- Nstates Scan Pipeline Finished Successfully ---\n');

    catch ME_nscan_main 
        fprintf(2,'\n--- ERROR in Nstates Scan Pipeline ---\n');fprintf(2,'ErrID:%s\nErrMsg:%s\n',ME_nscan_main.identifier,ME_nscan_main.message);rpt_str=getReport(ME_nscan_main,'extended','hyperlinks','off');fprintf(2,'StackTrace:\n%s\n',rpt_str);
        if exist('diary_filename_final','var')&&~isempty(diary_filename_final)&&strcmp(get(0,'Diary'),'on') 
            fprintf('DIARY_ERROR:Err in NScanPipe.\nID:%s\nMsg:%s\n',ME_nscan_main.identifier,ME_nscan_main.message);rpt_str_d=strrep(rpt_str,'%','%%');fprintf('DIARY_ERR:Stack:\n%s\n',rpt_str_d);
        end
        fprintf(2,'\n--- Nstates Scan Pipeline Terminated Due to Error ---\n');
    end
    
    if exist('cleanupObj_diary','var') && ~isempty(cleanupObj_diary) && isvalid(cleanupObj_diary)
        delete(cleanupObj_diary); 
    elseif strcmp(get(0,'Diary'),'on') 
        diary('off');
        fprintf('INFO: Diary logging stopped (fallback method).\n');
    end
    if exist('diary_filename_final', 'var') && ~isempty(diary_filename_final)
        fprintf('INFO: Diary file location: %s\n', diary_filename_final);
    end
end

% --- Helper function extract_trajectory_cell_scan ---
function final_traj_cell = extract_trajectory_cell_scan(input_var_name, loaded_var)
    final_traj_cell = {}; 
    if ~iscell(loaded_var), if ~isempty(loaded_var), fprintf('INFO (extract_helper): "%s" not cell.\n', input_var_name); end; return; end
    if isempty(loaded_var), return; end
    temp_data = loaded_var; nest_lvl = 0; max_nest = 3;
    while iscell(temp_data) && numel(temp_data)==1 && iscell(temp_data{1}) && nest_lvl < max_nest, temp_data=temp_data{1}; nest_lvl=nest_lvl+1; end
    if iscell(temp_data) && ~isempty(temp_data)
        peek = temp_data{1};
        if isnumeric(peek) && (isempty(peek) || (ismatrix(peek) && size(peek,2)>=2))
            processed_trajectories = cell(size(temp_data)); valid_c = 0;
            for k=1:numel(temp_data)
                traj_d = temp_data{k};
                if isnumeric(traj_d) && (isempty(traj_d) || (ismatrix(traj_d) && size(traj_d,2)>=2))
                    valid_c=valid_c+1;
                    if ~isempty(traj_d), processed_trajectories{valid_c} = traj_d(:,1:2); else, processed_trajectories{valid_c} = zeros(0,2); end
                elseif ~isempty(traj_d), warning('EXTRACT_H: El %d in "%s" not num NxM (M>=2). Skip.', k, input_var_name); end
            end
            final_traj_cell = processed_trajectories(1:valid_c);
        elseif iscell(peek)
            concatenated = {};
            for k_o=1:numel(temp_data)
                inner_c = temp_data{k_o};
                if iscell(inner_c)
                    for k_i=1:numel(inner_c)
                        traj_d_i = inner_c{k_i};
                        if isnumeric(traj_d_i) && (isempty(traj_d_i)||(ismatrix(traj_d_i)&&size(traj_d_i,2)>=2))
                            if ~isempty(traj_d_i), concatenated{end+1}=traj_d_i(:,1:2); else, concatenated{end+1}=zeros(0,2); end %#ok<AGROW>
                        elseif ~isempty(traj_d_i), warning('EXTRACT_H: InnerEl of "%s"(%d,%d) not num NxM. Skip.',input_var_name,k_o,k_i);end
                    end
                elseif ~isempty(inner_c), warning('EXTRACT_H: OuterEl %d of "%s" not cell. Skip batch.',k_o,input_var_name);end
            end
            final_traj_cell = concatenated;
        else 
             if ~isempty(temp_data) && ~isempty(peek), warning('EXTRACT_H: Format of "%s" unknown. FirstEl class: %s.',input_var_name,class(peek));end; final_traj_cell={};
        end
    else 
        if iscell(temp_data)&&isempty(temp_data),else warning('EXTRACT_H: Var "%s" not cell post-unwrap. Class:%s.',input_var_name,class(temp_data));end; final_traj_cell={};
    end
end