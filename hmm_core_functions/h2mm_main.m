% h2mm_main.m
% Core HMM learning function with data preprocessing, resume capability,
% and corrected iteration printing and nested function call.

function [LL, prior, transmat, obsmat, lastIteration, all_gamma] = ...
    h2mm_main(data_input_orig, prior_init, transmat_init, obsmat_init, max_iter_input, MaxRunTime_hours, OBSFix_flag, LastData_input, fileName_LastData, par_internal_flag)

% INPUTS:
% data_input_orig{ex}   : Cell array of trajectories. Each cell data{ex} is an Nx2 matrix:
%                         data{ex}(:,1) = photon arrival times (relative to burst start)
%                         data{ex}(:,2) = observation channel (discrete integer)
% prior_init(i)         : Initial prior probability for state i
% transmat_init(i,j)    : Initial transition probability from state i to j (for unit time step TAU)
% obsmat_init(i,o)      : Initial observation probability of symbol o given state i
% max_iter_input        : Maximum number of EM iterations
% MaxRunTime_hours      : Maximum wall clock time for this run (in hours)
% OBSFix_flag           : If 0, observation matrix is not updated. If 1, it is.
% LastData_input        : Struct containing data from a previous run for resuming.
% fileName_LastData     : Full path/name of the .mat file to save LastData for resuming.
% par_internal_flag     : If 1, use parfor for E-step over trajectories within compute_ess_dhmm.

%% Define some option variables
thresh = 1e-14;       
verbose = 1;          
obs_prior_weight = 0; 
adj_prior = 1;        
adj_trans = 1;        
adj_obs = OBSFix_flag;
backup_frequency = 100; 

%% Data Preprocessing: Inter-photon arrival times and Factorization
% This block processes the input 'data_input_orig' to create 'data_for_hmm_engine'
% with indexed inter-photon times, and also 'TotalArrivalDelta' and 'R'.

fprintf('DEBUG H2MM_MAIN: Entering corrected TotalArrivalDelta and data (formerly data_red) generation.\n');
TAU = 1; 
original_input_data_h2mm_main = data_input_orig; % Keep a clear name for original input

num_trajectories = length(original_input_data_h2mm_main);
if num_trajectories == 0
    error('H2MM_MAIN_ERROR: Input data is empty. No trajectories to process.');
end

all_inter_photon_diffs_list = {}; 
for traj_idx = 1:num_trajectories
    current_trajectory_times = original_input_data_h2mm_main{traj_idx}(:,1);
    if length(current_trajectory_times) > 1 
        inter_photon_diffs_this_traj = diff(current_trajectory_times);
        % Replace non-positive diffs with 1 (smallest unit time, assumes TAU=1 is in TotalArrivalDelta)
        inter_photon_diffs_this_traj(inter_photon_diffs_this_traj <= 0) = 1; 
        all_inter_photon_diffs_list{end+1} = inter_photon_diffs_this_traj; %#ok<AGROW>
    end
end

TotalArrivalDelta = []; % Initialize
if isempty(all_inter_photon_diffs_list)
    warning('H2MM_MAIN: No inter-photon time differences found (all trajectories might be too short). Defaulting TotalArrivalDelta to handle this.');
    TotalArrivalDelta = 1; % Default to TAU if no diffs, Factorization and downstream should handle this.
else
    all_diffs_collected_vector = vertcat(all_inter_photon_diffs_list{:});
    if isempty(all_diffs_collected_vector)
        warning('H2MM_MAIN: Collected diffs vector is empty after concatenation. Defaulting TotalArrivalDelta.');
        TotalArrivalDelta = 1;
    else
        unique_positive_diffs = unique(all_diffs_collected_vector(all_diffs_collected_vector > 0));
        if isempty(unique_positive_diffs)
            warning('H2MM_MAIN: No positive inter-photon time differences. Defaulting TotalArrivalDelta.');
            TotalArrivalDelta = 1;
        else
            TotalArrivalDelta = sort(unique_positive_diffs)'; % Sort and ensure ROW vector
        end
    end
end
fprintf('DEBUG H2MM_MAIN: Min TotalArrivalDelta = %g, Max TotalArrivalDelta = %g, NumUniqueDeltas = %d\n', ...
    min(TotalArrivalDelta), max(TotalArrivalDelta), length(TotalArrivalDelta));
if length(TotalArrivalDelta) > 0 && length(TotalArrivalDelta) < 20
    disp('DEBUG H2MM_MAIN: Few unique TotalArrivalDelta values:'); disp(TotalArrivalDelta);
end

R = Factorization(TotalArrivalDelta, TAU); % TAU=1
fprintf('DEBUG H2MM_MAIN: Factorization R computed using %d unique delta_t values.\n', length(TotalArrivalDelta));

data_processed_for_engine = cell(size(original_input_data_h2mm_main)); 
for traj_idx = 1:length(original_input_data_h2mm_main)
    current_original_trajectory = original_input_data_h2mm_main{traj_idx};
    if isempty(current_original_trajectory) || size(current_original_trajectory,1) < 1
        data_processed_for_engine{traj_idx} = zeros(0,2); 
        continue;
    end
    
    num_photons_this_traj = size(current_original_trajectory, 1);
    data_processed_for_engine{traj_idx} = zeros(num_photons_this_traj, 2); 
    data_processed_for_engine{traj_idx}(:,2) = current_original_trajectory(:,2); 
    data_processed_for_engine{traj_idx}(1,1) = 0; % First "time" is 0 (index or special value)
                                
    if num_photons_this_traj > 1
        inter_photon_times_this_traj = diff(current_original_trajectory(:,1));
        inter_photon_times_this_traj(inter_photon_times_this_traj <= 0) = 1; 
            
        for photon_event_idx = 1:length(inter_photon_times_this_traj)
            current_dt_value = inter_photon_times_this_traj(photon_event_idx);
            % Find index in master TotalArrivalDelta list
            found_master_idx_list = find(TotalArrivalDelta == current_dt_value, 1); 
            if isempty(found_master_idx_list)
                % This should ideally not happen if TotalArrivalDelta includes all unique positive diffs
                warning('H2MM_MAIN: CRITICAL - dt_value %g from traj %d not found in master TotalArrivalDelta. Using index 1. Check logic.', current_dt_value, traj_idx);
                found_master_idx_list = 1; % Fallback to first index
            end
            data_processed_for_engine{traj_idx}(photon_event_idx + 1, 1) = found_master_idx_list(1);
        end
    end
end
data = data_processed_for_engine; % 'data' now holds the indexed inter-photon times
clear original_input_data_h2mm_main data_processed_for_engine all_inter_photon_diffs_list all_diffs_collected_vector; 
fprintf('DEBUG H2MM_MAIN: Final data structure for HMM engine prepared.\n');
if ~isempty(data) && ~isempty(data{1}) && size(data{1},1) > 0
    fprintf('DEBUG H2MM_MAIN: Example final data{1}(1:min(5,end),1) (should be indices or 0):\n');
    disp(data{1}(1:min(5,size(data{1},1)),1));
else
    fprintf('DEBUG H2MM_MAIN: data{1} is empty or too short to display example.\n');
end
% --- END OF DATA PREPROCESSING BLOCK ---


%% Load parameters from LastData_input or use initial parameters
LastData = LastData_input; 
max_iter = max_iter_input; % Use renamed input from function signature

if LastData.Iteration == 0
    prior = prior_init; 
    transmat = transmat_init; 
    obsmat = obsmat_init;
    
    previous_loglik = -inf; 
    num_iter = 1;
    LL = -Inf(1, max_iter); 
    iteration_times = zeros(1, max_iter); % Corrected spelling from interation_times
    
    if ~isfield(LastData, 'PriorList'), LastData.PriorList = []; end
    if ~isfield(LastData, 'transmatList'), LastData.transmatList = []; end % Was TransmatList
    if ~isfield(LastData, 'ObsList'), LastData.ObsList = []; end
    if ~isfield(LastData, 'InteruptCount'), LastData.InteruptCount = 0; end
    
    LastData.LL = LL; 
    LastData.interation_times = iteration_times; % Corrected field name
else 
    prior = LastData.last_prior; 
    transmat = LastData.last_transmat; 
    obsmat = LastData.last_obsmat;
    num_iter = LastData.Iteration + 1; 
    previous_loglik = LastData.LL(LastData.Iteration);
    LL = LastData.LL; 
    iteration_times = LastData.interation_times; % Corrected field name
    
    if length(LL) < max_iter, LL(length(LL)+1 : max_iter) = -Inf; end
    if length(iteration_times) < max_iter, iteration_times(length(iteration_times)+1 : max_iter) = 0; end
    LastData.InteruptCount = LastData.InteruptCount + 1;
end
converged = 0;
first_iter_of_this_segment = true; % For refined printing

%% Main EM Loop
fprintf('H2MM_MAIN: Starting EM iterations (Max %d). Current iteration: %d.\n', max_iter, num_iter);
overall_start_time = tic; 

% Ensure data is cell, though data_processed_for_engine should already be
if ~iscell(data)
    data = num2cell(data, 2); 
end

while (num_iter <= max_iter) && ~converged
    iter_start_time_abs_loop = toc(overall_start_time); % Time at start of this loop iteration
    
    % E step
    % The nested compute_ess_dhmm will calculate transmat_t and Rho internally
    % using the current 'transmat', 'TotalArrivalDelta', and 'R'.
    [loglik, exp_num_trans, exp_num_visits1, exp_num_emit, ~, all_gamma_current_iter] = ...
        compute_ess_dhmm(prior, transmat, obsmat, data, obs_prior_weight, TotalArrivalDelta, R, par_internal_flag);
        % exp_num_visitsT output from compute_ess_dhmm is ignored with ~

    % M step
    if adj_prior
        if sum(exp_num_visits1(:)) > 1e-9 % Check sum before division
            prior = exp_num_visits1 / sum(exp_num_visits1);
        else
            warning('H2MM_MAIN: Sum of exp_num_visits1 near zero in iter %d. Prior not updated.', num_iter);
        end
        if ~isfield(LastData,'PriorList') || isempty(LastData.PriorList) LastData.PriorList=prior(:)'; else LastData.PriorList = [LastData.PriorList; prior(:)']; end
    end
    
    if adj_trans && ~isempty(exp_num_trans)
        transmat = mk_stochastic(exp_num_trans);
        if ~isfield(LastData,'transmatList') || isempty(LastData.transmatList) LastData.transmatList=transmat(:)'; else LastData.transmatList=[LastData.transmatList ; transmat(:)']; end
    end
    
    if adj_obs && ~isempty(exp_num_emit)
        obsmat = mk_stochastic(exp_num_emit);
        if ~isfield(LastData,'ObsList') || isempty(LastData.ObsList) LastData.ObsList=obsmat(:)'; else LastData.ObsList=[LastData.ObsList ; obsmat(:)']; end
    end
    
    time_for_this_iteration = toc(overall_start_time) - iter_start_time_abs_loop;
    
    if verbose
        print_frequency = 25; % Or make this an input option
        if mod(num_iter, print_frequency) == 0 || first_iter_of_this_segment
            fprintf(1, '%s - Iteration %d, LogLik = %.6f, IterTime = %.2f s\n', datestr(now, 'HH:MM:SS'), num_iter, loglik, time_for_this_iteration);
        end
        if first_iter_of_this_segment, first_iter_of_this_segment = false; end 
    end
    
    % Convergence Check
    % Only check convergence if not the very first iteration of a fresh run (where previous_loglik is -inf)
    if num_iter > 1 || (LastData.Iteration > 0 && num_iter == LastData.Iteration + 1) 
        converged = em_converged(loglik, previous_loglik, thresh);
        if converged, fprintf('H2MM_MAIN: Converged at iteration %d. Delta log-likelihood (%.2e) <= threshold (%.1e).\n', num_iter, loglik - previous_loglik, thresh); end
    end
    
    if ~converged && (toc(overall_start_time)/3600 > MaxRunTime_hours)
        converged = 1; fprintf('H2MM_MAIN: Max runtime (%.1f hours) exceeded at iteration %d. Stopping.\n', MaxRunTime_hours, num_iter);
    end
    
    previous_loglik = loglik;
    if num_iter <= length(LL), LL(num_iter) = loglik; else LL = [LL, loglik]; end % Append if needed
    if num_iter <= length(iteration_times), iteration_times(num_iter) = time_for_this_iteration; else iteration_times = [iteration_times, time_for_this_iteration]; end

    LastData.Iteration = num_iter; LastData.LL = LL;
    LastData.last_obsmat = obsmat; LastData.last_prior = prior; LastData.last_transmat = transmat;
    LastData.interation_times = iteration_times; % Ensure this field name is consistent
    
    if mod(num_iter, backup_frequency) == 0
        try 
            save(fileName_LastData, 'LastData'); 
        catch ME_save_resume
            warning('H2MM_MAIN_SAVE_RESUME_FAIL', 'Failed to save resume data at iter %d to %s. Error: %s', num_iter, fileName_LastData, ME_save_resume.message); 
        end
    end
    num_iter =  num_iter + 1;
end 

total_elapsed_time_loop = toc(overall_start_time);
fprintf('H2MM_MAIN: EM loop finished. Total time: %.2f seconds.\n', total_elapsed_time_loop);

if num_iter > max_iter && ~converged, fprintf('H2MM_MAIN: Reached max_iter (%d) threshold.\n', max_iter); end

lastIteration = num_iter - 1; 
if lastIteration > 0
    LL = LL(1:min(lastIteration, length(LL))); 
else % No iterations were performed (e.g. already converged or max_iter=0)
    LL = []; % Or handle as appropriate, e.g. initial loglik if calculated
end
all_gamma = all_gamma_current_iter; % Gamma from the last E-step
fprintf('Reached end of while loop (end of HMM optimization for this guess).\n');


% --- NESTED FUNCTION DEFINITION for compute_ess_dhmm ---
    function [loglik_nf, exp_num_trans_nf, exp_num_visits1_nf, exp_num_emit_nf, exp_num_visitsT_nf, all_gamma_nf] = ...
        compute_ess_dhmm(startprob_ess, transmat_ess, obsmat_ess, data_ess, dirichlet_ess, TotalArrivalDelta_ess, R_ess, par_ess)
    % COMPUTE_ESS_DHMM Compute the Expected Sufficient Statistics
    % This function is nested, so it uses variables from the h2mm_main workspace if not masked by inputs.
    % It calculates transmat_t and Rho internally using the current transmat_ess.
    
    % Ensure all input names are unique if they might clash with h2mm_main workspace vars
    % (e.g., _ess suffix added here for clarity if this function was separate).
    % Since it's nested, it will use its own local arguments first.

    % Renaming for clarity within this nested function, if desired, but not strictly necessary
    % as MATLAB handles scope. Using suffixes like _nf (nested function) for outputs.
    
    numex_nf = length(data_ess);
    [S_nf ,O_nf] = size(obsmat_ess); % S_nf = Nstates
    exp_num_trans_nf = zeros(S_nf, S_nf);
    exp_num_visits1_nf = zeros(S_nf, 1);
    exp_num_visitsT_nf = zeros(S_nf, 1); % Though not used in M-step of h2mm_main
    
    % Initialize exp_num_emits as a cell array, one per trajectory
    exp_num_emits_cell_nf = cell(1, numex_nf);
    for ex_nf = 1:numex_nf
        exp_num_emits_cell_nf{ex_nf} = dirichlet_ess * ones(S_nf, O_nf);
    end
    
    loglik_nf = 0;
    all_gamma_nf = cell(1, numex_nf); % Store gamma for each trajectory
    
    % These are critical and depend on the current iteration's transmat_ess
    % Ensure TotalArrivalDelta_ess is a row vector for Factorization if Factorization expects that.
    % R_ess was already computed based on TotalArrivalDelta (which doesn't change per iteration).
    % transmat_t and Rho need to be computed with the current transmat_ess.
    if isempty(TotalArrivalDelta_ess) || isempty(R_ess)
        % Handle scenario with no valid time deltas (e.g., all trajectories < 2 photons)
        % In this case, fwdback might only use prior*obslik, or need specific handling.
        % This assumes fwdback_photonByphoton_fast can handle empty/trivial transmat_t and Rho.
        transmat_t_calc_nf = zeros(S_nf,S_nf,0); % Empty, implies no time evolution beyond obslik
        Rho_calc_nf = {};
        warning('COMPUTE_ESS_DHMM: TotalArrivalDelta or R is empty. Time evolution in HMM might be trivial.');
    else
        transmat_t_calc_nf = CalculatePowerOfTransMatrices(R_ess, transmat_ess); % Uses current transmat
        Rho_calc_nf = Calc_Rho(transmat_t_calc_nf, R_ess);
    end

    if par_ess == 1 && ~isempty(gcp('nocreate')) % Check if parpool exists if using parfor
        % fprintf('DEBUG COMPUTE_ESS_DHMM: Using parfor for E-step.\n');
        parfor ex_nf = 1:numex_nf
            obs_nf = data_ess{ex_nf}(:,2);
            delta_t_nf = data_ess{ex_nf}(:,1);
            % delta_t_nf contains indices or 0. fwdback needs to handle delta_t_nf(1)=0,
            % typically by using index 1 for transmat_t or a special case for prior.
            % The original fwdback_photonByphoton_fast had delta_t(delta_t==0)=1;
            delta_t_nf(delta_t_nf == 0) = 1; % Ensure indices are at least 1
            
            T_nf = length(obs_nf);
            if T_nf == 0, continue; end % Skip empty trajectories
            
            obslik_nf = multinomial_prob(obs_nf, obsmat_ess);
            
            [~, ~, gamma_nf, current_ll_nf, xi_summed_nf] = ...
                fwdback_photonByphoton_fast(startprob_ess, transmat_ess, obslik_nf, delta_t_nf, Rho_calc_nf, transmat_t_calc_nf);
            
            % Accumulation inside parfor needs to be handled carefully.
            % Direct accumulation like loglik_nf = loglik_nf + current_ll_nf won't work.
            % Instead, collect results and sum after the loop.
            % For now, assuming fwdback returns components that can be summed/processed later.
            % However, the original nested compute_ess_dhmm did direct accumulation,
            % which suggests it was not intended for its internal loop to be the parfor target this way.
            % The parfor in original code was on 'ex'.
            
            % Let's match original structure which accumulated inside parfor using temporary variables
            % This requires a reduction variable for loglik_nf and careful handling of others.
            % For simplicity and directness, let's assume the original parfor accumulated correctly.
            % If not, this part needs a rewrite for proper parfor accumulation.
            % The simplest way if parfor is here, is to make outputs cell arrays.
            
            % Reverting to original accumulation style (will make parfor less efficient without proper reduction)
            % THIS IS WHERE A PARFOR WOULD NEED CAREFUL REDUCTION VARIABLES
            % OR OUTPUTS AS CELL ARRAYS AND SUMMING AFTERWARDS.
            % For now, mimicking the original structure that might not be parfor-optimal for these sums.
            
            % Store results for later aggregation if using parfor for 'ex_nf'
            loglik_temp_array(ex_nf) = current_ll_nf;
            exp_num_trans_temp_cell{ex_nf} = xi_summed_nf;
            exp_num_visits1_temp_cell{ex_nf} = gamma_nf(:,1);
            all_gamma_nf{ex_nf} = gamma_nf; % Store gamma for this trajectory

            % Handle exp_num_emits_cell_nf{ex_nf}
            current_exp_emit = dirichlet_ess * ones(S_nf, O_nf); % Re-initialize for each traj
            if T_nf < O_nf % If fewer time points than observation types
                for t_nf = 1:T_nf
                    current_exp_emit(:, obs_nf(t_nf)) = current_exp_emit(:, obs_nf(t_nf)) + gamma_nf(:, t_nf);
                end
            else % More time points, iterate over observation types for efficiency
                for o_nf = 1:O_nf
                    ndx_nf = find(obs_nf == o_nf);
                    if ~isempty(ndx_nf)
                        current_exp_emit(:, o_nf) = current_exp_emit(:, o_nf) + sum(gamma_nf(:, ndx_nf), 2);
                    end
                end
            end
            exp_num_emits_cell_nf{ex_nf} = current_exp_emit; % Store for this trajectory
        end % End parfor
        
        % Aggregate results from parfor
        loglik_nf = sum(loglik_temp_array);
        for ex_nf_agg = 1:numex_nf
            exp_num_trans_nf = exp_num_trans_nf + exp_num_trans_temp_cell{ex_nf_agg};
            exp_num_visits1_nf = exp_num_visits1_nf + exp_num_visits1_temp_cell{ex_nf_agg};
            % exp_num_emits_cell_nf already populated correctly per trajectory
        end

    else % Regular for loop (par_ess == 0 or no parpool)
        % fprintf('DEBUG COMPUTE_ESS_DHMM: Using serial for loop for E-step.\n');
        for ex_nf = 1:numex_nf
            obs_nf = data_ess{ex_nf}(:,2);
            delta_t_nf = data_ess{ex_nf}(:,1);
            delta_t_nf(delta_t_nf == 0) = 1; % Ensure indices are at least 1
            
            T_nf = length(obs_nf);
            if T_nf == 0, all_gamma_nf{ex_nf} = zeros(S_nf,0); continue; end % Skip empty trajectories
            
            obslik_nf = multinomial_prob(obs_nf, obsmat_ess);
            
            [~, ~, gamma_nf, current_ll_nf, xi_summed_nf] = ...
                fwdback_photonByphoton_fast(startprob_ess, transmat_ess, obslik_nf, delta_t_nf, Rho_calc_nf, transmat_t_calc_nf);
            
            loglik_nf = loglik_nf + current_ll_nf;
            exp_num_trans_nf = exp_num_trans_nf + xi_summed_nf;
            exp_num_visits1_nf = exp_num_visits1_nf + gamma_nf(:,1);
            all_gamma_nf{ex_nf} = gamma_nf; 
            
            % exp_num_emits_cell_nf{ex_nf} is initialized. Add gamma contributions.
            if T_nf < O_nf
                for t_nf = 1:T_nf
                    exp_num_emits_cell_nf{ex_nf}(:, obs_nf(t_nf)) = exp_num_emits_cell_nf{ex_nf}(:, obs_nf(t_nf)) + gamma_nf(:, t_nf);
                end
            else
                for o_nf = 1:O_nf
                    ndx_nf = find(obs_nf == o_nf);
                    if ~isempty(ndx_nf)
                        exp_num_emits_cell_nf{ex_nf}(:, o_nf) = exp_num_emits_cell_nf{ex_nf}(:, o_nf) + sum(gamma_nf(:, ndx_nf), 2);
                    end
                end
            end
        end % End for ex_nf
    end % End if par_ess

    % Summing emission probability measure from all trajectories
    exp_num_emit_nf = zeros(S_nf, O_nf); % Initialize accumulator
    if ~isempty(exp_num_emits_cell_nf) && ~isempty(exp_num_emits_cell_nf{1})
         exp_num_emit_nf = exp_num_emits_cell_nf{1}; % Start with the first one
        for ex_nf_sum = 2:numex_nf
            if ~isempty(exp_num_emits_cell_nf{ex_nf_sum})
                exp_num_emit_nf = exp_num_emit_nf + exp_num_emits_cell_nf{ex_nf_sum};
            end
        end
    else % if dirichlet_ess was 0 and no trajectories, exp_num_emit_nf could be all zeros
        exp_num_emit_nf = dirichlet_ess * ones(S_nf, O_nf) * numex_nf; % if dirichlet_ess > 0
    end

    end % End of NESTED function compute_ess_dhmm

end % End of MAIN function h2mm_main