% h2mm_wexac_ver1_local.m
% This script prepares to run the HMM analysis locally,
% mimicking the cluster environment for one job from an array.
% It expects 'varargin.mat' to be in its current directory (set by local_hmm_runner)
% and global variable LOCAL_JOB_INDEX_FOR_HMM_CORE to be set.
% It also uses LOCAL_PAR_CONTROL_FOR_HMM_CORE for parpool management.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Comments by Demian Liebermann 02.10.2021
% (Adapted for local execution)
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get variables for local execution and setup
% --- MODIFIED/ADDED FOR LOCAL EXECUTION ---
global LOCAL_JOB_INDEX_FOR_HMM_CORE; 
global LOCAL_PAR_CONTROL_FOR_HMM_CORE; 

if isempty(LOCAL_JOB_INDEX_FOR_HMM_CORE)
    error('HMM_CORE_ERROR: LOCAL_JOB_INDEX_FOR_HMM_CORE is not set. This script should be called via local_hmm_runner.');
end
ind = LOCAL_JOB_INDEX_FOR_HMM_CORE; % The specific index for this HMM run

warning off all; % Suppress some common warnings

% Paths to dependent HMM functions (e.g., h2mm_main, Factorization)
% are assumed to be managed by the main calling script (runHMM_local.m)
% by adding them to the MATLAB path before this script is run.
% --- END MODIFIED/ADDED FOR LOCAL EXECUTION ---

%% Debug: Check current directory and varargin.mat existence at script start
fprintf('HMM_CORE_DEBUG (Guess %d): Script started. Current directory is: %s\n', ind, pwd);
varargin_mat_filename = 'varargin.mat';
varargin_mat_full_path_check = fullfile(pwd, varargin_mat_filename);

if exist(varargin_mat_full_path_check, 'file')
    fprintf('HMM_CORE_DEBUG (Guess %d): %s FOUND in current directory at script start.\n', ind, varargin_mat_filename);
else
    fprintf('HMM_CORE_CRITICAL (Guess %d): %s NOT FOUND in current directory (%s) at script start!\n', ind, varargin_mat_filename, pwd);
    % Attempt to list files to see what is present
    fprintf('HMM_CORE_DEBUG (Guess %d): Contents of current directory (%s):\n', ind, pwd);
    ls;
    error('HMM_CORE_ERROR: Required varargin.mat not found in expected location.');
end

%% Loading varargin parameters and the initial guess array
try
    % Load the variable 'varargin_to_save' (which contains the cell array) from varargin.mat
    loaded_data_from_mat = load(varargin_mat_filename, 'varargin_to_save'); 
    if ~isfield(loaded_data_from_mat, 'varargin_to_save')
        error('HMM_CORE_ERROR: "varargin_to_save" variable not found inside %s.', varargin_mat_filename);
    end
    varargin_cell_loaded = loaded_data_from_mat.varargin_to_save; % This is the actual cell array {}
    fprintf('HMM_CORE_INFO (Guess %d): Successfully loaded %s and extracted "varargin_to_save".\n', ind, varargin_mat_filename);
catch ME_load_varg
    error('HMM_CORE_ERROR (Guess %d): Failed to load %s or extract "varargin_to_save" variable. Error: %s', ind, varargin_mat_filename, ME_load_varg.message);
end

% Initialize variables to be extracted from varargin_cell_loaded
max_iter_val_loaded = []; 
MaxRunTime_val_loaded = [];
OBSFix_val_loaded = 1; % Default, as in original cluster submission scripts
par_val_loaded = 0;    % Default for internal HMM parallelism (serial E-step)
fix_val_loaded = [];
range_FRET_val_loaded = [];
H_val_loaded = [];
InitialGues_loaded_struct = []; % Will hold the single struct for this 'ind'
data_loaded = [];

% Iterate through the loaded cell array 'varargin_cell_loaded'
for i_varg = 1:2:length(varargin_cell_loaded)
    param_name = varargin_cell_loaded{i_varg};
    param_value = varargin_cell_loaded{i_varg+1};
    try
        switch param_name
            case 'data'
                data_loaded = param_value;
            case 'number_of_initialguess'
                % Not directly used in this script, but good to know it's there
            case 'InitialGues'
                % param_value here is the *cell array* of all initial guesses
                if ind > 0 && ind <= length(param_value)
                    InitialGues_loaded_struct = param_value{ind}; % Select the one for current 'ind'
                else
                    error('HMM_CORE_ERROR (Guess %d): Invalid job index %d for InitialGues array of length %d.', ind, length(param_value));
                end
            case 'max_iter'
                max_iter_val_loaded = param_value;
            case 'MaxRunTime'
                MaxRunTime_val_loaded = param_value;
            case 'OBS'
                OBSFix_val_loaded = param_value;
            case 'par'
                par_val_loaded = param_value;
            case 'fixed_state'
                fix_val_loaded = param_value;
            case 'range_FRET'
                range_FRET_val_loaded = param_value;
            case 'H'
                H_val_loaded = param_value;
        end
    catch ME_varg_extract
        warning('HMM_CORE_WARNING (Guess %d): Error processing varargin field "%s". %s', ind, param_name, ME_varg_extract.message);
    end
end

% Assign to variables used by h2mm_main* functions for clarity
data_for_hmm = data_loaded;
InitialGuess_for_hmm = InitialGues_loaded_struct; 
max_iter_for_hmm = max_iter_val_loaded;
MaxRunTime_for_hmm_hours = MaxRunTime_val_loaded;
OBSFix_for_hmm = OBSFix_val_loaded;
par_internal_hmm_flag = par_val_loaded; 
fixed_state_param = fix_val_loaded;
range_FRET_param = range_FRET_val_loaded;
H_param = H_val_loaded;

% Basic validation of loaded parameters
if isempty(InitialGuess_for_hmm)
    error('HMM_CORE_ERROR (Guess %d): Initial guess for index %d was not loaded correctly from varargin.mat.', ind);
end
if isempty(data_for_hmm)
    error('HMM_CORE_ERROR (Guess %d): Trajectory data was not loaded correctly from varargin.mat.', ind);
end
fprintf('HMM_CORE_INFO (Guess %d): Parameters successfully extracted from varargin.mat.\n', ind);

%% Activate parpool (if par_internal_hmm_flag is 1 and not controlled externally)
manage_parpool_in_core_script = true; 
if ~isempty(LOCAL_PAR_CONTROL_FOR_HMM_CORE) && LOCAL_PAR_CONTROL_FOR_HMM_CORE == false
    manage_parpool_in_core_script = false; 
end

pool_was_started_by_this_script = false;
if par_internal_hmm_flag == 1 && manage_parpool_in_core_script
    current_pool_obj = gcp('nocreate'); 
    if isempty(current_pool_obj)
        fprintf('HMM_CORE_INFO (Guess %d, PID %d): Starting new parpool for internal HMM E-step.\n', ind, feature('getpid'));
        try
            parpool; 
            pool_was_started_by_this_script = true;
        catch ME_parpool
            warning('HMM_CORE_WARNING (Guess %d): Failed to start parpool. Running serially. Error: %s', ind, ME_parpool.message);
            par_internal_hmm_flag = 0; % Fallback to serial if pool fails
        end
    else
        fprintf('HMM_CORE_INFO (Guess %d, PID %d): Parpool already active. Using existing for internal HMM E-step.\n', ind, feature('getpid'));
    end
elseif par_internal_hmm_flag == 1 && ~manage_parpool_in_core_script
    fprintf('HMM_CORE_INFO (Guess %d, PID %d): Parpool management deferred to external control. Internal HMM E-step will use existing pool if available.\n', ind, feature('getpid'));
end

%% Load last iteration data, for resuming interrupted jobs (configurable for local runs)
RESUME_LOCAL_RUNS_ENABLED = false; % Configuration: true to enable resume, false for fresh runs.

fileName_LastData_mat = sprintf('LastData_job%d.mat', ind); % Use sprintf for consistency
LastData_struct = struct('Iteration', 0); % Initialize for a fresh run

if RESUME_LOCAL_RUNS_ENABLED && exist(fileName_LastData_mat, 'file')
    fprintf('HMM_CORE_INFO (Guess %d): Attempting to resume from %s.\n', ind, fileName_LastData_mat);
    try
        loaded_LastData = load(fileName_LastData_mat, 'LastData'); 
        LastData_struct = loaded_LastData.LastData; 
        fprintf('HMM_CORE_INFO (Guess %d): Successfully loaded LastData. Resuming from iteration %d.\n', ind, LastData_struct.Iteration + 1);
    catch ME_load_last_data
        warning('HMM_CORE_WARNING (Guess %d): Failed to load %s. Starting fresh. Error: %s', ind, fileName_LastData_mat, ME_load_last_data.message);
        LastData_struct.Iteration = 0; 
    end
else
    if RESUME_LOCAL_RUNS_ENABLED % File didn't exist
        fprintf('HMM_CORE_INFO (Guess %d): No resume file (%s) found. Starting fresh.\n', ind, fileName_LastData_mat);
    else % Resume was disabled
        fprintf('HMM_CORE_INFO (Guess %d): Resume disabled. Starting fresh.\n', ind);
    end
end

%% Assigns maximum run time and iteration numbers, if not provided in varargin
if isempty(max_iter_for_hmm)
    max_iter_for_hmm = 10000; 
end
if isempty(MaxRunTime_for_hmm_hours)
    MaxRunTime_for_hmm_hours = 24; 
end

%% Load initial guess parameters from InitialGuess_for_hmm or from LastData_struct
if LastData_struct.Iteration == 0 
    prior0_hmm = InitialGuess_for_hmm.prior0;
    transmat0_hmm = InitialGuess_for_hmm.transmat0;
    obsmat0_hmm = InitialGuess_for_hmm.obsmat0;

    LastData_struct.transmatList = []; 
    LastData_struct.ObsList = [];      
    LastData_struct.PriorList = [];    
    LastData_struct.last_prior = prior0_hmm;
    LastData_struct.last_transmat = transmat0_hmm;
    LastData_struct.last_obsmat = obsmat0_hmm;
    LastData_struct.LL = -Inf(1, max_iter_for_hmm);         
    LastData_struct.interation_times = zeros(1, max_iter_for_hmm); 
    LastData_struct.InteruptCount = 0; 
else 
    fprintf('HMM_CORE_INFO (Guess %d): Resuming. Last completed iteration was %d. Interruption count: %d.\n', ...
        ind, LastData_struct.Iteration, LastData_struct.InteruptCount); % Iteration is last completed
    LastData_struct.InteruptCount = LastData_struct.InteruptCount + 1; 
    prior0_hmm = LastData_struct.last_prior;
    transmat0_hmm = LastData_struct.last_transmat;
    obsmat0_hmm = LastData_struct.last_obsmat;
    
    current_ll_len = length(LastData_struct.LL);
    if current_ll_len < max_iter_for_hmm && current_ll_len > 0 && LastData_struct.LL(current_ll_len) == -Inf && LastData_struct.Iteration < current_ll_len
        % If last entry is -Inf and iteration is less, means it was preallocated beyond last run
        % No need to extend if already correctly preallocated by previous run.
    elseif current_ll_len < max_iter_for_hmm
        LastData_struct.LL(current_ll_len+1 : max_iter_for_hmm) = -Inf;
    end

    current_times_len = length(LastData_struct.interation_times);
    if current_times_len < max_iter_for_hmm
        LastData_struct.interation_times(current_times_len+1 : max_iter_for_hmm) = 0;
    end
end

%% Call the appropriate H2MM main function
LL_final_trace = []; 
prior_final = []; 
transmat_final = []; 
obsmat_final = []; 
lastIteration_completed = 0; 
all_gamma_final = {};

fprintf('HMM_CORE_INFO (Guess %d): Starting HMM optimization.\n', ind);
try
    if ~isempty(range_FRET_param) && ~isempty(fixed_state_param) 
        fprintf('HMM_CORE_INFO (Guess %d): Using h2mm_main_range.\n', ind);
        [LL_final_trace, prior_final, transmat_final, obsmat_final, lastIteration_completed, all_gamma_final] = ...
            h2mm_main_range(data_for_hmm, prior0_hmm, transmat0_hmm, obsmat0_hmm, max_iter_for_hmm, MaxRunTime_for_hmm_hours, OBSFix_for_hmm, LastData_struct, fileName_LastData_mat, par_internal_hmm_flag, fixed_state_param, range_FRET_param);
    elseif ~isempty(H_param)
        % Add specific cases for H_param if you have other h2mm_main variants
        switch H_param
            % case 1 
            %   [LL_final_trace, prior_final, transmat_final, obsmat_final, lastIteration_completed, all_gamma_final] = h2mm_main_Hagen(...);
            % case 3 
            %   [LL_final_trace, prior_final, transmat_final, obsmat_final, lastIteration_completed, all_gamma_final] = h2mm_main_blink(...);
            otherwise
                fprintf('HMM_CORE_WARNING (Guess %d): Unknown H_param value %d. Using default h2mm_main.\n', ind, H_param);
                [LL_final_trace, prior_final, transmat_final, obsmat_final, lastIteration_completed, all_gamma_final] = ...
                    h2mm_main(data_for_hmm, prior0_hmm, transmat0_hmm, obsmat0_hmm, max_iter_for_hmm, MaxRunTime_for_hmm_hours, OBSFix_for_hmm, LastData_struct, fileName_LastData_mat, par_internal_hmm_flag);
        end
    else 
        fprintf('HMM_CORE_INFO (Guess %d): Using default h2mm_main.\n', ind);
        [LL_final_trace, prior_final, transmat_final, obsmat_final, lastIteration_completed, all_gamma_final] = ...
            h2mm_main(data_for_hmm, prior0_hmm, transmat0_hmm, obsmat0_hmm, max_iter_for_hmm, MaxRunTime_for_hmm_hours, OBSFix_for_hmm, LastData_struct, fileName_LastData_mat, par_internal_hmm_flag);
    end
catch ME_hmm_main_call
    fprintf('HMM_CORE_ERROR (Guess %d): Error during HMM main function execution: %s\n', ind, ME_hmm_main_call.message);
    disp(ME_hmm_main_call.getReport()); 
    if exist(fileName_LastData_mat, 'file') 
        try
            loaded_LastData_on_error = load(fileName_LastData_mat, 'LastData');
            LastData_struct = loaded_LastData_on_error.LastData;
        catch
            % Use LastData_struct as it was before the h2mm_main call if reload fails
        end
    end
    % Populate output vars from LastData if they are empty after error
    if isempty(LL_final_trace) && isfield(LastData_struct, 'LL'), LL_final_trace = LastData_struct.LL; end
    if isempty(prior_final) && isfield(LastData_struct, 'last_prior'), prior_final = LastData_struct.last_prior; end
    if isempty(transmat_final) && isfield(LastData_struct, 'last_transmat'), transmat_final = LastData_struct.last_transmat; end
    if isempty(obsmat_final) && isfield(LastData_struct, 'last_obsmat'), obsmat_final = LastData_struct.last_obsmat; end
    if lastIteration_completed == 0 && isfield(LastData_struct, 'Iteration'), lastIteration_completed = LastData_struct.Iteration; end
end

fprintf('HMM_CORE_INFO (Guess %d): HMM optimization finished. Completed %d iterations.\n', ind, lastIteration_completed);

%% Save results of this individual HMM run
output_filename_for_this_hmm_run = sprintf('hmm_PhotonByPhotonJobArray%d.mat', ind);
fprintf('HMM_CORE_INFO (Guess %d): Saving results to %s.\n', ind, output_filename_for_this_hmm_run);

LL_out_save = LL_final_trace;
prior_out_save = prior_final;
transmat_out_save = transmat_final;
obsmat_out_save = obsmat_final;
lastIteration_out_save = lastIteration_completed;
all_gamma_out_save = all_gamma_final;
% LastData_struct now contains the final state of LastData from h2mm_main*
% InitialGuess_for_hmm is the initial guess structure used for this run

try
    save(output_filename_for_this_hmm_run, ...
        'LL_out_save', 'prior_out_save', 'transmat_out_save', 'obsmat_out_save', ...
        'lastIteration_out_save', 'all_gamma_out_save', ...
        'LastData_struct', 'InitialGuess_for_hmm', '-v7.3'); 
    fprintf('HMM_CORE_INFO (Guess %d): Results successfully saved.\n', ind);
catch ME_save
    fprintf('HMM_CORE_ERROR (Guess %d): Failed to save results to %s. Error: %s\n', ind, output_filename_for_this_hmm_run, ME_save.message);
end

%% Clean up parpool if it was started by this script
if pool_was_started_by_this_script
    fprintf('HMM_CORE_INFO (Guess %d, PID %d): Shutting down parpool started by this script.\n', ind, feature('getpid'));
    delete(gcp('nocreate'));
end

% The original "exit Matlab" command is REMOVED.
fprintf('HMM_CORE_INFO (Guess %d): Processing finished for this HMM run.\n', ind);

% It's generally better for the caller (local_hmm_runner) to clear the globals
% after it's done with this script, especially if run in a loop.
% clear global LOCAL_JOB_INDEX_FOR_HMM_CORE; 
% clear global LOCAL_PAR_CONTROL_FOR_HMM_CORE;