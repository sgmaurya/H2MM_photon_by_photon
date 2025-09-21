% viterbi_path_PhotnoByPhoton2_local.m
% Viterbi algorithm using a rate matrix approach for variable time steps,
% with a fallback to discrete power of transition matrix.

function path = viterbi_path_PhotnoByPhoton2_local(prior, model_transmat_unit_time, obs_likelihoods, actual_inter_photon_durations_sec, dt_analysis_sec)
% Inputs:
% prior: Initial state probabilities (Qx1 column vector).
% model_transmat_unit_time: Transition probability matrix P(t+dt_analysis | t) (QxQ).
% obs_likelihoods(state, photon_idx): P(Observation_k | State_s) for each photon (QxT).
% actual_inter_photon_durations_sec: Vector of actual time durations (in seconds)
%                                    between photons. Length is T-1 for T photons.
%                                    Can be empty if T <= 1.
% dt_analysis_sec: The unit time step (in seconds) for which model_transmat_unit_time is defined.

% --- Input Validation (Basic) ---
if dt_analysis_sec <= 0
    error('VITERBI_LOCAL_ERROR: dt_analysis_sec must be positive.');
end
if ~iscolumn(prior)
    prior = prior(:); % Ensure prior is a column vector
    if ~iscolumn(prior) % Check again after transpose if it was a row
         error('VITERBI_LOCAL_ERROR: prior must be a vector.');
    end
end

scaled = 1; % Use scaling to prevent underflow/overflow

T = size(obs_likelihoods, 2); % Number of photons (observations)
Q = length(prior);   % Number of hidden states

if Q ~= size(model_transmat_unit_time,1) || Q ~= size(model_transmat_unit_time,2) || Q ~= size(obs_likelihoods,1)
    error('VITERBI_LOCAL_ERROR: Dimension mismatch between prior, transmat, or obs_likelihoods regarding number of states (Q).');
end
if T > 0 && ~isempty(actual_inter_photon_durations_sec) && length(actual_inter_photon_durations_sec) ~= T-1
    error('VITERBI_LOCAL_ERROR: Length of actual_inter_photon_durations_sec must be T-1 if T > 1 and durations are provided.');
end


% Initialize Viterbi variables
delta = zeros(Q,T);  % delta(j,t) = max probability of any path ending at state j at time t
psi = zeros(Q,T);    % psi(j,t) = argmax state (predecessor) for state j at time t
path = zeros(1,T);   % Stores the Viterbi path

if T == 0
    path = []; % No photons, no path
    return;
end

% --- Convert P (Transition Probability Matrix) to K (Rate Matrix) ---
K_rate_matrix = zeros(Q,Q);
valid_K_matrix = false;
try
    % Ensure model_transmat_unit_time is suitable for logm
    temp_transmat = model_transmat_unit_time;
    min_diag_val = 1e-12; % Small positive number to avoid log(0) or log(negative)

    % Check for non-positive diagonal elements which make logm fail for TPMs
    if any(diag(temp_transmat) <= min_diag_val)
        % fprintf('VITERBI_LOCAL_INFO: Diagonal elements of transmat near/below zero found. Attempting clamp for logm.\n');
        for r_idx=1:Q
            if temp_transmat(r_idx,r_idx) <= min_diag_val
                temp_transmat(r_idx,r_idx) = min_diag_val;
            end
        end
        % Re-normalize rows if clamped
        temp_transmat = bsxfun(@rdivide, temp_transmat, sum(temp_transmat,2));
        temp_transmat(isnan(temp_transmat)) = 1/Q; % Handle rows that summed to 0 (should be rare)
    end

    K_rate_matrix_complex = logm(temp_transmat) / dt_analysis_sec;

    % Check for significant imaginary parts
    if max(abs(imag(K_rate_matrix_complex(:)))) > 1e-6 % Tolerance for numerical noise
        warning('VITERBI_LOCAL_WARN: logm resulted in a rate matrix with significant imaginary parts. This suggests the transition matrix might not be embeddable or dt_analysis_sec is too large. K matrix may be unreliable.');
    end
    K_rate_matrix = real(K_rate_matrix_complex);

    % Ensure K is a valid generator matrix (non-negative off-diagonals, rows sum to ~0)
    for r_idx = 1:Q
        for c_idx = 1:Q
            if r_idx ~= c_idx && K_rate_matrix(r_idx, c_idx) < 0
                K_rate_matrix(r_idx, c_idx) = 0; % Force non-negative off-diagonals
            end
        end
    end
    for r_idx = 1:Q
        K_rate_matrix(r_idx,r_idx) = 0; % Temporarily zero out diagonal
        K_rate_matrix(r_idx,r_idx) = -sum(K_rate_matrix(r_idx,:)); % Ensure row sum is 0
    end
    valid_K_matrix = true;
catch ME_logm
    warning('VITERBI_LOCAL_WARN: Could not compute valid rate matrix K from transmat using logm. Identifier: %s, Message: %s. Falling back to P^k for transitions.', ME_logm.identifier, ME_logm.message);
    valid_K_matrix = false;
end
% --- End K calculation ---

% Initialization (t=1, first photon)
delta(:,1) = prior .* obs_likelihoods(:,1);
if scaled
    [delta(:,1), norm_factor_init] = normalise(delta(:,1));
    if norm_factor_init == 0 && T > 0 % All initial paths are zero prob
         % This can happen if prior is e.g. [1;0] and obs_likelihoods(1,1) is 0.
         % fprintf('VITERBI_LOCAL_INFO: Initial delta calculation resulted in all zeros. Resetting to uniform.\n');
         delta(:,1) = ones(Q,1)/Q; % Reset to uniform to allow Viterbi to proceed
    end
end
psi(:,1) = 0; % Arbitrary, no predecessor for the first photon

% Recursion (t=2 to T)
for t_photon = 2:T
    delta_t_actual_sec_step = dt_analysis_sec; % Default to unit step
    if ~isempty(actual_inter_photon_durations_sec) && (t_photon-1) <= length(actual_inter_photon_durations_sec)
        delta_t_actual_sec_step = actual_inter_photon_durations_sec(t_photon-1);
        if delta_t_actual_sec_step <= 0 % Should have been handled by caller
            % warning('VITERBI_LOCAL_WARN: Inter-photon duration for step %d is non-positive (%g s). Using unit dt_analysis_sec.', t_photon, delta_t_actual_sec_step);
            delta_t_actual_sec_step = dt_analysis_sec;
        end
    else
        % This case should ideally not be reached if inputs are correct (T-1 durations for T photons)
        % warning('VITERBI_LOCAL_WARN: Missing inter-photon duration for step %d. Using unit dt_analysis_sec.', t_photon);
    end

    transmat_for_current_duration = [];
    if valid_K_matrix
        try
            transmat_for_current_duration = expm(K_rate_matrix * delta_t_actual_sec_step);
            % Ensure rows sum to 1 (expm can sometimes have small numerical deviations)
            row_sums = sum(transmat_for_current_duration,2);
            if any(abs(row_sums - 1) > 1e-6) % Check if sums are significantly off from 1
                 transmat_for_current_duration = bsxfun(@rdivide, transmat_for_current_duration, row_sums);
                 transmat_for_current_duration(isnan(transmat_for_current_duration) | isinf(transmat_for_current_duration)) = 1/Q; % Handle division by zero sum
            end
            if any(transmat_for_current_duration(:) < 0) % Probabilities cannot be negative
                % This can happen if K was problematic despite checks
                % warning('VITERBI_LOCAL_WARN: expm(K*t) resulted in negative probabilities at t_photon=%d. Fallback to P^k.', t_photon);
                valid_K_matrix = false; % Force fallback for this and subsequent steps
            end
        catch ME_expm
            warning('VITERBI_LOCAL_WARN: Error in expm(K*t) at t_photon=%d. Identifier: %s, Message: %s. Fallback to P^k.', t_photon, ME_expm.identifier, ME_expm.message);
            valid_K_matrix = false; % Force fallback
        end
    end

    if ~valid_K_matrix % Fallback for this step if K was invalid or expm failed
        num_unit_steps = round(delta_t_actual_sec_step / dt_analysis_sec);
        if num_unit_steps < 1, num_unit_steps = 1; end
        % fprintf('VITERBI_LOCAL_INFO: Using transmat^%d (fallback) for step %d.\n',num_unit_steps, t_photon);
        transmat_for_current_duration = model_transmat_unit_time ^ num_unit_steps;
    end

    for j_state = 1:Q % For each current state j at current photon t_photon
        % Maximize over previous states i: delta(i, t_photon-1) * P(State_j at t_photon | State_i at t_photon-1)
        % P(State_j at t_photon | State_i at t_photon-1) is transmat_for_current_duration(i, j_state)
        % So we need delta_prev_col .* transmat_col_for_j
        [max_val, best_prev_state_idx] = max(delta(:,t_photon-1) .* transmat_for_current_duration(:,j_state));
        delta(j_state,t_photon) = max_val;
        psi(j_state,t_photon) = best_prev_state_idx;
    end
    delta(:,t_photon) = delta(:,t_photon) .* obs_likelihoods(:,t_photon);

    if scaled
        [delta(:,t_photon), norm_factor] = normalise(delta(:,t_photon));
        if norm_factor == 0 && T > 1 % If all paths have zero probability
             % warning('VITERBI_LOCAL_WARN: All paths have zero probability at photon %d. Path may be arbitrary from here.', t_photon);
             delta(:,t_photon) = ones(Q,1)/Q; % Reset to uniform to allow continuation
        end
    end
end

% Termination: Find the best state for the last photon
if all(delta(:,T) == 0) && T > 0 % If all final states have zero probability
    % warning('VITERBI_LOCAL_WARN: All final states have zero probability. Assigning path arbitrarily.');
    % Try to make a reasonable guess based on the last observation and prior
    [~, path(T)] = max(prior .* obs_likelihoods(:,T));
    if isempty(path(T)) || path(T) == 0 % Should not happen if obs_likelihoods is valid
        path(T) = 1; % Absolute fallback
    end
else
    [~, path(T)] = max(delta(:,T));
end


% Path backtracking
for t_backtrack = T-1:-1:1
    % Ensure path(t_backtrack+1) is a valid index for psi
    if path(t_backtrack+1) >= 1 && path(t_backtrack+1) <= Q && psi(path(t_backtrack+1), t_backtrack+1) >=1 && psi(path(t_backtrack+1), t_backtrack+1) <=Q
        path(t_backtrack) = psi(path(t_backtrack+1), t_backtrack+1);
    else
        % warning('VITERBI_LOCAL_WARN: Invalid state index during backtracking at t=%d from state %d. Path may be corrupted. Using state 1.', t_backtrack+1, path(t_backtrack+1));
        path(t_backtrack) = 1; % Fallback, arbitrary choice
    end
end

end