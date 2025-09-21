% --- START OF FILE InitialGuess_PhotonByPhoton_v2.m ---
function [InitialGues] = InitialGuess_PhotonByPhoton(number_of_initialguess_per_job, Nstates, O_channels, dt_analysis_sec)
% This function constructs an initial guess array for HMM analysis.
% Enhanced version:
% - Accepts O_channels for obsmat0 generation.
% - Accepts dt_analysis_sec for rate-aware transmat0 generation.

    % --- Parameters for Rate-Aware transmat0 ---
    % These define a plausible range for physical kinetic rates (transitions per second)
    % Adjust these based on typical expectations for the system being analyzed.
    k_min_Hz = 1e2;   % Minimum rate (e.g., 100 Hz)
    k_max_Hz = 1e7;   % Maximum rate (e.g., 10 MHz)
    % Target sum of off-diagonal probabilities in a row for transmat0.
    % Should be < 1 to leave room for diagonal. A smaller value means more persistent states.
    target_sum_off_diagonal_P = 0.1; % e.g., sum of P_ij (j~=i) aims for 0.1

    InitialGues = cell(1, number_of_initialguess_per_job); % Pre-allocate



    for iii = 1:number_of_initialguess_per_job
        transmat_temp = zeros(Nstates, Nstates);
        if Nstates == 0 % Handle empty case
            InitialGues{iii}.transmat0 = [];
            InitialGues{iii}.prior0 = [];
            InitialGues{iii}.obsmat0 = [];
            InitialGues{iii}.Nstate = 0;
            continue;
        end

        if Nstates == 1
            transmat_temp = 1; % Only one state, must stay there
        else % Nstates > 1
            for r = 1:Nstates % Current row (from state r)
                row_off_diagonal_P = zeros(1, Nstates); % Store P_rc for this row r
                raw_rates_for_scaling = zeros(1, Nstates); % Store k_guess_Hz for scaling if needed

                for c = 1:Nstates % Current col (to state c)
                    if r == c, continue; end

                    log_k_min = log10(k_min_Hz);
                    log_k_max = log10(k_max_Hz);
                    k_guess_Hz = 10^(log_k_min + rand() * (log_k_max - log_k_min));
                    
                    raw_rates_for_scaling(c) = k_guess_Hz; % Store the rate itself
                    row_off_diagonal_P(c) = k_guess_Hz * dt_analysis_sec;
                end

                % --- Scaling logic for off-diagonal probabilities ---
                current_sum_off_diag_P = sum(row_off_diagonal_P);

                if abs(current_sum_off_diag_P) < 1e-9 % If sum of P_ij is effectively zero (all rates*dt were tiny)
                    % Distribute target_sum_off_diagonal_P equally among off-diagonal elements
                    % This happens if k_max_Hz * dt_analysis_sec is still very small for all sampled k.
                    num_off_diag_elements = Nstates - 1;
                    if num_off_diag_elements > 0
                        scaled_P_rc = target_sum_off_diagonal_P / num_off_diag_elements;
                        for c_assign = 1:Nstates
                            if r == c_assign, continue; end
                            transmat_temp(r, c_assign) = scaled_P_rc;
                        end
                    end
                else
                    % Scale based on the ratio of raw rates to ensure P_ij sums to target_sum_off_diagonal_P
                    % This preserves the relative magnitudes of the initially sampled rates.
                    sum_raw_rates_for_row = sum(raw_rates_for_scaling); % Sum of k_ij for j~=i
                    if abs(sum_raw_rates_for_row) < 1e-9 % Should not happen if k_min_Hz > 0
                         % Fallback: equal distribution if all raw rates were somehow zero
                        num_off_diag_elements = Nstates - 1;
                        scaled_P_rc = target_sum_off_diagonal_P / num_off_diag_elements;
                        for c_assign = 1:Nstates
                            if r == c_assign, continue; end
                            transmat_temp(r, c_assign) = scaled_P_rc;
                        end
                    else
                        for c_assign = 1:Nstates
                            if r == c_assign, continue; end
                            % P_rc = (k_rc / sum(k_rj)) * target_sum_off_diagonal_P
                            transmat_temp(r, c_assign) = (raw_rates_for_scaling(c_assign) / sum_raw_rates_for_row) * target_sum_off_diagonal_P;
                        end
                    end
                end
                % --- End Scaling Logic ---
                
                % Set diagonal element
                current_row_sum_off_diag = sum(transmat_temp(r, find((1:Nstates) ~= r) ));
                transmat_temp(r, r) = 1 - current_row_sum_off_diag;

                % Safety net for diagonal (should ideally not be needed with target_sum < 1)
                if transmat_temp(r,r) < 0
                    warning('InitialGuess_v2_refined: transmat0(%d,%d)=%.2e negative. Clamping and re-normalizing row %d.',r,r,transmat_temp(r,r),r);
                    % Reset off-diagonals for this row and make it strongly diagonal
                    transmat_temp(r,:) = 0;
                    transmat_temp(r,r) = 1.0;
                    if Nstates > 1 % Add a tiny probability to escape to one other state
                        other_state_idx = mod(r, Nstates) + 1; % pick another state
                        transmat_temp(r, other_state_idx) = 1e-5; % small escape probability
                        transmat_temp(r,r) = 1.0 - 1e-5;
                    end
                end
            end % End loop for row r
        end % End if Nstates == 1 else ...

        transi = transmat_temp ./ sum(transmat_temp, 2); % Final normalization for safety

        % --- Prior Vector (prior0) ---
        if Nstates > 0
            temp_prior = rand(1, Nstates) + 1e-9;
            priori = temp_prior / sum(temp_prior);
        else
            priori = [];
        end

        % --- Observation Matrix (obsmat0) ---
        if Nstates > 0 && O_channels > 0
            obsmat_temp = rand(Nstates, O_channels) + 1e-9;
            if Nstates <= O_channels
                for s_idx = 1:Nstates
                    preferred_channel = s_idx;
                    obsmat_temp(s_idx, preferred_channel) = obsmat_temp(s_idx, preferred_channel) * (O_channels * 2) + rand();
                end
            else
                 for s_idx = 1:Nstates
                    preferred_channel = mod(s_idx - 1, O_channels) + 1;
                    obsmat_temp(s_idx, preferred_channel) = obsmat_temp(s_idx, preferred_channel) * (Nstates) + rand();
                end
            end
            obsmati = obsmat_temp ./ sum(obsmat_temp, 2);
        else
            obsmati = [];
        end

        InitialGues{iii}.transmat0 = transi;
        InitialGues{iii}.prior0 = priori;
        InitialGues{iii}.obsmat0 = obsmati;
        InitialGues{iii}.Nstate = Nstates;
    end
end
% --- END OF FILE InitialGuess_PhotonByPhoton_v2.m (with refinement) ---% --- START OF FILE InitialGuess_PhotonByPhoton_v2.m ---
function [InitialGues] = InitialGuess_PhotonByPhoton_v2(number_of_initialguess_per_job, Nstates, O_channels, dt_analysis_sec)
% This function constructs an initial guess array for HMM analysis.
% Enhanced version:
% - Accepts O_channels for obsmat0 generation.
% - Accepts dt_analysis_sec for rate-aware transmat0 generation.

    % --- Parameters for Rate-Aware transmat0 ---
    % These define a plausible range for physical kinetic rates (transitions per second)
    % Adjust these based on typical expectations for the system being analyzed.
    k_min_Hz = 1e2;   % Minimum rate (e.g., 100 Hz)
    k_max_Hz = 1e7;   % Maximum rate (e.g., 10 MHz)
    % Target sum of off-diagonal probabilities in a row for transmat0.
    % Should be < 1 to leave room for diagonal. A smaller value means more persistent states.
    target_sum_off_diagonal_P = 0.1; % e.g., sum of P_ij (j~=i) aims for 0.1

    InitialGues = cell(1, number_of_initialguess_per_job); % Pre-allocate



    for iii = 1:number_of_initialguess_per_job
        transmat_temp = zeros(Nstates, Nstates);
        if Nstates == 0 % Handle empty case
            InitialGues{iii}.transmat0 = [];
            InitialGues{iii}.prior0 = [];
            InitialGues{iii}.obsmat0 = [];
            InitialGues{iii}.Nstate = 0;
            continue;
        end

        if Nstates == 1
            transmat_temp = 1; % Only one state, must stay there
        else % Nstates > 1
            for r = 1:Nstates % Current row (from state r)
                row_off_diagonal_P = zeros(1, Nstates); % Store P_rc for this row r
                raw_rates_for_scaling = zeros(1, Nstates); % Store k_guess_Hz for scaling if needed

                for c = 1:Nstates % Current col (to state c)
                    if r == c, continue; end

                    log_k_min = log10(k_min_Hz);
                    log_k_max = log10(k_max_Hz);
                    k_guess_Hz = 10^(log_k_min + rand() * (log_k_max - log_k_min));
                    
                    raw_rates_for_scaling(c) = k_guess_Hz; % Store the rate itself
                    row_off_diagonal_P(c) = k_guess_Hz * dt_analysis_sec;
                end

                % --- Scaling logic for off-diagonal probabilities ---
                current_sum_off_diag_P = sum(row_off_diagonal_P);

                if abs(current_sum_off_diag_P) < 1e-9 % If sum of P_ij is effectively zero (all rates*dt were tiny)
                    % Distribute target_sum_off_diagonal_P equally among off-diagonal elements
                    % This happens if k_max_Hz * dt_analysis_sec is still very small for all sampled k.
                    num_off_diag_elements = Nstates - 1;
                    if num_off_diag_elements > 0
                        scaled_P_rc = target_sum_off_diagonal_P / num_off_diag_elements;
                        for c_assign = 1:Nstates
                            if r == c_assign, continue; end
                            transmat_temp(r, c_assign) = scaled_P_rc;
                        end
                    end
                else
                    % Scale based on the ratio of raw rates to ensure P_ij sums to target_sum_off_diagonal_P
                    % This preserves the relative magnitudes of the initially sampled rates.
                    sum_raw_rates_for_row = sum(raw_rates_for_scaling); % Sum of k_ij for j~=i
                    if abs(sum_raw_rates_for_row) < 1e-9 % Should not happen if k_min_Hz > 0
                         % Fallback: equal distribution if all raw rates were somehow zero
                        num_off_diag_elements = Nstates - 1;
                        scaled_P_rc = target_sum_off_diagonal_P / num_off_diag_elements;
                        for c_assign = 1:Nstates
                            if r == c_assign, continue; end
                            transmat_temp(r, c_assign) = scaled_P_rc;
                        end
                    else
                        for c_assign = 1:Nstates
                            if r == c_assign, continue; end
                            % P_rc = (k_rc / sum(k_rj)) * target_sum_off_diagonal_P
                            transmat_temp(r, c_assign) = (raw_rates_for_scaling(c_assign) / sum_raw_rates_for_row) * target_sum_off_diagonal_P;
                        end
                    end
                end
                % --- End Scaling Logic ---
                
                % Set diagonal element
                current_row_sum_off_diag = sum(transmat_temp(r, find((1:Nstates) ~= r) ));
                transmat_temp(r, r) = 1 - current_row_sum_off_diag;

                % Safety net for diagonal (should ideally not be needed with target_sum < 1)
                if transmat_temp(r,r) < 0
                    warning('InitialGuess_v2_refined: transmat0(%d,%d)=%.2e negative. Clamping and re-normalizing row %d.',r,r,transmat_temp(r,r),r);
                    % Reset off-diagonals for this row and make it strongly diagonal
                    transmat_temp(r,:) = 0;
                    transmat_temp(r,r) = 1.0;
                    if Nstates > 1 % Add a tiny probability to escape to one other state
                        other_state_idx = mod(r, Nstates) + 1; % pick another state
                        transmat_temp(r, other_state_idx) = 1e-5; % small escape probability
                        transmat_temp(r,r) = 1.0 - 1e-5;
                    end
                end
            end % End loop for row r
        end % End if Nstates == 1 else ...

        transi = transmat_temp ./ sum(transmat_temp, 2); % Final normalization for safety

        % --- Prior Vector (prior0) ---
        if Nstates > 0
            temp_prior = rand(1, Nstates) + 1e-9;
            priori = temp_prior / sum(temp_prior);
        else
            priori = [];
        end

        % --- Observation Matrix (obsmat0) ---
        if Nstates > 0 && O_channels > 0
            obsmat_temp = rand(Nstates, O_channels) + 1e-9;
            if Nstates <= O_channels
                for s_idx = 1:Nstates
                    preferred_channel = s_idx;
                    obsmat_temp(s_idx, preferred_channel) = obsmat_temp(s_idx, preferred_channel) * (O_channels * 2) + rand();
                end
            else
                 for s_idx = 1:Nstates
                    preferred_channel = mod(s_idx - 1, O_channels) + 1;
                    obsmat_temp(s_idx, preferred_channel) = obsmat_temp(s_idx, preferred_channel) * (Nstates) + rand();
                end
            end
            obsmati = obsmat_temp ./ sum(obsmat_temp, 2);
        else
            obsmati = [];
        end

        InitialGues{iii}.transmat0 = transi;
        InitialGues{iii}.prior0 = priori;
        InitialGues{iii}.obsmat0 = obsmati;
        InitialGues{iii}.Nstate = Nstates;
    end
end
% --- END OF FILE InitialGuess_PhotonByPhoton_v2.m (with refinement) ---