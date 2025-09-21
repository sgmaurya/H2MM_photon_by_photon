% --- START OF FILE InitialGuess_FromPilot.m ---
function [InitialGuesMain] = InitialGuess_FromPilot(best_pilot_model_params, num_guesses_main, Nstates, O_channels, dt_analysis_sec, perturb_config)
    % Generates initial guesses for the main HMM run based on a pilot model,
    % by perturbing the pilot model's parameters.
    %
    % INPUTS:
    %   best_pilot_model_params: Struct with .prior, .transmat, .obsmat from pilot.
    %   num_guesses_main:        Number of initial guesses to generate for the main run.
    %   Nstates:                 Number of states.
    %   O_channels:              Number of observation channels.
    %   dt_analysis_sec:         Time resolution, for converting P to k if needed.
    %   perturb_config:          Struct with perturbation parameters, e.g.,
    %                            .rate_noise_factor (for transmat)
    %                            .obs_noise_level (for obsmat)
    %                            .prior_noise_level (for prior)
    %
    % OUTPUTS:
    %   InitialGuesMain:         Cell array of initial guess structs for the main run.

    InitialGuesMain = cell(1, num_guesses_main);

    if Nstates == 0
        warning('InitialGuess_FromPilot: Nstates is 0, cannot generate guesses.');
        for i = 1:num_guesses_main
            InitialGuesMain{i}.transmat0 = [];
            InitialGuesMain{i}.prior0 = [];
            InitialGuesMain{i}.obsmat0 = [];
            InitialGuesMain{i}.Nstate = 0;
        end
        return;
    end
    
    % --- Default Perturbation Config (if not fully provided) ---
    if ~isfield(perturb_config, 'rate_noise_factor'), perturb_config.rate_noise_factor = 0.25; end
    if ~isfield(perturb_config, 'obs_noise_level'), perturb_config.obs_noise_level = 0.1; end % Additive noise magnitude
    if ~isfield(perturb_config, 'prior_noise_level'), perturb_config.prior_noise_level = 0.1; end % Additive noise magnitude

    fprintf('INFO (InitialGuess_FromPilot): Generating %d guesses based on pilot model.\n', num_guesses_main);
    fprintf('  Perturbation: Rate Noise Factor=%.2f, Obs Noise=%.2f, Prior Noise=%.2f\n', ...
        perturb_config.rate_noise_factor, perturb_config.obs_noise_level, perturb_config.prior_noise_level);

    % --- 1. The pilot model itself is the first guess ---
    InitialGuesMain{1}.prior0 = best_pilot_model_params.prior;
    InitialGuesMain{1}.transmat0 = best_pilot_model_params.transmat;
    InitialGuesMain{1}.obsmat0 = best_pilot_model_params.obsmat;
    InitialGuesMain{1}.Nstate = Nstates;

    % --- 2. Generate perturbed guesses for the rest ---
    for i_guess = 2:num_guesses_main
        % --- Perturb Prior ---
        perturbed_prior = best_pilot_model_params.prior + randn(1, Nstates) * perturb_config.prior_noise_level;
        perturbed_prior(perturbed_prior < 1e-6) = 1e-6; % Ensure non-negative, floor at small value
        perturbed_prior = perturbed_prior / sum(perturbed_prior); % Normalize

        % --- Perturb Obsmat ---
        perturbed_obsmat = best_pilot_model_params.obsmat + randn(Nstates, O_channels) * perturb_config.obs_noise_level;
        perturbed_obsmat(perturbed_obsmat < 1e-6) = 1e-6; % Ensure non-negative
        perturbed_obsmat = perturbed_obsmat ./ sum(perturbed_obsmat, 2); % Normalize rows

        % --- Perturb Transmat (more involved: rates -> perturb rates -> probabilities) ---
        perturbed_transmat = best_pilot_model_params.transmat; % Start with pilot
        if Nstates > 1
            try
                % Convert P_ij from pilot to effective rates Q_ij (Q is rate matrix)
                % Q = logm(P) / dt. Off-diagonals of Q are k_ij.
                % Add a small epsilon to diagonal before logm to ensure P is non-singular if it's too identity-like
                P_pilot_reg = best_pilot_model_params.transmat;
                for r=1:Nstates, P_pilot_reg(r,r) = P_pilot_reg(r,r) + 1e-9 * rand(); end
                P_pilot_reg = P_pilot_reg ./ sum(P_pilot_reg,2);

                Q_matrix_pilot = logm(P_pilot_reg) / dt_analysis_sec;
                
                Q_perturbed = zeros(Nstates, Nstates);
                for r = 1:Nstates
                    sum_perturbed_k_row = 0;
                    for c = 1:Nstates
                        if r == c, continue; end
                        k_rc_pilot = Q_matrix_pilot(r,c);
                        if k_rc_pilot < 0, k_rc_pilot = 1e-9; end % logm can give small negatives for P_rc near 0

                        % Perturb multiplicatively in log-space for rates (ensures positive)
                        k_rc_perturbed = k_rc_pilot * exp(randn() * perturb_config.rate_noise_factor);
                        k_rc_perturbed = max(k_rc_perturbed, 1e-9); % Floor at a very small rate
                        
                        Q_perturbed(r,c) = k_rc_perturbed;
                        sum_perturbed_k_row = sum_perturbed_k_row + k_rc_perturbed;
                    end
                    Q_perturbed(r,r) = -sum_perturbed_k_row; % Diagonal of Q matrix
                end
                
                % Convert perturbed Q back to P: P = expm(Q*dt)
                perturbed_transmat_candidate = expm(Q_perturbed * dt_analysis_sec);
                
                % Ensure stochasticity (due to numerical precision of expm and logm)
                perturbed_transmat_candidate(perturbed_transmat_candidate < 0) = 0; % Remove negatives
                perturbed_transmat = perturbed_transmat_candidate ./ sum(perturbed_transmat_candidate, 2);

                % Final check for NaNs or Infs if expm/logm had issues
                if any(isnan(perturbed_transmat(:))) || any(isinf(perturbed_transmat(:)))
                    warning('InitialGuess_FromPilot: NaN/Inf in perturbed transmat for guess %d. Using pilot transmat.', i_guess);
                    perturbed_transmat = best_pilot_model_params.transmat;
                end

            catch ME_trans_perturb
                warning('InitialGuess_FromPilot: Error perturbing transmat for guess %d (msg: %s). Using pilot transmat.', i_guess, ME_trans_perturb.message);
                perturbed_transmat = best_pilot_model_params.transmat; % Fallback
            end
        end % if Nstates > 1

        InitialGuesMain{i_guess}.prior0 = perturbed_prior;
        InitialGuesMain{i_guess}.transmat0 = perturbed_transmat;
        InitialGuesMain{i_guess}.obsmat0 = perturbed_obsmat;
        InitialGuesMain{i_guess}.Nstate = Nstates;
    end
end
% --- END OF FILE InitialGuess_FromPilot.m ---