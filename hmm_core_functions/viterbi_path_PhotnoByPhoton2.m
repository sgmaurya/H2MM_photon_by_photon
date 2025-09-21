% viterbi_path_PhotnoByPhoton2.m
% Corrected to calculate transition matrices on-the-fly for each step.

function path = viterbi_path_PhotnoByPhoton2(prior, transmat, obslik, PhotonTimeSteps)
% VITERBI Find the most-probable (Viterbi) path through the HMM state trellis.
%
% Inputs:
% prior(i) = Pr(Q(1) = i)
% transmat(i,j) = Pr(Q(t+1)=j | Q(t)=i) FOR A UNIT TIME STEP (dt_analysis_sec)
% obslik(i,t) = Pr(y(t) | Q(t)=i) for the t-th photon in the current burst
% PhotonTimeSteps(t-1) = number of unit time steps between photon t-1 and photon t
%
% Outputs:
% path(t) = q(t), where q1 ... qT is the most probable state sequence for the photons.

% PhotonTimeSteps contains the number of fundamental time steps (dt_analysis_sec) 
% between photons. It can be empty if only one photon in the burst.

scaled = 1; % Use scaling to prevent underflow

T = size(obslik, 2); % Number of photons in the current burst
prior = prior(:);    % Ensure prior is a column vector
Q = length(prior);   % Number of hidden states

delta = zeros(Q,T);  % delta(j,t) = prob. of the best sequence up to photon t, ending in state j
psi = zeros(Q,T);    % psi(j,t) = the best predecessor state for state j at photon t
path = zeros(1,T);   % Stores the Viterbi path

if T == 0
    path = []; % No photons, no path
    return;
end

% Initialization (t=1, first photon)
delta(:,1) = prior .* obslik(:,1);
if scaled
    [delta(:,1), n] = normalise(delta(:,1)); % normalise is a BNT utility
    % scale(1) = 1/n; % Not strictly needed for Viterbi path, only for loglik
end
psi(:,1) = 0; % Arbitrary, no predecessor

% Recursion (t=2 to T)
for t = 2:T
    num_steps = PhotonTimeSteps(t-1); % Number of unit time steps between photon t-1 and t
                                      % This should be a positive integer (already handled in runViterbi_local)
    
    if num_steps < 1
        warning('Viterbi: PhotonTimeSteps contained non-positive value %d, using 1.', num_steps);
        num_steps = 1;
    end

    % Calculate transmat_for_duration = transmat ^ num_steps
    % MATLAB's A^k is efficient for integer k (uses exponentiation by squaring)
    transmat_for_duration = transmat ^ num_steps;
    
    for j = 1:Q % For each current state j at photon t
        % Find the best previous state i that maximizes the path to state j
        [max_val, best_prev_state_idx] = max(delta(:,t-1) .* transmat_for_duration(:,j));
        delta(j,t) = max_val;
        psi(j,t) = best_prev_state_idx;
    end
    delta(:,t) = delta(:,t) .* obslik(:,t); % Include observation likelihood for current photon
    
    if scaled
        [delta(:,t), n] = normalise(delta(:,t));
        % scale(t) = 1/n;
    end
end

% Termination: Find the best state for the last photon
[~, path(T)] = max(delta(:,T));

% Path backtracking
for t = T-1:-1:1
    path(t) = psi(path(t+1), t+1);
end

% loglik calculation (optional, not the primary output of Viterbi path)
% if scaled
%     loglik = -sum(log(scale(scale>0))); % Avoid log(0) if any scale factor was 0
% else
%     [p_final, ~] = max(delta(:,T)); % Probability of the Viterbi path
%     if p_final > 0
%         loglik = log(p_final);
%     else
%         loglik = -inf;
%     end
% end
% fprintf('Viterbi path log-likelihood (approx): %f\n', loglik);

end