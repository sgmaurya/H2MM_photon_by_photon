% --- START OF NEW FILE fwdback_photonByphoton_fast_simplified.m ---
% This is the new, high-performance fwdback script.
% It is used by the Direct Calculation pathway for dense high-flux data.
% It does NOT use the complex Rho matrix.

function [alpha, beta, gamma, loglik, xi_summed, gamma2] = fwdback_photonByphoton_fast_simplified(init_state_distrib, ...
    transmat_base, obslik, delta_t_indices, time_steps_in_units, transmat_t, varargin)
% INPUTS:
% transmat_base: The base transition matrix for a single time unit (not used directly here, but good practice)
% obslik(i,t): Observation likelihood P(y_t | Q_t=i)
% delta_t_indices(t): The index into the TotalArrivalDelta vector for time t.
% time_steps_in_units: The mapping from a delta_t_index to an integer number of time steps.
% transmat_t: The pre-calculated lookup table of transition matrices.

[Q, T] = size(obslik);
scaled = 1; % Always use scaling for numerical stability

scale = ones(1,T);
alpha = zeros(Q,T);
beta = zeros(Q,T);
gamma = zeros(Q,T);
xi_summed = zeros(Q,Q);
gamma2 = []; % Not implemented in this version

%% Forward pass
alpha(:,1) = init_state_distrib(:) .* obslik(:,1);
if scaled
    [alpha(:,1), scale(1)] = normalise(alpha(:,1));
end

for t=2:T
    % --- SIMPLIFIED TRANSITION MATRIX LOOKUP ---
    time_index = delta_t_indices(t);
    step_unit = time_steps_in_units(time_index);
    trans = transmat_t(:,:,step_unit);
    
    m = trans' * alpha(:,t-1);
    alpha(:,t) = m(:) .* obslik(:,t);
    
    if scaled
        [alpha(:,t), scale(t)] = normalise(alpha(:,t));
    end
end

if any(scale==0)
    loglik = -inf;
else
    loglik = sum(log(scale));
end

%% Backward pass
beta(:,T) = ones(Q,1);
gamma(:,T) = normalise(alpha(:,T) .* beta(:,T));

for t = T-1:-1:1
    b = beta(:,t+1) .* obslik(:,t+1);
    
    % --- SIMPLIFIED TRANSITION MATRIX LOOKUP ---
    time_index = delta_t_indices(t+1);
    step_unit = time_steps_in_units(time_index);
    trans = transmat_t(:,:,step_unit);
    
    beta(:,t) = trans * b;
    
    if scaled
        beta(:,t) = normalise(beta(:,t));
    end
    
    gamma(:,t) = normalise(alpha(:,t) .* beta(:,t));
    
    % --- SIMPLIFIED XI CALCULATION (NO RHO) ---
    xi_temp = normalise((trans .* (alpha(:,t) * b')));
    xi_summed = xi_summed + xi_temp;
end

end

function [M, z] = normalise(A, dim)
% Make the entries of a (multidimensional) array sum to 1
if nargin < 2
    z = sum(A(:));
    if z==0, z=1; end % handle the case where the vector is all zeros
    M = A/z;
elseif dim==1 % normalize each column
    z = sum(A,1);
    z(z==0) = 1;
    M = A ./ repmat(z, [size(A,1) 1]);
else % normalize each row
    z = sum(A,2);
    z(z==0) = 1;
    M = A ./ repmat(z, [1 size(A,2)]);
end
end
% --- END OF NEW FILE fwdback_photonByphoton_fast_simplified.m ---