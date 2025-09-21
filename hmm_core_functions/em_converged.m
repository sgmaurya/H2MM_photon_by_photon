function [converged] = em_converged(loglik, previous_loglik,thresh)
% EM_CONVERGED Has EM converged?
% [converged, decrease] = em_converged(loglik, previous_loglik, )
%
% We have converged if the slope of the log-likelihood is negative',
% i.e., loglik - previous_loglik < 0,

converged = 0;

if loglik - previous_loglik <= thresh
    converged = 1;
end
end

