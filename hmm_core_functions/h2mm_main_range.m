function [LL, prior, transmat, obsmat, lastIteration,all_gamma] = ...
    h2mm_main_range(data, prior, transmat, obsmat, max_iter,  MaxRunTime,OBSFix,LastData,fileName,par,fix,range_FRET)
%% Define some option variables
% See explenations above
thresh=1e-14; % Precision threshold
verbose=1; % Print output
obs_prior_weight=0;
adj_prior=1;
adj_trans=1;
adj_obs=OBSFix;

%% Create sorted list of inter-photon arrival times
TAU=1;
% data=varargin{1,2};
len=length(data);
lendata=0;

for ii=1:len
    lendata=lendata+size(data{ii},1);
end

TotalArrivalDelta=zeros(lendata,1);
templ=1;

for ii=1:len
    TotalArrivalDelta(templ:size(data{ii},1)+templ-2)=diff(data{ii}(:,1));
    templ=length(TotalArrivalDelta)+1;
end

TotalArrivalDelta=sort(unique(TotalArrivalDelta));
TotalArrivalDelta=TotalArrivalDelta(2:end)';

%% Perform factorization inter-photon arrival times 
R = Factorization(TotalArrivalDelta,1);

%% Per inter-photon arrival time in the trajectory, find the matching index in TotalArrivalDelta 
data_red{length(data)}=[];

for i=1:length(data)
    
    data_red{i}=data{i};
    temp=diff(data{i}(:,1));
    temp(temp==0)=1;
    temp_dt=temp;
    data_red{i}(1,1)=0;
    
    for t=2:length(data{i}(:,1))
        data_red{i}(t,1)=find(TotalArrivalDelta==temp_dt(t-1));
    end
end

clear data
data=data_red;
clear data_red

%% Collected info from last run

if LastData.Iteration==0
    
    previous_loglik = -inf; 
    num_iter = 1;
    LL = -Inf(1,max_iter); % This is recommanded for memory allocation
    interation_times=zeros(1,max_iter);
    
    LastData.LL=LL;
    LastData.InteruptCount=0;   
else
    num_iter=LastData.Iteration+1;% Notice the + 1
    previous_loglik = LastData.LL(LastData.Iteration);
    LL=LastData.LL;
    interation_times=LastData.interation_times;
    LastData.InteruptCount=LastData.InteruptCount+1;%How many time run was interupted

end
converged = 0;

%% Begin H2MM on current model


if ~iscell(data)
    data = num2cell(data, 2); % each row gets its own cell
end


tic;
while (num_iter <= max_iter) && ~converged
    
    current_time=toc;
%     prior=varargin{1,6}{1,1}.prior0;
%     transmat=varargin{1,6}{1,1}.transmat0;
%     obsmat=varargin{1,6}{1,1}.obsmat0;
%     par=0;
    %% E step
    [loglik, exp_num_trans, exp_num_visits1, exp_num_emit,exp_num_visitsT,all_gamma] = ...
        compute_ess_dhmm(prior, transmat, obsmat, data, obs_prior_weight,TotalArrivalDelta,R,par);
    
    %% M step
    
    if adj_prior
        prior = exp_num_visits1/sum(exp_num_visits1);
        LastData.PriorList=[LastData.PriorList ; prior];
    end
    
    if adj_trans && ~isempty(exp_num_trans)
        transmat = mk_stochastic(exp_num_trans);
        LastData.transmatList=[LastData.transmatList ; transmat(:)'];
    end
    
    if adj_obs
        obsmat = mk_stochastic(exp_num_emit);
        %bring the respective state back into the range
        if obsmat(fix,2)<range_FRET(1)
            obsmat(fix,2)=range_FRET(1);obsmat(fix,1)=1-range_FRET(1);
        end
        if obsmat(fix,2)>range_FRET(2)
            obsmat(fix,2)=range_FRET(2);obsmat(fix,1)=1-range_FRET(2);
        end
        LastData.ObsList=[LastData.ObsList ; obsmat(:)'];
    end
    if verbose
        t=toc-current_time;
        fprintf(1, '%s -iteration %d, loglik = %f. Took: %d seconds\n',datestr(now, 'HH:MM:SS'), num_iter, loglik,t)                      
    end
    
    converged = em_converged(loglik, previous_loglik,thresh);
    if converged 
        fprintf('Converged, delta log-likelihood exceeded threshold of %d\n',thresh)
    
    else
        interation_times(num_iter)=toc-current_time;
        current_time=toc;
        if (sum(interation_times)/3600>MaxRunTime)
            converged=1;
            fprintf('Converged, run time exceeded threshold of %d hours\n',MaxRunTime)
        end
    end
    
    previous_loglik = loglik;
    LL(num_iter) =loglik;
    
    % update LastData
    LastData.Iteration=num_iter;
    LastData.LL=LL;
    LastData.last_obsmat=obsmat;
    LastData.last_prior=prior;
    LastData.last_transmat=transmat;
    LastData.all_gamma=all_gamma;
    LastData.interation_times=interation_times;
    save(fileName,'LastData')   
    
    num_iter =  num_iter + 1;

end
if num_iter>max_iter
    fprintf('Converged, reached %d iteration threshold \n',max_iter)
end
lastIteration = num_iter - 1;
LL=LL(LL~=-inf);
fprintf('Reached end of while loop\n')


function [loglik, exp_num_trans, exp_num_visits1, exp_num_emit, exp_num_visitsT,all_gamma] = ...
    compute_ess_dhmm(startprob, transmat, obsmat, data, dirichlet,TotalArrivalDelta, R,par)
% COMPUTE_ESS_DHMM Compute the Expected Sufficient Statistics for an HMM with discrete outputs
% function [loglik, exp_num_trans, exp_num_visits1, exp_num_emit, exp_num_visitsT] = ...
%    compute_ess_dhmm(startprob, transmat, obsmat, data, dirichlet)
%
% INPUTS:
% startprob(i)
% transmat(i,j)
% obsmat(i,o)
% data{seq}(t)
% dirichlet - weighting term for uniform dirichlet prior on expected emissions
%
% OUTPUTS:
% exp_num_trans(i,j) = sum_l sum_{t=2}^T Pr(X(t-1) = i, X(t) = j| Obs(l))
% exp_num_visits1(i) = sum_l Pr(X(1)=i | Obs(l))
% exp_num_visitsT(i) = sum_l Pr(X(T)=i | Obs(l))
% exp_num_emit(i,o) = sum_l sum_{t=1}^T Pr(X(t) = i, O(t)=o| Obs(l))
% where Obs(l) = O_1 .. O_T for sequence l.


TotalArrivalDelta_max=length(TotalArrivalDelta);
numex = length(data);
[S ,O] = size(obsmat);
exp_num_trans = zeros(S,S);
exp_num_visits1 = zeros(S,1);
exp_num_visitsT = zeros(S,1);

for ex=1:numex
    exp_num_emits{ex} = dirichlet*ones(S,O);
end

loglik = 0;
alpha_free{TotalArrivalDelta_max}(1:S,1:S)=0;

for t=1:TotalArrivalDelta_max-1
    alpha_free{t}(1:S,1:S)=0;
end


transmat_t(1:S,1:S,1:TotalArrivalDelta_max)=0;
transmat_t(:,:,1)=transmat;
transmat_t = CalculatePowerOfTransMatrices(R,transmat);

[Rho] = Calc_Rho(transmat_t,R);

if par==1
%% Parallel for    
parfor ex=1:numex
    obs = data{ex}(:,2);
    delta_t = data{ex}(:,1);
    delta_t(delta_t==0)=1;% correction
    ex;
    T = length(obs);
    obslik = multinomial_prob(obs, obsmat);
    
    % Forward backwards
    [~,~, gamma, current_ll, xi_summed] = fwdback_photonByphoton_fast(startprob, transmat, obslik, delta_t, Rho, transmat_t);
    all_gamma{ex}=gamma;
    loglik = loglik +  current_ll;
    exp_num_trans = exp_num_trans + xi_summed;
    exp_num_visits1 = exp_num_visits1 + gamma(:,1);
        % exp_num_visitsT = exp_num_visitsT + gamma(:,T);   
    
    if T < O
        for t=1:T
            exp_num_emits{ex}(:,obs(t)) = exp_num_emits{ex}(:,obs(t)) + gamma(:,t);
        end
    else
        for o=1:O
            ndx = find(obs==o);
            if ~isempty(ndx)
                exp_num_emits{ex}(:,o) = exp_num_emits{ex}(:,o) + sum(gamma(:,ndx), 2);
            end
        end
    end
end
%% Regular loop
else
    for ex=1:numex
        ex=55;
    obs = data{ex}(:,2);
    delta_t = data{ex}(:,1);
    delta_t(delta_t==0)=1;% correction
    
    T = length(obs);
    obslik = multinomial_prob(obs, obsmat);
    
    % Forward backwards
    [~,~, gamma, current_ll, xi_summed] = fwdback_photonByphoton_fast(startprob, transmat, obslik, delta_t, Rho, transmat_t);
    
    loglik = loglik +  current_ll;
    exp_num_trans = exp_num_trans + xi_summed;
    exp_num_visits1 = exp_num_visits1 + gamma(:,1);
    all_gamma{ex}=gamma;
    % exp_num_visitsT = exp_num_visitsT + gamma(:,T);   
    
    if T < O
        for t=1:T
            exp_num_emits{ex}(:,obs(t)) = exp_num_emits{ex}(:,obs(t)) + gamma(:,t);
        end
    else
        for o=1:O
            ndx = find(obs==o);
            if ~isempty(ndx)
                exp_num_emits{ex}(:,o) = exp_num_emits{ex}(:,o) + sum(gamma(:,ndx), 2);
            end
        end
    end
    end
end

% Summing emmition probability measure
exp_num_emit=exp_num_emits{1};
for ex=2:numex
    exp_num_emit=exp_num_emit+exp_num_emits{ex};
end
