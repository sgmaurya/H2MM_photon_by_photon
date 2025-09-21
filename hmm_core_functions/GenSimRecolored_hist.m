function [S, photon_color,Recolored_FRET] = GenSimRecolored_hist(data,model,Realization_number)
% GenSimRecolored_hist uses a simple simulation of HMM model to generate a
% photon color realization similar to the one described in the reference
% below and generates a simulated fret histogram. Ref: Gopich, Irina V.,
% and Attila Szabo. "Decoding the pattern of photon colors in
% single-molecule FRET." The Journal of Physical Chemistry B 113.31 (2009):
% 10965-10973. inputs: model is a stucture containing all the HMM model
% parameters: model.Nstates is the number of states. model.obsmat is the
% model state emission matrix. model.prior is the model initial state
% probability vector. model.transmat is the model state to state transition
% probability matrix. data is a two coloumn array the first column is a
% list of photon arrival times and the second coloumn is the photon color.
% X is a vectore of histogram bin positions.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs: Recolored_hist is a vector of model simulated burst fret
% histogram

if ~exist('Realization_number')
    Realization_number = 1;
end

% Simulate photon trajectories
for nr=1:Realization_number
    
    n_temp = model.Nstates;
    obsmat_temp = model.obsmat;
    prior_temp = model.prior;
    transmat_temp = model.transmat;
    
    for k = 1:length(data)
        traj_temp = diff(data{k}(:,1));
        photon_color{nr}{k }= zeros([length(traj_temp)+1,1]);
        S_temp = OneMonteCarloChoise(prior_temp);
        photon_color{nr}{k}(1) = OneMonteCarloChoise(obsmat_temp(S_temp,:));
        S{nr}{k }(1)=S_temp;
        for t = 1:length(traj_temp)
            A_temp = transmat_temp^traj_temp(t);
            S_temp = OneMonteCarloChoise(A_temp(S_temp,:));
            photon_color{nr}{k}(t+1) = OneMonteCarloChoise(obsmat_temp(S_temp,:));
            S{nr}{k }(t+1)=S_temp;
        end
    end
    
    Recolored_FRET{nr}= arrayfun(@calcFRET,photon_color{nr}); % Get FRET from data
end

end
%% mini functions
function fret = calcFRET(data)
data=cell2mat(data);
fret = nnz(data==2)/(nnz(data==2)+nnz(data==1));
end