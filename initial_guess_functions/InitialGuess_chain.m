%% InitialGuess_chain
function [InitialGuess] = InitialGuess_chain(number_of_initialguess,Nstates)
% This function construct an initial guess array for HMM analysis The
% initial guess array contains for each one the HMM parameters: -
% Transition matrix - Pi vector - Observation matrix - and number of states

iii=0;
translist = [];% Used to plot matrix elements. Added by Demian

for state=Nstates

    for ii=1:number_of_initialguess
        iii=iii+1;
        % Added by Demian:24.05.20
        % Randomize a raw stochastic matrix: The expo variable determines
        % around which values the random off-diagonal probabilities will be
        % centered. The parameter expo determines the mean value. For
        % instance, for an experiment with resolution of 5e-8 sec, an expo
        % value of 5 will yield off diagonal probabilities, corresponding
        % to off diagonal kinetic rates around 10^5 Hz. Plot a log
        % histogram of the rates to convince yourself.

        temp=rand(1,state);
        priori=temp/sum(temp); % Pi vector

        temp=sort(rand(state,1)); % Sorted
        obsmati=[1-temp temp]; % Observation matrix

        F=rand(state,state);
        % default value for expo is 1e5
        expo = 1e5;

        % Here, you can set a range (small number corresponds to larger
        % kinetic rates). Adjust the values to your experiment and sample.
        expo = (10.^randi([3 6],state,1)); % random order of magnitude
        F=F-tril(F,-2)-triu(F,2); % making chain architecture
        F=F+diag(diag(F).*expo); % making diagonal elements large
        transi=diag(sum(F,2))^-1*(F); % row normalization



        % Add random guesses to array
        InitialGuess{iii}.transmat0=transi;
        InitialGuess{iii}.prior0=priori;
        InitialGuess{iii}.obsmat0=obsmati;
        InitialGuess{iii}.Nstate=state;

        translist = [translist, transi]; % Used to plot matrix elements. Added by Demian
    end
end
% Added by Demian. Use this section to plot the distribution of kinetic
% rates for your initial guesses.

% dt = 0.5*1e-7; % set dt resolution of your experiment
% histogram(log10(translist(:)/dt),50)
% xlabel('log_{10}(Rate/Hz)')
% ylabel('Counts')

end