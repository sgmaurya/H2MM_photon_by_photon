%% InitialGuess chain seperated
function [InitialGuess] = InitialGuess_chain_seperated(number_of_initialguess,Nstates)
% This function construct an initial guess array for HMM analysis The
% initial guess array contains for each one the HMM parameters: -
% Transition matrix - Pi vector - Observation matrix - and number of states

% For blink state, construct a model that accounts for 3 color FRET and
% also for Dark states
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

        % Each state has a bright and dark state, therefore, there the
        % number of state is always double. The dark states are always at
        % the end the arrays.

        % Pi vector :
        temp=rand(1,state+1);
        temp(end) =0;  %This is to ensure starting of bright state
        priori=temp/sum(temp);

        % Observation matrix:
        % The observation matrix has 3 colors. Column 1 for donor from
        % donor pulse, column 2 is for acceptor from donor pulse and column
        % 3 is for acceptor from acceptor pulse.

        temp = rand(state+1,3);
        [~,I]=sort(temp(:,2)./(temp(:,1)+temp(:,2)));
        temp=temp(I,:);
        temp(state+1:end,3)=0; % imposed acceptor dark state
        temp(state+1:end,1)=temp(state+1:end,1)+1000; % bias donor emission
        obsmati=temp./sum(temp,2); % normalize

        % Transition matrix:
        % One dark state for all conformational states
        template=rand(state+1,state+1);
        F=rand(state,state);
        % optional:
        F=F-tril(F,-2)-triu(F,2); % making chain architecture

        template(1:state,1:state)=F;

        % Here, you can set a range (small number corresponds to larger
        % kinetic rates). Adjust the values to your experiment and sample.
        expo = 1e5;% default value for expo is 1e5
        expo = (10.^randi([3 6],state+1,1)); % random order of magnitude

        template=template+diag(diag(template).*expo); % making diagonal elements large
        transi=diag(sum(template,2))^-1*(template); % row normalization

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
