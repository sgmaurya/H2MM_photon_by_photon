%% InitialGuess 4 colors Toy Model
function [InitialGuess] = InitialGuess_4colorsToyModel(number_of_initialguess,Nstates)
% This function construct an initial guess array for HMM analysis The
% initial guess array contains for each one the HMM parameters: -
% Transition matrix - Pi vector - Observation matrix - and number of states

% for 4 colors encoding for donor and acceptor from two subunits on a
% protein. Used for multicolor FRET.

iii=0;
translist = [];% Used to plot matrix elements. Added by Demian


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
    temp=rand(1,Nstates);
    priori=temp/sum(temp);

    % Observation matrix:
    % The observation matrix has 4 colors. For the toy model
    % OO-OC-CO-CC, shared FRET efficeinciey.

    FRETs = rand(2,2);
    FRET1 = [FRETs(1,1) FRETs(1,1) FRETs(1,2) FRETs(1,2)];
    FRET2 = [FRETs(2,1) FRETs(2,2) FRETs(2,1) FRETs(2,2)];

    obsmati = [ones(Nstates,1)-FRET1',FRET1', ones(Nstates,1)-FRET2', FRET2'];
    obsmati = obsmati./sum(obsmati,2);

    % Transition matrix:
    % Impose same transition rates between conformational states

    % default value for expo is 1e5
    expo = 1e5;
    expo = (10.^randi([2 5],Nstates,Nstates)); % random order of magnitude

    % Here, you can set a range (small number corresponds to larger
    % kinetic rates). Adjust the values to your experiment and sample.
    F = rand(Nstates,Nstates).*expo;
    F = F-diag(diag(F));
    F = F-diag(sum(F,2));
    F =F*5e-8 + eye(Nstates);
    transi = diag(sum(F,2))^-1*(F);



    % Add random guesses to array
    InitialGuess{iii}.transmat0=transi;
    InitialGuess{iii}.prior0=priori;
    InitialGuess{iii}.obsmat0=obsmati;
    InitialGuess{iii}.Nstate=Nstates;

    translist = [translist, transi]; % Used to plot matrix elements. Added by Demian
end

% Added by Demian. Use this section to plot the distribution of kinetic
% rates for your initial guesses.

% dt = 0.5*1e-7; % set dt resolution of your experiment
% histogram(log10(translist(:)/dt),50)
% xlabel('log_{10}(Rate/Hz)')
% ylabel('Counts')

end
