%% InitialGuess_block_chain
function [InitialGuess answer] = InitialGuess_chainBlock(number_of_initialguess,state)
% This function construct an initial guess array for HMM analysis The
% initial guess array contains for each one the HMM parameters: -
% Transition matrix - Pi vector - Observation matrix - and number of states


translist = [];% Used to plot matrix elements. Added by Demian

answer =[];
answer = input('Select intial guess mode: random [r], prior model [p]:\n','s');

while ~strcmp(answer,'r')&~strcmp(answer,'p')

    answer = input('Invalid input\nSelect intial guess mode: random [r], prior model [p]:\n','s');
end

if  strcmp(answer,'r')

    blockNum = input('Select number of models to integrate. Choose integer number between 2-4\n');
    while ~any(blockNum==[2:4])
        blockNum = input('Invalid input\nSelect number of models to integrate. Choose integer number between 2-4\n');
    end


    for i=1:number_of_initialguess

        % Added by Demian:24.05.20
        % Randomize a raw stochastic matrix: The expo variable determines
        % around which values the random off-diagonal probabilities will be
        % centered. The parameter expo determines the mean value. For
        % instance, for an experiment with resolution of 5e-8 sec, an expo
        % value of 5 will yield off diagonal probabilities, corresponding
        % to off diagonal kinetic rates around 10^5 Hz. Plot a log
        % histogram of the rates to convince yourself.


        temp = rand(1,state*blockNum);
        priori = temp/sum(temp); % Pi vector

        temp = sort(rand(state,1)); % Sorted
        obsmati = repmat([1-temp temp],3,1); % Observation matrix

        %block diagonals slow to fast and slow
        % first block (start from 100 Hz)

        K = rand(state,state)*10.^(2);
        K = K-tril(K,-2)-triu(K,2)-diag(diag(K)); % making chain architecture
        K = K - diag(sum(K,2));

        % default time unit of 50 ns
        P = eye(state)+K*5e-8;
        transi = P;

        % add blocks
        for b = 2:blockNum

            K = rand(state,state)*10.^(b+2);
            K = K-tril(K,-2)-triu(K,2)-diag(diag(K)); % making chain architecture
            K = K - diag(sum(K,2));

            % default time unit of 50 ns
            P = eye(state)+K*5e-8;
            transi = blkdiag(transi,P);

        end



        % confrim normalizaion
        transi = transi./sum(transi,2);

        % Add random guesses to array
        InitialGuess{i}.transmat0=transi;
        InitialGuess{i}.prior0=priori;
        InitialGuess{i}.obsmat0=obsmati;
        InitialGuess{i}.Nstate=state;

        %         translist = [translist, transi]; % Used to plot matrix elements. Added by Demian
    end
elseif strcmp(answer,'p')

    number_of_initialguess = 1;
    [files path] = uigetfile('MultiSelect','on');


    fileTemp = load([path files{1}],'model');
    obsmati = fileTemp.model.obsmat;
    priori =  fileTemp.model.prior';
    transi = fileTemp.model.transmat;

    for i=2:numel(files)
        fileTemp = load([path files{i}],'model');
        obsmati = [obsmati;fileTemp.model.obsmat];
        priori = [priori fileTemp.model.prior'];
        transi = blkdiag(transi,fileTemp.model.transmat);
    end
    InitialGuess{1}.transmat0=transi;
    InitialGuess{1}.prior0=priori;
    InitialGuess{1}.obsmat0=obsmati;
    InitialGuess{1}.Nstate=state;
    InitialGuess{1}.modelOrder = files;

end

end