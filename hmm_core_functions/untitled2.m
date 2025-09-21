% --- START OF FILE test_InitialGuessFunctions.m ---
clearvars; close all; clc;

fprintf('--- Testing Initial Guess Functions ---\n');

% --- Define Test Parameters ---
num_guesses_to_generate = 5;
test_Nstates = 3;
test_O_channels = 2; % For obsmat comparison (v2) and to match v1's implicit 2 channels
test_dt_analysis_sec = 100e-9; % 100 ns, example instrument resolution

fprintf('Test Parameters:\n');
fprintf('  Number of Guesses: %d\n', num_guesses_to_generate);
fprintf('  Number of States (Nstates): %d\n', test_Nstates);
fprintf('  Number of Observation Channels (O_channels for v2): %d\n', test_O_channels);
fprintf('  Time Resolution (dt_analysis_sec for v2): %.2e s\n\n', test_dt_analysis_sec);

% --- Call Original InitialGuess_PhotonByPhoton ---
% Note: The original function has a loop `for state=Nstates`.
% If Nstates is a scalar, it runs once for that Nstates.
% We pass test_Nstates as a scalar.
fprintf('--- Generating Guesses with ORIGINAL InitialGuess_PhotonByPhoton.m ---\n');
try
    guesses_v1 = InitialGuess_PhotonByPhoton(num_guesses_to_generate, test_Nstates);
    
    fprintf('Successfully generated %d guesses using V1.\n', length(guesses_v1));
    for i = 1:min(2, length(guesses_v1)) % Display first few
        fprintf('\nGuess V1 - %d:\n', i);
        disp('  Nstate:'); disp(guesses_v1{i}.Nstate);
        disp('  prior0:'); disp(guesses_v1{i}.prior0);
        disp('  transmat0:'); disp(guesses_v1{i}.transmat0);
        disp('  obsmat0:'); disp(guesses_v1{i}.obsmat0);
        
        % Check properties
        if ~isfield(guesses_v1{i},'Nstate') || guesses_v1{i}.Nstate ~= test_Nstates
            warning('TestV1: Nstate mismatch for guess %d.', i);
        end
        if abs(sum(guesses_v1{i}.prior0) - 1) > 1e-9
            warning('TestV1: Prior for guess %d does not sum to 1.', i);
        end
        if any(abs(sum(guesses_v1{i}.transmat0, 2) - 1) > 1e-9)
            warning('TestV1: Transmat rows for guess %d do not sum to 1.', i);
        end
        if size(guesses_v1{i}.obsmat0, 2) ~= 2 % Original hardcodes 2 channels
            warning('TestV1: Obsmat for guess %d does not have 2 columns.',i);
        end
        if any(abs(sum(guesses_v1{i}.obsmat0, 2) - 1) > 1e-9)
            warning('TestV1: Obsmat rows for guess %d do not sum to 1.', i);
        end
    end
catch ME_v1
    fprintf('ERROR generating guesses with V1: %s\n', ME_v1.message);
    disp(ME_v1.getReport('extended','hyperlinks','off'));
    guesses_v1 = {};
end

% --- Call Modified InitialGuess_PhotonByPhoton_v2 ---
fprintf('\n--- Generating Guesses with MODIFIED InitialGuess_PhotonByPhoton_v2.m ---\n');
try
    guesses_v2 = InitialGuess_PhotonByPhoton_v2(num_guesses_to_generate, test_Nstates, test_O_channels, test_dt_analysis_sec);
    
    fprintf('Successfully generated %d guesses using V2.\n', length(guesses_v2));
    for i = 1:min(2, length(guesses_v2)) % Display first few
        fprintf('\nGuess V2 - %d:\n', i);
        disp('  Nstate:'); disp(guesses_v2{i}.Nstate);
        disp('  prior0:'); disp(guesses_v2{i}.prior0);
        disp('  transmat0:'); disp(guesses_v2{i}.transmat0);
        disp('  obsmat0:'); disp(guesses_v2{i}.obsmat0);

        % Check properties
        if ~isfield(guesses_v2{i},'Nstate') || guesses_v2{i}.Nstate ~= test_Nstates
            warning('TestV2: Nstate mismatch for guess %d.', i);
        end
        if abs(sum(guesses_v2{i}.prior0) - 1) > 1e-9
            warning('TestV2: Prior for guess %d does not sum to 1.', i);
        end
        if any(abs(sum(guesses_v2{i}.transmat0, 2) - 1) > 1e-9)
            warning('TestV2: Transmat rows for guess %d do not sum to 1.', i);
        end
        if size(guesses_v2{i}.obsmat0, 1) ~= test_Nstates || size(guesses_v2{i}.obsmat0, 2) ~= test_O_channels
             warning('TestV2: Obsmat for guess %d has wrong dimensions. Expected %dx%d, got %dx%d.', i, test_Nstates, test_O_channels, size(guesses_v2{i}.obsmat0,1), size(guesses_v2{i}.obsmat0,2));
        end
        if any(abs(sum(guesses_v2{i}.obsmat0, 2) - 1) > 1e-9)
            warning('TestV2: Obsmat rows for guess %d do not sum to 1.', i);
        end

        % Check rate awareness (qualitative for now)
        % Off-diagonal elements of transmat0 should be roughly k_eff * dt_analysis_sec
        % where k_eff is within [k_min_Hz, k_max_Hz] but scaled by target_sum_off_diagonal_P
        % So, P_ij (off-diag) should be <= target_sum_off_diagonal_P
        % And sum(P_ij for j~=i) should be close to target_sum_off_diagonal_P
        fprintf('  Transmat0 V2 - %d, off-diagonal sums:\n', i);
        for r_idx = 1:test_Nstates
            off_diag_sum_v2 = sum(guesses_v2{i}.transmat0(r_idx, find((1:test_Nstates) ~= r_idx) ));
            fprintf('    Row %d: sum(off-diag P) = %.4f (target_sum_off_diagonal_P was ~0.1)\n', r_idx, off_diag_sum_v2);
            if off_diag_sum_v2 > 0.5 % A loose check
                warning('TestV2: Off-diagonal sum for transmat0 row %d seems high: %.3f', r_idx, off_diag_sum_v2);
            end
        end
    end

    % Example: Test with different O_channels for V2
    if test_Nstates == 2 && test_O_channels == 2 % For a fair comparison with v1 first
        test_O_channels_alt = 3;
        fprintf('\n--- Generating Guesses with V2 and O_channels = %d ---\n', test_O_channels_alt);
        guesses_v2_alt_O = InitialGuess_PhotonByPhoton_v2(1, test_Nstates, test_O_channels_alt, test_dt_analysis_sec);
        if ~isempty(guesses_v2_alt_O)
            disp('Guess V2 (alt O_channels) - 1:');
            disp('  obsmat0:'); disp(guesses_v2_alt_O{1}.obsmat0);
             if size(guesses_v2_alt_O{1}.obsmat0, 2) ~= test_O_channels_alt
                warning('TestV2_alt_O: Obsmat for guess 1 does not have %d columns.', test_O_channels_alt);
            end
        end
    end


catch ME_v2
    fprintf('ERROR generating guesses with V2: %s\n', ME_v2.message);
    disp(ME_v2.getReport('extended','hyperlinks','off'));
    guesses_v2 = {};
end

fprintf('\n--- Test Complete ---\n');
% --- END OF FILE test_InitialGuessFunctions.m ---