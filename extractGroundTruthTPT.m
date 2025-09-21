function extractGroundTruthTPT()
    % This script extracts the ground truth transition path times (TPTs) from
    % large molecular dynamics trajectory files without loading them into memory.
    % v2: Corrects the state-machine logic to handle fast transitions where
    %     the particle "skips" over the defined barrier region in a single timestep.

    clearvars; close all; clc;
    fprintf('--- Starting Ground Truth TPT Extraction (v2 - Robust Logic) ---\n');

    %% 1. Define Constants and State Boundaries
    base_path = 'D:\Energy_barrier_random_walk_sim\Asymmetric_barrier';
    WELL_1_BOUNDARY_MAX = 5.7;
    WELL_2_BOUNDARY_MIN = 6.3;
    STATE_WELL_1 = 1;
    STATE_BARRIER = 0;
    STATE_WELL_2 = 2;

    %% 2. Find Simulation Directories and Loop Through Them
    fprintf('Searching for simulation directories in: %s\n', base_path);
    dir_list = dir(fullfile(base_path, 'D_*_Asymmetric'));
    folder_names = {dir_list.name};
    if isempty(folder_names), error('No simulation directories found. Check the base_path.'); end
    
    all_results = cell(1, length(folder_names));

    for i_dir = 1:length(folder_names)
        current_folder = folder_names{i_dir};
        fprintf('\n======================================================\n');
        fprintf('Processing Dataset: %s\n', current_folder);
        
        % Pre-allocate for speed, will trim later
        tpts_W1_to_W2 = zeros(1, 500000); 
        tpts_W2_to_W1 = zeros(1, 500000);
        count12 = 0;
        count21 = 0;
        
        for i_sim = 1:10
            % >>>>>>>>>>>>>>>>>>>> THIS IS THE CORRECTED LINE <<<<<<<<<<<<<<<<<<<<<<
            file_path = fullfile(base_path, current_folder, ['Simulation_', num2str(i_sim)], 'trajectory.txt');
            % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END CORRECTION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            if ~exist(file_path, 'file'), warning('File not found, skipping: %s', file_path); continue; end
            fprintf('  Processing file: %s ...\n', file_path);

            %% 3. Memory-Efficient File Processing (Line-by-Line) - CORRECTED LOGIC
            
            fileID = fopen(file_path, 'r');
            if fileID == -1, error('Cannot open file: %s', file_path); end

            % --- NEW, more robust state variables ---
            last_well_occupied = NaN;   % Tracks the last stable state (1 or 2)
            exit_time_from_last_well = NaN; % Time of the most recent exit from ANY well

            fgetl(fileID); % Skip header
            
            % Initialize with the first line of data
            line_text = fgetl(fileID);
            if ~ischar(line_text), fclose(fileID); continue; end % Skip empty files
            data = sscanf(line_text, '%f\t%f');
            previous_position = data(2);
            previous_state = classify_position(previous_position, WELL_1_BOUNDARY_MAX, WELL_2_BOUNDARY_MIN, STATE_WELL_1, STATE_BARRIER, STATE_WELL_2);
            if previous_state ~= STATE_BARRIER
                last_well_occupied = previous_state;
            end
            
            % Loop through the rest of the file
            while ~feof(fileID)
                line_text = fgetl(fileID);
                if ~ischar(line_text), break; end
                
                data = sscanf(line_text, '%f\t%f');
                time_ns = data(1);
                position_nm = data(2);
                
                current_state = classify_position(position_nm, WELL_1_BOUNDARY_MAX, WELL_2_BOUNDARY_MIN, STATE_WELL_1, STATE_BARRIER, STATE_WELL_2);

                % --- NEW State Transition Logic ---
                if current_state ~= previous_state
                    
                    % Event 1: Particle just EXITED a well.
                    % This is our potential start of a TPT. Record the time.
                    if (previous_state == STATE_WELL_1 || previous_state == STATE_WELL_2)
                        exit_time_from_last_well = time_ns;
                    end
                    
                    % Event 2: Particle just ENTERED a well.
                    % This is the end of a potential TPT. Now we analyze it.
                    if (current_state == STATE_WELL_1 || current_state == STATE_WELL_2)
                        % Check if this is a COMMITTED transition (e.g., last well was 1, current is 2)
                        if ~isnan(last_well_occupied) && current_state ~= last_well_occupied
                            % We have a valid, committed transition. Calculate its TPT.
                            tpt = time_ns - exit_time_from_last_well;
                            
                            if last_well_occupied == STATE_WELL_1 % It was a 1 -> 2 transition
                                count12 = count12 + 1;
                                tpts_W1_to_W2(count12) = tpt;
                            else % It was a 2 -> 1 transition
                                count21 = count21 + 1;
                                tpts_W2_to_W1(count21) = tpt;
                            end
                        end
                        % Crucially, update the last occupied well, regardless of transition type
                        last_well_occupied = current_state;
                    end
                end
                previous_state = current_state;
            end
            
            fclose(fileID);
        end % End of simulation loop
        
        % Trim the pre-allocated arrays to the actual number of transitions found
        tpts_W1_to_W2 = tpts_W1_to_W2(1:count12);
        tpts_W2_to_W1 = tpts_W2_to_W1(1:count21);

        fprintf('    ... Done. Found %d (W1->W2) and %d (W2->W1) total transitions.\n', count12, count21);

        % Store results for this D value
        all_results{i_dir}.D_value_folder = current_folder;
        all_results{i_dir}.TPT_ns_W1_to_W2 = tpts_W1_to_W2;
        all_results{i_dir}.TPT_ns_W2_to_W1 = tpts_W2_to_W1;
        
    end % End of directory loop
    
    %% 4. Save and Plot Final Results (Unchanged)
    fprintf('\n--- All processing complete. Saving results and generating plots. ---\n');
    
    output_filename = fullfile(base_path, 'GroundTruth_TPT_Results.mat');
    save(output_filename, 'all_results', 'WELL_1_BOUNDARY_MAX', 'WELL_2_BOUNDARY_MIN');
    fprintf('Results saved to: %s\n', output_filename);
    
    for i_dir = 1:length(all_results)
        fig = figure('Name', all_results{i_dir}.D_value_folder, 'Position', [100, 100, 800, 600]);
        
        tpts_us_12 = all_results{i_dir}.TPT_ns_W1_to_W2 / 1000;
        tpts_us_21 = all_results{i_dir}.TPT_ns_W2_to_W1 / 1000;
        
        mean_12 = mean(tpts_us_12);
        mean_21 = mean(tpts_us_21);
        
        histogram(tpts_us_12, 'DisplayName', sprintf('Well 1 -> 2 (N=%d, Mean=%.3f us)', length(tpts_us_12), mean_12), 'Normalization', 'pdf', 'BinWidth', 0.1);
        hold on;
        histogram(tpts_us_21, 'DisplayName', sprintf('Well 2 -> 1 (N=%d, Mean=%.3f us)', length(tpts_us_21), mean_21), 'Normalization', 'pdf', 'BinWidth', 0.1);
        
        xlabel('Ground Truth Transition Path Time (\mus)');
        ylabel('Probability Density');
        title({'Ground Truth TPT Distribution', all_results{i_dir}.D_value_folder}, 'Interpreter', 'none');
        legend('show');
        grid on;
        
        plot_filename = fullfile(base_path, [all_results{i_dir}.D_value_folder, '_GroundTruth_TPT_Hist.png']);
        saveas(fig, plot_filename);
        close(fig);
    end
    
    fprintf('Summary plots saved in the base directory.\n');
end

function state = classify_position(pos, max_w1, min_w2, s1, sb, s2)
    if pos < max_w1
        state = s1;
    elseif pos > min_w2
        state = s2;
    else
        state = sb;
    end
end