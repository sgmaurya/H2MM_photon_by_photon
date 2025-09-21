% Traj_Vet_acceptorPulse_local.m
% Visualizes trajectories with Viterbi decoded paths,
% with potential for indicating acceptor pulse regions.

% NOTE: This version assumes acceptor pulse timing information would be passed
% as an additional argument if this functionality is fully implemented.

function Traj_Vet_acceptorPulse_local(data_tv, Q0_tv, ind_trans_tv, trans_num_tv, AA_tv, dt_tv, Nstates_tv, acceptor_pulse_info_tv)
    fprintf('--- Running Traj_Vet_acceptorPulse_local (Visualization Function) ---\n');

    % --- Input Validations ---
    if nargin < 7 % Basic arguments
        error('Traj_Vet_acceptorPulse_local: Not enough input arguments. Expected at least 7.');
    end
    if nargin < 8 || isempty(acceptor_pulse_info_tv)
        acceptor_pulse_info_tv = []; % Default to no pulse info if not provided
        fprintf('INFO: No acceptor_pulse_info_tv provided. Acceptor pulse will not be specifically indicated.\n');
    elseif length(acceptor_pulse_info_tv) ~= length(data_tv) && ~isempty(acceptor_pulse_info_tv)
        warning('Traj_Vet_acceptorPulse_local: Length of acceptor_pulse_info_tv does not match data_tv. Pulse info might be incorrect.');
    end

    if isempty(data_tv)
        warning('Traj_Vet_acceptorPulse_local: Input "data_tv" is empty. No trajectories to plot.');
        return;
    end
    if isempty(Q0_tv) || length(Q0_tv) ~= length(data_tv)
        warning('Traj_Vet_acceptorPulse_local: Input "Q0_tv" is empty or its length does not match "data_tv". Skipping plotting.');
        return;
    end
    if isempty(AA_tv) || size(AA_tv,1) ~= 2 || size(AA_tv,2) ~= Nstates_tv
        error('Traj_Vet_acceptorPulse_local: AA_tv matrix is invalid. It should be 2xNstates_tv.');
    end
     if ~isequal(AA_tv(1,:), 1:Nstates_tv)
        warning('Traj_Vet_acceptorPulse_local: AA_tv(1,:) is not [1:Nstates_tv]. Ensure AA_tv(2, state_label) correctly corresponds to HMM state_label for plotting.');
    end


    % --- Get user input for plotting parameters (same as Traj_Vet_local) ---
    prompt_tv = {'Enter time resolution for binning traces (us):', ...
                 'Enter maximum number of trajectories to display:', ...
                 'Enter sorting option (1=Transitions, 2=Flux, 3=None, 4=Length, 5=Specific State):'};
    dlgtitle_tv = 'Trajectory Visualization Input (Acceptor Pulse Mode)';
    dims_tv = [1 60];
    definput_tv = {'50', '10', '1'};

    answer_tv = inputdlg(prompt_tv, dlgtitle_tv, dims_tv, definput_tv);
    if isempty(answer_tv)
        disp('Trajectory visualization (acceptor pulse mode) cancelled by user.');
        return;
    end

    try
        Time_Res_us = str2double(answer_tv{1});
        num_traj_to_display_max = str2double(answer_tv{2});
        sort_option = str2double(answer_tv{3});
        if isnan(Time_Res_us) || Time_Res_us <= 0 || isnan(num_traj_to_display_max) || num_traj_to_display_max < 1 || isnan(sort_option) || ~ismember(sort_option, 1:5)
            error('Invalid input for visualization parameters. Check values and types.');
        end
    catch ME_tv_input
        error('Error parsing visualization parameters: %s', ME_tv_input.message);
    end

    % --- Arrange data based on sorting option (identical to Traj_Vet_local) ---
    fprintf('Arranging data based on sort option %d...\n', sort_option);
    num_total_trajectories = length(data_tv);
    Trans_list_indices = 1:num_total_trajectories;

    switch sort_option
        case 1 % Sort by number of transitions
            if ~isempty(trans_num_tv) && length(trans_num_tv) == num_total_trajectories
                [~, sorted_indices] = sort(trans_num_tv, 'descend');
                Trans_list_indices = sorted_indices;
                fprintf('Sorting trajectories by number of transitions (descending).\n');
            else
                warning('Traj_Vet_acceptorPulse_local: trans_num_tv is empty or mismatched. Using default order.');
            end
        case 2 % Sort by photon flux
            Av_rate_of_photons = zeros(1, num_total_trajectories);
            for i_flux = 1:num_total_trajectories
                if ~isempty(data_tv{i_flux}) && size(data_tv{i_flux},1) > 0
                    num_photons = size(data_tv{i_flux},1);
                    burst_duration_sec = data_tv{i_flux}(end,1) - data_tv{i_flux}(1,1) + dt_tv;
                    if burst_duration_sec > 1e-9, Av_rate_of_photons(i_flux) = num_photons / burst_duration_sec; end
                end
            end
            if all(Av_rate_of_photons==0), warning('Traj_Vet_acceptorPulse_local: Photon flux is zero. Using default order.');
            else, [~, sorted_indices] = sort(Av_rate_of_photons, 'descend'); Trans_list_indices = sorted_indices; fprintf('Sorting by photon flux.\n'); end
        case 3 % No sorting
            fprintf('No sorting applied to trajectories.\n');
        case 4 % Sort by burst length (duration)
            burst_length_sec = zeros(1, num_total_trajectories);
            for i_len = 1:num_total_trajectories
                if ~isempty(data_tv{i_len}) && size(data_tv{i_len},1) > 0, burst_length_sec(i_len) = data_tv{i_len}(end,1) - data_tv{i_len}(1,1) + dt_tv; end
            end
            if all(burst_length_sec==0), warning('Traj_Vet_acceptorPulse_local: Burst length is zero. Using default order.');
            else, [~, sorted_indices] = sort(burst_length_sec, 'descend'); Trans_list_indices = sorted_indices; fprintf('Sorting by burst length.\n'); end
        case 5 % Select trajectories containing specific state(s)
            st_ans = inputdlg('Enter HMM state(s) to select (comma-separated, e.g., 1 or 1,2):', 'Select State(s)', [1 50], {'1'});
            if isempty(st_ans) || isempty(st_ans{1}), disp('State selection cancelled.'); return; end
            try
                states_to_find = str2num(st_ans{1}); %#ok<ST2NM>
                if isempty(states_to_find) || any(states_to_find < 1) || any(states_to_find > Nstates_tv), error('Invalid state numbers.'); end
            catch ME_st, error('Invalid input for states: %s', ME_st.message); end
            selected_indices = [];
            for i_findst = 1:num_total_trajectories
                if ~isempty(Q0_tv{i_findst}) && any(ismember(unique(Q0_tv{i_findst}), states_to_find)), selected_indices(end+1) = i_findst; %#ok<AGROW>
                end
            end
            Trans_list_indices = selected_indices;
            if isempty(Trans_list_indices), fprintf('No trajectories found for state(s): %s\n', num2str(states_to_find)); return; end
            fprintf('Selected %d trajectories containing specified state(s).\n', length(Trans_list_indices));
        otherwise, error('Invalid sorting option.');
    end

    num_traj_to_plot_actual = min(num_traj_to_display_max, length(Trans_list_indices));
    if num_traj_to_plot_actual == 0, fprintf('No trajectories to plot based on selection criteria.\n'); return; end
    fprintf('Will display up to %d trajectories.\n', num_traj_to_plot_actual);

    fig_traj_vet_ap = figure('Name', 'Trajectory Visualization (Acceptor Pulse Mode)', 'NumberTitle', 'off', 'Color', 'w');
    for i_plot_loop = 1:num_traj_to_plot_actual
        idx_in_data = Trans_list_indices(i_plot_loop);
        current_trajectory_photons = data_tv{idx_in_data};
        current_q0_path = Q0_tv{idx_in_data};

        if isempty(current_trajectory_photons) || size(current_trajectory_photons,1) < 1
            fprintf('Skipping plot for original trajectory index %d - empty photon data.\n', idx_in_data);
            continue;
        end
        % (Further checks for Q0 length mismatch as in Traj_Vet_local can be added here)

        Time_axis_photons_sec = current_trajectory_photons(:,1);
        Time_axis_photons_usec = Time_axis_photons_sec * 1e6;
        Photon_Symbols = current_trajectory_photons(:,2);

        acceptor_photon_times_usec = Time_axis_photons_usec(Photon_Symbols == 2);
        donor_photon_times_usec = Time_axis_photons_usec(Photon_Symbols == 1);

        min_time_usec_traj = Time_axis_photons_usec(1);
        max_time_usec_traj = Time_axis_photons_usec(end);
        if max_time_usec_traj < min_time_usec_traj + Time_Res_us, max_time_usec_traj = min_time_usec_traj + Time_Res_us; end

        X_bins_usec = min_time_usec_traj : Time_Res_us : max_time_usec_traj;
        if isempty(X_bins_usec) || length(X_bins_usec) < 2, X_bins_usec = [min_time_usec_traj, min_time_usec_traj + Time_Res_us]; end
        if X_bins_usec(end) < max_time_usec_traj, X_bins_usec(end+1) = X_bins_usec(end) + Time_Res_us; end
        X_bins_plot_centers = X_bins_usec(1:end-1) + Time_Res_us/2;

        Donor_binned = histcounts(donor_photon_times_usec, X_bins_usec);
        Acceptor_binned = histcounts(acceptor_photon_times_usec, X_bins_usec);
        Total_photons_binned = Donor_binned + Acceptor_binned;
        FRET_binned = Acceptor_binned ./ Total_photons_binned;
        FRET_binned(isnan(FRET_binned) | isinf(FRET_binned)) = 0;

        % --- Top Subplot: FRET and Viterbi Path ---
        subplot(2,1,1); cla;
        if ~isempty(X_bins_plot_centers) && ~isempty(FRET_binned)
            plot(X_bins_plot_centers, FRET_binned, '-b', 'LineWidth', 1.5, 'DisplayName', 'Binned FRET (D/A Symbols)');
        else, plot(0,0,'-b'); end
        hold on;

        % Plot Viterbi path
        if ~isempty(current_q0_path)
            viterbi_plot_values = AA_tv(2, current_q0_path);
            stairs(Time_axis_photons_usec, viterbi_plot_values, '-k', 'LineWidth', 2.5, 'DisplayName', 'Viterbi Path (P(A|State))');
        end

        % Plot transitions
        if length(ind_trans_tv) >= idx_in_data && ~isempty(ind_trans_tv{idx_in_data})
            transition_phot_indices = ind_trans_tv{idx_in_data};
            valid_trans_indices = transition_phot_indices(transition_phot_indices <= length(Time_axis_photons_usec));
            if ~isempty(valid_trans_indices)
                transition_times_usec = Time_axis_photons_usec(valid_trans_indices);
                for k_trans = 1:length(transition_times_usec)
                     plot([transition_times_usec(k_trans), transition_times_usec(k_trans)], [-0.1, 1.1], '--m', 'LineWidth', 1.2, 'HandleVisibility','off');
                end
            end
        end
        
        % --- START: Acceptor Pulse Indication (Conceptual) ---
        if ~isempty(acceptor_pulse_info_tv) && length(acceptor_pulse_info_tv) >= idx_in_data && ...
           ~isempty(acceptor_pulse_info_tv{idx_in_data}) && ...
           isnumeric(acceptor_pulse_info_tv{idx_in_data}) && length(acceptor_pulse_info_tv{idx_in_data}) == 2
            
            pulse_start_sec = acceptor_pulse_info_tv{idx_in_data}(1);
            pulse_end_sec   = acceptor_pulse_info_tv{idx_in_data}(2);
            pulse_start_usec = pulse_start_sec * 1e6;
            pulse_end_usec   = pulse_end_sec * 1e6;

            % Draw a shaded region for the pulse
            yl = ylim;
            patch([pulse_start_usec, pulse_end_usec, pulse_end_usec, pulse_start_usec], ...
                  [yl(1), yl(1), yl(2), yl(2)], ...
                  'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Acceptor Pulse');
            fprintf('INFO: Indicated acceptor pulse from %.2f us to %.2f us for trajectory %d.\n', pulse_start_usec, pulse_end_usec, idx_in_data);
        end
        % --- END: Acceptor Pulse Indication ---

        xlim_plot = [min_time_usec_traj, max_time_usec_traj];
        if diff(xlim_plot) < 1e-6, xlim_plot(2) = xlim_plot(1) + Time_Res_us; end
        xlim(xlim_plot);
        ylim([-0.1, 1.1]);
        ylabel('FRET Efficiency / P(A|State)');
        title(sprintf('Trajectory #%d (Original Index: %d) - Acceptor Pulse Mode', i_plot_loop, idx_in_data));
        legend('show', 'Location', 'best'); grid on; hold off;

        % --- Bottom Subplot: Photon Counts ---
        subplot(2,1,2); cla;
        if ~isempty(X_bins_plot_centers)
            plot(X_bins_plot_centers, Donor_binned, '-g', 'LineWidth',1.2, 'Marker','.', 'MarkerSize',8, 'DisplayName', 'Donor (Symbol 1)');
            hold on;
            plot(X_bins_plot_centers, Acceptor_binned, '-r', 'LineWidth',1.2, 'Marker','.', 'MarkerSize',8, 'DisplayName', 'Acceptor (Symbol 2)');
            plot(X_bins_plot_centers, Total_photons_binned, '--c', 'LineWidth',1.2, 'Marker','.', 'MarkerSize',8, 'DisplayName', 'Total (D+A)');
        else, plot(0,0,'-k'); end
        
        % --- START: Acceptor Pulse Indication (Bottom Plot) ---
        if ~isempty(acceptor_pulse_info_tv) && length(acceptor_pulse_info_tv) >= idx_in_data && ...
           ~isempty(acceptor_pulse_info_tv{idx_in_data}) && ...
           isnumeric(acceptor_pulse_info_tv{idx_in_data}) && length(acceptor_pulse_info_tv{idx_in_data}) == 2
            
            pulse_start_usec = acceptor_pulse_info_tv{idx_in_data}(1) * 1e6;
            pulse_end_usec   = acceptor_pulse_info_tv{idx_in_data}(2) * 1e6;
            yl_photons = ylim;
            patch([pulse_start_usec, pulse_end_usec, pulse_end_usec, pulse_start_usec], ...
                  [yl_photons(1), yl_photons(1), yl_photons(2), yl_photons(2)], ...
                  'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off'); % No separate legend entry
        end
        % --- END: Acceptor Pulse Indication ---
        
        ylabel(sprintf('Photons / %g us', Time_Res_us));
        xlabel('Time (\mus)');
        xlim(xlim_plot);
        legend('show','Location','best'); grid on; hold off;

        sgtitle(fig_traj_vet_ap, sprintf('Displaying Trajectory %d of %d (Acceptor Pulse Mode)', i_plot_loop, num_traj_to_plot_actual), 'FontSize', 10, 'FontWeight', 'bold');

        if i_plot_loop < num_traj_to_plot_actual
            fprintf('Displaying trajectory %d of %d. Press any key in Command Window to continue (or Ctrl+C to stop)...\n', i_plot_loop, num_traj_to_plot_actual);
            pause; % Waits for user to press a key
        else
            fprintf('Displayed last selected trajectory (%d).\n', i_plot_loop);
            % pause(0.1);
        end
    end

    if exist('fig_traj_vet_ap','var') && ishandle(fig_traj_vet_ap)
        try saveas(fig_traj_vet_ap, 'Last_Viewed_Viterbi_Trajectories_AcceptorPulse.png');
            fprintf('Saved plot to Last_Viewed_Viterbi_Trajectories_AcceptorPulse.png\n');
                catch ME_save
            warning(ME_save.identifier, 'Traj_Vet_acceptorPulse_local: Could not save figure. Error: %s', ME_save.message);
        end
    end
    fprintf('--- Traj_Vet_acceptorPulse_local Finished ---\n');

end % --- END OF FUNCTION Traj_Vet_acceptorPulse_local ---