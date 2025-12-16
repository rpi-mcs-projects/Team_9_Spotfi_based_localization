function h_fig = plot_localization_scenario(tx_locs, ap_locs)
% PLOT_LOCALIZATION_SCENARIO Visualizes Tx and AP locations.
%
% Inputs:
%   tx_locs: (N x 2) matrix of Transmitter [x, y] coordinates.
%   ap_locs: (A x B x 2) matrix of AP coordinates.
%            A = Number of APs
%            B = Number of Antennas per AP
%            Dimension 3 = [x, y] coordinates

    % Get the number of APs from the first dimension
    num_aps = size(ap_locs, 1);

    % Create figure
    h_fig = figure('Color', 'w');
    
    % --- 1. Plot Transmitters (Tx) ---
    scatter(tx_locs(:,1), tx_locs(:,2), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
    
    hold on;
    
    % --- 2. Plot Access Points (AP) ---
    % Extract X and Y coordinates for all antennas
    ap_x_all = ap_locs(:,:,1); 
    ap_y_all = ap_locs(:,:,2);
    
    % Plot all AP antennas (flattened)
    scatter(ap_x_all(:), ap_y_all(:), 60, 'r', 's', 'filled', 'MarkerEdgeColor', 'k');
    
    % --- 3. Label APs ---
    for i = 1:num_aps
        % Calculate centroid of the antennas for this specific AP
        centroid_x = mean(ap_x_all(i, :));
        centroid_y = mean(ap_y_all(i, :));
        
        % Add label
        text(centroid_x, centroid_y, ['  AP' num2str(i)], ...
             'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k', ...
             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    end
    
    % --- 4. Styling ---
    grid on;
    box on;
    axis equal;
    xlabel('X Coordinate (m)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y Coordinate (m)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Ground Truth Device Locations', 'FontSize', 14);
    legend('Ground Truth Tx', 'Access Points', 'Location', 'best');
    
    % --- 5. Increase Tick Size and Make Labels Bold ---
    ax = gca; % Get current axes
    ax.FontSize = 14; % Increase tick size
    ax.XAxis.FontWeight = 'bold'; % Make x-axis labels bold
    ax.YAxis.FontWeight = 'bold'; % Make y-axis labels bold
    
    hold off;
end