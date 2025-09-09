function visualize_final(cfg, doPlot, xy1, W, Nt, x_min, x_max, y_min, y_max, Dx, Dy)
    %VISUALIZE_FINAL Create visualization of final simulation results
    %
    % This function creates a plot of the final vorticity field if plotting is enabled.
    %
    % INPUTS:
    %   cfg    - Configuration structure
    %   doPlot - Boolean indicating if plotting should be enabled
    %   xy1    - Complete velocity grid coordinates
    %   W      - Velocity solution matrix [U; V] for all time steps
    %   Nt     - Number of time steps
    %   x_min  - Minimum x-coordinate of domain
    %   x_max  - Maximum x-coordinate of domain
    %   y_min  - Minimum y-coordinate of domain
    %   y_max  - Maximum y-coordinate of domain
    %   Dx     - RBF-FD differentiation matrix for d/dx operator
    %   Dy     - RBF-FD differentiation matrix for d/dy operator

    % Only plot if plotting is enabled
    if ~doPlot
        return
    end

    % Extract coordinates
    x1 = xy1(:, 1); % x-coordinates of velocity nodes
    y1 = xy1(:, 2); % y-coordinates of velocity nodes

    % Create figure for vorticity visualization
    figure('Name', 'Vorticity Field (1/Re = 1e-2)');
    colormap(jet);

    % Extract final time step velocity components
    j = Nt;
    U = W(1:length(xy1), (j - 1) * 1 + 1); % Final u-velocity
    V = W(length(xy1) + 1:end, (j - 1) * 1 + 1); % Final v-velocity

    % Compute vorticity: omega = dv/dx - du/dy
    dvdx = Dx * V; % dv/dx using RBF-FD operator
    dudy = Dy * U; % du/dy using RBF-FD operator
    vorticity = dvdx - dudy; % Vorticity field

    % Plot vorticity field
    scatter(x1, y1, cfg.visualization.scatter_size * ones(length(xy1), 1), vorticity, '.');
    axis equal;
    axis tight;
    hold on;
    shading interp;

    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    yticks(cfg.visualization.plot_tick_y);
    xticks(cfg.visualization.plot_tick_x);
    title('Vorticity (omega = dv/dx - du/dy)');

    xlabel('x');
    ylabel('y');
    colorbar; % Add colorbar for vorticity

    drawnow; % Update display immediately

end
