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
    % W is structured as W(:, time_step) where each column is a time step
    % The velocity vector is [U; V] stacked vertically
    n_nodes = length(xy1);
    U_final = W(1:n_nodes, end); % Final u-velocity (last time step)
    V_final = W(n_nodes + 1:end, end); % Final v-velocity (last time step)

    % Compute vorticity: omega = dv/dx - du/dy using RBF-FD operators
    % These operators are designed for scattered nodes and use local stencils
    dvdx = Dx * V_final; % dv/dx using RBF-FD operator on scattered nodes
    dudy = Dy * U_final; % du/dy using RBF-FD operator on scattered nodes
    vorticity = dvdx - dudy; % Vorticity field: ω = ∂v/∂x - ∂u/∂y

    % Validate vorticity computation
    fprintf('Vorticity field statistics:\n');
    fprintf('  Min vorticity: %.6f\n', min(vorticity));
    fprintf('  Max vorticity: %.6f\n', max(vorticity));
    fprintf('  Mean vorticity: %.6f\n', mean(vorticity));
    fprintf('  RMS vorticity: %.6f\n', sqrt(mean(vorticity.^2)));

    % Plot vorticity field with improved visualization
    scatter(x1, y1, cfg.visualization.scatter_size * ones(length(xy1), 1), vorticity, 'filled');
    axis equal;
    axis tight;
    hold on;

    % Set reasonable color limits to avoid extreme values dominating the plot
    vort_std = std(vorticity);
    vort_mean = mean(vorticity);
    color_limit = max(abs(vort_mean) + 3 * vort_std, max(abs(vorticity)) * 0.8);

    % Handle case where vorticity is all zeros or NaN
    if color_limit == 0 || ~isfinite(color_limit)
        color_limit = 1.0; % Default color range
        fprintf('Warning: Vorticity field is zero or contains NaN values. Using default color range.\n');
    end
    caxis([-color_limit, color_limit]);

    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    yticks(cfg.visualization.plot_tick_y);
    xticks(cfg.visualization.plot_tick_x);
    title(sprintf('Vorticity Field: ω = ∂v/∂x - ∂u/∂y (Re = %.0f)', 1 / cfg.simulation.viscosity));

    xlabel('x');
    ylabel('y');
    colorbar; % Add colorbar for vorticity

    % Add grid for better visualization
    grid on;
    grid minor;

    drawnow; % Update display immediately

    %% Plot u-velocity component
    figure('Name', 'U-Velocity Component');
    scatter(x1, y1, cfg.visualization.scatter_size * ones(length(xy1), 1), U_final, 'filled');
    axis equal;
    axis tight;
    hold on;

    % Set reasonable color limits for u-velocity
    u_std = std(U_final);
    u_mean = mean(U_final);
    u_color_limit = max(abs(u_mean) + 3 * u_std, max(abs(U_final)) * 0.8);

    % Handle case where u-velocity is all zeros or NaN
    if u_color_limit == 0 || ~isfinite(u_color_limit)
        u_color_limit = 1.0; % Default color range
        fprintf('Warning: U-velocity field is zero or contains NaN values. Using default color range.\n');
    end
    caxis([-u_color_limit, u_color_limit]);

    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    yticks(cfg.visualization.plot_tick_y);
    xticks(cfg.visualization.plot_tick_x);
    title(sprintf('U-Velocity Component (Re = %.0f)', 1 / cfg.simulation.viscosity));

    xlabel('x');
    ylabel('y');
    colorbar;
    colormap(jet);
    grid on;
    grid minor;

    fprintf('U-velocity field statistics:\n');
    fprintf('  Min u-velocity: %.6f\n', min(U_final));
    fprintf('  Max u-velocity: %.6f\n', max(U_final));
    fprintf('  Mean u-velocity: %.6f\n', mean(U_final));
    fprintf('  RMS u-velocity: %.6f\n', sqrt(mean(U_final.^2)));

    drawnow;

    %% Plot v-velocity component
    figure('Name', 'V-Velocity Component');
    scatter(x1, y1, cfg.visualization.scatter_size * ones(length(xy1), 1), V_final, 'filled');
    axis equal;
    axis tight;
    hold on;

    % Set reasonable color limits for v-velocity
    v_std = std(V_final);
    v_mean = mean(V_final);
    v_color_limit = max(abs(v_mean) + 3 * v_std, max(abs(V_final)) * 0.8);

    % Handle case where v-velocity is all zeros or NaN
    if v_color_limit == 0 || ~isfinite(v_color_limit)
        v_color_limit = 1.0; % Default color range
        fprintf('Warning: V-velocity field is zero or contains NaN values. Using default color range.\n');
    end
    caxis([-v_color_limit, v_color_limit]);

    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    yticks(cfg.visualization.plot_tick_y);
    xticks(cfg.visualization.plot_tick_x);
    title(sprintf('V-Velocity Component (Re = %.0f)', 1 / cfg.simulation.viscosity));

    xlabel('x');
    ylabel('y');
    colorbar;
    colormap(jet);
    grid on;
    grid minor;

    fprintf('V-velocity field statistics:\n');
    fprintf('  Min v-velocity: %.6f\n', min(V_final));
    fprintf('  Max v-velocity: %.6f\n', max(V_final));
    fprintf('  Mean v-velocity: %.6f\n', mean(V_final));
    fprintf('  RMS v-velocity: %.6f\n', sqrt(mean(V_final.^2)));

    drawnow;

    %% Plot velocity magnitude
    velocity_magnitude = sqrt(U_final.^2 + V_final.^2);

    figure('Name', 'Velocity Magnitude');
    scatter(x1, y1, cfg.visualization.scatter_size * ones(length(xy1), 1), velocity_magnitude, 'filled');
    axis equal;
    axis tight;
    hold on;

    % Color limits for velocity magnitude (always positive)
    mag_max = max(velocity_magnitude);

    % Handle case where velocity magnitude is all zeros or NaN
    if mag_max == 0 || ~isfinite(mag_max)
        mag_max = 1.0; % Default color range
        fprintf('Warning: Velocity magnitude is zero or contains NaN values. Using default color range.\n');
    end
    caxis([0, mag_max]);

    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    yticks(cfg.visualization.plot_tick_y);
    xticks(cfg.visualization.plot_tick_x);
    title(sprintf('Velocity Magnitude |V| (Re = %.0f)', 1 / cfg.simulation.viscosity));

    xlabel('x');
    ylabel('y');
    colorbar;
    colormap(jet);
    grid on;
    grid minor;

    fprintf('Velocity magnitude statistics:\n');
    fprintf('  Min |V|: %.6f\n', min(velocity_magnitude));
    fprintf('  Max |V|: %.6f\n', max(velocity_magnitude));
    fprintf('  Mean |V|: %.6f\n', mean(velocity_magnitude));
    fprintf('  RMS |V|: %.6f\n', sqrt(mean(velocity_magnitude.^2)));

    drawnow;

end
