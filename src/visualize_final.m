function visualize_final(cfg, doPlot, xy1, W, Nt, x_min, x_max, y_min, y_max)
%VISUALIZE_FINAL Create visualization of final simulation results
%
% This function creates plots of the final velocity field if plotting is enabled.
% It shows both u-velocity perturbation (u-1) and v-velocity components.
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

% Only plot if plotting is enabled
if ~doPlot
    return;
end

% Extract coordinates
x1 = xy1(:,1);   % x-coordinates of velocity nodes
y1 = xy1(:,2);   % y-coordinates of velocity nodes

% Create figure for results
figure('Name','1/Re = 1e-2');
colormap(jet)

% Extract final time step velocity components  
j = Nt;
U = W(1:length(xy1), (j-1)*1+1);              % Final u-velocity
V = W(length(xy1)+1:end, (j-1)*1+1);          % Final v-velocity

% Plot u-velocity field (perturbation from uniform flow)
subplot(2,1,1);
% Plot u-1 to show deviation from uniform flow (u=1)
scatter(x1, y1, cfg.visualization.scatter_size*ones(length(xy1),1), 1*(U-1), '.');
axis equal, axis tight, hold on;
xlim([x_min x_max]);
ylim([y_min y_max]);
yticks(cfg.visualization.plot_tick_y)
xticks(cfg.visualization.plot_tick_x)
title('u-velocity perturbation (u-1)');
ylabel('y');
xlabel('x');
shading interp;
caxis([-cfg.visualization.color_axis_range cfg.visualization.color_axis_range]);

% Plot v-velocity field (should be zero in uniform flow)  
subplot(2,1,2);
scatter(x1, y1, cfg.visualization.scatter_size*ones(length(xy1),1), 1*V, '.');
axis equal, axis tight, hold on;
shading interp;

xlim([x_min x_max]);
ylim([y_min y_max]);
yticks(cfg.visualization.plot_tick_y)
xticks(cfg.visualization.plot_tick_x)
title('v-velocity');

xlabel('x');
ylabel('y');
set(gca,'Ytick',[]);  % Remove y-tick labels for cleaner appearance
caxis([-cfg.visualization.color_axis_range cfg.visualization.color_axis_range]);

drawnow;  % Update display immediately

end
