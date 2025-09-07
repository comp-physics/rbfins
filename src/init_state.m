function [W, p0] = init_state(xy1, xy1_s, boundary_obs, Nt)
%INIT_STATE Initialize simulation state variables
%
% This function initializes the velocity and pressure fields for the simulation.
% The initial condition is uniform flow (u=1, v=0) with no-slip boundary
% conditions on the obstacle.
%
% INPUTS:
%   xy1         - Complete velocity grid (interior + boundary nodes)
%   xy1_s       - Complete pressure grid (interior + boundary nodes)
%   boundary_obs - Obstacle boundary nodes
%   Nt          - Number of time steps
%
% OUTPUTS:
%   W  - Velocity storage matrix [U; V] for all time steps
%   p0 - Initial pressure field (zero everywhere)

% Initial velocity field: uniform flow (u=1, v=0) everywhere except obstacle
V0 = zeros(length(xy1), 1); % Initial v-velocity (zero everywhere)
U0 = ones(length(xy1), 1); % Initial u-velocity (unit flow)
U0(end-length(boundary_obs)+1:end) = zeros(length(boundary_obs), 1); % No-slip on obstacle
W0 = [U0; V0]; % Combined velocity vector [U; V]

% Storage for velocity history (needed for multi-step time integration)
W = zeros(length(xy1)*2, Nt+1);
W(:, 1) = W0; % Store initial condition

% Initial pressure field (zero everywhere)
p0 = zeros(length(xy1_s), 1); % Pressure field on pressure grid

end
