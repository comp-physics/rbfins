function [doPlot, isCI, isTest, Nt] = setup_environment(cfg, scriptDir)
%SETUP_ENVIRONMENT Set up environment and dependencies for simulation
%
% This function handles:
% - CI/test environment detection
% - Plotting configuration
% - Time step configuration (reduced for CI/testing)
% - Random seed setup for reproducibility
% - DistMesh library availability
%
% INPUTS:
%   cfg       - Configuration structure from config()
%   scriptDir - Directory where main script is located
%
% OUTPUTS:
%   doPlot - Boolean indicating if plotting should be enabled
%   isCI   - Boolean indicating if running in CI environment
%   isTest - Boolean indicating if running in test environment
%   Nt     - Number of time steps (reduced for CI/testing)

% Check if running in CI environment or headless mode
isCI = strcmpi(getenv('CI'), 'true');
isTest = strcmpi(getenv('MATLAB_TEST'), 'true');

% Completely disable plotting for CI/test environments
doPlot = ~isCI && ~isTest;

if ~doPlot
    % Aggressively disable all graphics for tests
    set(0, 'DefaultFigureVisible', 'off');
    fprintf('Plotting completely disabled for testing/CI\n');
    hasDisplay = false;
else
    % Only check display if we might want to plot
    try
        % Try to create a figure to test if display is available
        f = figure('Visible', 'off');
        close(f);
        hasDisplay = true;
    catch
        hasDisplay = false;
        doPlot = false;
    end
end

% Set random seed for DistMesh reproducibility (DistMesh uses rand() for rejection method)
if isfield(cfg.simulation, 'random_seed') && ~isempty(cfg.simulation.random_seed)
    rng(cfg.simulation.random_seed);
    fprintf('Random seed set to %d for DistMesh reproducibility\n', cfg.simulation.random_seed);
end

% Ensure DistMesh library is available
if ~exist('distmesh2d', 'file')
    % Try to add distmesh from lib directory
    distmeshPath = fullfile(scriptDir, 'lib', 'distmesh');
    if exist(distmeshPath, 'dir')
        addpath(distmeshPath);
    else
        error('distmesh library not found. Please run setup_paths() or add lib/distmesh to your path manually.');
    end
end

% Configure number of time steps based on environment
% Use reduced time steps for CI environment to keep tests fast
Nt = cfg.simulation.num_time_steps;
if isCI || isTest
    Nt = cfg.simulation.num_time_steps_ci; % Reduced steps for fast CI testing
end

end
