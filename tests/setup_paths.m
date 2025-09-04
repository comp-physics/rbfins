function setup_paths()
% SETUP_PATHS - Ensure all required paths are added to MATLAB path
%
% This function is used to set up the correct paths for the tests,
% especially when running in CI environments where the working directory
% might be different.

% Get the directory where this script is located
scriptDir = fileparts(mfilename('fullpath'));

% Add the parent directory to the path (where cylinder.m and config.m are)
addpath(fullfile(scriptDir, '..'));

% Add the distmesh directory to the path
addpath(fullfile(scriptDir, '..', 'distmesh'));

% Print the current path for debugging
fprintf('MATLAB path set up for tests\n');

end