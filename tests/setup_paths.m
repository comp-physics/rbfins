function setup_paths()
% SETUP_PATHS - Ensure all required paths are added to MATLAB path
%
% This function is used to set up the correct paths for the tests,
% especially when running in CI environments where the working directory
% might be different.

% Disable all figure visibility for tests
set(0, 'DefaultFigureVisible', 'off');

% Get the directory where this script is located
scriptDir = fileparts(mfilename('fullpath'));

% Add the src directory to the path (where cylinder.m and config.m are)
addpath(fullfile(scriptDir, '..', 'src'));

% Add the distmesh directory to the path
distmeshPath = fullfile(scriptDir, '..', 'src', 'distmesh');
if exist(distmeshPath, 'dir')
    addpath(distmeshPath);
else
    fprintf('Warning: distmesh directory not found at %s\n', distmeshPath);
end

% Print the current path for debugging
fprintf('MATLAB path set up for tests\n');

end