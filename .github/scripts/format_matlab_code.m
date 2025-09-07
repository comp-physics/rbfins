function format_matlab_code()
% FORMAT_MATLAB_CODE Format MATLAB code using MBeautifier
%
% This function formats all MATLAB files in the src/ and tests/ directories
% using MBeautifier. It requires MBeautifier to be installed in ci/MBeautifier.
%
% To use:
% 1. Clone MBeautifier: git clone https://github.com/davidvarga/MBeautifier.git .github/scripts/MBeautifier
% 2. Run this script: matlab -batch "addpath('.github/scripts'); format_matlab_code"

% Add MBeautifier to path
addpath(fullfile('.github', 'scripts', 'MBeautifier'));

% Initialize MBeautifier
addpath(genpath(fullfile('.github', 'scripts', 'MBeautifier')));

% Get all MATLAB files
srcFiles = dir('src/**/*.m');
testFiles = dir('tests/**/*.m');
allFiles = [srcFiles; testFiles];

% Format files
fprintf('Formatting %d MATLAB files...\n', length(allFiles));
for i = 1:length(allFiles)
    file = fullfile(allFiles(i).folder, allFiles(i).name);
    fprintf('Formatting %s\n', file);
    
    % Format the file
    MBeautify.formatFile(file, file);
end

fprintf('\nDone! All files formatted.\n');
end
