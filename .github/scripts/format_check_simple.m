function format_check_simple()
% FORMAT_CHECK_SIMPLE - Simple formatting checks that work in CI
%
% This function performs basic formatting checks that don't require
% the MATLAB Editor, making it suitable for CI environments.

fprintf('🔍 Running simple MATLAB formatting checks...\n');

% Get all MATLAB files
srcFiles = dir('src/**/*.m');
testFiles = dir('tests/**/*.m');
allFiles = [srcFiles; testFiles];

fprintf('Found %d MATLAB files to check\n', length(allFiles));

hasIssues = false;
issueCount = 0;

for i = 1:length(allFiles)
    file = fullfile(allFiles(i).folder, allFiles(i).name);
    fprintf('Checking %s... ', file);
    
    % Read the file
    fid = fopen(file, 'r');
    if fid == -1
        fprintf('❌ Cannot read file\n');
        hasIssues = true;
        issueCount = issueCount + 1;
        continue;
    end
    
    lines = {};
    lineNum = 0;
    fileIssues = 0;
    
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line)
            lineNum = lineNum + 1;
            lines{end+1} = line;
            
            % Check for trailing whitespace
            if ~isempty(line) && (line(end) == ' ' || line(end) == char(9))
                fprintf('\n  Line %d: Trailing whitespace', lineNum);
                fileIssues = fileIssues + 1;
            end
            
            % Check for inconsistent indentation (tabs mixed with spaces)
            if length(line) > 0
                leadingChars = line(1:find(line ~= ' ' & line ~= char(9), 1, 'first')-1);
                if ~isempty(leadingChars)
                    if any(leadingChars == ' ') && any(leadingChars == char(9))
                        fprintf('\n  Line %d: Mixed tabs and spaces', lineNum);
                        fileIssues = fileIssues + 1;
                    end
                end
            end
            
            % Check for very long lines (>150 characters - reasonable for MATLAB)
            if length(line) > 150
                fprintf('\n  Line %d: Line too long (%d chars)', lineNum, length(line));
                fileIssues = fileIssues + 1;
            end
        end
    end
    fclose(fid);
    
    % Check for multiple consecutive empty lines
    emptyCount = 0;
    for j = 1:length(lines)
        if isempty(strtrim(lines{j}))
            emptyCount = emptyCount + 1;
            if emptyCount > 2
                fprintf('\n  Line %d: More than 2 consecutive empty lines', j);
                fileIssues = fileIssues + 1;
                emptyCount = 0; % Reset to avoid multiple warnings for same block
            end
        else
            emptyCount = 0;
        end
    end
    
    % Note: MBeautifier removes trailing empty lines, so we don't check for them
    
    if fileIssues > 0
        fprintf('\n  ❌ %d formatting issues found\n', fileIssues);
        hasIssues = true;
        issueCount = issueCount + fileIssues;
    else
        fprintf('✅ OK\n');
    end
end

fprintf('\n📊 Formatting Check Summary:\n');
fprintf('Files checked: %d\n', length(allFiles));

if hasIssues
    fprintf('❌ Found %d formatting issues\n', issueCount);
    fprintf('\n💡 To fix formatting issues:\n');
    fprintf('   1. Run locally: ./format_matlab.sh\n');
    fprintf('   2. Or manually: matlab -batch "addpath(''.github/scripts''); format_matlab_code"\n');
    fprintf('\n⚠️  Formatting issues found - failing CI.\n');
    
    % Fail CI for formatting issues
    error('Formatting issues found. Please run formatter locally and commit the changes.');
else
    fprintf('✅ All files are properly formatted!\n');
end

end
