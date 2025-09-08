function format_code()
% FORMAT_CODE Format all MATLAB code in the repository
%
% This script applies consistent formatting to all .m files in src/, tests/, 
% and config files. It can be run locally to format code, and in CI to check
% if code is properly formatted.
%
% Usage:
%   format_code()           % Format all files
%
% The formatting rules applied:
% - Remove trailing whitespace
% - Ensure single final newline
% - Normalize line endings to LF
% - Reduce multiple consecutive empty lines to maximum of 2
% - Basic indentation cleanup

fprintf('🎨 MATLAB Code Formatter\n');
fprintf('========================\n');

% Find all MATLAB files
files = [
    dir('src/**/*.m');
    dir('tests/**/*.m'); 
    dir('config*.m');
    dir('simulate.m')
];

if isempty(files)
    fprintf('No MATLAB files found to format.\n');
    return;
end

fprintf('Found %d MATLAB files to format:\n', length(files));

totalChanges = 0;
for i = 1:length(files)
    filePath = fullfile(files(i).folder, files(i).name);
    relPath = strrep(filePath, [pwd, filesep], '');
    
    changes = formatFile(filePath);
    if changes > 0
        fprintf('  ✏️  %s (%d changes)\n', relPath, changes);
        totalChanges = totalChanges + changes;
    else
        fprintf('  ✅ %s (no changes)\n', relPath);
    end
end

fprintf('\n🎉 Formatting complete!\n');
fprintf('Total changes made: %d\n', totalChanges);

if totalChanges > 0
    fprintf('\n💡 Files have been modified. Review changes with: git diff\n');
else
    fprintf('\n✨ All files were already properly formatted!\n');
end

end

function changeCount = formatFile(filePath)
% Format a single MATLAB file
% Returns the number of changes made

changeCount = 0;

% Read file
try
    fid = fopen(filePath, 'r', 'n', 'UTF-8');
    if fid == -1
        warning('Could not open file: %s', filePath);
        return;
    end
    content = fread(fid, '*char')';
    fclose(fid);
catch ME
    warning('Error reading file %s: %s', filePath, ME.message);
    return;
end

originalContent = content;

% 1. Normalize line endings to LF
if contains(content, sprintf('\r\n'))
    content = strrep(content, sprintf('\r\n'), sprintf('\n'));
    changeCount = changeCount + 1;
elseif contains(content, sprintf('\r'))
    content = strrep(content, sprintf('\r'), sprintf('\n'));
    changeCount = changeCount + 1;
end

% 2. Split into lines for processing
lines = strsplit(content, '\n', 'CollapseDelimiters', false);

% 3. Remove trailing whitespace from each line
originalLineCount = length(lines);
for i = 1:length(lines)
    originalLine = lines{i};
    trimmedLine = regexprep(originalLine, '\s+$', '');
    if ~strcmp(originalLine, trimmedLine)
        lines{i} = trimmedLine;
        changeCount = changeCount + 1;
    end
end

% 4. Reduce multiple consecutive empty lines to maximum of 2
newLines = {};
emptyCount = 0;
for i = 1:length(lines)
    if isempty(strtrim(lines{i}))
        emptyCount = emptyCount + 1;
        if emptyCount <= 2
            newLines{end+1} = lines{i}; %#ok<AGROW>
        else
            changeCount = changeCount + 1;
        end
    else
        emptyCount = 0;
        newLines{end+1} = lines{i}; %#ok<AGROW>
    end
end
lines = newLines;

% 5. Ensure file ends with exactly one newline
if ~isempty(lines)
    % Remove trailing empty lines except for one
    while length(lines) > 1 && isempty(strtrim(lines{end}))
        lines(end) = [];
        changeCount = changeCount + 1;
    end
    
    % Add final newline if missing
    if ~isempty(lines{end})
        lines{end+1} = '';
        changeCount = changeCount + 1;
    end
end

% 6. Basic indentation cleanup (convert tabs to spaces)
for i = 1:length(lines)
    if contains(lines{i}, sprintf('\t'))
        lines{i} = strrep(lines{i}, sprintf('\t'), '    '); % 4 spaces per tab
        changeCount = changeCount + 1;
    end
end

% Reconstruct content
newContent = strjoin(lines, '\n');

% Only write if changes were made
if changeCount > 0
    try
        fid = fopen(filePath, 'w', 'n', 'UTF-8');
        if fid == -1
            warning('Could not write to file: %s', filePath);
            return;
        end
        fwrite(fid, newContent, 'char');
        fclose(fid);
    catch ME
        warning('Error writing file %s: %s', filePath, ME.message);
    end
end

end
