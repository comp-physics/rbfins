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
% - Convert tabs to 4 spaces
% - Fix indentation for nested blocks (if/for/while/function/etc.)

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
if contains(content, [char(13) newline])
    content = strrep(content, [char(13) newline], newline);
    changeCount = changeCount + 1;
elseif contains(content, char(13))
    content = strrep(content, char(13), newline);
    changeCount = changeCount + 1;
end

% 2. Split into lines for processing
lines = strsplit(content, '\n', 'CollapseDelimiters', false);

% 3. Remove trailing whitespace from each line
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

% 6. Fix indentation (convert tabs to spaces and ensure consistent nesting)
originalLines = lines;
lines = fixIndentation(lines);
if ~isequal(lines, originalLines)
    changeCount = changeCount + 1;
end

% Reconstruct content
newContent = strjoin(lines, '\n');

% Only write if content actually changed
if ~strcmp(originalContent, newContent)
    try
        fid = fopen(filePath, 'w', 'n', 'UTF-8');
        if fid == -1
            warning('Could not write to file: %s', filePath);
            changeCount = 0;
            return;
        end
        fwrite(fid, newContent, 'char');
        fclose(fid);
    catch ME
        warning('Error writing file %s: %s', filePath, ME.message);
        changeCount = 0;
    end
else
    changeCount = 0;  % No actual changes made
end

end

function lines = fixIndentation(lines)
% Fix indentation for MATLAB code
% Converts tabs to spaces and ensures consistent 4-space indentation for nested blocks

indentSize = 4; % MATLAB standard
currentIndent = 0;

% Keywords that increase indentation
increaseKeywords = {
    'if', 'for', 'while', 'switch', 'try', 'parfor', 'spmd', ...
    'methods', 'properties', 'events', 'enumeration'
};

% Keywords that decrease indentation  
decreaseKeywords = {
    'end', 'else', 'elseif', 'case', 'otherwise', 'catch'
};

% Keywords that both decrease and increase (like else, elseif, case, otherwise, catch)
neutralKeywords = {'else', 'elseif', 'case', 'otherwise', 'catch'};

% Add these to the increase keywords since they do increase after being processed
increaseKeywords = [increaseKeywords, neutralKeywords];

for i = 1:length(lines)
    line = lines{i};
    
    % Convert tabs to spaces first
    line = strrep(line, sprintf('\t'), repmat(' ', 1, indentSize));
    
    % Skip empty lines and comments-only lines for indentation logic
    trimmedLine = strtrim(line);
    if isempty(trimmedLine)
        lines{i} = ''; % Keep empty lines empty
        continue;
    end
    
    % Check if this line should decrease indentation before applying
    shouldDecreaseBefore = false;
    for j = 1:length(decreaseKeywords)
        keyword = decreaseKeywords{j};
        % Check if line starts with the keyword (after whitespace)
        if startsWith(trimmedLine, keyword) && ...
           (length(trimmedLine) == length(keyword) || ...
            ~isstrprop(trimmedLine(length(keyword)+1), 'alphanum'))
            shouldDecreaseBefore = true;
            break;
        end
    end
    
    % Decrease indentation before applying to this line
    if shouldDecreaseBefore
        currentIndent = max(0, currentIndent - 1);
    end
    
    % Apply current indentation to the line
    if ~isempty(trimmedLine)
        properIndent = repmat(' ', 1, currentIndent * indentSize);
        lines{i} = [properIndent, trimmedLine];
    end
    
    % Check if this line should increase indentation after
    shouldIncreaseAfter = false;
    for j = 1:length(increaseKeywords)
        keyword = increaseKeywords{j};
        % Check if line starts with the keyword (after whitespace)
        if startsWith(trimmedLine, keyword) && ...
           (length(trimmedLine) == length(keyword) || ...
            ~isstrprop(trimmedLine(length(keyword)+1), 'alphanum'))
            shouldIncreaseAfter = true;
            break;
        end
    end
    
    % Handle indentation increase after this line
    if shouldIncreaseAfter
        % Check if this is a neutral keyword (decrease before, increase after)
        isNeutral = false;
        for j = 1:length(neutralKeywords)
            if startsWith(trimmedLine, neutralKeywords{j})
                isNeutral = true;
                break;
            end
        end
        
        % Only increase if not neutral, or if it's a block-starting neutral keyword
        if ~isNeutral || any(strcmp(trimmedLine(1:min(4,length(trimmedLine))), {'else', 'elif', 'case', 'othe', 'catc'}))
            currentIndent = currentIndent + 1;
        end
    end
    
    % Special case: single-line if statements don't increase indentation
    if startsWith(trimmedLine, 'if ') && contains(trimmedLine, ';') && ~contains(trimmedLine, 'end')
        % This is likely a single-line if, don't increase indentation
        if shouldIncreaseAfter
            currentIndent = currentIndent - 1;
        end
    end
end

end
