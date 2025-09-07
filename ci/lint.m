function lint
% Lint MATLAB code using Code Analyzer.
folders = {'src','tests'};
% Ignore common issues that are less critical for CI
% CABE:    Cyclomatic complexity warnings
% ACABE:   Anonymous function cyclomatic complexity
% NASGU:   Value assigned but not used
% ASGLU:   Value assigned but not used (variable)
% INUSD:   Input argument might be unused
% NCHKE:   Use NARGOUTCHK without ERROR
% NCHKN:   Use NARGINCHK without ERROR
% CAXIS:   caxis is not recommended, use clim instead
% DATST:   datestr is not recommended
% TNOW1:   now is not recommended
% FNCOLND: Consider explicitly defining the array
% USENS:   Explicitly initialize this variable
ignoreIDs = {
    'CABE';
    'ACABE';
    'NASGU';
    'ASGLU';
    'INUSD';
    'NCHKE';
    'NCHKN';
    'CAXIS';
    'DATST';
    'TNOW1';
    'FNCOLND';
    'USENS'
};
logFile = fullfile('ci','lint.log');
if ~exist('ci','dir')
    mkdir('ci');
end
fid = fopen(logFile,'w');

hasIssues = false;
for f = 1:numel(folders)
    if ~isfolder(folders{f})
        continue;
    end
    S = dir(fullfile(folders{f},'**','*.m'));
    for k = 1:numel(S)
        file = fullfile(S(k).folder,S(k).name);
        % Add more flags if desired
        msgs = checkcode(file,'-cyc','-id','-struct');
        if ~isempty(ignoreIDs)
            msgs = msgs(~ismember({msgs.id},ignoreIDs));
        end
        for m = 1:numel(msgs)
            hasIssues = true;
            line = sprintf('%s:%d:%d: %s [%s]\n', file, msgs(m).line, msgs(m).column, msgs(m).message, msgs(m).id);
            fprintf('%s', line);
            fprintf(fid, '%s', line);
        end
    end
end

fclose(fid);
if hasIssues
    error('Code Analyzer found issues. See ci/lint.log.');
end
end
