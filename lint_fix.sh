#!/bin/bash
# MATLAB Code Auto-Fixer
# This script automatically fixes MATLAB code issues where possible
#
# Usage: ./lint_fix.sh

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}[LINT-FIX] MATLAB Code Auto-Fixer${NC}"
echo "=================================="

# Function to run MATLAB auto-fixes
run_matlab_fixes() {
    echo -e "${YELLOW}Running MATLAB automatic fixes...${NC}"
    
    # Create a temporary MATLAB script for auto-fixing
    cat > temp_fix_script.m << 'EOF'
try
    % Analyze code issues in src and tests directories
    fprintf('Analyzing code issues...\n');
    issues = codeIssues({'src', 'tests'}, 'IncludeSubfolders', true);
    
    if height(issues.Issues) == 0
        fprintf('No fixable issues found.\n');
    else
        fprintf('Found %d total issues, attempting to fix...\n', height(issues.Issues));
        
        % Try to fix issues automatically
        % Note: fix() may not work in all MATLAB versions or with all issue types
        try
            fixedIssues = fix(issues);
        catch fixError
            fprintf('Auto-fix not available or failed: %s\n', fixError.message);
            fprintf('Manual fixes may be required for some issues.\n');
            fixedIssues = [];
        end
        
        if ~isempty(fixedIssues) && height(fixedIssues.Issues) > 0
            fprintf('Successfully fixed %d issues:\n', height(fixedIssues.Issues));
            for i = 1:height(fixedIssues.Issues)
                issue = fixedIssues.Issues(i,:);
                fprintf('  - %s:%d - %s [%s]\n', issue.FullFilename{1}, issue.LineNumber, issue.Description{1}, issue.CheckID{1});
            end
        else
            fprintf('No issues could be automatically fixed.\n');
        end
    end
    
    fprintf('MATLAB auto-fix completed.\n');
catch ME
    fprintf('Error during MATLAB auto-fix: %s\n', ME.message);
end

% Clean up
if exist('temp_fix_script.m', 'file')
    delete('temp_fix_script.m');
end
EOF
    
    # Run the MATLAB script
    if matlab -batch "run('temp_fix_script.m')" -nodesktop -nosplash; then
        echo -e "${GREEN}[SUCCESS] MATLAB auto-fix completed!${NC}"
        # Clean up the temp script if MATLAB didn't
        [ -f "temp_fix_script.m" ] && rm -f temp_fix_script.m
        return 0
    else
        echo -e "${RED}[FAILED] MATLAB auto-fix encountered errors${NC}"
        # Clean up the temp script if MATLAB didn't
        [ -f "temp_fix_script.m" ] && rm -f temp_fix_script.m
        return 1
    fi
}


# Main execution
echo "Attempting to automatically fix MATLAB code issues..."
echo ""

# Run MATLAB auto-fixes
if run_matlab_fixes; then
    echo ""
    echo "=================================="
    echo -e "${GREEN}[SUCCESS] MATLAB auto-fixes completed! ✓${NC}"
    echo -e "${YELLOW}[INFO] Run ./lint.sh to check if all MATLAB issues are resolved${NC}"
    echo -e "${YELLOW}[INFO] Run ./format.sh to check/fix any style issues${NC}"
    exit 0
else
    echo ""
    echo "=================================="
    echo -e "${RED}[FAILED] MATLAB auto-fixes failed! ✗${NC}"
    echo -e "${YELLOW}[INFO] Some issues may require manual fixing${NC}"
    echo -e "${YELLOW}[INFO] Run ./lint.sh to see remaining MATLAB issues${NC}"
    exit 1
fi
