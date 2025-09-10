#!/bin/bash
# MATLAB Code Analyzer - Pre-commit Lint Check
# This script runs the same MATLAB linting as the CI lint job
#
# Usage: ./lint.sh

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}[LINT] MATLAB Code Analyzer${NC}"
echo "=========================="

# Function to run MATLAB linting
run_matlab_lint() {
    # Run the same lint function as CI
    if matlab -batch "addpath('.github/scripts'); lint" -nodesktop -nosplash -nojvm; then
        echo -e "${GREEN}[SUCCESS] MATLAB Code Analyzer passed!${NC}"
        return 0
    else
        echo -e "${RED}[FAILED] MATLAB Code Analyzer found issues${NC}"
        if [ -f ".github/scripts/lint.log" ]; then
            echo -e "${YELLOW}Issues found:${NC}"
            cat .github/scripts/lint.log
        fi
        return 1
    fi
}

# Main execution
echo -e "${YELLOW}Running MATLAB Code Analyzer (same as CI linter)...${NC}"

# Run MATLAB linting
if run_matlab_lint; then
    echo ""
    echo "=========================="
    echo -e "${GREEN}[SUCCESS] MATLAB Code Analyzer passed!${NC}"
    echo -e "${GREEN}Your code will pass the CI lint job.${NC}"
    exit 0
else
    echo ""
    echo "=========================="
    echo -e "${RED}[FAILED] MATLAB Code Analyzer found issues!${NC}"
    echo -e "${YELLOW}Please fix the issues above before committing.${NC}"
    echo ""
    echo -e "${BLUE}Additional checks:${NC}"
    echo -e "  • For style issues: ${YELLOW}./format.sh${NC} (check/fix formatting)"
    echo -e "  • For auto-fixes: ${YELLOW}./fix.sh${NC} (attempt automatic fixes)"
    exit 1
fi
