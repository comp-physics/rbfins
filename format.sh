#!/bin/bash
# MATLAB Code Formatter using MISS_HIT
# This script formats MATLAB code using the professional MISS_HIT formatter
#
# Usage: ./format.sh

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}[FORMAT] MATLAB Code Formatter (MISS_HIT)${NC}"
echo "=================================="

# Check if we're in a virtual environment or if MISS_HIT is available
if ! command -v mh_style &> /dev/null; then
    if [ -d "venv" ]; then
        echo "Activating existing virtual environment..."
        source venv/bin/activate
    else
        echo -e "${YELLOW}Creating virtual environment...${NC}"
        python3 -m venv venv
        if [ $? -ne 0 ]; then
            echo -e "${RED}[ERROR] Failed to create virtual environment. Please ensure Python 3 is installed.${NC}"
            exit 1
        fi
        
        echo "Activating virtual environment..."
        source venv/bin/activate
        
        echo -e "${YELLOW}Installing MISS_HIT...${NC}"
        pip install --upgrade pip
        pip install miss_hit
        if [ $? -ne 0 ]; then
            echo -e "${RED}[ERROR] Failed to install MISS_HIT${NC}"
            exit 1
        fi
        
        echo -e "${GREEN}[SUCCESS] Virtual environment created and MISS_HIT installed!${NC}"
    fi
fi

# Check if MISS_HIT is now available
if ! command -v mh_style &> /dev/null; then
    echo -e "${RED}[ERROR] MISS_HIT still not found after setup${NC}"
    echo -e "${YELLOW}[INFO] Try installing manually: pip install miss_hit${NC}"
    exit 1
fi

# Format all MATLAB files
echo "Formatting MATLAB files..."

# Apply automatic fixes
echo -e "${YELLOW}Applying automatic fixes...${NC}"
mh_style --line_length 140 --fix src/ tests/ config.m simulate.m --brief

# Check final status
echo -e "${YELLOW}Final formatting check...${NC}"
if mh_style --line_length 140 src/ tests/ config.m simulate.m --brief; then
    echo -e "${GREEN}[SUCCESS] All files are properly formatted!${NC}"
else
    echo -e "${YELLOW}[WARNING] Some style issues remain (may require manual fixing)${NC}"
    echo -e "${YELLOW}[INFO] Check the output above for details${NC}"
    exit 1
fi
