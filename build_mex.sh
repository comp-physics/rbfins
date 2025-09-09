#!/bin/bash
# Build script for MEX files
# This script compiles the C++ MEX implementations for performance-critical functions

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}[BUILD] Compiling MEX files${NC}"
echo "=========================="

# Change to geometry directory
cd src/geometry

echo -e "${YELLOW}Compiling dairfoil.cpp...${NC}"
if matlab -batch "mex dairfoil.cpp" -nodesktop -nosplash; then
    echo -e "${GREEN}[SUCCESS] dairfoil MEX compiled successfully!${NC}"
else
    echo -e "${RED}[FAILED] dairfoil MEX compilation failed!${NC}"
    exit 1
fi


# Return to root directory
cd ../..

echo ""
echo "=========================="
echo -e "${GREEN}[SUCCESS] dairfoil MEX compiled! ✓${NC}"
echo -e "${YELLOW}[INFO] dairfoil MEX provides significant speedup for airfoil mesh generation${NC}"
echo -e "${YELLOW}[INFO] DistMesh calls the distance function thousands of times during mesh generation${NC}"
