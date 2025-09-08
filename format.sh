#!/bin/bash
# Format all MATLAB files in the workspace using MBeautifier
#
# Usage: ./format_matlab.sh
#
# This script will:
# 1. Check if MBeautifier is installed
# 2. Install it if needed
# 3. Run the formatter on all .m files

set -e  # Exit on any error

echo "🎨 MATLAB Code Formatter"
echo "======================="

# Check if MBeautifier is installed
if [ ! -d ".github/scripts/MBeautifier" ]; then
    echo "📦 Installing MBeautifier..."
    git clone https://github.com/davidvarga/MBeautifier.git .github/scripts/MBeautifier
    echo "✅ MBeautifier installed!"
else
    echo "✅ MBeautifier already installed"
fi

# Count .m files
m_files=$(find src tests -name "*.m" 2>/dev/null | wc -l)
echo "📁 Found $m_files MATLAB files to format"

# Run the formatter
echo "🚀 Running MBeautifier..."
matlab -batch "addpath('.github/scripts'); format_matlab_code"

echo "🎉 Formatting complete!"
echo ""
echo "💡 To see what changed, run: git diff --stat"
