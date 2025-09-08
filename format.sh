#!/bin/bash
# Format MATLAB code in the repository
#
# Usage: ./format.sh

set -e

echo "🎨 Running MATLAB Code Formatter..."
matlab -batch "format_code"
echo "✅ Formatting complete!"
