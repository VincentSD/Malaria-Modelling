#!/bin/bash

# Malaria Model Simulation - R Example Runner
# ==========================================
#
# This script runs the complete R implementation of the malaria model simulation project.

echo "🦟 Malaria Model Simulation - R Implementation"
echo "=============================================="
echo ""

# Check if R is installed
if ! command -v Rscript &> /dev/null; then
    echo "❌ Error: R is not installed or not in PATH"
    echo "Please install R from https://www.r-project.org/"
    exit 1
fi

echo "✅ R found: $(Rscript --version | head -n1)"
echo ""

# Check if required packages are installed
echo "🔍 Checking required R packages..."
Rscript -e "
required_packages <- c('deSolve', 'ggplot2', 'dplyr', 'gridExtra', 'reshape2', 'viridis')
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
    cat('❌ Missing packages:', paste(missing_packages, collapse = ', '), '\n')
    cat('Installing missing packages...\n')
    install.packages(missing_packages, repos = 'https://cran.r-project.org/')
} else {
    cat('✅ All required packages are installed\n')
}
"
echo ""

# Create results directory if it doesn't exist
if [ ! -d "results" ]; then
    echo "📁 Creating results directory..."
    mkdir -p results/plots results/analysis
fi

# Run the main example
echo "🚀 Running malaria model simulation example..."
echo "This may take a few minutes..."
echo ""

Rscript run_example.R

# Check if the example ran successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Example completed successfully!"
    echo ""
    echo "📊 Generated files:"
    echo "  - results/plots/: Model visualizations"
    echo "  - results/analysis/: Analysis results"
    echo ""
    echo "🔍 To view results:"
    echo "  - Open the generated PNG files in results/plots/"
    echo "  - Check CSV files in results/analysis/ for numerical results"
    echo ""
    echo "📚 Next steps:"
    echo "  - Try individual models: source('models/basic_sir.R')"
    echo "  - Run parameter analysis: source('analysis/parameter_analysis.R')"
    echo "  - Compare with Python implementation: ../python/run_example.py"
    echo ""
    echo "🌐 For web interface, use Python version: ../start_web_interface.sh"
else
    echo ""
    echo "❌ Example failed to run"
    echo "Check the error messages above for details"
    exit 1
fi


