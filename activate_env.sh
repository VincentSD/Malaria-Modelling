#!/bin/bash
# Malaria Model Simulation Project - Environment Activation Script
# =============================================================

echo "🦟 Malaria Model Simulation Project"
echo "===================================="
echo ""

# Check if virtual environment exists
if [ ! -d "malaria_env" ]; then
    echo "❌ Virtual environment not found!"
    echo "Please run: python3 -m venv malaria_env"
    echo "Then run: source malaria_env/bin/activate"
    echo "Finally run: pip install -r python/requirements.txt"
    exit 1
fi

# Activate virtual environment
echo "🔧 Activating virtual environment..."
source malaria_env/bin/activate

echo "✅ Environment activated successfully!"
echo ""
echo "📚 Available commands:"
echo "  ./start_web_interface.sh                 # 🌐 Start web dashboard (RECOMMENDED)"
echo "  python run_example.py                    # Run complete example"
echo "  python python/models/basic_sir.py        # Run basic SIR model"
echo "  python python/models/seir_model.py      # Run SEIR model"
echo "  python python/models/age_structured_sir.py  # Run age-structured model"
echo "  python python/models/vector_dynamics.py # Run vector dynamics model"
echo ""
echo "📊 Analysis tools:"
echo "  python python/analysis/visualization.py     # Visualization tools"
echo "  python python/analysis/parameter_analysis.py # Parameter analysis"
echo ""
echo "🔬 R models:"
echo "  Rscript r/models/basic_sir.R"
echo "  Rscript r/models/age_structured_sir.R"
echo ""
echo "📖 Documentation:"
echo "  cat README.md"
echo "  cat EXPERIMENTAL_SETUPS.md              # 🧪 Experimental scenarios"
echo ""
echo "Happy modeling! 🧬"
