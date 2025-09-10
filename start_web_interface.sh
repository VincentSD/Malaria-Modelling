#!/bin/bash
# Malaria Model Simulation - Web Interface Startup Script
# =====================================================

echo "🦟 Malaria Model Simulation - Web Interface"
echo "==========================================="
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

# Check if Streamlit is installed
if ! python -c "import streamlit" 2>/dev/null; then
    echo "❌ Streamlit not found! Installing..."
    pip install streamlit
fi

echo "✅ Environment ready!"
echo ""
echo "🌐 Starting web interface..."
echo "📱 The dashboard will open in your default web browser"
echo "🔗 If it doesn't open automatically, go to: http://localhost:8501"
echo ""
echo "💡 Tips:"
echo "  - Use the sidebar to adjust model parameters"
echo "  - Click 'Run Simulation' to see results"
echo "  - Try different model types for comparison"
echo "  - Use the sensitivity analysis features"
echo ""
echo "🛑 Press Ctrl+C to stop the web interface"
echo ""

# Start Streamlit
streamlit run web_interface.py --server.port 8501 --server.address localhost


