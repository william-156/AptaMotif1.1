#!/bin/bash

# AptaMotif Analyzer Launch Script
echo "========================================="
echo "  AptaMotif Analyzer - Starting Up"
echo "========================================="
echo ""

# Check if streamlit is installed
if ! command -v streamlit &> /dev/null; then
    echo "❌ Streamlit not found. Installing required packages..."
    pip install streamlit biopython logomaker plotly scipy pandas numpy matplotlib seaborn statsmodels --break-system-packages
    echo "✅ Installation complete!"
    echo ""
fi

# Check for ViennaRNA (optional)
python3 -c "import RNA" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "⚠️  ViennaRNA Python package not found."
    echo "    Structure prediction will use fallback algorithm."
    echo "    For full ViennaRNA support, install with:"
    echo "    conda install -c bioconda viennarna"
    echo ""
fi

echo "🚀 Launching AptaMotif Analyzer..."
echo ""
echo "📖 Usage Tips:"
echo "  1. Configure your pool in the sidebar"
echo "  2. Input sequences (FASTA or plain text)"
echo "  3. Run motif and structure analysis"
echo "  4. Export results"
echo ""
echo "📁 Example sequences available in: example_sequences.fa"
echo ""
echo "Press Ctrl+C to stop the server"
echo "========================================="
echo ""

# Launch streamlit
streamlit run aptamer_motif_analyzer.py
