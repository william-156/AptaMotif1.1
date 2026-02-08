# 🧬 AptaMotif

**Aptamer Motif Discovery & Enrichment Analysis Tool**

A comprehensive web-based application for discovering enriched sequence motifs and analyzing secondary structures in aptamer libraries from SELEX experiments.

## Features

- ✨ **Motif Discovery**: Identifies enriched k-mers (configurable 4-25 bp) shared across sequences
- 📊 **Statistical Analysis**: Binomial test with FDR (Benjamini-Hochberg) correction
- 🧬 **Structure Prediction**: RNA/DNA secondary structure analysis using ViennaRNA
- 📈 **Interactive Visualizations**: Heatmaps, volcano plots, sequence logos, and more
- ⚙️ **Flexible Configuration**: Customizable primers, parameters, and molecule types
- 💾 **Export Results**: Download motif tables, structures, and visualizations

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager
- ViennaRNA package (for structure prediction)

### Step 1: Install ViennaRNA

**On Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install vienna-rna
```

**On macOS (using Homebrew):**
```bash
brew install viennarna
```

**On Windows:**
Download and install from: https://www.tbi.univie.ac.at/RNA/

### Step 2: Install Python Dependencies

```bash
# Navigate to the aptamotif directory
cd aptamotif

# Install required packages
pip install -r requirements.txt
```

## Quick Start

### Launch the Web Application

```bash
streamlit run app.py
```

This will open AptaMotif in your default web browser (typically at http://localhost:8501).

### Basic Usage

1. **Configure Settings** (left sidebar):
   - Select molecule type (RNA or ssDNA)
   - Verify/edit primer sequences
   - Adjust analysis parameters

2. **Input Sequences**:
   - Paste sequences directly or upload a file
   - Supports FASTA or plain text format
   - Sequences should include primer regions

3. **Run Analysis**:
   - Click "🔬 Run Analysis"
   - Wait for processing (typically 30 seconds - 2 minutes)

4. **Explore Results**:
   - View enriched motifs table
   - Explore interactive visualizations
   - Download results as CSV

## Usage Examples

### Example 1: Basic Motif Discovery

```
Input sequences (FASTA format):
>Clone1
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTGGATCCGGATCCGGATCCAGATAGTAAGTGCAATCT

>Clone2
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTGGATCCGGATCCGGATCCAGATAGTAAGTGCAATCT

>Clone3
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTCCAATTCCAATTCCAATTCCAGATAGTAAGTGCAATCT
```

Expected output:
- Motif "GGATCC" found in 2/3 sequences (enriched)
- Motif "CCAATT" found in 1/3 sequences
- Statistical significance calculated for each

### Example 2: Custom Primers

For a different SELEX experiment with custom primers:

1. Open "Edit Primers" in sidebar
2. Enter your forward primer (5' → 3')
3. Enter your reverse primer region
4. Specify expected N-region length
5. Run analysis as usual

## Configuration Options

### Analysis Parameters

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| Min Motif Length | 5 bp | 4-20 bp | Minimum length of motifs to search |
| Max Motif Length | 15 bp | 5-25 bp | Maximum length of motifs to search |
| Min Sequences | 2 | 2-20 | Min. sequences motif must appear in |
| FDR Threshold | 0.05 | 0.01-0.10 | Significance threshold after correction |

### Primer Configuration

**Default N71 Pool Configuration:**
- Forward Primer: `TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATT`
- Reverse Region: `AGATAGTAAGTGCAATCT`
- Expected N-length: 71 bp

These can be modified in the sidebar for different SELEX libraries.

## Understanding Results

### Motif Table Columns

- **Motif**: The sequence motif (k-mer)
- **Length**: Length of motif in base pairs
- **Sequences**: Number of sequences containing the motif
- **Frequency**: Proportion of sequences with motif
- **Enrichment**: Fold enrichment over expected (observed/expected)
- **p-value**: Raw p-value from binomial test
- **Adj. p-value**: FDR-corrected p-value (Benjamini-Hochberg)
- **GC%**: GC content of the motif

### Statistical Interpretation

**p-value < 0.05 (after FDR correction)**: Motif is significantly enriched

**Fold Enrichment > 2**: Motif appears at least twice as often as expected by chance

### Structure Analysis

When enabled, provides:
- Individual secondary structures (dot-bracket notation)
- Minimum free energy (MFE) values
- Consensus structural properties
- Common hairpin loop motifs

## Output Files

### Downloadable Results

1. **motif_enrichment.csv**: Complete table of motifs with statistics
2. **structures.txt**: Secondary structure predictions for all sequences

### File Formats

**Motif Table (CSV):**
```csv
Motif,Length,Sequences,Frequency,Enrichment,p-value,Adj. p-value,GC%
GGATCC,6,15,75.0%,4.50x,1.23e-05,3.45e-04,66.7%
```

**Structures (Text):**
```
>Seq1
Sequence:  GGAUCCGGAUCCGGAUCC...
Structure: (((....)))........
MFE:       -12.34 kcal/mol
```

## Troubleshooting

### Common Issues

**Problem**: "ViennaRNA not found"
- **Solution**: Install ViennaRNA package (see Installation section)

**Problem**: "No motifs found"
- **Solution**: Try decreasing min_sequences or adjusting motif length range

**Problem**: Analysis takes too long
- **Solution**: Disable structure prediction for faster results, or reduce sequence count

**Problem**: Primers not detected
- **Solution**: Check primer sequences match your data; tool allows up to 2 mismatches

### Performance Tips

- For >50 sequences, consider disabling structure prediction initially
- Use stricter parameters (higher min_sequences) for large datasets
- Structure analysis time scales with number of sequences × sequence length

## Algorithm Details

### Motif Discovery

1. **k-mer Extraction**: Identifies all k-mers of specified lengths in random regions
2. **Frequency Counting**: Counts occurrences across all sequences
3. **Filtering**: Retains motifs appearing in ≥ min_sequences

### Statistical Testing

1. **Null Model**: Assumes equal base probabilities (25% each)
2. **Binomial Test**: Tests if observed frequency exceeds expected
   - P(motif at position) = (0.25)^k
   - P(motif in sequence) = 1 - (1 - p_position)^n_positions
   - Test: Does motif appear in more sequences than expected?
3. **FDR Correction**: Benjamini-Hochberg procedure to control false discoveries

### Structure Prediction

Uses ViennaRNA's RNA.fold() with default energy parameters:
- Temperature: 37°C
- Algorithm: Minimum Free Energy (MFE) folding
- Output: Dot-bracket notation + MFE value

## Technical Details

### Dependencies

- **streamlit**: Web application framework
- **biopython**: Sequence manipulation and analysis
- **pandas**: Data handling and export
- **numpy**: Numerical computations
- **scipy**: Statistical tests
- **matplotlib/seaborn**: Static visualizations
- **plotly**: Interactive plots
- **logomaker**: Sequence logo generation
- **ViennaRNA**: RNA structure prediction

### Performance

- Typical analysis time: 30-120 seconds
- Scales linearly with number of sequences
- Memory usage: ~100-500 MB for typical datasets

## Citation

If you use AptaMotif in your research, please cite:

**ViennaRNA Package:**
Lorenz, R., et al. (2011). ViennaRNA Package 2.0. Algorithms Mol Biol 6:26

**Biopython:**
Cock, P.J.A., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25(11):1422-1423

## License

This tool is provided for research use. Please check individual package licenses for dependencies.

## Support

For questions, issues, or feature requests:
1. Check the "Help" tab in the application
2. Review this README
3. Contact your bioinformatics core facility

## Version History

**v1.0** (February 2026)
- Initial release
- Core motif discovery and enrichment analysis
- ViennaRNA integration for structure prediction
- Interactive web interface
- Configurable primers and parameters

---

**Built with** ❤️ **using Python, Streamlit, and ViennaRNA**
