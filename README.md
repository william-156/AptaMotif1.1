# AptaMotif Analyzer

A comprehensive web-based tool for RNA aptamer motif discovery and secondary structure analysis.

## Features

- **Motif Discovery**: Identifies enriched sequence motifs (≥5 bp) shared across multiple aptamer clones
- **Statistical Analysis**: Binomial tests with FDR correction for significance
- **Secondary Structure Prediction**: ViennaRNA integration for RNA folding
- **Interactive Visualizations**: Heatmaps, sequence logos, enrichment plots, and volcano plots
- **Configurable Pools**: Save and manage multiple primer configurations
- **Export Capabilities**: CSV tables, structure files, and plots

## Installation

### Requirements
- Python 3.8+
- pip

### Quick Start

1. Install required packages:
```bash
pip install streamlit biopython logomaker plotly scipy pandas numpy matplotlib seaborn statsmodels ViennaRNA --break-system-packages
```

2. Run the application:
```bash
streamlit run aptamer_motif_analyzer.py
```

3. Open your browser to the URL shown (typically http://localhost:8501)

## Usage

### 1. Configure Your Pool

- Select or create a pool configuration in the sidebar
- Default configuration is for N71 pools with standard primers
- Add custom configurations for different aptamer pools
- Export/import configurations for team sharing

### 2. Input Sequences

- Paste sequences in FASTA format or plain text (one per line)
- Or upload a file (.fasta, .fa, .txt)
- DNA sequences will be automatically converted to RNA for structure analysis

### 3. Run Analysis

**Motif Analysis**:
- Set parameters: minimum/maximum motif length, minimum occurrences, FDR threshold
- Click "Run Motif Analysis"
- View enriched motifs with statistics

**Structure Analysis**:
- Enable structure prediction in sidebar
- Set temperature for folding
- Click "Run Structure Prediction"
- View predicted structures with MFE values

### 4. Export Results

- Download motif enrichment table (CSV)
- Download structure predictions (TXT)
- Export visualizations

## Input Format

### FASTA Format
```
>Clone_1
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTACGTACGTACGTAGATAGTAAGTGCAATCT
>Clone_2
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTGGCTAGCTAGCTAGATAGTAAGTGCAATCT
```

### Plain Text (one sequence per line)
```
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTACGTACGTACGTAGATAGTAAGTGCAATCT
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTGGCTAGCTAGCTAGATAGTAAGTGCAATCT
```

## Pool Configuration

The default N71 pool configuration:
- **Forward Primer**: 5' -TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATT- 3'
- **Random Region**: N71
- **Reverse Complement**: AGATAGTAAGTGCAATCT

To add custom pools:
1. Click "Add/Edit Pool Configuration" in sidebar
2. Enter primer sequences and parameters
3. Save configuration
4. Share via Export/Import feature

## Analysis Methods

### Motif Discovery
- Extracts all k-mers (k = 5 to 15 bp by default) from random regions
- Identifies motifs present in ≥2 sequences (configurable)
- Removes redundant sub-motifs

### Statistical Testing
- **Binomial Test**: Tests if motif frequency exceeds random expectation
- **Expected Probability**: Based on (0.25)^k for k-mer length
- **FDR Correction**: Benjamini-Hochberg procedure for multiple testing
- **Fold Enrichment**: Observed count / Expected count

### Structure Prediction
- Uses ViennaRNA RNA.fold() for minimum free energy (MFE) structures
- Converts DNA to RNA (T→U) automatically
- Identifies common structural motifs (stems, loops, bulges)
- Calculates consensus structures

## Visualizations

1. **Heatmap**: Motif presence/absence across sequences
2. **Sequence Logos**: Position-specific nucleotide frequency
3. **Enrichment Bar Plot**: Fold enrichment of top motifs
4. **Volcano Plot**: Enrichment vs. significance
5. **Length Distribution**: Distribution of motif sizes

## Output Files

### Motif Enrichment CSV
Columns:
- Motif: Sequence motif
- Length: Motif length (bp)
- Count: Number of sequences containing motif
- Expected_Count: Expected occurrences by chance
- Fold_Enrichment: Observed/Expected ratio
- Frequency: Proportion of sequences
- P_value: Raw p-value from binomial test
- FDR: False discovery rate (corrected p-value)
- Significant: Boolean (True if FDR < threshold)
- Sequences: List of sequence IDs containing motif

### Structure Predictions TXT
Format:
```
>Clone_1
Sequence: AUGCUAGCUAGC...
Structure: (((...)))....
MFE: -12.34 kcal/mol
```

## Tips for Best Results

1. **Quality Control**: Ensure sequences are clean and properly trimmed
2. **Primer Verification**: Double-check primer sequences match your pool
3. **Parameter Tuning**: 
   - Lower FDR threshold (e.g., 0.01) for more stringent results
   - Increase min_occurrences for highly conserved motifs
   - Adjust motif length range based on expected binding sites
4. **Structure Analysis**: Use physiologically relevant temperature
5. **Interpretation**: Consider both statistical significance and biological relevance

## Troubleshooting

**Issue**: "Forward primer not found"
- Solution: Verify primer sequence matches exactly
- Check for sequencing quality issues

**Issue**: "No significant motifs found"
- Solution: Lower FDR threshold or min_occurrences
- Check if sample size is sufficient (need ≥10 sequences for good power)

**Issue**: ViennaRNA errors
- Solution: Tool will fall back to simple structure prediction
- For full ViennaRNA: `conda install -c bioconda viennarna`

## Citation

If you use this tool in your research, please cite:
- ViennaRNA Package: Lorenz et al. (2011) "ViennaRNA Package 2.0" Algorithms Mol Biol 6:26

## Support

For questions, issues, or feature requests:
- Check the documentation in each tab
- Review parameter tooltips
- Export configurations for reproducibility

## Version

Version 1.0
Last Updated: February 2026

## License

This tool is provided for research use in aptamer discovery and analysis.
