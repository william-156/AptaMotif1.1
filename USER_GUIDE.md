# AptaMotif Analyzer - Complete User Guide

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Getting Started](#getting-started)
4. [Detailed Workflow](#detailed-workflow)
5. [Understanding Results](#understanding-results)
6. [Advanced Features](#advanced-features)
7. [Troubleshooting](#troubleshooting)
8. [FAQs](#faqs)

---

## Introduction

AptaMotif Analyzer is a comprehensive tool designed for RNA aptamer researchers to:
- Discover enriched sequence motifs in aptamer pools
- Perform statistical analysis with proper multiple testing correction
- Predict secondary structures using ViennaRNA
- Generate publication-quality visualizations
- Manage multiple pool configurations for different experiments

### When to Use This Tool

Use AptaMotif Analyzer when you have:
- ✅ Sequenced aptamer clones from a selection experiment
- ✅ Clean, quality-controlled sequences (QC passed, good chromatograms)
- ✅ At least 10-20 sequences for meaningful statistical analysis
- ✅ Known primer sequences flanking your random region

---

## Installation

### Step 1: Install Python Dependencies

```bash
# Option A: Using requirements.txt (recommended)
pip install -r requirements.txt --break-system-packages

# Option B: Manual installation
pip install streamlit biopython logomaker plotly scipy pandas numpy matplotlib seaborn statsmodels ViennaRNA --break-system-packages
```

### Step 2: Verify Installation

```bash
python3 test_modules.py
```

You should see:
```
✅ All tests passed!
Ready to launch AptaMotif Analyzer!
```

### Step 3: Launch the Application

```bash
# Option A: Using launch script
./launch.sh

# Option B: Direct streamlit command
streamlit run aptamer_motif_analyzer.py
```

The application will open in your web browser at `http://localhost:8501`

---

## Getting Started

### Quick Start (5 minutes)

1. **Test with Example Data**
   - Launch the application
   - Go to "Input Sequences" tab
   - Click "Upload File" and select `example_sequences.fa`
   - Click "Process Uploaded File"

2. **Run Analysis**
   - Go to "Motif Analysis" tab
   - Click "Run Motif Analysis"
   - View the enriched motifs table and visualizations

3. **View Structures**
   - Go to "Structure Analysis" tab
   - Click "Run Structure Prediction"
   - Explore predicted RNA structures

### Your First Real Analysis

1. **Prepare Your Sequences**
   - Ensure sequences include primer regions
   - Format as FASTA or plain text
   - Verify sequence quality

2. **Configure Pool Settings**
   - In sidebar, click "Add/Edit Pool Configuration"
   - Enter your forward primer sequence
   - Enter the reverse complement sequence after the N region
   - Set the random region length
   - Save configuration

3. **Input Your Sequences**
   - Paste or upload your sequences
   - Process sequences
   - Verify correct number loaded

4. **Run Analysis**
   - Adjust parameters as needed
   - Run motif analysis
   - Run structure prediction
   - Export results

---

## Detailed Workflow

### Step 1: Pool Configuration

#### Understanding Pool Structure

Your aptamer library has this structure:
```
5' - [Forward Primer] - [Random Region] - [Reverse Complement] - 3'
```

Example (N71 pool):
```
5' - TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATT - N71 - AGATAGTAAGTGCAATCT - 3'
                    (Forward Primer)                          (Rev Comp)
```

#### Adding a New Configuration

1. In sidebar, expand "Add/Edit Pool Configuration"
2. Enter configuration name (e.g., "My_N50_Pool")
3. Paste forward primer sequence (5'→3')
4. Paste reverse complement sequence
5. Set random region length
6. Add description
7. Click "Save Configuration"

#### Sharing Configurations

**Export**: 
- Expand "Export/Import Configurations"
- Click "Download Configurations"
- Share JSON file with team

**Import**:
- Click "Upload Configuration File"
- Select JSON file
- Configurations automatically merge

### Step 2: Sequence Input

#### Format Requirements

**FASTA Format** (recommended):
```
>Clone_001
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTACGTACGTAGATAGTAAGTGCAATCT
>Clone_002
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTGGCTAGCTAGATAGTAAGTGCAATCT
```

**Plain Text** (one per line):
```
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTACGTACGTAGATAGTAAGTGCAATCT
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTGGCTAGCTAGATAGTAAGTGCAATCT
```

#### Input Methods

**Method 1: Paste**
1. Go to "Input Sequences" tab
2. Select "Paste Sequences"
3. Copy/paste your sequences
4. Click "Process Sequences"

**Method 2: Upload**
1. Go to "Input Sequences" tab
2. Select "Upload File"
3. Choose your .fasta, .fa, or .txt file
4. Click "Process Uploaded File"

#### Quality Checks

The tool will:
- ✅ Parse sequence format
- ✅ Extract random regions
- ✅ Warn if primers not found
- ✅ Display sequence summary

### Step 3: Motif Discovery

#### Parameter Selection

**Minimum Motif Length** (default: 5 bp)
- Shorter motifs: More hits, less specific
- Longer motifs: Fewer hits, more specific
- Recommendation: Start with 5 bp

**Maximum Motif Length** (default: 15 bp)
- Consider expected binding site size
- Larger ranges = more computation time
- Recommendation: 5-15 bp covers most motifs

**Minimum Occurrences** (default: 2)
- Higher values: Focus on highly conserved motifs
- Lower values: Discover rare but enriched motifs
- Recommendation: 2 for discovery, 3+ for validation

**FDR Threshold** (default: 0.05)
- 0.05 = 5% false discovery rate
- Lower values (0.01): More stringent
- Higher values (0.10): More permissive
- Recommendation: 0.05 for initial analysis

#### Running Analysis

1. Set parameters in sidebar
2. Click "Run Motif Analysis"
3. Wait for processing (usually <30 seconds)
4. View results table and visualizations

### Step 4: Structure Prediction

#### Settings

**Enable Structure Analysis**: Check box in sidebar
**Temperature**: Set folding temperature (default: 37°C)
- Use your selection temperature
- Affects MFE calculations

#### Running Prediction

1. Enable structure analysis
2. Click "Run Structure Prediction"
3. Processing time: ~1-5 seconds per sequence
4. View individual structures and MFE values

#### ViennaRNA Status

The tool will indicate:
- ✅ "Using ViennaRNA" - Full algorithm
- ⚠️ "Using fallback algorithm" - Simplified prediction

For full ViennaRNA:
```bash
conda install -c bioconda viennarna
```

---

## Understanding Results

### Motif Enrichment Table

**Columns Explained**:

- **Motif**: The sequence motif (e.g., "GGCTAG")
- **Length**: Number of bases
- **Count**: How many sequences contain it
- **Expected_Count**: Expected by random chance
- **Fold_Enrichment**: Count / Expected_Count
  - >1 means enriched
  - >2 is moderate enrichment
  - >5 is strong enrichment
- **Frequency**: Proportion of sequences (0-1)
- **P_value**: Raw statistical significance
- **FDR**: Corrected p-value (controls false discoveries)
- **Significant**: TRUE if FDR < threshold
- **Sequences**: Which clones contain the motif

### Interpreting Statistics

#### P-values
- Tests null hypothesis: "Motif occurs by random chance"
- Small p-value (< 0.05) = unlikely to be random
- **Important**: Don't use raw p-values for interpretation!

#### FDR (False Discovery Rate)
- Controls proportion of false positives
- FDR < 0.05 means < 5% of "significant" motifs are false positives
- **Use FDR, not p-values, for interpretation**

#### Fold Enrichment
- Practical measure of enrichment strength
- Consider both FDR and fold enrichment
- Example: FDR=0.001, Fold=10x → Strong candidate

### Visualizations

#### 1. Heatmap
- **Rows**: Individual sequences
- **Columns**: Motifs
- **Colors**: Blue = present, Gray = absent
- **Use**: See which sequences share motifs

#### 2. Sequence Logos
- **Height**: Information content (conservation)
- **Position**: Relative to motif
- **Colors**: A=green, C=blue, G=orange, T=red
- **Yellow highlight**: Core motif region
- **Use**: See conservation and flanking preferences

#### 3. Enrichment Bar Plot
- **X-axis**: Fold enrichment
- **Y-axis**: Motifs (sorted by significance)
- **Colors**: Blue=significant, Purple=not significant
- **Use**: Quick overview of top motifs

#### 4. Volcano Plot
- **X-axis**: Fold enrichment
- **Y-axis**: -log₁₀(FDR)
- **Red line**: FDR threshold
- **Use**: See enrichment vs. significance trade-off

### Structure Results

**For Each Sequence**:
- **Sequence**: RNA sequence (DNA converted to RNA)
- **Structure**: Dot-bracket notation
  - `(` and `)` = base pairs
  - `.` = unpaired bases
- **MFE**: Minimum free energy (kcal/mol)
  - More negative = more stable
  - Typical range: -5 to -30 kcal/mol

**Example Structure**:
```
Sequence:  GGCUAGCUAGC
Structure: (((...)))..
MFE:       -8.6 kcal/mol
```

This shows:
- First 3 bases paired with last 3 bases (stem)
- Middle 3 bases form a loop
- Last 2 bases unpaired

---

## Advanced Features

### Multiple Pool Management

**Use Case**: Different aptamer libraries in your lab

**Workflow**:
1. Create configuration for each pool
2. Export configurations to JSON
3. Share with lab members
4. Switch between configurations as needed

### Batch Analysis

**For Multiple Experiments**:
1. Create separate folders for each experiment
2. Keep pool configurations
3. Process each experiment separately
4. Compare results across experiments

### Custom Analysis Parameters

**GC Content Adjustment**:
- If your sequences have biased GC content
- Advanced users can modify `statistics_module.py`
- Use `adjust_for_gc_bias()` method

**Permutation Testing**:
- More robust but slower
- Use `permutation_test()` in statistics module
- Good for small sample sizes (<10 sequences)

### Exporting for Publications

**Motif Table**:
- CSV format, Excel-compatible
- Include in supplementary materials
- Reference in methods section

**Visualizations**:
- High-resolution figures
- Right-click and save images
- Adjust figure sizes in code if needed

**Structure Files**:
- Standard Vienna format
- Compatible with structure visualization tools
- Can be imported into RNA structure programs

---

## Troubleshooting

### Common Issues

#### 1. "Forward primer not found in sequence"

**Causes**:
- Primer sequence mismatch
- Sequences already trimmed
- Sequencing errors at primer region

**Solutions**:
- Double-check primer sequence in configuration
- Verify sequences include full primer
- Check for lowercase/uppercase mismatches
- Try allowing 1-2 mismatches (requires code modification)

#### 2. "No significant motifs found"

**Causes**:
- Sample size too small (< 10 sequences)
- No real enrichment in your selection
- FDR threshold too stringent
- Random region truly random (no selection pressure)

**Solutions**:
- Increase sample size (sequence more clones)
- Lower FDR threshold to 0.10 for exploration
- Check fold enrichment (may be enriched but not significant)
- Verify selection worked (do you have binding clones?)

#### 3. Structure prediction errors

**Causes**:
- ViennaRNA not installed
- Invalid RNA sequence
- Sequence too long

**Solutions**:
- Tool will use fallback algorithm
- For full ViennaRNA: `conda install -c bioconda viennarna`
- Check for invalid characters in sequence
- Split very long sequences

#### 4. Slow performance

**Causes**:
- Large number of sequences (> 100)
- Wide motif length range
- Many unique sequences

**Solutions**:
- Start with smaller motif length range (5-10)
- Process in batches
- Use more powerful computer
- Reduce max_motif_length

### Getting Help

1. **Check the README**: Basic information and setup
2. **Review this guide**: Detailed workflows and explanations
3. **Test with examples**: Use example_sequences.fa
4. **Check parameters**: Verify all settings are appropriate
5. **Review error messages**: Often indicate the exact problem

---

## FAQs

### General Questions

**Q: How many sequences do I need?**
A: Minimum 10, recommended 20-50 for good statistical power. More sequences = better detection of enriched motifs.

**Q: Can I analyze DNA aptamers?**
A: Currently optimized for RNA. DNA aptamers can be analyzed for motifs, but structure prediction assumes RNA folding. Modify temperature and energy parameters for DNA.

**Q: What file formats are supported?**
A: FASTA (.fasta, .fa) and plain text (.txt). One sequence per line in plain text, or standard FASTA format with headers.

**Q: Can I run this on a server?**
A: Yes! Deploy with: `streamlit run aptamer_motif_analyzer.py --server.port 8501 --server.address 0.0.0.0`

### Analysis Questions

**Q: What's a good FDR threshold?**
A: 0.05 (5%) is standard. Use 0.01 for high confidence, 0.10 for exploratory analysis.

**Q: Should I filter by fold enrichment?**
A: Yes! Consider both FDR < 0.05 AND Fold_Enrichment > 2. Statistical significance alone doesn't mean biological relevance.

**Q: How do I know if a motif is real?**
A: Check:
1. Low FDR (< 0.05)
2. High fold enrichment (> 2)
3. Present in multiple sequences
4. Makes biological sense
5. Validation experiments

**Q: Can overlapping motifs be meaningful?**
A: Yes! The tool removes exact redundancy, but overlapping motifs may represent:
- Extended binding sites
- Variations of the same motif
- Different parts of the same structural element

### Technical Questions

**Q: Why use Binomial test instead of Fisher's exact?**
A: Binomial test is more appropriate for this design - we're testing if a motif appears more often than expected by random chance in independent sequences.

**Q: What's the difference between p-value and FDR?**
A: P-value tests one motif. FDR corrects for testing many motifs simultaneously (multiple testing correction using Benjamini-Hochberg procedure).

**Q: Why does the tool remove redundant motifs?**
A: If "GGCTAG" and "GGCTAGC" appear in the exact same sequences, the longer one is kept. This reduces noise while preserving information.

**Q: How accurate is structure prediction?**
A: ViennaRNA is industry-standard for MFE prediction (~70-80% accurate). Remember:
- Prediction is for isolated sequences
- In vivo structures may differ
- MFE is one possible structure
- Validate important structures experimentally

### Workflow Questions

**Q: Should I analyze before or after clustering?**
A: Before! Clustering removes information. This tool handles redundancy appropriately.

**Q: Can I compare enrichment between rounds?**
A: Yes, but analyze each round separately. Look for motifs that increase in frequency and significance across rounds.

**Q: What if my primers have degeneracies?**
A: Use the most common variant as your primer sequence. You may need to manually trim sequences if variation is high.

**Q: Can I use this for other applications?**
A: Yes! Any experiment with:
- Sequences with flanking regions
- Need to find enriched motifs
- Statistical significance testing
Examples: DNA sequencing, peptide display, etc.

---

## Citation and Acknowledgments

### Software Citation

If you use AptaMotif Analyzer in your research, please cite:

**ViennaRNA Package**:
- Lorenz, R., et al. (2011). "ViennaRNA Package 2.0." Algorithms for Molecular Biology, 6:26.

**Key Libraries**:
- Streamlit: https://streamlit.io
- Biopython: https://biopython.org
- logomaker: https://logomaker.readthedocs.io

### Acknowledgments

This tool was developed to support RNA aptamer research and combines best practices from:
- Bioinformatics pipeline development
- Statistical genomics
- RNA structure prediction
- Data visualization

---

## Appendix: Example Analysis

### Complete Workflow Example

**Scenario**: You have 20 sequences from Round 5 of SELEX against a protein target.

**Step 1**: Prepare sequences
```
>R5_Clone_01
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTGGCTAGCTAGCGGGCCCAGATAGTAAGTGCAATCT
>R5_Clone_02
...
```

**Step 2**: Configure pool
- Forward: TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATT
- Rev Comp: AGATAGTAAGTGCAATCT
- N region: 71 bp

**Step 3**: Run analysis
- Parameters: 5-15 bp motifs, ≥2 occurrences, FDR < 0.05
- Results: 12 motifs found, 5 significant

**Step 4**: Interpret results
Top motif: GGCTAGC
- FDR: 0.001
- Fold enrichment: 8.2x
- Present in: 14/20 sequences (70%)
- **Conclusion**: Strong candidate for binding site

**Step 5**: Validate
- Design mutations of GGCTAGC
- Test binding affinity
- Confirm importance

---

## Version History

**Version 1.0** (February 2026)
- Initial release
- Motif discovery with statistical analysis
- ViennaRNA structure prediction
- Interactive visualizations
- Pool configuration management
- Export functionality

---

## Support and Updates

For questions, issues, or feature requests:
- Review this guide thoroughly
- Test with example sequences
- Check parameter settings
- Verify input format

**Future Development**:
- Batch processing multiple experiments
- Advanced structural motif analysis
- Integration with binding data
- Comparison tools for multiple rounds
- Extended statistical tests

---

*This tool is provided for research purposes. Always validate computational results with experimental data.*
