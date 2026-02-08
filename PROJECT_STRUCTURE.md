# AptaMotif Project Structure

## Directory Layout

```
aptamotif/
├── app.py                      # Main Streamlit web application
├── sequence_parser.py          # Sequence parsing and primer trimming
├── motif_finder.py             # Motif discovery engine
├── statistics.py               # Statistical analysis (binomial test, FDR)
├── structure_analyzer.py       # ViennaRNA integration for structure prediction
├── visualizer.py               # Visualization tools (heatmaps, logos, plots)
├── requirements.txt            # Python package dependencies
├── README.md                   # Comprehensive documentation
├── QUICKSTART.md              # Quick start guide
├── test_installation.py       # Installation verification script
├── example_sequences.fasta    # Example data for testing
└── PROJECT_STRUCTURE.md       # This file
```

## Module Descriptions

### Core Analysis Modules

#### `sequence_parser.py`
**Purpose**: Handle input sequences and extract random regions

**Key Classes**:
- `SequenceParser`: Parse FASTA/text, trim primers, prepare for analysis

**Main Functions**:
- `parse_sequences()`: Parse input text
- `extract_random_region()`: Identify and extract N-region between primers
- `prepare_for_folding()`: Convert DNA→RNA for structure prediction
- `get_stats()`: Calculate sequence statistics

**Usage**:
```python
parser = SequenceParser(
    forward_primer="ACGT...",
    reverse_primer_region="TGCA...",
    molecule_type="RNA"
)
sequences = parser.parse_sequences(input_text)
sequences = parser.extract_random_region(sequences)
```

---

#### `motif_finder.py`
**Purpose**: Discover enriched k-mer motifs

**Key Classes**:
- `MotifFinder`: Find and analyze motif patterns

**Main Functions**:
- `find_all_kmers()`: Extract all k-mers of specified length
- `find_enriched_motifs()`: Identify motifs appearing in ≥min_sequences
- `create_motif_matrix()`: Generate presence/absence matrix
- `filter_redundant_motifs()`: Remove overlapping/redundant patterns

**Usage**:
```python
finder = MotifFinder(min_length=5, max_length=15, min_sequences=2)
motifs = finder.find_enriched_motifs(sequences)
```

---

#### `statistics.py`
**Purpose**: Calculate statistical significance

**Key Classes**:
- `MotifStatistics`: Statistical tests and corrections

**Main Functions**:
- `calculate_binomial_pvalue()`: Test motif enrichment
- `apply_fdr_correction()`: Benjamini-Hochberg FDR correction
- `calculate_enrichment_score()`: Fold enrichment over expected
- `calculate_base_composition()`: Actual base frequencies

**Statistical Methods**:
1. **Binomial Test**: P(motif) = (0.25)^k per position
2. **FDR Correction**: Controls false discovery rate
3. **Enrichment**: Observed / Expected ratio

**Usage**:
```python
stats = MotifStatistics()
motifs = stats.add_pvalues_to_motifs(motifs, mean_length)
motifs = stats.apply_fdr_correction(motifs, method='fdr_bh')
```

---

#### `structure_analyzer.py`
**Purpose**: Predict RNA/DNA secondary structures

**Key Classes**:
- `StructureAnalyzer`: ViennaRNA wrapper

**Main Functions**:
- `fold_sequence()`: Predict MFE structure for one sequence
- `fold_all_sequences()`: Batch structure prediction
- `find_consensus_structure()`: Identify common structural features
- `identify_structural_motifs()`: Find recurring hairpins/loops

**Output Format**:
- Dot-bracket notation: `(((...)))`
- Minimum Free Energy (MFE) in kcal/mol
- Structural element statistics

**Usage**:
```python
analyzer = StructureAnalyzer(molecule_type="RNA", temperature=37.0)
sequences = analyzer.fold_all_sequences(sequences)
consensus = analyzer.find_consensus_structure(sequences)
```

---

#### `visualizer.py`
**Purpose**: Create visualizations

**Key Classes**:
- `MotifVisualizer`: Generate plots and heatmaps

**Main Functions**:
- `create_motif_heatmap()`: Presence/absence heatmap (Plotly)
- `create_sequence_logo()`: Sequence logo (Matplotlib + Logomaker)
- `create_enrichment_barplot()`: Top enriched motifs
- `create_pvalue_volcano_plot()`: Enrichment vs. significance
- `create_structure_comparison()`: MFE vs. % paired

**Visualization Types**:
- Interactive Plotly figures (web app)
- Static Matplotlib figures (export)
- Sequence logos (information content)
- Statistical plots (volcano, distributions)

**Usage**:
```python
viz = MotifVisualizer()
fig = viz.create_motif_heatmap(matrix, motif_names, seq_names)
```

---

### Application Layer

#### `app.py`
**Purpose**: Streamlit web interface

**Structure**:
1. **Configuration Sidebar**:
   - Molecule type selection
   - Primer editing
   - Analysis parameters
   - Structure options

2. **Main Tabs**:
   - **Input & Analysis**: Sequence input, run analysis
   - **Results**: Tables, visualizations, downloads
   - **Help**: Documentation and guidance

3. **Session State**:
   - Stores analysis results between interactions
   - Prevents recomputation on UI updates

**Key Features**:
- Drag-and-drop file upload
- Real-time parameter adjustment
- Interactive visualizations
- CSV/text export

---

## Data Flow

```
User Input
    ↓
[sequence_parser.py]
    ↓ Parsed sequences with extracted N-regions
[motif_finder.py]
    ↓ Discovered motifs with counts
[statistics.py]
    ↓ Motifs with p-values and enrichment scores
[structure_analyzer.py] (optional)
    ↓ Sequences with predicted structures
[visualizer.py]
    ↓ Interactive plots and heatmaps
[app.py]
    ↓
Display Results
```

## Dependencies

### Required Python Packages

| Package | Version | Purpose |
|---------|---------|---------|
| streamlit | 1.31.0 | Web application framework |
| biopython | 1.83 | Sequence analysis |
| pandas | 2.2.0 | Data manipulation |
| numpy | 1.26.3 | Numerical computing |
| scipy | 1.12.0 | Statistical tests |
| matplotlib | 3.8.2 | Static plots |
| seaborn | 0.13.2 | Statistical visualization |
| logomaker | 0.8 | Sequence logos |
| plotly | 5.18.0 | Interactive plots |
| viennarna | 2.6.4 | RNA structure prediction |
| statsmodels | 0.14.1 | Additional statistics |

### External Dependencies

- **ViennaRNA Package**: System-level installation required
  - Ubuntu/Debian: `apt-get install vienna-rna`
  - macOS: `brew install viennarna`

## Configuration

### Default Settings

**Primers (N71 pool)**:
```python
forward_primer = "TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATT"
reverse_primer_region = "AGATAGTAAGTGCAATCT"
n_region_length = 71
```

**Analysis Parameters**:
```python
min_motif_length = 5
max_motif_length = 15
min_sequences = 2
fdr_threshold = 0.05
```

**Structure Settings**:
```python
molecule_type = "RNA"
temperature = 37.0
```

### Customization

All settings can be modified through the Streamlit sidebar without code changes.

## Testing

### Quick Installation Test
```bash
python test_installation.py
```

Checks:
- ✓ All Python packages installed
- ✓ Modules can be imported
- ✓ ViennaRNA is functional
- ✓ Example data exists
- ✓ Quick analysis runs successfully

### Example Data

`example_sequences.fasta`: 15 sequences with known enriched motifs
- GGATCC (appears in ~47% - highly significant)
- CCAATT (appears in ~40% - significant)
- GGTTAA (appears in ~27% - significant)

## Performance Characteristics

### Time Complexity

- **Motif Discovery**: O(n * L * k_max)
  - n = number of sequences
  - L = average sequence length
  - k_max = maximum motif length

- **Statistical Testing**: O(m)
  - m = number of unique motifs

- **Structure Prediction**: O(n * L^3)
  - Cubic in sequence length (RNA folding algorithm)

### Typical Performance

| Dataset Size | Motif Discovery | Structure Prediction | Total |
|--------------|----------------|---------------------|--------|
| 10 sequences | 5-10 seconds | 10-20 seconds | ~30 sec |
| 50 sequences | 15-30 seconds | 45-90 seconds | ~2 min |
| 100 sequences | 30-60 seconds | 90-180 seconds | ~4 min |

### Memory Usage

- Typical: 100-500 MB
- Scales with: number of sequences × sequence length × motif count

## Error Handling

### Common Errors

1. **Primers Not Found**:
   - Allows up to 2 mismatches
   - Falls back to full sequence if primers missing

2. **No Motifs Found**:
   - Check min_sequences threshold
   - Verify sequence quality

3. **ViennaRNA Errors**:
   - Invalid sequence characters
   - Sequence too long (>10,000 nt)

### Logging

- Errors displayed in Streamlit UI
- Warnings shown as info boxes
- Progress indicators for long operations

## Extension Points

### Adding New Features

1. **New Motif Types**:
   - Edit `motif_finder.py`
   - Add new discovery algorithms

2. **Additional Statistics**:
   - Edit `statistics.py`
   - Implement new tests

3. **Custom Visualizations**:
   - Edit `visualizer.py`
   - Add plotting functions

4. **Alternative Structure Predictors**:
   - Edit `structure_analyzer.py`
   - Wrap other tools (RNAfold, mfold)

## Best Practices

### For Users

1. Quality-check sequences before analysis
2. Start with default parameters
3. Review primer sequences carefully
4. Export results regularly
5. Use example data to test

### For Developers

1. Keep modules independent
2. Use type hints
3. Document functions thoroughly
4. Test with edge cases
5. Profile performance for optimizations

## Troubleshooting

See README.md "Troubleshooting" section for:
- Installation issues
- Runtime errors
- Performance problems
- Data format questions

## Future Enhancements

Potential additions:
- [ ] Batch processing multiple experiments
- [ ] Comparison across SELEX rounds
- [ ] Position-specific motif analysis
- [ ] Advanced structure clustering
- [ ] Integration with sequence databases
- [ ] API endpoint for programmatic access

---

**Version**: 1.0
**Last Updated**: February 2026
