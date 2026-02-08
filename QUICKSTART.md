# Quick Start Guide for AptaMotif

## Installation (5 minutes)

### Step 1: Install ViennaRNA
```bash
# Ubuntu/Debian
sudo apt-get update && sudo apt-get install vienna-rna

# macOS
brew install viennarna
```

### Step 2: Install Python packages
```bash
cd aptamotif
pip install -r requirements.txt
```

## Running the Application

### Launch the web app:
```bash
streamlit run app.py
```

The app will open in your browser at: http://localhost:8501

## First Analysis (2 minutes)

### 1. Load Example Data
- Click "Upload File" in the Input & Analysis tab
- Select `example_sequences.fasta`
- Or copy-paste the contents

### 2. Run Analysis
- Keep default settings
- Click "🔬 Run Analysis"
- Wait ~30 seconds

### 3. View Results
- Switch to "Results" tab
- Explore the motif table
- Check visualizations

## Expected Results from Example Data

The example file contains 15 sequences with three enriched motifs:

1. **GGATCC** - appears in ~47% of sequences (highly significant)
2. **CCAATT** - appears in ~40% of sequences (significant)  
3. **GGTTAA** - appears in ~27% of sequences (significant)
4. **ACGTGC** - appears in ~13% of sequences (may be significant)

You should see:
- 50+ total motifs discovered
- 3-5 significantly enriched motifs (adj. p-value < 0.05)
- Fold enrichments between 2x - 10x
- Clear patterns in the heatmap

## Customizing for Your Data

### For a Different SELEX Library:

1. **Update Primers** (sidebar → Edit Primers):
   ```
   Forward Primer: YOUR_FWD_PRIMER
   Reverse Region: YOUR_REV_REGION
   N-Region Length: YOUR_LENGTH
   ```

2. **Adjust Parameters**:
   - Increase "Min Sequences" for stricter filtering
   - Adjust motif length range based on expected motifs
   - Change FDR threshold for more/fewer significant calls

3. **Input Your Sequences**:
   - FASTA format (with >headers) OR
   - Plain text (one sequence per line)
   - Include full sequences with primers

## Troubleshooting

**No motifs found?**
- Decrease "Min Sequences" to 2
- Check that primers are correctly specified
- Verify sequences include primer regions

**Analysis too slow?**
- Uncheck "Perform Secondary Structure Analysis"
- This speeds up analysis 5-10x

**Primers not detected?**
- Tool allows up to 2 mismatches
- Check for typos in primer sequences
- Ensure sequences are in correct orientation

## Key Features to Try

1. **Interactive Heatmap**: Hover over cells to see details
2. **Volcano Plot**: Visualize enrichment vs. significance
3. **Sequence Logos**: See base preferences in motifs
4. **Structure Analysis**: View predicted RNA structures
5. **CSV Export**: Download results for further analysis

## Next Steps

After your first analysis:

1. Try different parameter combinations
2. Compare results across SELEX rounds
3. Look for structural motifs in enriched sequences
4. Export data for downstream analysis

## Getting Help

- Check the **Help** tab in the application
- Read the full README.md
- Review example outputs

---

**Ready to discover your aptamer motifs!** 🧬
