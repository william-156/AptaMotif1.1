# 🧬 AptaMotif Analyzer - Quick Start Guide

## What You Have

A complete web-based tool for RNA aptamer motif analysis with:
- ✅ Motif discovery and statistical enrichment analysis
- ✅ ViennaRNA secondary structure prediction
- ✅ Interactive visualizations (heatmaps, logos, plots)
- ✅ Configurable pool management
- ✅ CSV/file export capabilities

---

## 🚀 Installation & Launch (5 Minutes)

### Step 1: Install Python Packages

Open a terminal and run:

```bash
pip install -r requirements.txt --break-system-packages
```

**Or install individually:**
```bash
pip install streamlit biopython logomaker plotly scipy pandas numpy matplotlib seaborn statsmodels ViennaRNA --break-system-packages
```

### Step 2: Verify Installation

```bash
python3 test_modules.py
```

✅ You should see: "All tests passed!"

### Step 3: Launch the Application

**Option A - Using the launch script:**
```bash
./launch.sh
```

**Option B - Direct command:**
```bash
streamlit run aptamer_motif_analyzer.py
```

### Step 4: Open in Browser

The app will automatically open at: **http://localhost:8501**

If it doesn't open automatically, copy the URL from the terminal.

---

## 📝 First Analysis (5 Minutes)

### Try with Example Data

1. **Load Sequences**
   - Go to "📝 Input Sequences" tab
   - Click "Upload File"
   - Select `example_sequences.fa`
   - Click "🚀 Process Uploaded File"
   - ✅ You should see "Processed 12 sequences"

2. **Run Motif Analysis**
   - Go to "🔍 Motif Analysis" tab
   - Click "🔬 Run Motif Analysis"
   - Wait ~10 seconds
   - ✅ View the enriched motifs table

3. **View Visualizations**
   - Scroll down to see:
     - Heatmap of motif occurrences
     - Sequence logos for top motifs

4. **Run Structure Prediction**
   - Go to "🧬 Structure Analysis" tab
   - Make sure "Run ViennaRNA Structure Prediction" is checked (sidebar)
   - Click "🧬 Run Structure Prediction"
   - ✅ View predicted structures with MFE values

5. **Export Results**
   - Go to "📊 Results & Export" tab
   - Click "📥 Download Motif Table (CSV)"
   - Click "📥 Download Structures (TXT)"

---

## 🔧 Your First Real Analysis

### 1. Configure Your Pool

**In the sidebar:**
- Expand "➕ Add/Edit Pool Configuration"
- Enter your pool details:
  - Configuration name (e.g., "Lab_N50_Pool")
  - Forward primer (5'→3')
  - Reverse complement sequence
  - Random region length
  - Description
- Click "💾 Save Configuration"

### 2. Prepare Your Sequences

**Format as FASTA:**
```
>Clone_1
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTACGTACGTAGATAGTAAGTGCAATCT
>Clone_2
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTGGCTAGCTAGATAGTAAGTGCAATCT
```

**Or plain text (one per line):**
```
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTACGTACGTAGATAGTAAGTGCAATCT
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTGGCTAGCTAGATAGTAAGTGCAATCT
```

### 3. Input Sequences

- Go to "📝 Input Sequences" tab
- Either paste sequences OR upload file
- Click "🚀 Process Sequences"

### 4. Adjust Parameters (Sidebar)

**Analysis Parameters:**
- Minimum Motif Length: 5 bp (default)
- Maximum Motif Length: 15 bp (default)
- Minimum Occurrences: 2 sequences (default)
- FDR Threshold: 0.05 (default)

**Structure Analysis:**
- Enable: ✅ Checked
- Temperature: 37°C (default)

### 5. Run Analysis

- Go to "🔍 Motif Analysis" tab → Click "Run Motif Analysis"
- Go to "🧬 Structure Analysis" tab → Click "Run Structure Prediction"

### 6. Interpret Results

**In the motif table, look for:**
- FDR < 0.05 (statistically significant)
- Fold_Enrichment > 2 (biologically meaningful)
- High frequency across sequences

**Check visualizations:**
- Heatmap shows which sequences share motifs
- Sequence logos show conservation patterns
- Structure predictions reveal folding

### 7. Export Everything

- Go to "📊 Results & Export"
- Download CSV table
- Download structures
- Save visualizations (right-click images)

---

## 📚 Documentation

### Quick Reference
- **README.md** - Basic info and features

### Complete Manual
- **USER_GUIDE.md** - 20 KB comprehensive guide with:
  - Detailed workflows
  - Result interpretation
  - Troubleshooting
  - FAQs
  - Best practices

### Technical Details
- **PROJECT_SUMMARY.md** - Complete technical overview

---

## ⚙️ Key Parameters Explained

### Motif Analysis

**Minimum Motif Length (5 bp)**
- Shorter = more hits, less specific
- Longer = fewer hits, more specific

**Minimum Occurrences (2)**
- Lower = find rare but enriched motifs
- Higher = focus on highly conserved motifs

**FDR Threshold (0.05)**
- 0.05 = 5% false discovery rate
- Lower (0.01) = more stringent
- Higher (0.10) = more permissive

### Structure Analysis

**Temperature (37°C)**
- Use your experimental selection temperature
- Affects folding thermodynamics

---

## 🎯 Understanding Results

### Motif Enrichment Table

| Column | Meaning |
|--------|---------|
| Motif | The sequence motif (e.g., "GGCTAG") |
| Count | Number of sequences containing it |
| Expected_Count | How many by random chance |
| Fold_Enrichment | Count / Expected (>2 is good) |
| FDR | Corrected p-value (use this, not p-value) |
| Significant | TRUE if FDR < threshold |

### What to Look For

✅ **Strong Candidates:**
- FDR < 0.05
- Fold_Enrichment > 2
- Present in many sequences

⚠️ **Weak Candidates:**
- FDR > 0.05
- Low fold enrichment
- Present in only 2-3 sequences

---

## 🐛 Troubleshooting

### "Forward primer not found"
→ Check that primer sequence in config matches your sequences

### "No significant motifs found"  
→ Try lowering FDR threshold or min_occurrences  
→ Check that you have enough sequences (need ≥10)

### ViennaRNA warnings
→ Tool will use fallback algorithm  
→ For full ViennaRNA: `conda install -c bioconda viennarna`

### Slow performance
→ Reduce motif length range (try 5-10 instead of 5-15)  
→ Reduce number of sequences in single analysis

---

## 💡 Pro Tips

1. **Start with defaults** - They work well for most cases
2. **Sample size matters** - Need ≥20 sequences for good power
3. **Look at both metrics** - FDR AND fold enrichment
4. **Validate experimentally** - Confirm top motifs with binding assays
5. **Compare rounds** - Analyze each SELEX round separately
6. **Save configurations** - Export pool settings for reproducibility

---

## 📊 What Makes a Good Motif?

**Statistical Significance:**
- FDR < 0.05 ✅
- FDR < 0.01 ✅✅ (even better)

**Biological Relevance:**
- Fold enrichment > 2 ✅
- Fold enrichment > 5 ✅✅ (strong)
- Fold enrichment > 10 ✅✅✅ (very strong)

**Prevalence:**
- Present in >50% of sequences ✅
- Present in >75% of sequences ✅✅

**Example of Excellent Motif:**
```
Motif: GGCTAGC
Count: 18/20 (90%)
Fold_Enrichment: 12.5
FDR: 1.2e-08
→ Strong binding site candidate!
```

---

## 🔬 Typical Workflows

### Workflow 1: Basic SELEX Analysis
1. Sequence final round clones (20-50)
2. Input sequences
3. Run motif analysis
4. Identify top 3-5 motifs (FDR < 0.05, Fold > 2)
5. Design validation experiments

### Workflow 2: Round Progression
1. Analyze each round separately
2. Export motif tables for each
3. Compare frequency and FDR across rounds
4. Track emergence of consensus sequences

### Workflow 3: Structure-Function
1. Run structure prediction for all clones
2. Group by similar structures
3. Correlate with binding data
4. Identify structure-binding relationships

---

## 🎓 Learning Resources

### Understand the Statistics
- **Binomial test**: Tests if motif frequency exceeds random expectation
- **FDR**: Controls false discovery rate (use this, not raw p-values!)
- **Fold enrichment**: Practical measure of enrichment strength

### Interpret Visualizations
- **Heatmap**: Which sequences share which motifs
- **Sequence logos**: Base preferences at each position
- **Volcano plot**: See enrichment vs. significance tradeoff

### Structure Analysis
- **Dot-bracket notation**: ( ) = paired, . = unpaired
- **MFE**: Minimum free energy (more negative = more stable)
- **Typical range**: -5 to -30 kcal/mol

---

## ✉️ Need More Help?

1. **Read USER_GUIDE.md** - Comprehensive 20 KB manual
2. **Check troubleshooting section** - Common issues and solutions
3. **Test with examples** - Use example_sequences.fa
4. **Review tooltips** - Hover over parameters in the app

---

## 📦 Files Included

```
✅ aptamer_motif_analyzer.py    - Main web app
✅ motif_analysis.py            - Motif discovery
✅ statistics_module.py         - Statistical tests
✅ structure_analysis.py        - RNA folding
✅ visualizations.py            - Plotting
✅ README.md                    - Quick reference
✅ USER_GUIDE.md               - Full manual
✅ PROJECT_SUMMARY.md          - Technical specs
✅ requirements.txt            - Dependencies
✅ launch.sh                   - Easy startup
✅ test_modules.py             - Verification
✅ example_sequences.fa        - Test data
```

---

## 🏁 Ready to Start?

```bash
# 1. Install
pip install -r requirements.txt --break-system-packages

# 2. Test
python3 test_modules.py

# 3. Launch
./launch.sh

# 4. Open browser to http://localhost:8501

# 5. Load example_sequences.fa and explore!
```

---

**Version 1.0 | February 2026 | Ready for Production ✅**

*Happy Analyzing! 🧬*
