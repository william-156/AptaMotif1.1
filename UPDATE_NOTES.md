# AptaMotif v1.1 - Update Notes

## Major Improvements Based on Feedback

### 1. ✅ Corrected Forward Primer
**Change**: Updated forward primer from full T7 promoter sequence to just the binding region
- **Old**: `TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATT`
- **New**: `GGGAGATACCAGCTTATTCAATT`

### 2. ✅ Full Sequence Folding
**Change**: Structure prediction now uses the FULL ~112 nt sequence, not just the N71 random region

### 3. ✅ SVG Structure Diagrams
**New Feature**: Generate graphical structure diagrams using ViennaRNA's RNAplot

### 4. ✅ No Truncation - Full Results Display
All motifs, structures, and logos now accessible (not limited to top 100 or 10)

### 5. ✅ Sequence Identity Detection
Automatically detects duplicate and near-duplicate sequences

### 6. ✅ Enhanced Redundancy Filtering
Identifies and merges overlapping short motifs (e.g., CCTAT + TATGG → CCTATGG)

### 7. ✅ Motif Table - Sequence IDs Column
Shows which specific sequences contain each motif

### 8. ✅ Pairwise Similarity Heatmap
New visualization showing shared significant motifs between sequence pairs

### 9. ✅ Enrichment Calculation Documentation
Complete formula with worked example in Help tab

## All Requested Features Implemented! 🎉
