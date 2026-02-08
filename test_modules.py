#!/usr/bin/env python3
"""
Test script for AptaMotif Analyzer modules
"""

print("Testing AptaMotif Analyzer modules...")
print("=" * 50)

# Test imports
print("\n1. Testing imports...")
try:
    from motif_analysis import MotifAnalyzer
    from statistics_module import StatisticalAnalyzer
    from visualizations import MotifVisualizer
    from structure_analysis import StructureAnalyzer
    print("✅ All modules imported successfully")
except ImportError as e:
    print(f"❌ Import error: {e}")
    exit(1)

# Test with example sequences
print("\n2. Testing motif discovery...")
test_sequences = {
    'Seq1': 'GGCTAGCTAGCTAGCGGGCCCAAAT',
    'Seq2': 'ACGTGGCTAGCTAGCTAGCAAAGGGCCC',
    'Seq3': 'TGCATGGCTAGCTAGCTAGCGGGCCCATC'
}

try:
    analyzer = MotifAnalyzer(
        test_sequences,
        min_length=5,
        max_length=10,
        min_occurrences=2
    )
    motifs = analyzer.find_motifs()
    print(f"✅ Found {len(motifs)} motifs")
    print(f"   Top motif: {motifs.iloc[0]['Motif'] if len(motifs) > 0 else 'None'}")
except Exception as e:
    print(f"❌ Motif analysis error: {e}")
    exit(1)

# Test statistical analysis
print("\n3. Testing statistical analysis...")
try:
    stats_analyzer = StatisticalAnalyzer(
        motifs,
        test_sequences,
        random_region_length=25
    )
    results = stats_analyzer.calculate_enrichment(fdr_threshold=0.05)
    print(f"✅ Statistical analysis complete")
    if len(results) > 0:
        print(f"   Most significant motif FDR: {results.iloc[0]['FDR']:.2e}")
except Exception as e:
    print(f"❌ Statistical analysis error: {e}")
    exit(1)

# Test structure prediction
print("\n4. Testing structure prediction...")
try:
    structure_analyzer = StructureAnalyzer(temperature=37)
    rna_seqs = {'Test': 'GGCUAGCUAGCUAGCGGGCCCAAAU'}
    structures = structure_analyzer.predict_structures(rna_seqs)
    print(f"✅ Structure prediction complete")
    if 'Test' in structures and structures['Test']['structure']:
        print(f"   Structure: {structures['Test']['structure'][:30]}...")
        print(f"   MFE: {structures['Test']['mfe']:.2f} kcal/mol")
    if not structure_analyzer.vienna_available:
        print("   ⚠️  Using fallback algorithm (ViennaRNA not available)")
except Exception as e:
    print(f"❌ Structure prediction error: {e}")
    exit(1)

# Test visualization
print("\n5. Testing visualization...")
try:
    visualizer = MotifVisualizer()
    # Just check that visualizer can be created
    print("✅ Visualizer initialized successfully")
except Exception as e:
    print(f"❌ Visualization error: {e}")
    exit(1)

print("\n" + "=" * 50)
print("✅ All tests passed!")
print("\nReady to launch AptaMotif Analyzer!")
print("Run: streamlit run aptamer_motif_analyzer.py")
print("Or use: ./launch.sh")
