#!/usr/bin/env python3
"""
Test script to verify AptaMotif installation
Run this before launching the web app to check dependencies
"""

import sys

def test_imports():
    """Test if all required packages can be imported"""
    print("Testing package imports...")
    
    packages = [
        ('streamlit', 'Streamlit'),
        ('Bio', 'Biopython'),
        ('pandas', 'Pandas'),
        ('numpy', 'NumPy'),
        ('scipy', 'SciPy'),
        ('matplotlib', 'Matplotlib'),
        ('seaborn', 'Seaborn'),
        ('logomaker', 'Logomaker'),
        ('plotly', 'Plotly'),
        ('RNA', 'ViennaRNA')
    ]
    
    failed = []
    
    for package, name in packages:
        try:
            __import__(package)
            print(f"  ✓ {name}")
        except ImportError as e:
            print(f"  ✗ {name} - FAILED")
            failed.append((name, str(e)))
    
    return failed

def test_modules():
    """Test if AptaMotif modules can be imported"""
    print("\nTesting AptaMotif modules...")
    
    modules = [
        'sequence_parser',
        'motif_finder',
        'statistics',
        'structure_analyzer',
        'visualizer'
    ]
    
    failed = []
    
    for module in modules:
        try:
            __import__(module)
            print(f"  ✓ {module}.py")
        except ImportError as e:
            print(f"  ✗ {module}.py - FAILED")
            failed.append((module, str(e)))
    
    return failed

def test_viennarna():
    """Test ViennaRNA functionality"""
    print("\nTesting ViennaRNA functionality...")
    
    try:
        import RNA
        
        # Test basic folding
        test_seq = "GGGAAACCCUUUGGG"
        fc = RNA.fold_compound(test_seq)
        structure, mfe = fc.mfe()
        
        print(f"  ✓ ViennaRNA is working")
        print(f"    Test sequence: {test_seq}")
        print(f"    Structure: {structure}")
        print(f"    MFE: {mfe:.2f} kcal/mol")
        
        return []
        
    except Exception as e:
        print(f"  ✗ ViennaRNA test failed: {str(e)}")
        return [('ViennaRNA', str(e))]

def test_example_file():
    """Check if example file exists"""
    print("\nChecking example data...")
    
    import os
    
    if os.path.exists('example_sequences.fasta'):
        print("  ✓ example_sequences.fasta found")
        
        # Count sequences
        with open('example_sequences.fasta', 'r') as f:
            content = f.read()
            seq_count = content.count('>')
        
        print(f"    Contains {seq_count} sequences")
        return []
    else:
        print("  ✗ example_sequences.fasta not found")
        return [('Example file', 'File not found')]

def run_quick_test():
    """Run a quick analysis test"""
    print("\nRunning quick analysis test...")
    
    try:
        from sequence_parser import SequenceParser
        from motif_finder import MotifFinder
        
        # Test data
        test_sequences = """
>Seq1
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTGGATCCGGATCCAGATAGTAAGTGCAATCT
>Seq2
TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATTGGATCCGGATCCAGATAGTAAGTGCAATCT
""".strip()
        
        # Parse
        parser = SequenceParser(
            forward_primer="TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATT",
            reverse_primer_region="AGATAGTAAGTGCAATCT",
            molecule_type="RNA"
        )
        
        sequences = parser.parse_sequences(test_sequences)
        sequences = parser.extract_random_region(sequences)
        
        print(f"  ✓ Parsed {len(sequences)} sequences")
        
        # Find motifs
        finder = MotifFinder(min_length=5, max_length=10, min_sequences=2)
        motifs = finder.find_enriched_motifs(sequences)
        
        print(f"  ✓ Found {len(motifs)} motifs")
        
        if motifs:
            top_motif = motifs[0]
            print(f"    Top motif: {top_motif['motif']} "
                  f"(appears in {top_motif['num_sequences']} sequences)")
        
        return []
        
    except Exception as e:
        print(f"  ✗ Quick test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return [('Quick test', str(e))]

def main():
    """Run all tests"""
    print("="*60)
    print("AptaMotif Installation Test")
    print("="*60)
    
    all_failures = []
    
    # Test imports
    failures = test_imports()
    all_failures.extend(failures)
    
    # Test modules
    failures = test_modules()
    all_failures.extend(failures)
    
    # Test ViennaRNA
    failures = test_viennarna()
    all_failures.extend(failures)
    
    # Check example file
    failures = test_example_file()
    all_failures.extend(failures)
    
    # Run quick test
    failures = run_quick_test()
    all_failures.extend(failures)
    
    # Summary
    print("\n" + "="*60)
    if not all_failures:
        print("✓ ALL TESTS PASSED!")
        print("\nYou can now run the application with:")
        print("  streamlit run app.py")
    else:
        print("✗ SOME TESTS FAILED")
        print("\nFailed components:")
        for name, error in all_failures:
            print(f"  - {name}: {error}")
        print("\nPlease fix these issues before running the application.")
        print("See README.md for installation instructions.")
    print("="*60)
    
    return 0 if not all_failures else 1

if __name__ == "__main__":
    sys.exit(main())
