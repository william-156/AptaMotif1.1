"""
Structure Analysis Module
Predicts secondary structures using ViennaRNA Python bindings
"""

import subprocess
import tempfile
import os
from collections import defaultdict

# Try to import ViennaRNA
try:
    import RNA
    VIENNA_AVAILABLE = True
except ImportError:
    VIENNA_AVAILABLE = False
    print("Warning: ViennaRNA Python package not available.")
    print("Structure prediction will use a simple algorithm.")


class StructureAnalyzer:
    """
    Predicts and analyzes RNA secondary structures using ViennaRNA.
    Falls back to simple algorithm if ViennaRNA is not available.
    """
    
    def __init__(self, temperature=37):
        """
        Initialize the structure analyzer.
        
        Parameters:
        -----------
        temperature : float
            Temperature in Celsius for folding (default: 37)
        """
        self.temperature = temperature
        self.vienna_available = VIENNA_AVAILABLE
    
    def predict_structure(self, rna_sequence, seq_id="sequence"):
        """
        Predict secondary structure for a single RNA sequence.
        
        Parameters:
        -----------
        rna_sequence : str
            RNA sequence (with U, not T)
        seq_id : str
            Sequence identifier
            
        Returns:
        --------
        dict : Dictionary with structure, MFE, and other information
        """
        if self.vienna_available:
            return self._predict_with_vienna(rna_sequence, seq_id)
        else:
            return self._predict_simple(rna_sequence, seq_id)
    
    def _predict_with_vienna(self, rna_sequence, seq_id):
        """Use ViennaRNA for structure prediction."""
        try:
            # Set temperature
            RNA.cvar.temperature = self.temperature
            
            # Fold the sequence
            fc = RNA.fold_compound(rna_sequence)
            structure, mfe = fc.mfe()
            
            return {
                'sequence': rna_sequence,
                'structure': structure,
                'mfe': mfe,
                'seq_id': seq_id
            }
        except Exception as e:
            print(f"Error with ViennaRNA prediction: {str(e)}")
            return self._predict_simple(rna_sequence, seq_id)
    
    def _predict_simple(self, rna_sequence, seq_id):
        """
        Simple structure prediction using Nussinov-like algorithm.
        This is a fallback when ViennaRNA is not available.
        """
        n = len(rna_sequence)
        
        # Initialize DP table
        dp = [[0] * n for _ in range(n)]
        trace = [[None] * n for _ in range(n)]
        
        # Base pairing rules
        pairs = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'G': 'U', 'U': 'G'}
        
        # Fill DP table
        for length in range(2, n + 1):
            for i in range(n - length + 1):
                j = i + length - 1
                
                # Case 1: j unpaired
                dp[i][j] = dp[i][j-1]
                trace[i][j] = (i, j-1, 'unpaired')
                
                # Case 2: j pairs with some k
                for k in range(i, j):
                    if (rna_sequence[k] in pairs and 
                        pairs[rna_sequence[k]] == rna_sequence[j]):
                        score = 1
                        if k > i:
                            score += dp[i][k-1]
                        if k + 1 < j:
                            score += dp[k+1][j-1]
                        
                        if score > dp[i][j]:
                            dp[i][j] = score
                            trace[i][j] = (k, j, 'paired')
        
        # Traceback to get structure
        structure = ['.'] * n
        self._traceback_simple(trace, structure, 0, n-1, rna_sequence, pairs)
        
        # Estimate MFE (very rough approximation)
        mfe = -dp[0][n-1] * 2.5  # Rough estimate
        
        return {
            'sequence': rna_sequence,
            'structure': ''.join(structure),
            'mfe': mfe,
            'seq_id': seq_id
        }
    
    def _traceback_simple(self, trace, structure, i, j, seq, pairs):
        """Traceback for simple structure prediction."""
        if i >= j or trace[i][j] is None:
            return
        
        k, j_end, pair_type = trace[i][j]
        
        if pair_type == 'unpaired':
            self._traceback_simple(trace, structure, i, j-1, seq, pairs)
        elif pair_type == 'paired':
            structure[k] = '('
            structure[j] = ')'
            if k > i:
                self._traceback_simple(trace, structure, i, k-1, seq, pairs)
            if k + 1 < j:
                self._traceback_simple(trace, structure, k+1, j-1, seq, pairs)
    
    def predict_structures(self, rna_sequences):
        """
        Predict structures for multiple RNA sequences.
        
        Parameters:
        -----------
        rna_sequences : dict
            Dictionary of {sequence_id: rna_sequence}
            
        Returns:
        --------
        dict : Dictionary of {sequence_id: structure_info}
        """
        structures = {}
        
        for seq_id, rna_seq in rna_sequences.items():
            try:
                structure_info = self.predict_structure(rna_seq, seq_id)
                structures[seq_id] = structure_info
            except Exception as e:
                print(f"Error predicting structure for {seq_id}: {str(e)}")
                structures[seq_id] = {
                    'sequence': rna_seq,
                    'structure': None,
                    'mfe': None,
                    'seq_id': seq_id,
                    'error': str(e)
                }
        
        return structures
    
    def find_structural_motifs(self, structures, min_occurrences=2):
        """
        Find common structural motifs across sequences.
        
        Parameters:
        -----------
        structures : dict
            Dictionary of structure information from predict_structures
        min_occurrences : int
            Minimum number of sequences that must share a motif
            
        Returns:
        --------
        dict : Dictionary of structural motifs and their occurrences
        """
        # Extract structural elements
        structural_elements = defaultdict(list)
        
        for seq_id, struct_info in structures.items():
            if struct_info['structure'] is None:
                continue
            
            structure = struct_info['structure']
            
            # Find stems (consecutive base pairs)
            stems = self._extract_stems(structure)
            for stem in stems:
                structural_elements[f"stem_{stem}"].append(seq_id)
            
            # Find loops
            loops = self._extract_loops(structure)
            for loop in loops:
                structural_elements[f"loop_{loop}"].append(seq_id)
            
            # Find bulges
            bulges = self._extract_bulges(structure)
            for bulge in bulges:
                structural_elements[f"bulge_{bulge}"].append(seq_id)
        
        # Filter by minimum occurrences
        common_motifs = {
            motif: seq_list 
            for motif, seq_list in structural_elements.items()
            if len(seq_list) >= min_occurrences
        }
        
        return common_motifs
    
    def _extract_stems(self, structure, min_length=3):
        """Extract stem regions (consecutive base pairs)."""
        stems = []
        current_stem = 0
        
        for char in structure:
            if char == '(':
                current_stem += 1
            elif char == ')':
                if current_stem >= min_length:
                    stems.append(current_stem)
                current_stem = 0
            else:
                if current_stem >= min_length:
                    stems.append(current_stem)
                current_stem = 0
        
        if current_stem >= min_length:
            stems.append(current_stem)
        
        return stems
    
    def _extract_loops(self, structure):
        """Extract loop regions (unpaired bases within structure)."""
        loops = []
        in_loop = False
        loop_size = 0
        depth = 0
        
        for char in structure:
            if char == '(':
                depth += 1
                if in_loop and loop_size > 0:
                    loops.append(loop_size)
                in_loop = False
                loop_size = 0
            elif char == ')':
                depth -= 1
                if in_loop and loop_size > 0:
                    loops.append(loop_size)
                in_loop = False
                loop_size = 0
            elif char == '.':
                if depth > 0:  # Internal loop
                    in_loop = True
                    loop_size += 1
        
        return loops
    
    def _extract_bulges(self, structure):
        """Extract bulge regions (unpaired bases on one side of stem)."""
        # Simplified bulge detection
        bulges = []
        unpaired = 0
        
        for i, char in enumerate(structure):
            if char == '.':
                unpaired += 1
            else:
                if unpaired > 0 and unpaired < 5:  # Small unpaired regions are bulges
                    bulges.append(unpaired)
                unpaired = 0
        
        return bulges
    
    def calculate_structure_similarity(self, structure1, structure2):
        """
        Calculate similarity between two structures using base pair distance.
        
        Parameters:
        -----------
        structure1 : str
            First structure in dot-bracket notation
        structure2 : str
            Second structure in dot-bracket notation
            
        Returns:
        --------
        float : Similarity score (0-1)
        """
        if len(structure1) != len(structure2):
            # Align to shorter length
            min_len = min(len(structure1), len(structure2))
            structure1 = structure1[:min_len]
            structure2 = structure2[:min_len]
        
        # Count matching positions
        matches = sum(1 for a, b in zip(structure1, structure2) if a == b)
        similarity = matches / len(structure1) if len(structure1) > 0 else 0
        
        return similarity
    
    def get_consensus_structure(self, structures):
        """
        Get consensus structure from multiple sequences.
        
        Parameters:
        -----------
        structures : dict
            Dictionary of structure information
            
        Returns:
        --------
        str : Consensus structure in dot-bracket notation
        """
        # Get all structures
        struct_list = [s['structure'] for s in structures.values() 
                      if s['structure'] is not None]
        
        if not struct_list:
            return None
        
        # Find maximum length
        max_len = max(len(s) for s in struct_list)
        
        # Pad structures to same length
        padded = [s + '.' * (max_len - len(s)) for s in struct_list]
        
        # Consensus by majority vote at each position
        consensus = []
        for pos in range(max_len):
            chars = [s[pos] for s in padded]
            # Most common character at this position
            most_common = max(set(chars), key=chars.count)
            consensus.append(most_common)
        
        return ''.join(consensus)
    
    def visualize_structure(self, structure_info, output_file=None):
        """
        Create a visualization of the RNA structure.
        
        Parameters:
        -----------
        structure_info : dict
            Structure information from predict_structure
        output_file : str
            Output file path for the visualization (optional)
            
        Returns:
        --------
        str : Path to the visualization file
        """
        seq_id = structure_info['seq_id']
        sequence = structure_info['sequence']
        structure = structure_info['structure']
        
        if structure is None:
            return None
        
        # Create temporary fasta file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fa') as f:
            f.write(f">{seq_id}\n{sequence}\n{structure}\n")
            temp_input = f.name
        
        # Determine output file
        if output_file is None:
            output_file = f"{seq_id}_structure.svg"
        
        try:
            # Run RNAplot to create visualization
            subprocess.run(
                ['RNAplot', '--output-format=svg', f'--outfile={output_file}'],
                stdin=open(temp_input),
                capture_output=True
            )
            
            if os.path.exists(output_file):
                return output_file
            else:
                return None
                
        except Exception as e:
            print(f"Error visualizing structure: {str(e)}")
            return None
        finally:
            if os.path.exists(temp_input):
                os.remove(temp_input)
