"""
Motif Analysis Module
Finds shared k-mers across sequences
"""

from collections import defaultdict, Counter
import pandas as pd
from itertools import combinations


class MotifAnalyzer:
    """
    Identifies motifs (k-mers) shared between multiple sequences.
    """
    
    def __init__(self, sequences, min_length=5, max_length=15, min_occurrences=2):
        """
        Initialize the motif analyzer.
        
        Parameters:
        -----------
        sequences : dict
            Dictionary of {sequence_id: sequence}
        min_length : int
            Minimum motif length (default: 5)
        max_length : int
            Maximum motif length (default: 15)
        min_occurrences : int
            Minimum number of sequences that must share a motif (default: 2)
        """
        self.sequences = sequences
        self.min_length = min_length
        self.max_length = max_length
        self.min_occurrences = min_occurrences
        self.motifs = None
    
    def extract_kmers(self, sequence, k):
        """
        Extract all k-mers from a sequence.
        
        Parameters:
        -----------
        sequence : str
            DNA sequence
        k : int
            k-mer length
            
        Returns:
        --------
        set : Set of k-mers found in the sequence
        """
        kmers = set()
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            # Only include k-mers with valid DNA bases
            if all(base in 'ATCG' for base in kmer):
                kmers.add(kmer)
        return kmers
    
    def find_motifs(self):
        """
        Find all motifs shared between sequences.
        
        Returns:
        --------
        pandas.DataFrame : DataFrame with motif information
        """
        # Dictionary to store which sequences contain each motif
        motif_occurrences = defaultdict(set)
        
        # For each k-mer length
        for k in range(self.min_length, self.max_length + 1):
            # Extract k-mers from each sequence
            for seq_id, sequence in self.sequences.items():
                kmers = self.extract_kmers(sequence, k)
                
                # Record which sequences contain each k-mer
                for kmer in kmers:
                    motif_occurrences[kmer].add(seq_id)
        
        # Filter motifs by minimum occurrences
        filtered_motifs = {
            motif: seq_ids 
            for motif, seq_ids in motif_occurrences.items()
            if len(seq_ids) >= self.min_occurrences
        }
        
        # Remove redundant motifs (if a motif is contained in a longer motif with same occurrence)
        filtered_motifs = self._remove_redundant_motifs(filtered_motifs)
        
        # Create DataFrame
        motif_data = []
        for motif, seq_ids in filtered_motifs.items():
            motif_data.append({
                'Motif': motif,
                'Length': len(motif),
                'Count': len(seq_ids),
                'Frequency': len(seq_ids) / len(self.sequences),
                'Sequences': ','.join(sorted(seq_ids))
            })
        
        self.motifs = pd.DataFrame(motif_data)
        
        # Sort by count (descending) and then by length (descending)
        if len(self.motifs) > 0:
            self.motifs = self.motifs.sort_values(
                by=['Count', 'Length'], 
                ascending=[False, False]
            ).reset_index(drop=True)
        
        return self.motifs
    
    def _remove_redundant_motifs(self, motif_dict):
        """
        Remove motifs that are substrings of longer motifs with identical occurrence.
        
        This helps reduce redundancy while keeping the most informative motifs.
        """
        motifs_by_length = defaultdict(list)
        for motif in motif_dict.keys():
            motifs_by_length[len(motif)].append(motif)
        
        non_redundant = {}
        
        # Process from longest to shortest
        for length in sorted(motifs_by_length.keys(), reverse=True):
            for motif in motifs_by_length[length]:
                seq_set = motif_dict[motif]
                
                # Check if this motif is redundant with any already added longer motif
                is_redundant = False
                for existing_motif, existing_set in non_redundant.items():
                    # If this motif is a substring of an existing one
                    # and they occur in the exact same sequences, it's redundant
                    if motif in existing_motif and seq_set == existing_set:
                        is_redundant = True
                        break
                
                if not is_redundant:
                    non_redundant[motif] = seq_set
        
        return non_redundant
    
    def get_motif_positions(self, motif):
        """
        Get the positions of a motif in each sequence.
        
        Parameters:
        -----------
        motif : str
            The motif sequence to search for
            
        Returns:
        --------
        dict : Dictionary mapping sequence_id to list of positions
        """
        positions = defaultdict(list)
        
        for seq_id, sequence in self.sequences.items():
            # Find all occurrences
            start = 0
            while True:
                pos = sequence.find(motif, start)
                if pos == -1:
                    break
                positions[seq_id].append(pos)
                start = pos + 1
        
        return positions
    
    def create_motif_matrix(self):
        """
        Create a binary matrix of motif presence/absence.
        
        Returns:
        --------
        pandas.DataFrame : Binary matrix (sequences x motifs)
        """
        if self.motifs is None:
            raise ValueError("Must run find_motifs() first")
        
        # Create matrix
        seq_ids = sorted(self.sequences.keys())
        motif_list = self.motifs['Motif'].tolist()
        
        matrix = pd.DataFrame(0, index=seq_ids, columns=motif_list)
        
        for motif in motif_list:
            for seq_id, sequence in self.sequences.items():
                if motif in sequence:
                    matrix.loc[seq_id, motif] = 1
        
        return matrix
    
    def get_consensus_sequence(self, sequences_subset):
        """
        Generate a consensus sequence from a subset of sequences.
        
        Parameters:
        -----------
        sequences_subset : list
            List of sequences to generate consensus from
            
        Returns:
        --------
        str : Consensus sequence
        """
        if not sequences_subset:
            return ""
        
        # Align sequences (simple approach - assumes similar length)
        max_len = max(len(seq) for seq in sequences_subset)
        
        consensus = []
        for pos in range(max_len):
            bases = []
            for seq in sequences_subset:
                if pos < len(seq):
                    bases.append(seq[pos])
            
            if bases:
                # Most common base at this position
                counter = Counter(bases)
                consensus.append(counter.most_common(1)[0][0])
        
        return ''.join(consensus)
