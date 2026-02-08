"""
Statistical Analysis Module
Calculates p-values and FDR for motif enrichment
"""

import numpy as np
import pandas as pd
from scipy.stats import binom
from statsmodels.stats.multitest import multipletests


class StatisticalAnalyzer:
    """
    Performs statistical analysis on motif enrichment.
    """
    
    def __init__(self, motif_df, sequences, random_region_length):
        """
        Initialize the statistical analyzer.
        
        Parameters:
        -----------
        motif_df : pandas.DataFrame
            DataFrame with motif information from MotifAnalyzer
        sequences : dict
            Dictionary of {sequence_id: sequence}
        random_region_length : int
            Expected length of random region
        """
        self.motif_df = motif_df.copy() if motif_df is not None else pd.DataFrame()
        self.sequences = sequences
        self.n_sequences = len(sequences)
        self.random_region_length = random_region_length
        
        # Base probabilities (assume equal)
        self.base_prob = 0.25
    
    def calculate_motif_probability(self, motif_length):
        """
        Calculate the probability of a specific motif occurring by chance.
        
        For a motif of length k, probability = (0.25)^k
        
        Parameters:
        -----------
        motif_length : int
            Length of the motif
            
        Returns:
        --------
        float : Probability of the motif occurring at any position
        """
        return self.base_prob ** motif_length
    
    def calculate_positions_per_sequence(self, motif_length):
        """
        Calculate number of possible positions for a motif in a sequence.
        
        Parameters:
        -----------
        motif_length : int
            Length of the motif
            
        Returns:
        --------
        int : Number of possible starting positions
        """
        return max(1, self.random_region_length - motif_length + 1)
    
    def binomial_test(self, motif_length, observed_count):
        """
        Perform binomial test for motif enrichment.
        
        Tests the null hypothesis that the motif occurs at the expected 
        frequency by chance across all sequences.
        
        Parameters:
        -----------
        motif_length : int
            Length of the motif
        observed_count : int
            Number of sequences containing the motif
            
        Returns:
        --------
        float : p-value
        """
        # Probability that a specific motif appears at any position
        p_motif = self.calculate_motif_probability(motif_length)
        
        # Number of positions to check per sequence
        n_positions = self.calculate_positions_per_sequence(motif_length)
        
        # Probability that the motif appears at least once in a sequence
        # P(at least once) = 1 - P(never) = 1 - (1-p)^n
        p_in_sequence = 1 - (1 - p_motif) ** n_positions
        
        # Binomial test: out of n_sequences, how many contain the motif?
        # P(X >= observed_count) where X ~ Binomial(n_sequences, p_in_sequence)
        p_value = binom.sf(observed_count - 1, self.n_sequences, p_in_sequence)
        
        return p_value
    
    def calculate_enrichment(self, fdr_threshold=0.05):
        """
        Calculate enrichment statistics for all motifs.
        
        Adds p-values and FDR-corrected values to the motif DataFrame.
        
        Parameters:
        -----------
        fdr_threshold : float
            FDR threshold for significance (default: 0.05)
            
        Returns:
        --------
        pandas.DataFrame : Updated DataFrame with statistical values
        """
        if len(self.motif_df) == 0:
            return self.motif_df
        
        # Calculate p-values for each motif
        p_values = []
        expected_counts = []
        fold_enrichments = []
        
        for idx, row in self.motif_df.iterrows():
            motif_length = row['Length']
            observed_count = row['Count']
            
            # P-value from binomial test
            p_val = self.binomial_test(motif_length, observed_count)
            p_values.append(p_val)
            
            # Calculate expected count
            p_motif = self.calculate_motif_probability(motif_length)
            n_positions = self.calculate_positions_per_sequence(motif_length)
            p_in_sequence = 1 - (1 - p_motif) ** n_positions
            expected = p_in_sequence * self.n_sequences
            expected_counts.append(expected)
            
            # Fold enrichment
            if expected > 0:
                fold_enrichment = observed_count / expected
            else:
                fold_enrichment = np.inf
            fold_enrichments.append(fold_enrichment)
        
        # Add to DataFrame
        self.motif_df['Expected_Count'] = expected_counts
        self.motif_df['Fold_Enrichment'] = fold_enrichments
        self.motif_df['P_value'] = p_values
        
        # FDR correction (Benjamini-Hochberg)
        if len(p_values) > 0:
            reject, fdr_values, _, _ = multipletests(
                p_values, 
                alpha=fdr_threshold, 
                method='fdr_bh'
            )
            self.motif_df['FDR'] = fdr_values
            self.motif_df['Significant'] = reject
        else:
            self.motif_df['FDR'] = []
            self.motif_df['Significant'] = []
        
        # Sort by FDR
        self.motif_df = self.motif_df.sort_values('FDR').reset_index(drop=True)
        
        # Reorder columns for better readability
        column_order = [
            'Motif', 'Length', 'Count', 'Expected_Count', 
            'Fold_Enrichment', 'Frequency', 'P_value', 'FDR', 
            'Significant', 'Sequences'
        ]
        self.motif_df = self.motif_df[column_order]
        
        return self.motif_df
    
    def permutation_test(self, motif, n_permutations=1000):
        """
        Perform permutation test for a specific motif.
        
        This is more computationally intensive but makes fewer assumptions.
        
        Parameters:
        -----------
        motif : str
            The motif sequence to test
        n_permutations : int
            Number of permutations (default: 1000)
            
        Returns:
        --------
        float : Empirical p-value
        """
        # Count observed occurrences
        observed = sum(1 for seq in self.sequences.values() if motif in seq)
        
        # Generate null distribution
        null_counts = []
        for _ in range(n_permutations):
            # Shuffle each sequence
            permuted_count = 0
            for seq in self.sequences.values():
                # Shuffle sequence
                seq_list = list(seq)
                np.random.shuffle(seq_list)
                permuted_seq = ''.join(seq_list)
                
                # Check for motif
                if motif in permuted_seq:
                    permuted_count += 1
            
            null_counts.append(permuted_count)
        
        # Calculate p-value
        p_value = (np.sum(np.array(null_counts) >= observed) + 1) / (n_permutations + 1)
        
        return p_value
    
    def calculate_gc_content(self, sequence):
        """
        Calculate GC content of a sequence.
        
        Parameters:
        -----------
        sequence : str
            DNA sequence
            
        Returns:
        --------
        float : GC content (0-1)
        """
        if len(sequence) == 0:
            return 0.0
        
        gc_count = sequence.count('G') + sequence.count('C')
        return gc_count / len(sequence)
    
    def adjust_for_gc_bias(self):
        """
        Adjust expected counts based on actual GC content of sequences.
        
        This provides a more accurate null model if sequences have 
        non-uniform base composition.
        """
        # Calculate average GC content
        gc_contents = [self.calculate_gc_content(seq) 
                      for seq in self.sequences.values()]
        avg_gc = np.mean(gc_contents)
        
        # Adjust probabilities
        # If GC content is higher, G and C are more likely
        p_gc = avg_gc / 2  # probability of G or C
        p_at = (1 - avg_gc) / 2  # probability of A or T
        
        # Recalculate p-values with adjusted probabilities
        # This is a simplified adjustment - more sophisticated methods exist
        adjusted_p_values = []
        
        for idx, row in self.motif_df.iterrows():
            motif = row['Motif']
            observed_count = row['Count']
            
            # Calculate probability of this specific motif
            p_motif = 1.0
            for base in motif:
                if base in 'GC':
                    p_motif *= p_gc
                else:
                    p_motif *= p_at
            
            # Rest of calculation is same as binomial_test
            n_positions = self.calculate_positions_per_sequence(len(motif))
            p_in_sequence = 1 - (1 - p_motif) ** n_positions
            p_value = binom.sf(observed_count - 1, self.n_sequences, p_in_sequence)
            
            adjusted_p_values.append(p_value)
        
        self.motif_df['P_value_GC_adjusted'] = adjusted_p_values
        
        # Recalculate FDR
        if len(adjusted_p_values) > 0:
            reject, fdr_values, _, _ = multipletests(
                adjusted_p_values, 
                alpha=0.05, 
                method='fdr_bh'
            )
            self.motif_df['FDR_GC_adjusted'] = fdr_values
        
        return self.motif_df
