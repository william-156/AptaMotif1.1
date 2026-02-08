"""
Visualization Module
Creates heatmaps, sequence logos, and other visualizations
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as npx
import pandas as pd
import logomaker


class MotifVisualizer:
    """
    Creates visualizations for motif analysis results.
    """
    
    def __init__(self):
        """Initialize the visualizer with default settings."""
        # Set style
        sns.set_style("whitegrid")
        plt.rcParams['figure.dpi'] = 100
    
    def create_heatmap(self, motif_df, sequences, top_n=20):
        """
        Create a heatmap showing motif presence/absence across sequences.
        
        Parameters:
        -----------
        motif_df : pandas.DataFrame
            DataFrame with motif information
        sequences : dict
            Dictionary of {sequence_id: sequence}
        top_n : int
            Number of top motifs to display (default: 20)
            
        Returns:
        --------
        matplotlib.figure.Figure : Heatmap figure
        """
        if len(motif_df) == 0:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, 'No motifs found', 
                   ha='center', va='center', fontsize=16)
            ax.axis('off')
            return fig
        
        # Select top motifs
        top_motifs = motif_df.nsmallest(min(top_n, len(motif_df)), 'FDR')
        motif_list = top_motifs['Motif'].tolist()
        
        # Create binary matrix
        seq_ids = sorted(sequences.keys())
        matrix = pd.DataFrame(0, index=seq_ids, columns=motif_list)
        
        for motif in motif_list:
            for seq_id, sequence in sequences.items():
                if motif in sequence:
                    matrix.loc[seq_id, motif] = 1
        
        # Create figure
        fig, ax = plt.subplots(figsize=(max(12, len(motif_list) * 0.5), 
                                        max(8, len(seq_ids) * 0.3)))
        
        # Create heatmap
        sns.heatmap(
            matrix,
            cmap=['#f0f0f0', '#2E86AB'],  # Light gray for 0, blue for 1
            cbar_kws={'label': 'Motif Present', 'ticks': [0, 1]},
            linewidths=0.5,
            linecolor='white',
            square=False,
            ax=ax
        )
        
        # Customize
        ax.set_xlabel('Motif', fontsize=12, fontweight='bold')
        ax.set_ylabel('Sequence ID', fontsize=12, fontweight='bold')
        ax.set_title('Motif Occurrence Heatmap', fontsize=14, fontweight='bold', pad=20)
        
        # Rotate x labels
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
        plt.setp(ax.get_yticklabels(), rotation=0, fontsize=9)
        
        plt.tight_layout()
        
        return fig
    
    def create_sequence_logo(self, motif, sequences, flanking=2):
        """
        Create a sequence logo for a motif and its flanking regions.
        
        Parameters:
        -----------
        motif : str
            The motif sequence
        sequences : dict
            Dictionary of {sequence_id: sequence}
        flanking : int
            Number of bases to include on each side (default: 2)
            
        Returns:
        --------
        matplotlib.figure.Figure : Logo figure
        """
        # Find all occurrences of the motif with flanking regions
        aligned_seqs = []
        
        for seq_id, sequence in sequences.items():
            # Find all positions of motif
            start = 0
            while True:
                pos = sequence.find(motif, start)
                if pos == -1:
                    break
                
                # Extract with flanking
                left_flank = max(0, pos - flanking)
                right_flank = min(len(sequence), pos + len(motif) + flanking)
                
                segment = sequence[left_flank:right_flank]
                
                # Pad if at edges
                if pos < flanking:
                    segment = 'N' * (flanking - pos) + segment
                if pos + len(motif) + flanking > len(sequence):
                    segment = segment + 'N' * (pos + len(motif) + flanking - len(sequence))
                
                aligned_seqs.append(segment)
                start = pos + 1
        
        if len(aligned_seqs) == 0:
            fig, ax = plt.subplots(figsize=(10, 3))
            ax.text(0.5, 0.5, f'Motif "{motif}" not found in sequences', 
                   ha='center', va='center')
            ax.axis('off')
            return fig
        
        # Create position frequency matrix
        max_len = max(len(seq) for seq in aligned_seqs)
        counts_matrix = pd.DataFrame(0, 
                                     index=['A', 'C', 'G', 'T'], 
                                     columns=range(max_len))
        
        for seq in aligned_seqs:
            for pos, base in enumerate(seq):
                if base in 'ACGT':
                    counts_matrix.loc[base, pos] += 1
        
        # Convert to frequency matrix
        freq_matrix = counts_matrix / counts_matrix.sum(axis=0)
        freq_matrix = freq_matrix.T  # Transpose for logomaker
        
        # Create logo
        fig, ax = plt.subplots(figsize=(max(8, max_len * 0.5), 3))
        
        logo = logomaker.Logo(freq_matrix, ax=ax, color_scheme='classic')
        
        # Highlight the motif region
        motif_start = flanking
        motif_end = flanking + len(motif)
        ax.axvspan(motif_start - 0.5, motif_end - 0.5, 
                  alpha=0.2, color='yellow', zorder=-1)
        
        # Customize
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_xlabel('Position', fontsize=11)
        ax.set_title(f'Sequence Logo: {motif} (n={len(aligned_seqs)} occurrences)', 
                    fontsize=12, fontweight='bold')
        
        # Add position labels
        positions = list(range(-flanking, len(motif) + flanking))
        ax.set_xticks(range(len(positions)))
        ax.set_xticklabels(positions)
        
        plt.tight_layout()
        
        return fig
    
    def create_enrichment_plot(self, motif_df, top_n=15):
        """
        Create a bar plot of fold enrichment for top motifs.
        
        Parameters:
        -----------
        motif_df : pandas.DataFrame
            DataFrame with motif enrichment information
        top_n : int
            Number of top motifs to display (default: 15)
            
        Returns:
        --------
        matplotlib.figure.Figure : Bar plot figure
        """
        if len(motif_df) == 0:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, 'No motifs found', 
                   ha='center', va='center', fontsize=16)
            ax.axis('off')
            return fig
        
        # Select top motifs by FDR
        top_motifs = motif_df.nsmallest(min(top_n, len(motif_df)), 'FDR')
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, max(6, len(top_motifs) * 0.4)))
        
        # Create color map based on significance
        colors = ['#2E86AB' if sig else '#A23B72' 
                 for sig in top_motifs['Significant']]
        
        # Create horizontal bar plot
        y_pos = np.arange(len(top_motifs))
        ax.barh(y_pos, top_motifs['Fold_Enrichment'], color=colors, alpha=0.8)
        
        # Customize
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_motifs['Motif'])
        ax.invert_yaxis()  # Top to bottom
        ax.set_xlabel('Fold Enrichment', fontsize=12, fontweight='bold')
        ax.set_ylabel('Motif', fontsize=12, fontweight='bold')
        ax.set_title('Top Enriched Motifs', fontsize=14, fontweight='bold', pad=20)
        
        # Add FDR values as text
        for i, (idx, row) in enumerate(top_motifs.iterrows()):
            ax.text(row['Fold_Enrichment'] + 0.1, i, 
                   f"FDR={row['FDR']:.2e}", 
                   va='center', fontsize=9)
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#2E86AB', alpha=0.8, label='Significant'),
            Patch(facecolor='#A23B72', alpha=0.8, label='Not Significant')
        ]
        ax.legend(handles=legend_elements, loc='lower right')
        
        plt.tight_layout()
        
        return fig
    
    def create_length_distribution(self, motif_df):
        """
        Create a histogram of motif lengths.
        
        Parameters:
        -----------
        motif_df : pandas.DataFrame
            DataFrame with motif information
            
        Returns:
        --------
        matplotlib.figure.Figure : Histogram figure
        """
        if len(motif_df) == 0:
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.text(0.5, 0.5, 'No motifs found', 
                   ha='center', va='center', fontsize=16)
            ax.axis('off')
            return fig
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Separate significant and non-significant
        significant = motif_df[motif_df['Significant']]
        non_significant = motif_df[~motif_df['Significant']]
        
        # Create histogram
        bins = np.arange(motif_df['Length'].min(), motif_df['Length'].max() + 2) - 0.5
        
        ax.hist(non_significant['Length'], bins=bins, 
               alpha=0.7, label='Not Significant', color='#A23B72')
        ax.hist(significant['Length'], bins=bins, 
               alpha=0.7, label='Significant', color='#2E86AB')
        
        # Customize
        ax.set_xlabel('Motif Length (bp)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Count', fontsize=12, fontweight='bold')
        ax.set_title('Distribution of Motif Lengths', fontsize=14, fontweight='bold', pad=20)
        ax.legend()
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        
        return fig
    
    def create_volcano_plot(self, motif_df):
        """
        Create a volcano plot of fold enrichment vs. -log10(FDR).
        
        Parameters:
        -----------
        motif_df : pandas.DataFrame
            DataFrame with motif enrichment information
            
        Returns:
        --------
        matplotlib.figure.Figure : Volcano plot figure
        """
        if len(motif_df) == 0:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, 'No motifs found', 
                   ha='center', va='center', fontsize=16)
            ax.axis('off')
            return fig
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Calculate -log10(FDR)
        motif_df['neg_log_FDR'] = -np.log10(motif_df['FDR'] + 1e-300)  # Add small value to avoid log(0)
        
        # Separate significant and non-significant
        significant = motif_df[motif_df['Significant']]
        non_significant = motif_df[~motif_df['Significant']]
        
        # Plot
        ax.scatter(non_significant['Fold_Enrichment'], non_significant['neg_log_FDR'],
                  alpha=0.5, s=50, c='#A23B72', label='Not Significant')
        ax.scatter(significant['Fold_Enrichment'], significant['neg_log_FDR'],
                  alpha=0.7, s=50, c='#2E86AB', label='Significant')
        
        # Add horizontal line for FDR threshold
        fdr_threshold = 0.05
        ax.axhline(-np.log10(fdr_threshold), color='red', linestyle='--', 
                  linewidth=1, alpha=0.7, label=f'FDR = {fdr_threshold}')
        
        # Label top motifs
        top_motifs = motif_df.nsmallest(5, 'FDR')
        for idx, row in top_motifs.iterrows():
            ax.annotate(row['Motif'], 
                       xy=(row['Fold_Enrichment'], row['neg_log_FDR']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.8)
        
        # Customize
        ax.set_xlabel('Fold Enrichment', fontsize=12, fontweight='bold')
        ax.set_ylabel('-log₁₀(FDR)', fontsize=12, fontweight='bold')
        ax.set_title('Volcano Plot: Motif Enrichment', fontsize=14, fontweight='bold', pad=20)
        ax.legend()
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        
        return fig
