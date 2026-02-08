"""
AptaMotif Analyzer - RNA Aptamer Motif Discovery and Analysis Tool
Author: Lab Bioinformatics Pipeline
Version: 1.0
"""

import streamlit as st
import pandas as pd
import numpy as np
from io import StringIO
import json
from datetime import datetime

# Import custom modules
from motif_analysis import MotifAnalyzer
from structure_analysis import StructureAnalyzer
from visualizations import MotifVisualizer
from statistics_module import StatisticalAnalyzer

# Page configuration
st.set_page_config(
    page_title="AptaMotif Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
    <style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        text-align: center;
        margin-bottom: 2rem;
    }
    </style>
""", unsafe_allow_html=True)

# Initialize session state
if 'pool_configs' not in st.session_state:
    st.session_state.pool_configs = {
        'Default_N71': {
            'forward_primer': 'TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATT',
            'reverse_complement': 'AGATAGTAAGTGCAATCT',
            'random_region_length': 71,
            'description': 'Standard N71 pool'
        }
    }

if 'current_config' not in st.session_state:
    st.session_state.current_config = 'Default_N71'

if 'sequences' not in st.session_state:
    st.session_state.sequences = None

if 'results' not in st.session_state:
    st.session_state.results = None


def main():
    # Header
    st.markdown('<h1 class="main-header">🧬 AptaMotif Analyzer</h1>', unsafe_allow_html=True)
    st.markdown('<p class="sub-header">RNA Aptamer Motif Discovery & Secondary Structure Analysis</p>', 
                unsafe_allow_html=True)
    
    # Sidebar for configuration
    with st.sidebar:
        st.header("⚙️ Configuration")
        
        # Pool configuration selection
        st.subheader("Pool Configuration")
        config_names = list(st.session_state.pool_configs.keys())
        selected_config = st.selectbox(
            "Select Pool Configuration",
            config_names,
            index=config_names.index(st.session_state.current_config)
        )
        st.session_state.current_config = selected_config
        
        # Display current configuration
        current = st.session_state.pool_configs[selected_config]
        st.info(f"**Description:** {current['description']}\n\n"
                f"**Forward Primer:** {current['forward_primer'][:20]}...\n\n"
                f"**Random Region:** N{current['random_region_length']}")
        
        # Add/Edit configuration
        with st.expander("➕ Add/Edit Pool Configuration"):
            new_config_name = st.text_input("Configuration Name", value="New_Pool")
            new_forward = st.text_area("Forward Primer (5'→3')", 
                                       value=current['forward_primer'])
            new_reverse_comp = st.text_area("Reverse Complement (after N region)", 
                                            value=current['reverse_complement'])
            new_n_length = st.number_input("Random Region Length", 
                                          min_value=10, max_value=200, 
                                          value=current['random_region_length'])
            new_description = st.text_input("Description", value=current['description'])
            
            if st.button("💾 Save Configuration"):
                st.session_state.pool_configs[new_config_name] = {
                    'forward_primer': new_forward.strip().upper(),
                    'reverse_complement': new_reverse_comp.strip().upper(),
                    'random_region_length': new_n_length,
                    'description': new_description
                }
                st.session_state.current_config = new_config_name
                st.success(f"Configuration '{new_config_name}' saved!")
                st.rerun()
        
        # Export/Import configurations
        with st.expander("💾 Export/Import Configurations"):
            # Export
            config_json = json.dumps(st.session_state.pool_configs, indent=2)
            st.download_button(
                label="📥 Download Configurations",
                data=config_json,
                file_name="aptamer_pool_configs.json",
                mime="application/json"
            )
            
            # Import
            uploaded_config = st.file_uploader("📤 Upload Configuration File", type=['json'])
            if uploaded_config is not None:
                try:
                    imported_configs = json.load(uploaded_config)
                    st.session_state.pool_configs.update(imported_configs)
                    st.success("Configurations imported successfully!")
                    st.rerun()
                except Exception as e:
                    st.error(f"Error importing configurations: {str(e)}")
        
        st.divider()
        
        # Analysis parameters
        st.subheader("Analysis Parameters")
        min_motif_length = st.slider("Minimum Motif Length", 5, 15, 5)
        max_motif_length = st.slider("Maximum Motif Length", 5, 20, 15)
        min_occurrences = st.number_input("Minimum Sequences Sharing Motif", 
                                         min_value=2, max_value=100, value=2)
        fdr_threshold = st.number_input("FDR Threshold", 
                                       min_value=0.001, max_value=0.5, 
                                       value=0.05, step=0.01)
        
        st.divider()
        
        # Structure analysis options
        st.subheader("Structure Analysis")
        run_structure = st.checkbox("Run ViennaRNA Structure Prediction", value=True)
        if run_structure:
            temperature = st.slider("Temperature (°C)", 20, 50, 37)
    
    # Main content area
    tab1, tab2, tab3, tab4 = st.tabs(["📝 Input Sequences", "🔍 Motif Analysis", 
                                        "🧬 Structure Analysis", "📊 Results & Export"])
    
    with tab1:
        st.header("Input DNA Sequences")
        st.info("💡 Enter your DNA sequences (they will be converted to RNA for structure analysis). "
                "The tool will automatically extract the random region based on your pool configuration.")
        
        # Input method selection
        input_method = st.radio("Input Method", ["Paste Sequences", "Upload File"])
        
        if input_method == "Paste Sequences":
            st.subheader("Paste Sequences")
            st.caption("Format: FASTA format or one sequence per line")
            
            sequence_input = st.text_area(
                "Sequences",
                height=300,
                placeholder=">Clone_1\nTTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATT...\n>Clone_2\nTTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATT..."
            )
            
            if st.button("🚀 Process Sequences", type="primary"):
                if sequence_input:
                    sequences = parse_sequences(sequence_input)
                    if sequences:
                        st.session_state.sequences = sequences
                        st.success(f"✅ Processed {len(sequences)} sequences successfully!")
                        display_sequence_summary(sequences)
                else:
                    st.warning("Please enter sequences first.")
        
        else:  # Upload File
            st.subheader("Upload Sequence File")
            uploaded_file = st.file_uploader(
                "Choose a file (FASTA or text format)",
                type=['fasta', 'fa', 'txt', 'seq']
            )
            
            if uploaded_file is not None:
                sequence_input = StringIO(uploaded_file.getvalue().decode("utf-8")).read()
                
                if st.button("🚀 Process Uploaded File", type="primary"):
                    sequences = parse_sequences(sequence_input)
                    if sequences:
                        st.session_state.sequences = sequences
                        st.success(f"✅ Processed {len(sequences)} sequences successfully!")
                        display_sequence_summary(sequences)
        
        # Display loaded sequences
        if st.session_state.sequences:
            with st.expander("👀 View Loaded Sequences"):
                df_seqs = pd.DataFrame([
                    {'ID': seq_id, 'Length': len(seq), 'Sequence': seq[:50] + '...' if len(seq) > 50 else seq}
                    for seq_id, seq in st.session_state.sequences.items()
                ])
                st.dataframe(df_seqs, use_container_width=True)
    
    with tab2:
        st.header("Motif Discovery & Enrichment Analysis")
        
        if st.session_state.sequences is None:
            st.warning("⚠️ Please input sequences in the 'Input Sequences' tab first.")
        else:
            if st.button("🔬 Run Motif Analysis", type="primary"):
                with st.spinner("Analyzing motifs..."):
                    # Extract random regions
                    config = st.session_state.pool_configs[st.session_state.current_config]
                    random_regions = extract_random_regions(
                        st.session_state.sequences,
                        config['forward_primer'],
                        config['reverse_complement']
                    )
                    
                    # Run motif analysis
                    analyzer = MotifAnalyzer(
                        random_regions,
                        min_length=min_motif_length,
                        max_length=max_motif_length,
                        min_occurrences=min_occurrences
                    )
                    
                    motifs = analyzer.find_motifs()
                    
                    # Statistical analysis
                    stats_analyzer = StatisticalAnalyzer(
                        motifs,
                        random_regions,
                        config['random_region_length']
                    )
                    
                    results = stats_analyzer.calculate_enrichment(fdr_threshold=fdr_threshold)
                    
                    st.session_state.results = {
                        'motifs': results,
                        'random_regions': random_regions,
                        'config': config
                    }
                    
                    st.success("✅ Motif analysis complete!")
            
            # Display results
            if st.session_state.results:
                display_motif_results(st.session_state.results['motifs'])
    
    with tab3:
        st.header("Secondary Structure Analysis")
        
        if st.session_state.sequences is None:
            st.warning("⚠️ Please input sequences in the 'Input Sequences' tab first.")
        elif not run_structure:
            st.info("ℹ️ Structure analysis is disabled. Enable it in the sidebar to run ViennaRNA prediction.")
        else:
            if st.button("🧬 Run Structure Prediction", type="primary"):
                with st.spinner("Predicting secondary structures with ViennaRNA..."):
                    config = st.session_state.pool_configs[st.session_state.current_config]
                    random_regions = extract_random_regions(
                        st.session_state.sequences,
                        config['forward_primer'],
                        config['reverse_complement']
                    )
                    
                    # Convert DNA to RNA
                    rna_sequences = {seq_id: seq.replace('T', 'U') 
                                    for seq_id, seq in random_regions.items()}
                    
                    # Run structure prediction
                    structure_analyzer = StructureAnalyzer(temperature=temperature)
                    structures = structure_analyzer.predict_structures(rna_sequences)
                    
                    if 'results' not in st.session_state.results or st.session_state.results is None:
                        st.session_state.results = {}
                    
                    st.session_state.results['structures'] = structures
                    
                    st.success("✅ Structure prediction complete!")
            
            # Display structures
            if st.session_state.results and 'structures' in st.session_state.results:
                display_structure_results(st.session_state.results['structures'])
    
    with tab4:
        st.header("Results Summary & Export")
        
        if st.session_state.results is None:
            st.warning("⚠️ No results available. Please run motif analysis first.")
        else:
            st.subheader("📊 Export Options")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                if 'motifs' in st.session_state.results:
                    motif_df = st.session_state.results['motifs']
                    csv = motif_df.to_csv(index=False)
                    st.download_button(
                        label="📥 Download Motif Table (CSV)",
                        data=csv,
                        file_name=f"motif_enrichment_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                        mime="text/csv"
                    )
            
            with col2:
                if 'structures' in st.session_state.results:
                    structure_text = format_structures_for_export(
                        st.session_state.results['structures']
                    )
                    st.download_button(
                        label="📥 Download Structures (TXT)",
                        data=structure_text,
                        file_name=f"structures_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                        mime="text/plain"
                    )
            
            with col3:
                # Generate full report
                if st.button("📄 Generate Full Report"):
                    st.info("Full report generation coming soon!")
            
            # Display summary statistics
            st.subheader("📈 Summary Statistics")
            
            if 'motifs' in st.session_state.results:
                motif_df = st.session_state.results['motifs']
                significant_motifs = motif_df[motif_df['FDR'] < fdr_threshold]
                
                col1, col2, col3, col4 = st.columns(4)
                col1.metric("Total Sequences", len(st.session_state.sequences))
                col2.metric("Motifs Found", len(motif_df))
                col3.metric("Significant Motifs", len(significant_motifs))
                col4.metric("FDR Threshold", f"{fdr_threshold:.3f}")


def parse_sequences(text):
    """Parse sequences from FASTA or plain text format."""
    sequences = {}
    current_id = None
    current_seq = []
    
    lines = text.strip().split('\n')
    
    for i, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        
        if line.startswith('>'):
            # Save previous sequence
            if current_id:
                sequences[current_id] = ''.join(current_seq).upper()
            # Start new sequence
            current_id = line[1:].strip() or f"Sequence_{i+1}"
            current_seq = []
        else:
            # Handle both FASTA and plain text
            if current_id is None:
                current_id = f"Sequence_{len(sequences)+1}"
            # Remove any whitespace and non-DNA characters
            cleaned = ''.join(c for c in line.upper() if c in 'ATCG')
            current_seq.append(cleaned)
    
    # Save last sequence
    if current_id:
        sequences[current_id] = ''.join(current_seq).upper()
    
    return sequences


def extract_random_regions(sequences, forward_primer, reverse_complement):
    """Extract the random region from full sequences."""
    random_regions = {}
    
    for seq_id, seq in sequences.items():
        # Find forward primer
        fwd_pos = seq.find(forward_primer)
        if fwd_pos == -1:
            st.warning(f"⚠️ Forward primer not found in {seq_id}")
            continue
        
        # Find reverse complement
        rev_pos = seq.find(reverse_complement, fwd_pos + len(forward_primer))
        if rev_pos == -1:
            st.warning(f"⚠️ Reverse complement not found in {seq_id}")
            continue
        
        # Extract random region
        start = fwd_pos + len(forward_primer)
        random_region = seq[start:rev_pos]
        random_regions[seq_id] = random_region
    
    return random_regions


def display_sequence_summary(sequences):
    """Display summary of loaded sequences."""
    st.subheader("Sequence Summary")
    
    lengths = [len(seq) for seq in sequences.values()]
    
    col1, col2, col3 = st.columns(3)
    col1.metric("Total Sequences", len(sequences))
    col2.metric("Avg Length", f"{np.mean(lengths):.0f} bp")
    col3.metric("Length Range", f"{min(lengths)}-{max(lengths)} bp")


def display_motif_results(motif_df):
    """Display motif analysis results."""
    st.subheader("Enriched Motifs")
    
    # Filter significant motifs
    significant = motif_df[motif_df['FDR'] < 0.05]
    
    if len(significant) == 0:
        st.info("No statistically significant motifs found at FDR < 0.05")
    else:
        st.success(f"Found {len(significant)} significant motifs!")
    
    # Display table
    st.dataframe(
        motif_df.style.background_gradient(subset=['FDR'], cmap='RdYlGn_r'),
        use_container_width=True
    )
    
    # Visualizations
    if len(motif_df) > 0:
        visualizer = MotifVisualizer()
        
        # Create heatmap
        st.subheader("Motif Occurrence Heatmap")
        fig_heatmap = visualizer.create_heatmap(
            st.session_state.results['motifs'],
            st.session_state.results['random_regions']
        )
        st.pyplot(fig_heatmap)
        
        # Create sequence logos for top motifs
        st.subheader("Sequence Logos (Top 5 Motifs)")
        top_motifs = motif_df.nsmallest(5, 'FDR')
        
        for idx, row in top_motifs.iterrows():
            with st.expander(f"Motif: {row['Motif']} (FDR: {row['FDR']:.2e})"):
                fig_logo = visualizer.create_sequence_logo(
                    row['Motif'],
                    st.session_state.results['random_regions']
                )
                st.pyplot(fig_logo)


def display_structure_results(structures):
    """Display secondary structure results."""
    st.subheader("Predicted Secondary Structures")
    
    for seq_id, struct_data in structures.items():
        with st.expander(f"{seq_id} (MFE: {struct_data['mfe']:.2f} kcal/mol)"):
            col1, col2 = st.columns([1, 1])
            
            with col1:
                st.text("Sequence:")
                st.code(struct_data['sequence'])
                st.text("Structure (dot-bracket):")
                st.code(struct_data['structure'])
            
            with col2:
                if 'visualization' in struct_data:
                    st.image(struct_data['visualization'])


def format_structures_for_export(structures):
    """Format structures for text export."""
    lines = []
    for seq_id, struct_data in structures.items():
        lines.append(f">{seq_id}")
        lines.append(f"Sequence: {struct_data['sequence']}")
        lines.append(f"Structure: {struct_data['structure']}")
        lines.append(f"MFE: {struct_data['mfe']:.2f} kcal/mol")
        lines.append("")
    
    return "\n".join(lines)


if __name__ == "__main__":
    main()
