import streamlit as st
import pandas as pd
import numpy as np
import requests
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# PAGE CONFIGURATION
# ==========================================
st.set_page_config(
    page_title="Gene_Explora",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Gene_Explora")
st.markdown("### Population Genomics & Exon-Level $F_{ST}$ Visualization Tool")

# Standard 1000 Genomes continental order (Non-Alphabetical Geographic Grouping)
CONTINENT_ORDER = {
    'AFR': ['YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB'],
    'AMR': ['MXL', 'PUR', 'CLM', 'PEL'],
    'EAS': ['CHB', 'JPT', 'CHS', 'CDX', 'KHV'],
    'EUR': ['CEU', 'TSI', 'FIN', 'GBR', 'IBS'],
    'SAS': ['GIH', 'PJL', 'BEB', 'STU', 'ITU']
}

POP_TO_CONTINENT = {pop: cont for cont, pops in CONTINENT_ORDER.items() for pop in pops}
SORTED_POPS = [pop for cont in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS'] for pop in CONTINENT_ORDER[cont]]

# ==========================================
# DATA LOADING & ENSEMBL API FUNCTIONS
# ==========================================

@st.cache_data
def load_gene_table(file_path_or_buffer):
    """Dynamically reads gene coordinates file."""
    df = pd.read_csv(file_path_or_buffer, sep='\t')
    valid_df = df[df['status'] == 'ok'].copy()
    valid_df['label'] = valid_df.apply(
        lambda r: f"{r['input_symbol']} ({r['gene_name']}) | Ensembl: {r['ensembl_id']} | Chr{r['chromosome']}", axis=1
    )
    return valid_df

@st.cache_data
def fetch_ensembl_exons(ensembl_id):
    """Fetches exon coordinates directly from Ensembl REST API."""
    url = f"https://rest.ensembl.org/lookup/id/{ensembl_id}?expand=1"
    headers = {"Content-Type": "application/json"}
    
    try:
        response = requests.get(url, headers=headers, timeout=10)
        if not response.ok:
            st.error(f"Ensembl API error: {response.status_code}")
            return []
        
        data = response.json()
        transcripts = data.get('Transcript', [])
        
        chosen_tx = None
        for tx in transcripts:
            if tx.get('is_canonical') == 1:
                chosen_tx = tx
                break
        if not chosen_tx and transcripts:
            chosen_tx = transcripts[0]
            
        exons = []
        if chosen_tx:
            for ex in chosen_tx.get('Exon', []):
                exons.append({
                    'exon_id': ex.get('id'),
                    'start': int(ex.get('start')),
                    'end': int(ex.get('end')),
                    'strand': ex.get('strand'),
                    'length': int(ex.get('end')) - int(ex.get('start')) + 1
                })
        exons.sort(key=lambda x: x['start'])
        return exons
    except Exception as e:
        st.warning(f"Failed to retrieve Ensembl exons for {ensembl_id}: {e}")
        return []

def generate_sample_fst_data(gene_row, num_snps=80):
    """Generates synthetic population FST values across geographic regions."""
    start = int(gene_row['start'])
    end = int(gene_row['end'])
    positions = np.sort(np.random.randint(start - 2000, end + 2000, size=num_snps))
    
    data = []
    np.random.seed(42)
    for pos in positions:
        row = {'chromosome': str(gene_row['chromosome']), 'position': pos}
        base_afr = np.random.beta(1, 5)
        base_eur = np.random.beta(2, 4)
        base_eas = np.random.beta(3, 3)
        
        for pop in SORTED_POPS:
            cont = POP_TO_CONTINENT[pop]
            if cont == 'AFR':
                val = np.clip(base_afr + np.random.normal(0, 0.05), 0, 1)
            elif cont in ['EUR', 'AMR']:
                val = np.clip(base_eur + np.random.normal(0, 0.05), 0, 1)
            else:
                val = np.clip(base_eas + np.random.normal(0, 0.05), 0, 1)
            row[pop] = round(val, 4)
        data.append(row)
        
    return pd.DataFrame(data)

def classify_coordinate(pos, exons, gene_start, gene_end, flank_bp=100):
    """Compares SNP coordinate against exon boundaries, 100bp flanking regions, and introns."""
    for idx, ex in enumerate(exons, 1):
        if ex['start'] <= pos <= ex['end']:
            return f"Exon_{idx} ({ex['exon_id']})"
        if (ex['start'] - flank_bp <= pos < ex['start']) or (ex['end'] < pos <= ex['end'] + flank_bp):
            return f"Exon_{idx} Flank (±{flank_bp}bp)"
            
    if gene_start <= pos <= gene_end:
        return "Intron / Non-Coding Region"
    elif pos < gene_start:
        return f"Upstream ({gene_start - pos} bp)"
    else:
        return f"Downstream ({pos - gene_end} bp)"

# ==========================================
# SIDEBAR CONTROLS
# ==========================================
st.sidebar.header("⚙️ Project Controls")

uploaded_genes = st.sidebar.file_uploader("Upload gene_coordinates.tsv", type=['tsv', 'csv'])
gene_file = uploaded_genes if uploaded_genes is not None else "gene_coordinates.tsv"

try:
    genes_df = load_gene_table(gene_file)
except Exception as e:
    st.error(f"Could not load gene coordinates file: {e}")
    st.stop()

selected_label = st.sidebar.selectbox(
    "Select Blood Gene:",
    options=genes_df['label'].tolist()
)

selected_gene = genes_df[genes_df['label'] == selected_label].iloc[0]

st.sidebar.subheader("Data Customization")
fst_file = st.sidebar.file_uploader("Optional: Upload Custom FST File (TSV/CSV)", type=['tsv', 'csv'])
window_size = st.sidebar.slider("FST Window Size (bp for averaging)", min_value=1, max_value=5000, value=500, step=100)

# ==========================================
# MAIN APPLICATION LOGIC
# ==========================================

st.subheader(f"Selected Gene: **{selected_gene['input_symbol']}** ({selected_gene['gene_name']})")

col1, col2, col3, col4 = st.columns(4)
col1.metric("Ensembl ID", str(selected_gene['ensembl_id']))
col2.metric("Chromosome", f"Chr {selected_gene['chromosome']}")
col3.metric("Gene Span", f"{int(selected_gene['end'] - selected_gene['start']):,} bp")
col4.metric("Strand", str(selected_gene['strand']))

with st.spinner("Fetching exon coordinates from Ensembl API..."):
    exons = fetch_ensembl_exons(selected_gene['ensembl_id'])

st.markdown("#### 📐 Exon Coordinates & 100bp Boundary Lengths")
if exons:
    exon_data = []
    total_coding_len = 0
    total_flank_100_len = 0
    
    for idx, ex in enumerate(exons, 1):
        ex_len = ex['length']
        flank_100_len = ex_len + 200
        total_coding_len += ex_len
        total_flank_100_len += flank_100_len
        
        exon_data.append({
            'Exon Number': f"Exon {idx}",
            'Ensembl Exon ID': ex['exon_id'],
            'Start': ex['start'],
            'End': ex['end'],
            'Exon Length (bp)': ex_len,
            'Length + 100bp Flank (bp)': flank_100_len
        })
        
    exon_df = pd.DataFrame(exon_data)
    st.dataframe(exon_df, use_container_width=True)
    st.caption(f"**Total Coding Exon Length:** {total_coding_len:,} bp | **Total Exon + 100bp Flank Length:** {total_flank_100_len:,} bp")
else:
    st.info("No detailed exon structure returned from Ensembl. Using full gene span.")

if fst_file is not None:
    fst_df = pd.read_csv(fst_file, sep='\t' if fst_file.name.endswith('.tsv') else ',')
else:
    fst_df = generate_sample_fst_data(selected_gene)

gene_start = int(selected_gene['start'])
gene_end = int(selected_gene['end'])
fst_df['Genomic_Region'] = fst_df['position'].apply(
    lambda pos: classify_coordinate(pos, exons, gene_start, gene_end, flank_bp=100)
)

pop_cols = [c for c in fst_df.columns if c in POP_TO_CONTINENT or c in SORTED_POPS]
if not pop_cols:
    pop_cols = [c for c in fst_df.columns if c not in ['chromosome', 'position', 'Genomic_Region']]

ordered_pop_cols = [p for p in SORTED_POPS if p in pop_cols] + [p for p in pop_cols if p not in SORTED_POPS]

with st.expander("🔍 View Per-SNP $F_{ST}$ Coordinate & Exon Mapping Table"):
    st.dataframe(fst_df[['chromosome', 'position', 'Genomic_Region'] + ordered_pop_cols], use_container_width=True)

# ==========================================
# VISUALIZATIONS
# ==========================================

st.markdown("---")
st.subheader("📊 Visualizations")

tab1, tab2 = st.tabs(["🔥 Population/Continent $F_{ST}$ Heatmap", "🧬 Genomic Structure & $F_{ST}$ Track Plot"])

with tab1:
    st.markdown("#### $F_{ST}$ Pairwise Distance Matrix Grouped by Geographic Continent")
    st.caption("Populations are ordered geographically by continent (AFR → AMR → EAS → EUR → SAS).")
    
    available_pops = [p for p in SORTED_POPS if p in fst_df.columns]
    
    if available_pops:
        pop_data = fst_df[available_pops]
        pairwise_matrix = 1 - pop_data.corr()
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        sns.heatmap(
            pairwise_matrix,
            cmap="YlOrRd",
            cbar_kws={'label': 'Fst Distance'},
            xticklabels=available_pops,
            yticklabels=available_pops,
            ax=ax
        )
        
        gene_name = selected_gene['input_symbol']
        ax.set_title(f"Pairwise $F_{{ST}}$ Distance Heatmap for {gene_name}", fontsize=13, fontweight='bold')
        
        plt.xticks(rotation=45, ha='right', fontsize=8)
        plt.yticks(rotation=0, fontsize=8)
        plt.tight_layout()
        st.pyplot(fig)
    else:
        st.warning("No population columns found in the dataset matching standard population codes.")

with tab2:
    st.markdown("#### Genomic Map: Exon Architecture with Overlaid $F_{ST}$ Values")
    
    fig, (ax_gene, ax_fst) = plt.subplots(2, 1, figsize=(14, 8), sharex=True, gridspec_kw={'height_ratios': [1, 2.5]})
    
    # 1. Gene Architecture Track
    ax_gene.hlines(1, gene_start, gene_end, color='black', linewidth=3, label='Intron / Gene Body')
    
    for idx, ex in enumerate(exons, 1):
        ax_gene.add_patch(
            plt.Rectangle(
                (ex['start'], 0.7), ex['length'], 0.6,
                color='crimson', alpha=0.8, zorder=3
            )
        )
        ax_gene.text(
            ex['start'] + ex['length'] / 2, 1.35, f"E{idx}",
            ha='center', va='bottom', fontsize=8, fontweight='bold', color='darkred'
        )
        ax_gene.add_patch(
            plt.Rectangle(
                (ex['start'] - 100, 0.8), ex['length'] + 200, 0.4,
                color='orange', alpha=0.3, zorder=2
            )
        )
        
    ax_gene.set_ylim(0, 1.8)
    ax_gene.set_yticks([])
    ax_gene.set_title(f"Genomic Structure & $F_{{ST}}$ Signals for {selected_gene['input_symbol']} (Chr{selected_gene['chromosome']})", fontsize=14, fontweight='bold')
    ax_gene.legend(["Gene Body", "Exon", "100bp Boundary"], loc='upper right', fontsize=8)
    
    # 2. FST Signals Track
    pop_continents = [POP_TO_CONTINENT.get(p, 'Other') for p in ordered_pop_cols]
    colors = sns.color_palette("Set2", len(set(pop_continents)))
    cont_color_map = dict(zip(set(pop_continents), colors))
    
    for cont in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        cont_pops = [p for p in ordered_pop_cols if POP_TO_CONTINENT.get(p) == cont]
        if cont_pops:
            avg_fst = fst_df[cont_pops].mean(axis=1)
            ax_fst.plot(
                fst_df['position'], avg_fst,
                marker='o', linestyle='-', linewidth=1.5, markersize=4,
                label=f"Continent: {cont}", color=cont_color_map.get(cont)
            )
            
    ax_fst.set_xlabel("Genomic Coordinate (Chr Position)", fontsize=11)
    ax_fst.set_ylabel("$F_{{ST}}$ Value", fontsize=11)
    ax_fst.set_ylim(-0.05, 1.05)
    ax_fst.grid(True, linestyle='--', alpha=0.5)
    ax_fst.legend(loc='upper right', title="Continental Groups")
    
    plt.tight_layout()
    st.pyplot(fig)