🧬 Gene_Exploracdcdc
Is an interactive, Streamlit-based bioinformatics dashboard for exploring population genomics, exon-level genetic architecture, and fixated index ($F_{ST}$) divergence across global populations.The tool dynamically retrieves exon boundaries from the Ensembl REST API, maps genetic variants to coding and flanking boundaries, and provides interactive visualizations categorized by continental geographic region rather than alphabetical order.


Key FeaturesDynamic Ensembl REST API Integration
Automatically fetches transcript and exon coordinate boundaries using Ensembl Gene IDs (e.g., ENSG00000244734).Exon & 100bp Flank Boundary Mapping: Calculates exact coding lengths as well as extended boundaries ($\pm 100\text{ bp}$) to classify variants into Exons, Flanks, Introns, or Intergenic regions.Geographic Continent-Grouped Heatmaps: Displays population pairwise $F_{ST}$ distance matrices ordered geographically ($\text{AFR} \rightarrow \text{AMR} \rightarrow \text{EAS} \rightarrow \text{EUR} \rightarrow \text{SAS}$).Dual-Track Genomic Structure Plot: Simultaneously plots transcript exon architecture alongside per-SNP or per-window $F_{ST}$ signal lines.Fully Flexible: Zero hardcoded paths or fixed gene lists—dynamically loads any input gene_coordinates.tsv file.


🛠️ Installation & Setup 
Ensure you have Python 3.9+ installed. If running on Ubuntu or WSL (Windows Subsystem for Linux), set up a virtual environment to avoid PEP 668 system package restrictions.

