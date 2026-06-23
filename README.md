#========================================#
# CpG Stability Analysis - Pangenome Graph
#========================================#
#
# Goal:
#   Measure CpG variation across pangenome haplotypes.
#
# Main idea:
#   Reference CpGs are not invariant.
#   SNPs may create or destroy CpG sites across haplotypes.
#
# Key hypothesis being tested:
#   CpG sites are not invariant genomic features.
#   Sequence variation can create or disrupt CpGs,
#   causing their presence to vary among haplotypes.
#
# Workflow:
#  1. Download HPRC GBZ graph
#  2. Verify graph contents
#  3. Enumerate embedded haplotypes
#  4. Extract chr6 subgraph
#  5. Export paths and graph structure
#  6. Identify reference CpG sites
#  7. Measure CpG retention across haplotypes
#
# Biological next step (not in this script):
#  1. Extract reference CCGG motifs
#  2. Retrieve ±15 bp flanks
#  3. Search haplotype paths
#  4. Evaluate CpG state
#  5. Classify loci: Stable / Variable / Lost / Novel
#
# Data products:
#
# GBZ
# ├── extract reference paths
# │      ├── GRCh38_chr6.fa
# │      └── CHM13_chr6.fa (NOT USED)
# ├── extract chr6 chunk
# │      └── chunk_0_chr6.gbz
# │             └── chr6_paths.txt
# └── scan CCGG motifs (in R)
#
#========================================#
