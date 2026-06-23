#!/usr/bin/env bash

set -euo pipefail

#-------[ Required tools ]--------#
brew install vg

echo
echo "=== Environment ==="
vg version

#-------[ Directories ]--------#
# Input directory
WD_DIR="/Users/beuser/Desktop/vg-files" #FIXME

cd "${WD_DIR}"

#-------[ Input files ]--------#
# Required input files are available at:
#'https://github.com/human-pangenomics/hpp_pangenome_resources'

URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz"

GBZ="$(basename "${URL}")"

CHUNK="chunk_0_chr6.gbz"

CHUNK_PATHS="chr6_paths.txt"

if [ ! -f "${GBZ}" ]; then
echo
echo "Downloading HPRC graph..."
wget -c "${URL}"
fi

echo
echo "Graph file:"
ls -lh "${GBZ}"

# Number of samples
#vg gbwt -Z "${GBZ}" -S

# Number of haplotypes
#vg gbwt -Z "${GBZ}" -H

# GBWT listing
#vg gbwt -Z "${GBZ}" -T

#-------[ Extract reference genome ]--------#
echo
echo "Extracting reference genome..."

vg paths \
   -x "${GBZ}" \
   -F \
   -Q GRCh38#0#chr6 \
   > GRCh38_chr6.fa

#vg paths \
#   -x "${GBZ}" \
#   -F \
#   -Q CHM13#0#chr6 \
#   > CHM13_chr6.fa

#-------[ Build chunk graph ]--------#
echo
echo "Building chunk graph..."

vg chunk \
   -x "${GBZ}" \
   --gbz \
   --contig chr6 \
   > "${CHUNK}"

#-------[ Export paths in chunk  ]--------#
echo
echo "Exporting paths in chunk ..."

vg paths \
   -x "${CHUNK}" \
   -L \
   > "${CHUNK_PATHS}"

echo
echo "===========End of script==========="