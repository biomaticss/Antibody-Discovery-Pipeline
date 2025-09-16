#!/bin/bash

# MIXCR.sh (Parameterized for Snakemake)
# Usage: bash MIXCR.sh <sample_id> <r1.fq.gz> <r2.fq.gz> <output_dir>

SAMPLE_ID=$1
R1=$2
R2=$3
OUTPUT_DIR=$4
OUTPUT_VDJCA="${OUTPUT_DIR}/${SAMPLE_ID}.vdjca"
OUTPUT_CLNS="${OUTPUT_DIR}/${SAMPLE_ID}_assemble.clns"
OUTPUT_TSV_IGH="${OUTPUT_DIR}/${SAMPLE_ID}_VH_assemble_Clones_IGH.tsv"
OUTPUT_TSV_IGK="${OUTPUT_DIR}/${SAMPLE_ID}_VK_assemble_Clones_IGK.tsv"

echo "Processing sample: $SAMPLE_ID"

# Align
mixcr -Xmx500g align \
    --preset generic-amplicon \
    --species hsa \
    --rna \
    --assemble-clonotypes-by CDR3 \
    --rigid-left-alignment-boundary \
    --rigid-right-alignment-boundary C \
    -OmergerParameters.minimalOverlap=10 \
    -OmergerParameters.minimalIdentity=0.9 \
    "$R1" "$R2" "$OUTPUT_VDJCA"

# Assemble
mixcr -Xmx500g assemble \
    --assemble-clonotypes-by "{FR1Begin:FR4End}" \
    -OseparateByC=true \
    -ObadQualityThreshold=0 \
    "$OUTPUT_VDJCA" "$OUTPUT_CLNS"

# Export IGH
mixcr -Xmx500g exportClones \
    "$OUTPUT_CLNS" \
    "$OUTPUT_TSV_IGH"

# Export IGK (Light chain; adjust if IGL)
mixcr -Xmx500g exportClonesLightChains \
    "$OUTPUT_CLNS" \
    "$OUTPUT_TSV_IGK"

# Cleanup
rm -f "$OUTPUT_VDJCA" "$OUTPUT_CLNS"

if [[ -s "$OUTPUT_TSV_IGH" && -s "$OUTPUT_TSV_IGK" ]]; then
    echo "Completed: $SAMPLE_ID"
else
    echo "Error: Empty output for $SAMPLE_ID"
    exit 1
fi
