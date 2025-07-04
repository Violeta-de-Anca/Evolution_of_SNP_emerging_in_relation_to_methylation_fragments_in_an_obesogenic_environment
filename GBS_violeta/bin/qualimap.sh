#!/bin/bash

# load modules
module load bioinfo-tools
module load QualiMap/2.2.1

prefix="${1/.bam/}"

qualimap bamqc \
-bam $1 \
-outdir  "$prefix" \
-outfile "$prefix" \
-outformat "PDF"  >& "${prefix}-qualimap.log"
