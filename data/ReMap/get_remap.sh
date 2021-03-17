#!/usr/bin/env bash

# Get ReMap 2018 peaks 
REMAP="remap2018_all_macs2_hg38_v1_2.bed.gz"
if [ ! -f ${REMAP} ]; then
    wget http://remap.univ-amu.fr/storage/remap2018/hg38/MACS/${REMAP}
fi

# Extract peak summits and sort
SUMMITS="summits.bed"
if [ ! -f ${SUMMITS} ]; then
    zless ${REMAP} | cut -f 1,4,7,8 | awk '{print($1"\t"$3"\t"$4"\t"$2);}' | \
    LC_ALL=C sort --parallel=8 -T ./ -k1,1 -k2,2n > ${SUMMITS}
fi

# Extract ENCODE assay accessions (ENCSRs) and sample names
ENC2SAM="encsr2sample.txt"
if [ ! -f ${ENC2SAM} ]; then
    cut -f 4 ${SUMMITS} | grep "^ENCSR" | cut -d "." -f 1,3 | sort | uniq > ${ENC2SAM}
fi