#!/bin/bash

# Get ReMap 2020 peaks
REMAP="remap2020_all_macs2_hg38_v1_0.bed.gz"
if [ ! -f ${REMAP} ]; then
    wget http://remap.univ-amu.fr/storage/remap2020/hg38/MACS2/${REMAP}
fi

# Extract datasets
DATASETS="remap2020_datasets.txt"
if [ ! -f ${DATASETS} ]; then
    zless remap2020_all_macs2_hg38_v1_0.bed.gz | cut -f 4 | sort | uniq -c > ${DATASETS}
fi
