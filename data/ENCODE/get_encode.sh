#!/usr/bin/env bash

for GENOME in hg38 mm10; do
    ./get_encode.py --genome ${GENOME} --feature accessibility --out-dir ./${GENOME}/
    ./get_encode.py --genome ${GENOME} --feature tf --out-dir ./${GENOME}/
#    # Extract uniform DNase-seq regions of length 150 bp
#    if [ ! -f ./${G}/DNase-seq.150bp.bed ]; then
#        cut -f 1,2,4,6 ./${G}/DNase-seq.bed | awk '{print($1"\t"$2+$4-75"\t"$2+$4+75"\t"$3);}' | \
#        LC_ALL=C sort --parallel=8 -T ./ -k1,1 -k2,2n > ./${G}/DNase-seq.150bp.bed
#    fi
done
