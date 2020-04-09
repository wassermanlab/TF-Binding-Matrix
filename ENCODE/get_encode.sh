#!/bin/bash

for G in hg38 mm10
do
    python get_encode.py --genome ${G} --feature accessibility --out-dir ./${G}/
    python get_encode.py --genome ${G} --feature tf --out-dir ./${G}/

    # Extract uniform DNase-seq regions of length 150 bp
    if [ ! -f ./${G}/summits.bed ]; then

        cut -f 1,2,4,6 ./${G}/DNase-seq.bed | awk '{print($1"\t"$2+$4-75"\t"$2+$4+75"\t"$3);}' | \
        LC_ALL=C sort --parallel=8 -T ./ -k1,1 -k2,2n > ./${G}/DNase-seq.150bp.bed
    fi
done