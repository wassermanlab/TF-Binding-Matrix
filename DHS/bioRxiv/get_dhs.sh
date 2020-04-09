#!/bin/bash

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# From https://www.biorxiv.org/content/10.1101/822510v2
# Index and biological spectrum of accessible DNA elements in the human genome
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# DNase I Hypersensitive Sites (DHSs) defined using a consensus-based approach 
# across 733 human biosamples. Each DHS is annotated with a signal level which
# can be used as a confidence score, as well as estimates for its summit posi-
# tion and summit variability (as captured by the 'core' coordinates).
ENCODE="ENCFF503GCK"
if [ ! -f ${ENCODE}.tsv ]; then
    wget https://www.encodeproject.org/files/${ENCODE}/@@download/${ENCODE}.tsv
fi

# Get chrom sizes (i.e. option "-g") for bedtools slop
GENOME="hg38.chrom.sizes"
if [ ! -f ${GENOME} ]; then
    wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/${GENOME}
fi

# Extract DHS summits and sort
SUMMITS="summits.bed"
if [ ! -f ${SUMMITS} ]; then
    grep "^chr" ${ENCODE}.tsv | cut -f 1,7 | awk '{print($1"\t"$2"\t"$2);}' | \
    LC_ALL=C sort --parallel=8 -T ./ -k1,1 -k2,2n > ${SUMMITS}
fi

# L = 200-1000 bp
# Increase in steps of 100 bp
for L in 200 300 400 500 600 700 800 900 1000
do
    # Extend summits L/2 bp in each direction
    if [ ! -f DHS.${L}bp.bed ]; then
        bedtools slop -i ${SUMMITS} -g ${GENOME} -b $(expr ${L} / 2) > \
        DHS.${L}bp.bed
    fi
    # i.e. bedtools slop might produce overlapping regions
    if [ ! -f DHS.${L}bp_nov.bed ]; then
        less DHS.${L}bp.bed | \
        perl -e '
            $pc="";
            while(<>){
                chomp;
                ($c,$s,$e)=split(/\t/,$_);
                if($c ne $pc){$pc=$c;}
                elsif($s<$pe){next;}
                print("$_\n");
                $pe=$e;
            }' > DHS.${L}bp_nov.bed
    fi
done
