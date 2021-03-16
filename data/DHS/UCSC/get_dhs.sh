#!/usr/bin/env bash

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# From https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=wgEncodeRegDnase
# DNase I Hypersensitivity in 95 cell types from ENCODE
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CLUSTS="wgEncodeRegDnaseClustered.txt.gz"
if [ ! -f ${CLUSTS} ]; then
    wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/${CLUSTS}
fi
INPUTS="wgEncodeRegDnaseClusteredInputs.txt.gz"
if [ ! -f ${INPUTS} ]; then
    wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/${INPUTS}
fi
SOURCS="wgEncodeRegDnaseClusteredSources.txt.gz"
if [ ! -f ${SOURCS} ]; then
    wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/${SOURCS}
fi

# Get chrom sizes (i.e. option "-g") for bedtools slop
GENOME="hg38.chrom.sizes"
if [ ! -f ${GENOME} ]; then
    wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/${GENOME}
fi

# Extract DHS centers
CENTERS="centers.bed"
if [ ! -f ${CENTERS} ]; then
    zless ${CLUSTS} | cut -f 2-4,8 | grep -v "_" | \
    awk '{print($1"\t"$2+(($3-$2)/2)"\t"$2+(($3-$2)/2)"\t"$4);}' > ${CENTERS}
fi

# L = 200-1000 bp
# Increase in steps of 100 bp
# for L in 200 300 400 500 600 700 800 900 1000
for L in 200 500 1000
do
    # Extend summits L/2 bp in each direction
    if [ ! -f DHS.${L}bp.bed ]; then
        bedtools slop -i ${CENTERS} -g ${GENOME} -b $(expr ${L} / 2) > \
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
