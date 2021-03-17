#!/usr/bin/env bash

# Get UniBind TFBSs for prediction model PWM and peak-caller MACS
UNIBIND="pwm_tfbs.tar.gz"
if [ ! -f ${UNIBIND} ]; then
    wget https://unibind.uio.no/static/data/bulk/${UNIBIND}
    tar xvfz ${UNIBIND}
fi

# Collapse all TFBSs into a single file
ALLTFBS="tfbs.bed"
if [ ! -f ${ALLTFBS} ]; then
    ls *.pwm.bed | \
    perl -e '
        while($bed=<>){
            chomp $bed;
            open($fh,$bed);
            while($line=<$fh>){
                chomp $line;
                @l=split(/\t/,$line);
                print($l[0]."\t".$l[1]."\t".$l[2]."\t".$bed."\n");
            }
            close($f);
        }' | \
        LC_ALL=C sort --parallel=8 -T ./ -k1,1 -k2,2n > ${ALLTFBS}
fi

# Extract ENCODE assay accessions (ENCSRs) and sample names
ENC2SAM="encsr2sample.txt"
if [ ! -f ${ENC2SAM} ]; then
    ls ENCSR*.pwm.bed | cut -d "." -f 1,2 | sort | uniq > ${ENC2SAM}
fi

# Extract transcription factors
TF="tfs.txt"
if [ ! -f ${TF} ]; then
    ls *.pwm.bed | cut -d "." -f 3 | sort | uniq > ${TF}
fi