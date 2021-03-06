#!/usr/bin/env bash

if [ ! -f hg38.fa ]; then
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
    tar xvfz hg38.chromFa.tar.gz 
    find chroms -type f | grep -v "_" | perl -e 'while(<>){chomp;system("mv $_ .");}'
    rm -r chroms/ hg38.chromFa.tar.gz
    cat *.fa > hg38.fa
fi