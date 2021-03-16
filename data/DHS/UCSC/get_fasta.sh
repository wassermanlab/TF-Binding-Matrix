#!/usr/bin/env bash

bedtools getfasta -fi ../../Genomes/hg38/hg38.fa -bed DHS.1000bp.bed -fo DHS.1000bp.fa
bedtools getfasta -fi ../../Genomes/hg38/hg38.fa -bed DHS.500bp.bed -fo DHS.500bp.fa
bedtools getfasta -fi ../../Genomes/hg38/hg38.fa -bed DHS.200bp.bed -fo DHS.200bp.fa