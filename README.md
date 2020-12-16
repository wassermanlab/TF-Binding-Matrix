# TF binding matrix
A sparse 3D matrix of ~~1,817,918~~ 2,503,732 bound and open regions across ~~163~~ 175 transcription factors and ~~52~~ 70 cell and tissue types

![alt text](https://github.com/wassermanlab/TF-Binding-Matrix/blob/master/matrix.png?raw=true)

## News
01/09/2020 We have expanded the matrix using recent data from [ENCODE](https://www.encodeproject.org/files/ENCFF503GCK/)

## Content
* The `examples` folder contains the sequences of two transcription factors (TFs) and one protein that is not a transcription factor, such as the human serine/threonine-protein kinase [mTOR](https://www.uniprot.org/uniprot/P42345)
* The `files` folder contains the output of the script [`get_files.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/files/get_files.py), which downloads TF sequences from [UniProt](https://www.uniprot.org/), DNA-binding domains (DBDs) from [Pfam](https://pfam.xfam.org/), and retrieves cut-offs on the DBD percentage of sequence identity from [Cis-BP](http://cisbp.ccbr.utoronto.ca/), etc.
* The `models` folder contains the similarity regression models created by calling the script [`pairwise.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/models/pairwise.py) followed by [`regression.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/models/regression.py)
* The script [`infer_profile.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/infer_profile.py) takes as input the folders `files` and `models`, plus one or more proteic sequences in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) (_e.g._ a proteome), and infers DNA-binding profiles from JASPAR 
* The file [`environment.yml`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/environment.yml) contains the conda environment used to develop the profile inference tool for JASPAR 2020 (see installation)

The original scripts used for the publication of [JASPAR 2016](https://doi.org/10.1093/nar/gkv1176) have been placed in the folder [`version-1.0`](https://github.com/wassermanlab/JASPAR-profile-inference/tree/master/version-1.0).

## Dependencies
* [GNU core utilities](https://www.gnu.org/software/coreutils/) with [Wget](https://www.gnu.org/software/wget/)
* [Python 3.7](https://www.python.org/download/releases/3.7/) with the following libraries: [Biopython](http://biopython.org) (<1.74), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [pybedtools](https://daler.github.io/pybedtools/), and [sparse](https://sparse.pydata.org/en/stable/) 

All dependencies can be easily installed through the [conda](https://docs.conda.io/en/latest/) package manager:
```
conda create -c bioconda -c conda-forge python=3.7 biopython coreutils numpy \
    pandas pybedtools sparse wget
```

## Steps
The following steps were followed to generate the TF binding matrices for the transfer learning manuscript.
### 1. Data
#### 1.1 DNase I hypersensitive sites
Download clustered DHS data in 95 cell and tissue types from the [ENCODE DHS peak clusters at the UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=wgEncodeRegDnase). Then, extract the center of each cluster and expand it 100 bp in each direction using bedtools slop for a final length of 200 bp. In addition, download information about the names of the clustered cells and tissues and their correspondance with the different cluster IDs.
```
cd ./DHS/UCSC/
./get_dhs.sh
```
#### 1.2 ENCODE accessibility and TF binding
Download all human DNase-seq and TF ChIP-seq data from [ENCODE](https://www.encodeproject.org/matrix/?type=Experiment&status=released&perturbed=false&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assembly=GRCh38&assay_title=TF+ChIP-seq&assay_title=DNase-seq). Then, resize all the DNase-seq data using bedtools slop for a final length of 150 bp and store it in a single file. Finally, store all the ChIP-seq data for each TF in a separate file.
```
cd ./ENCODE/
./get_encode.sh
```
#### 1.3 Hg38 genome sequence
Download the FASTA sequence of the build 38 of the Genome Reference Consortium human genome (*i.e.* [hg38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/)). Discard any non-standard chromosomes.
```
cd ./Genomes/hg38/
./get_hg38.sh
```
#### 1.4 ReMap TF binding
Download all human TF ChIP-seq peaks from [ReMap 2018](http://remap.univ-amu.fr/download_page#remap2018tab). Then, extract the peak summits and the sample names given to the different ENCODE experiments (*i.e.* files whose name starts with _ENCSR_).
```
cd ./ReMap/
./get_remap.sh
```
#### 1.5 UniBind TFBSs
Download all human PWM-based TFBS predictions from [UniBind](https://unibind.uio.no/downloads/). Then, collapse all TFBSs into a single file, and extract the names of the different TFs as well as the sample names given to the different ENCODE experiments (*i.e.* files whose name starts with _ENCSR_).
```
cd ./UniBind/
./get_unibind.sh 
```
