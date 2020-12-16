# TF binding matrix
A sparse 3D matrix of ~~1,817,918~~ 2,503,732 bound and open regions across ~~163~~ 175 transcription factors and ~~52~~ 70 cell and tissue types
![alt text](https://github.com/wassermanlab/TF-Binding-Matrix/blob/master/matrix.png?raw=true)

## News
01/09/2020 We have expanded the matrix using recent data from [ENCODE](https://www.encodeproject.org/files/ENCFF503GCK/)

## Content
* The `data` folder contains scripts to download all the data necessary to build the matrices
* The `lib` folder contains global functions to be used by all Python scripts
* The `matrix` folder contains the Python scripts to build the matrices
* The file [`environment.yml`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/environment.yml) contains the conda environment used to build the matrices (see dependencies)

## Dependencies
* [GNU core utilities](https://www.gnu.org/software/coreutils/) with [Wget](https://www.gnu.org/software/wget/)
* [Python 3.7](https://www.python.org/download/releases/3.7/) with the following libraries: [Biopython](http://biopython.org) (<1.74), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [pybedtools](https://daler.github.io/pybedtools/), and [sparse](https://sparse.pydata.org/en/stable/) 

All dependencies can be easily installed through the [conda](https://docs.conda.io/en/latest/) package manager:
```
conda create -n TfBindingMatrix -c bioconda -c conda-forge python=3.7 biopython \
    coreutils numpy pandas pybedtools sparse wget
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
### 2. Matrices
Build two TF binding matrices (i.e. data structures containing information about TF binding events, not motif models), one more sparse and the other less sparse. The matrices aggregate binding data, both from ChIP-seq experiments and TFBS predictions, of 163 TFs to 1,817,918 accessible genomic regions (*i.e.* DHSs) in 52 cell and tissue types. The matrices are saved as 2D numpy arrays, with rows and columns being individual TFs and DHS regions, respectively.
```
cd ./matrix/UCSC/
./matrix.py --dhs-file ../../data/DHS/UCSC/DHS.200bp.bed \
            --encode-dir ../../data/ENCODE/hg38/ \
            --fasta-file ../../data/Genomes/hg38/hg38.fa \
            --remap-dir ../../data/ReMap/ \
            --unibind-dir ../../data/UniBind/
```
The final matrices can be found under the `200bp` folder with the names [matrix2d.ReMap+UniBind.sparse.npz](https://github.com/wassermanlab/TF-Binding-Matrix/blob/master/matrix/UCSC/200bp/matrix2d.ReMap%2BUniBind.sparse.npz) and [matrix2d.ReMap+UniBind.less-sparse.npz](https://github.com/wassermanlab/TF-Binding-Matrix/blob/master/matrix/UCSC/200bp/matrix2d.ReMap%2BUniBind.less-sparse.npz).