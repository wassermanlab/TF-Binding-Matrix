# Transcription factor binding matrix
A sparse 3D matrix of 2,503,732 bound and open regions across 175 transcription factors and 70 cells/tissues

```
git clone https://github.com/wassermanlab/TF-Binding-Matrix.git
```

## News
01/09/2019 We have improved the profile inference tool by implementing our own [similarity regression](https://www.nature.com/articles/s41588-019-0411-1) method.

## Content
* The `examples` folder contains the sequences of two transcription factors (TFs) and one protein that is not a transcription factor, such as the human serine/threonine-protein kinase [mTOR](https://www.uniprot.org/uniprot/P42345)
* The `files` folder contains the output of the script [`get_files.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/files/get_files.py), which downloads TF sequences from [UniProt](https://www.uniprot.org/), DNA-binding domains (DBDs) from [Pfam](https://pfam.xfam.org/), and retrieves cut-offs on the DBD percentage of sequence identity from [Cis-BP](http://cisbp.ccbr.utoronto.ca/), etc.
* The `models` folder contains the similarity regression models created by calling the script [`pairwise.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/models/pairwise.py) followed by [`regression.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/models/regression.py)
* The script [`infer_profile.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/infer_profile.py) takes as input the folders `files` and `models`, plus one or more proteic sequences in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) (_e.g._ a proteome), and infers DNA-binding profiles from JASPAR 
* The file [`environment.yml`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/environment.yml) contains the conda environment used to develop the profile inference tool for JASPAR 2020 (see installation)

The original scripts used for the publication of [JASPAR 2016](https://doi.org/10.1093/nar/gkv1176) have been placed in the folder [`version-1.0`](https://github.com/wassermanlab/JASPAR-profile-inference/tree/master/version-1.0).

## Dependencies
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [HMMER](http://hmmer.org/) (version ≥3.0)
* [Python 3.7](https://www.python.org/download/releases/3.7/) with the following libraries: [Biopython](http://biopython.org) (<1.74), [CoreAPI](http://www.coreapi.org), [glmnet](https://github.com/civisanalytics/python-glmnet), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [ProDy](http://prody.csb.pitt.edu/), [SciPy](https://www.scipy.org/), [scikit-learn](https://scikit-learn.org/stable/) and [tqdm](https://tqdm.github.io) 
* [Tomtom](http://meme-suite.org/doc/tomtom.html) as distributed in the [MEME](http://meme-suite.org/index.html) suite (version ≥5.0)

Note that for running `infer_profile.py`, neither Tomtom nor the CoreAPI, glmnet, ProDy and scikit-learn python packages are required.

## Installation
All dependencies can be installed through the [conda](https://docs.conda.io/en/latest/) package manager:
```
conda create -c bioconda -c conda-forge python=3.7 biopython coreutils pandas \
    pybedtools pyliftover sparse wget
```

## Steps
The following steps were followed to generate the TF binding matrices for the transfer learning manuscript.
### 1. Data download
#### 1.1 Hg38
#### 1.2 ENCODE
#### 1.3 UCSC
#### 1.4 UCSC
#### 1.5 UCSC
