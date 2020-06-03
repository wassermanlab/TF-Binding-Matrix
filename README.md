A sparse 3D matrix of 2,503,732 bound and open regions across 175 transcription factors and 70 cells/tissues

```
git lfs clone https://github.com/wassermanlab/TF-Binding-Matrix.git
```

## Requirements
* [GNU Core Utilities](https://www.gnu.org/software/coreutils/), including [Wget](https://www.gnu.org/software/wget/)
* [Python ≥3.7](https://www.python.org) with the following libraries: [Biopython](https://biopython.org), [pandas](https://pandas.pydata.org/), [pybedtools](https://daler.github.io/pybedtools/), [pyliftover](https://github.com/konstantint/pyliftover), [sparse](https://sparse.pydata.org/en/stable/)

## Configuration

```
conda create -n matrix -c bioconda -c conda-forge python=3.7 biopython coreutils git-lfs \
    pandas pybedtools pyliftover sparse wget
```
