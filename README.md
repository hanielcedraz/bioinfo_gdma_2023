# bioinfo_gdma_2023


# Introduction to bioinformatics

Open the file <a href="https://github.com/hanielcedraz/bioinfo_gdma_2023/blob/main/rnaSeqAnalysis.ipynb" target="_blank">rnaSeqAnalysis.ipynb</a> and click on "Open In Colab" button to follow the tutorial in Google Colab.



<!---
```
https://github.com/hanielcedraz/RNA-Seq_Course/blob/main/rnaSeqAnalysis.ipynb">rnaSeqAnalysis.ipynb
```
-->



## Downloading fastq files
Download the fastq files in <a href="https://github.com/hanielcedraz/bioinfo_gdma_2023/tree/main/00-Fastq" target="_blank">00-Fastq</a> 


## Loading conda env

#### Download latest version of Miniconda installer
access https://docs.conda.io/en/latest/miniconda.html and choose Miniconda3 Linux 64-bit

```shell
wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.11.0-Linux-x86_64.sh
```

#### Install miniconda
```
bash Miniconda3-py37_4.11.0-Linux-x86_64.sh
```

#### Activate conda base env
```
conda activate base
```

#### Update conda to the latest version if needed
```
conda update -n base -c defaults conda
```

#### Download the yaml file that contains the conda env: <a href="https://raw.githubusercontent.com/hanielcedraz/bioinfo_gdma_2023/main/curso_RNA-Seq.yaml" target="_blank">curso_RNA-Seq.yaml</a>  

#### Create curso_RNA-Seq env from yaml
```
conda env create -f curso_RNA-Seq.yaml

```

