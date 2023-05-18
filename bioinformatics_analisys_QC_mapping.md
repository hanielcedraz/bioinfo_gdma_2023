Bioinformatics Analisys
================
2023-05-18

## *Setting up the environment*

#### <br>

<br>

### *Creating the conda env*

``` bash
## We will use tmux for keeping the terminal active. https://tmuxcheatsheet.com
## tmux cheat sheet: https://tmuxcheatsheet.com/
$ tmux new -s curso_bioinfo

# access https://docs.conda.io/en/latest/miniconda.html and choose Miniconda3 Linux 64-bit

$ wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.11.0-Linux-x86_64.sh

# Install miniconda
$ bash Miniconda3-py37_4.11.0-Linux-x86_64.sh

# Activate conda base env
$ conda activate base

# Update conda to the latest version if needed
$ conda update -n base -c defaults conda


# Download the yaml file that contains the conda env: curso_RNA-Seq.yaml
# Create curso_RNA-Seq env from yaml
$ conda env create -f curso_RNA-Seq.yaml


#Activating the new environment
$ conda activate curso_RNA-Seq
```

#### <br>

<br>

### *Accessing reads quality before QC*

``` bash
# Using fastqc
$ mkdir fastqc_raw/
$ fastqc 00-Fastq/* -o fastqc_raw/


# Using multiqc
$ multiqc_raw/
$ multiqc fastqc_raw/* -o multiqc_raw/
```

#### <br>

<br>

<div id="htmlwidget-ad30b8e9d15d9733d3de" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-ad30b8e9d15d9733d3de">{"x":{"filter":"none","vertical":false,"caption":"<caption>Note that these samples are different from the files provided. The files you downloaded were supseted from T16A_ACAGTG and T17N_CAGATC<\/caption>","data":[["1","2","3","4"],["T16A_ACAGTG_L001_R1_001","T16A_ACAGTG_L001_R2_001","T17N_CAGATC_L001_R1_001","T17N_CAGATC_L001_R2_001"],[12424755,12424755,18742925,18742925],[101,101,101,101]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Sample<\/th>\n      <th>Total Sequences<\/th>\n      <th>avg_sequence_length<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":20,"lengthMenu":[5,10,15,20],"columnDefs":[{"className":"dt-right","targets":[2,3]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

## *Quality Control*

### *Running Quality Control with Trimmomatic*

#### *Trimmomatic: A flexible read trimming tool for Illumina NGS data*

*Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible
trimmer for Illumina Sequence Data. Bioinformatics, btu170*

[GitHub - Trimmomatic](https://github.com/usadellab/Trimmomatic)

[Trimmomatic homepage](http://www.usadellab.org/cms/?page=trimmomatic)

<br>

<br>

*This is how a basic trimmomatic command line looks like*

``` text
trimmomatic PE input_forward.fq.gz input_reverse.fq.gz \
  output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
   output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
   ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \ 
   LEADING:3 TRAILING:3 MINLEN:36
```

### *This will perform the following:*

- *Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)*

- *Remove leading low quality or N bases (below quality 3) (LEADING:3)*

- *Remove trailing low quality or N bases (below quality 3)
  (TRAILING:3)*

- *Scan the read with a 4-base wide sliding window, cutting when the
  average \\*

*quality per base drops below 15 (SLIDINGWINDOW:4:15)*

- *Drop reads below the 36 bases long (MINLEN:36)*

*Step options:*

- *ILLUMINACLIP:fastaWithAdaptersEtc:seedMismatches:palindromeClipThreshold:simpleClipThreshold:2:keepBothReads*

  - *fastaWithAdaptersEtc: specifies the path to a fasta file containing
    all the adapters, PCR sequences etc. The naming of the various
    sequences within this file determines how they are used. See below.*

  - *seedMismatches: specifies the maximum mismatch count which will
    still allow a full match to be performed*

  - *palindromeClipThreshold: specifies how accurate the match between
    the two ‘adapter ligated’ reads must be for PE palindrome read
    alignment.*

  - *simpleClipThreshold: specifies how accurate the match between any
    adapter etc. sequence must be against a read.*

  - *2 is the minimum adapter length in palindrome mode.*

  - *keepBothReads can be useful when working with paired end data, you
    will keep even redunfant information but this likely makes your
    pipelines more manageable.*

- *SLIDINGWINDOW:windowSize:requiredQuality*

  - *windowSize: specifies the number of bases to average across*

  - *requiredQuality: specifies the average quality required.*

- *LEADING:quality*

  - *quality: Specifies the minimum quality required to keep a base.*

- *TRAILING:quality*

  - *quality: Specifies the minimum quality required to keep a base.*

- *MINLEN:length*

  - *length: Specifies the minimum length of reads to be kept.*

- -threads

  - *Number of processors to use*

``` bash
# Running Quality control with Trimmomatic

# Creating a folder to store results
$ mkdir 01-CleanedReads # This will be used to store the results from the QC


# Run sample T16A-500K
$ trimmomatic PE  \
00-Fastq/T16A-500K_L001_R1_001.fastq.gz  \
00-Fastq/T16A-500K_L001_R2_001.fastq.gz  \
01-CleanedReads/T16A-500K_L001_PE1.fastq.gz  \
01-CleanedReads/T16A-500K_L001_SE1.fastq.gz  \
01-CleanedReads/T16A-500K_L001_PE2.fastq.gz  \
01-CleanedReads/T16A-500K_L001_SE2.fastq.gz  \
ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads  \
LEADING:3 TRAILING:3 -threads 50

# Run sample T17N-500K
$ trimmomatic PE  \
00-Fastq/T17N-500K_L001_R1_001.fastq.gz  \
00-Fastq/T17N-500K_L001_R2_001.fastq.gz  \
01-CleanedReads/T17N-500K_L001_PE1.fastq.gz  \
01-CleanedReads/T17N-500K_L001_SE1.fastq.gz  \
01-CleanedReads/T17N-500K_L001_PE2.fastq.gz  \
01-CleanedReads/T17N-500K_L001_SE2.fastq.gz  \
ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads  \
LEADING:3 TRAILING:3 -threads 50
```

### *Accessing reads quality after QC*

``` bash

# Using fastqc
$ mkdir fastqc_cleaned/
$ fastqc 01-CleanedReads/* -o fastqc_cleaned/


# Using multiqc
$ multiqc fastqc_cleaned/* -o multiqc_cleaned/
```

<div id="htmlwidget-59990cfb86ec7a7908f4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-59990cfb86ec7a7908f4">{"x":{"filter":"none","vertical":false,"caption":"<caption>Note that these samples are different from the files provided. The files you downloaded were supseted from T16A_ACAGTG and T17N_CAGATC<\/caption>","data":[["1","2","3","4"],["T16A_ACAGTG_L001_PE1","T16A_ACAGTG_L001_PE2","T17N_CAGATC_L001_PE1","T17N_CAGATC_L001_PE2"],[12308124,12308124,18549802,18549802],["2-101","2-101","2-101","2-101"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Sample<\/th>\n      <th>total_sequences<\/th>\n      <th>sequence_length<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":20,"lengthMenu":[5,10,15,20],"columnDefs":[{"className":"dt-right","targets":2},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
<div id="htmlwidget-f61207ee51a389681a60" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f61207ee51a389681a60">{"x":{"filter":"none","vertical":false,"caption":"<caption>Note that these samples are different from the files provided. The files you downloaded were supseted from T16A_ACAGTG and T17N_CAGATC<\/caption>","data":[["1","2","3","4"],["T16A_ACAGTG_L001_R1_001","T16A_ACAGTG_L001_R2_001","T17N_CAGATC_L001_R1_001","T17N_CAGATC_L001_R2_001"],[12424755,12424755,18742925,18742925],[12308124,12308124,18549802,18549802],[101,101,101,101],[116631,116631,193123,193123]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Sample<\/th>\n      <th>Total_Seq_before_qc<\/th>\n      <th>Total_Seq_after_qc<\/th>\n      <th>Avg_seq_len<\/th>\n      <th>removed_seq<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":20,"lengthMenu":[5,10,15,20],"columnDefs":[{"className":"dt-right","targets":[2,3,4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

#### <br>

<br>

## *Mapping*

For aligning the cleaned reads agains the regerence genome you will need
an aligner.  
For RNA-Seq data you can use
[STAR](https://academic.oup.com/bioinformatics/article/29/1/15/272537).

We will need two required files which need to be downloaded. Go to
[Ensembl ftp](https://www.ensembl.org/info/data/ftp/index.html) and
search for the specie (Chicken for this case).

Download the data:

- Reference fasta file:
  [Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz](https://ftp.ensembl.org/pub/release-109/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz)

- Annotation gtf file:
  [Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf.gz](https://ftp.ensembl.org/pub/release-109/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf.gz)

``` bash
$ mkdir reference_genome

$ wget https://ftp.ensembl.org/pub/release-109/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz -O reference_genome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz

$ wget https://ftp.ensembl.org/pub/release-109/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf.gz -O reference_genome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf.gz



$ cd reference_genome/
$ gunzip Gallus_gallus*
```

<br>

<br>

### STAR

<https://doi.org/10.1093/bioinformatics/bts635>

#### Running STAR for genome indexation:

    STAR --runMode genomeGenerate \
    --runThreadN 10 \
    --genomeDir index_Folder \
    --genomeFastaFiles mappingTarget \
    --sjdbGTFfile gtfTarget \
    --sjdbOverhang 99

- *`--runMode` genomeGenerate*

  - Directs STAR to run genome indices generation job.

- *`--runThreadN`*

  - Defines the number of threads to be used for genome generation, it
    has to be set to the number of available cores on the server node

- *`--genomeDir`*

  - specifies path to the directory (henceforth called “genome
    directory” where the genome indices are stored. This directory has
    to be created (with mkdir) before STAR run and needs to have writing
    permissions. The file system needs to have at least 100GB of disk
    space available for a typical mammalian genome. It is recommended to
    remove all files from the genome directory before running the genome
    generation step. This directory path will have to be supplied at the
    mapping step to identify the reference genome.

- *`--genomeFastaFiles`*

  - One or more FASTA files with the genome reference sequences.
    Multiple reference sequences (henceforth called “chromosomes”) are
    allowed for each fasta file. You can rename the chromosomes’ names
    in the chrName.txt keeping the order of the chromosomes in the file:
    the names from this file will be used in all output alignment files
    (such as .sam). The tabs are not allowed in chromosomes’ names, and
    spaces are not recommended.

- *`--sjdbGTFfile`*

  - The path to the file with annotated transcripts in the standard GTF
    format. STAR will extract splice junctions from this file and use
    them to greatly improve accuracy of the mapping. While this is
    optional, and STAR can be run without annotations, using annotations
    is highly recommended whenever they are available.

- *`--sjdbOverhang`*

  - The length of the genomic sequence around the annotated junction to
    be used in constructing the splice junctions database. Ideally, this
    length should be equal to the ReadLength-1, where ReadLength is the
    length of the reads. For instance, for Illumina 2x100b paired-end
    reads, the ideal value is 100-1=99. In case of reads of varying
    length, the ideal value is max(ReadLength)-1. In most cases, the
    default value of 100 will work as well as the ideal value.

``` bash
# create the index folder inside of the reference_genome directory

$ mkdir reference_genome/index_STAR

# running the indexation
# It will take a while

# To find out the number of available CPUs run
$ lscpu | grep "^CPU(s):"

$ STAR --runMode genomeGenerate --runThreadN 80  \
--genomeDir reference_genome/index_STAR  \
--genomeFastaFiles reference_genome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa  \
--sjdbGTFfile reference_genome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf --sjdbOverhang 99
```

<br>

<br>

#### Running STAR for mapping

The basic parameters for mapping samples using STAR. For more, read the
[documentation](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).

    STAR --genomeDir reference_genome/index_STAR/ 
    --readFilesIn inFile_PE1.fastq.gz \
    inFile_PE2.fastq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbGTFfile annotation_file.gtf \
    --quantMode GeneCounts \
    --outFileNamePrefix 02-MappedReads/outFile_

- *`--readFilesIn`*

  - Name(s) (with path) of the files containing the sequences to be
    mapped. If using Illumina paired-end reads, the read1 and read2
    files have to be supplied. STAR can process both FASTA and FASTQ
    files. If the read files are compressed, use the
    `--readFilesCommand UncompressionCommand option`, where
    UncompressionCommand is the un-compression command that takes the
    file name as input parameter, and sends the uncompressed output to
    stdout. For example, for gzipped files (\*.gz) use
    --readFilesCommand zcat.

- *--readFilesCommand*

  - string(s): command line to execute for each of the input file. This
    command should generate FASTA or FASTQ text and send it to stdout
    For example: zcat - to uncompress .gz files, bzcat - to uncompress
    .bz2 files, etc.

- *`--outSAMtype`*

  - BAM Unsorted

    - output unsorted Aligned.out.bam file. The paired ends of an
      alignment are always adjacent, and multiple alignments of a read
      are adjacent as well. This “unsorted” file can be directly used
      with downstream software such as HTseq, without the need of name
      sorting.

  - BAM SortedByCoordinate

    - output sorted by coordinate Aligned.sortedByCoord.out.bam file,
      similar to samtools sort command.

  - BAM Unsorted SortedByCoordinate

    - output both unsorted and sorted files.

- *`--quantMode`*

  - GeneCounts

    - STAR will count number reads per gene while mapping. A read is
      counted if it overlaps (1nt or more) one and only one gene. Both
      ends of the paired-end read are checked for overlaps. The counts
      coincide with those produced by htseq-count with default
      parameters. This option requires annotations in GTF format
      (i.e. gene id tag for each exon) specified in --sjdbGTFfile at the
      genome generation step or at the mapping step provided in option.

    <br>

    STAR outputs read counts per gene into ReadsPerGene.out.tab file
    with 4 columns which correspond to different strandedness options:

    - column 1: gene ID

    - column 2: counts for unstranded RNA-seq

    - column 3: counts for the 1st read strand aligned with RNA
      (htseq-count option -s yes)

    - column 4: counts for the 2nd read strand aligned with RNA
      (htseq-count option -s reverse)

``` bash
# Running mapping
# create a folder to store the mapped and unmapped reads

$ mkdir 02-MappedReads


# Mapping sample T16A-500K
$ STAR --genomeDir reference_genome/index_STAR \
--readFilesCommand zcat \
--readFilesIn 01-CleanedReads/T16A-500K_L001_PE1.fastq.gz \
01-CleanedReads/T16A-500K_L001_PE2.fastq.gz \
--outSAMtype BAM Unsorted SortedByCoordinate \
--sjdbGTFfile reference_genome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf  \
--quantMode GeneCounts \
--outFileNamePrefix 02-MappedReads/T16A-500K_


# Mapping sample T17N-500K
$ STAR --genomeDir reference_genome/index_STAR \
--readFilesCommand zcat \
--readFilesIn 01-CleanedReads/T17N-500K_L001_PE1.fastq.gz \
01-CleanedReads/T17N-500K_L001_PE2.fastq.gz \
--outSAMtype BAM Unsorted SortedByCoordinate \
--sjdbGTFfile reference_genome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf  \
--quantMode GeneCounts \
--outFileNamePrefix 02-MappedReads/T17N-500K_
```

<br>

<br>

------------------------------------------------------------------------

#### Understanding STAR output

Accessing the log files

``` bash
$ less T16A-500K_Log.final.out

$ less T16A-500K_Log.out
```

<div id="htmlwidget-e07921792bf772a403b4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-e07921792bf772a403b4">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7"],["Query name","FLAG","Reference name (chr)","Mapping quality","Reference sequence name of the primary alignment of the mate: '*' for no mate and '=' for Same chromosome)","Sequence","Quality"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Field<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":20,"lengthMenu":[5,10,15,20],"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
<div id="htmlwidget-72458451827010178e23" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-72458451827010178e23">{"x":{"filter":"none","vertical":false,"data":[["1"],["HWI-1KL182:99:C3GGUACXX:1:1102:18257:68144"],[163],[1],[255],["="],["ATCCCATCAACATACAGTCCCTCCCCCCCTTGCGCGTGCGTAGCAGGTCGTGATGATGAAGGCTCGTAGTCTTCCCCGCTGCAATCTTCTGACCCTGAAGG"],["00BBBFFFFFFFFFFIIFIIIIIIIIIIIFIFFBFB&lt;B&lt;B&lt;BBFBFF&lt;B&lt;7B&lt;&lt;BBBBBFBBBBBB7&lt;B7&lt;0BBBFF&lt;B&lt;0&lt;BBBB&lt;BBB0&lt;BBB&lt;B0&lt;&lt;B"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Query name<\/th>\n      <th>FLAG<\/th>\n      <th>Reference name (chr)<\/th>\n      <th>Mapping quality<\/th>\n      <th>Reference sequence name of the primary alignment of the mate<\/th>\n      <th>Sequence<\/th>\n      <th>Quality<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"lengthMenu":[[10,25,50,-1],["10","25","50","All"]],"dom":"lrtip","search":{"regex":true,"caseInsensitive":false},"columnDefs":[{"className":"nowrap","targets":"_all"},{"orderable":false,"targets":0}],"scrollX":true,"autoWidth":true,"order":[],"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

<br>

[SAM Format](https://www.samformat.info/sam-format-flag)

![Sam format. Available in
<https://www.samformat.info/>.](analysis_linux/bioinfo_gdma_2023/sam_format.png)

<br>

For visualizing the bam file we will need to use samtools view

``` bash
$ samtools view T16A-500K_Aligned.sortedByCoord.out.bam  

$ samtools view T17N-500K_Aligned.sortedByCoord.out.bam
```

Accessing the counting reads file

``` text
$ less -S T16A-500K_ReadsPerGene.out.tab

$ less -S T17N-500K_ReadsPerGene.out.tab
```

<div id="htmlwidget-b38651889980e332cfb3" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-b38651889980e332cfb3">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100"],["N_unmapped","N_multimapping","N_noFeature","N_ambiguous","ENSGALG00010000711","ENSGALG00010000715","ENSGALG00010000717","ENSGALG00010000720","ENSGALG00010000724","ENSGALG00010000726","ENSGALG00010000729","ENSGALG00010000734","ENSGALG00010000757","ENSGALG00010000762","ENSGALG00010000777","ENSGALG00010000784","ENSGALG00010000785","ENSGALG00010000788","ENSGALG00010000791","ENSGALG00010000793","ENSGALG00010000795","ENSGALG00010000796","ENSGALG00010000800","ENSGALG00010000804","ENSGALG00010000810","ENSGALG00010000812","ENSGALG00010000816","ENSGALG00010000818","ENSGALG00010000820","ENSGALG00010000823","ENSGALG00010000843","ENSGALG00010000859","ENSGALG00010000864","ENSGALG00010000870","ENSGALG00010000876","ENSGALG00010000901","ENSGALG00010000902","ENSGALG00010000903","ENSGALG00010000905","ENSGALG00010000906","ENSGALG00010000907","ENSGALG00010000908","ENSGALG00010000910","ENSGALG00010000912","ENSGALG00010000913","ENSGALG00010000914","ENSGALG00010000916","ENSGALG00010000917","ENSGALG00010000918","ENSGALG00010000920","ENSGALG00010000921","ENSGALG00010000922","ENSGALG00010000923","ENSGALG00010000925","ENSGALG00010000926","ENSGALG00010000927","ENSGALG00010000928","ENSGALG00010000930","ENSGALG00010000931","ENSGALG00010000932","ENSGALG00010000933","ENSGALG00010000935","ENSGALG00010000936","ENSGALG00010000938","ENSGALG00010000939","ENSGALG00010000941","ENSGALG00010000942","ENSGALG00010000943","ENSGALG00010000944","ENSGALG00010000947","ENSGALG00010000948","ENSGALG00010000951","ENSGALG00010000953","ENSGALG00010000954","ENSGALG00010000956","ENSGALG00010000957","ENSGALG00010000958","ENSGALG00010000960","ENSGALG00010000961","ENSGALG00010000963","ENSGALG00010000964","ENSGALG00010000967","ENSGALG00010000968","ENSGALG00010000969","ENSGALG00010000970","ENSGALG00010000972","ENSGALG00010000973","ENSGALG00010000974","ENSGALG00010000975","ENSGALG00010000976","ENSGALG00010000977","ENSGALG00010000979","ENSGALG00010000981","ENSGALG00010000982","ENSGALG00010000983","ENSGALG00010000984","ENSGALG00010000986","ENSGALG00010000987","ENSGALG00010000988","ENSGALG00010000989"],[11097,2846,29593,8734,3,0,30,29,5,0,0,0,1,1,44,6,9,9,1,0,1,1,38,0,14,12,13,3,1,2,2,2,3,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[11097,2846,128501,605,2,0,15,18,5,0,0,0,0,0,27,3,5,6,0,0,1,0,24,0,6,5,2,1,1,1,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[11097,2846,127295,652,1,0,15,11,0,0,0,0,1,1,17,3,4,4,2,0,0,1,14,0,8,7,11,2,0,1,0,1,2,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>gene ID<\/th>\n      <th>unstranded<\/th>\n      <th>1st read strand<\/th>\n      <th>2nd read strand<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
