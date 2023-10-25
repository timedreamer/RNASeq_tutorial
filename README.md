# Tutorial on RNA-Seq analysis

Author: Ji Huang


## 1. Environment preparation

### Terminal

You will need a terminal to interact with the server. Many cases, you will also need a FTP client to upload/download files to the server.

On windows, I like: Windows Terminal, Putty, MobaXterm, WinSCP, and FileZilla.

[tmux](https://github.com/tmux/tmux) is very good helper.

### Text editor

You will need a text editor to create/modify files on the server.

I like: Vim, VScode and Sublime Text.

### R

Once you generate the expression count matrix, you will mostly work in R for the RNA-Seq analysis.

RStudio, RStudio server, Jupyter Notebook.

## 2. Pipeline

The Basic pipeline that works with most of the RNA-seq data without UMI.

![](content/2023-10-20-16-37-23.png)


## 3. Slurm

Some slurm command that I use:

```shell
1. srun --mem 4GB --cpus-per-task 1 -t02:00:00 --pty /bin/bash
2. sacct -j <jobid> --format=JobID,JobName,state,exitcode,derivedexitcode,MaxRSS,Elapsed,MaxVMSize,MaxVMSizeNode,ReqMem
3. scontrol show jobid -dd <jobid> # for job details
4. seff <jobid> # efficiency of resource usage by the completed job
```
[Convenient SLURM Commands](https://docs.rc.fas.harvard.edu/kb/convenient-slurm-commands/)

A simple slurm job header:

```shell
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1  # modify this
#SBATCH --time=1:00:00     # modify this
#SBATCH --mem=2GB          # modify this
#SBATCH --job-name=test    # modify this
#SBATCH --mail-type=NONE
#SBATCH --mail-user=whoever@nyu.edu
#SBATCH --output=slurm_%j.out

module purge

```

## 4. Reads cleaning, aligning and counting

All for **single-end reads**

```shell
# cleaning
fastp -l 20 --thread 1 -y -t 1 -x -a AGATCGGAAGAGC -f 2 \
    -i input.fq.gz -o output2.fq.gz;

# aligning
STAR --genomeDir $STARREF --readFilesCommand zcat \
    --runThreadN 8 --readFilesIn output2.fq.gz \
    --outFilterType BySJout --outFilterMultimapNmax 20 \
    --outSAMattributes NH HI NM MD \
    --outSAMtype BAM SortedByCoordinate

# counting
featureCounts -s 1 -T 2 -a $GTF -o final.feacureCounts output2.bam
```

## 5. Files

1. `fastq` format. [Wiki](https://en.wikipedia.org/wiki/FASTQ_format)

> A FASTQ file has four line-separated fields per sequence:
Field 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
Field 2 is the raw sequence letters.
Field 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
Field 4 encodes the quality values for the sequence in Field 2, and must contain the same number of symbols as letters in the sequence.

>@SRR8699958.1 1/1
CGGGACTATACATTTACAACAAAAAGAAACAAATCTTGTGGTCAAAGTTTCCATACGTAGCTTCTCTTCTCTACAC
>+
AAA/A/EAE6AAAE/EEEA/EEEEE6AEEEEAE6/E/////EEEEE6E/E/EA/EA//E/6/6EEEEE/E<EAE/6



## 6. Hands-on

```shell
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117857
## Download fastq files
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR869/000/SRR8699970/SRR8699970.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR869/008/SRR8699958/SRR8699958.fastq.gz

## Change name
mv SRR8699970.fastq.gz DIV1_rep1.fastq.gz
mv SRR8699958.fastq.gz EV_rep1.fastq.gz

## Only use 10000 reads
zcat DIV1_rep1.fastq.gz |head -n 40000 > DIV1_p1.fq.gz
zcat EV_rep1.fastq.gz |head -n 40000 > EV_p1.fq.gz

## Clean reads
module load fastp/intel/0.20.1
fastp -l 20 --thread 1 -y -t 1 -x -a AGATCGGAAGAGC -f 2 -i DIV1_p1.fq.gz -o DIV1_p1_clean.fq.gz;
fastp -l 20 --thread 1 -y -t 1 -x -a AGATCGGAAGAGC -f 2 -i EV_p1.fq.gz -o EV_p1_clean.fq.gz;

## Continue with HISAT2-build
```
## 4. From count matrix to results

Once you have the count matrix, you can use R to analyze it. DESeq2, limma, edgeR are the most popular packages for Differential Expression analysis.


