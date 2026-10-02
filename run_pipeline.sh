# Miniconda (if not already installed)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc

# Accept Terms of Service for the default channels
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Core tools environment
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n chipseq -c bioconda -c conda-forge fastqc trim-galore bowtie2 samtools macs2 seqtk -y
conda activate chipseq

# Data download
mkdir -p ~/chipseq_exercise/raw_fastq
mv /path/to/ENCFF001HTO.fastq.gz ~/chipseq_exercise/raw_fastq/ctcf.fastq.gz
mv /path/to/ENCFF001HTT.fastq.gz ~/chipseq_exercise/raw_fastq/control.fastq.gz

# Quality control (raw data)
mkdir -p ~/chipseq_exercise/qc_raw
fastqc raw_fastq/ctcf.fastq.gz raw_fastq/control.fastq.gz -o qc_raw/

# Genome indexing (chr21, not the full genome)
mkdir -p ~/chipseq_exercise/genome
cd ~/chipseq_exercise/genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz
gunzip chr21.fa.gz
bowtie2-build chr21.fa chr21_index

# Alignment
mkdir -p ~/chipseq_exercise/aligned
bowtie2 -x genome/chr21_index -U raw_fastq/ctcf.fastq.gz -S aligned/ctcf.sam --no-unal -p 4
bowtie2 -x genome/chr21_index -U raw_fastq/control.fastq.gz -S aligned/control.sam --no-unal -p 4

# FIltering and duplication
cd ~/chipseq_exercise/aligned

# MAPQ >= 30 (uniquely/high-confidence aligned reads only), then sort
samtools view -bq 30 ctcf.sam | samtools sort -o ctcf.sorted.bam -
samtools index ctcf.sorted.bam
samtools view -bq 30 control.sam | samtools sort -o control.sorted.bam -
samtools index control.sorted.bam

# PCR duplicate removal (after alignment, not before — requires genomic position)
samtools markdup -r ctcf.sorted.bam ctcf.dedup.bam
samtools markdup -r control.sorted.bam control.dedup.bam
samtools index ctcf.dedup.bam
samtools index control.dedup.bam

# Check
samtools flagstat ctcf.dedup.bam
samtools flagstat control.dedup.bam

# Peak calling (MACS2)
mkdir -p ~/chipseq_exercise/peaks
macs2 callpeak \
  -t aligned/ctcf.dedup.bam \
  -c aligned/control.dedup.bam \
  -f BAM -g hs -n ctcf_chr21 \
  --outdir peaks -q 0.05

# Output and results
wc -l peaks/ctcf_chr21_peaks.narrowPeak
head -5 peaks/ctcf_chr21_peaks.narrowPeak
