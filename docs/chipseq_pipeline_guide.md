# ChIP-seq Analysis — Practical Guide (CTCF, K562, chr21)

Complete pipeline from raw FASTQ to peak calling, run on real ENCODE data.
Dataset: `ENCSR000DWE` (CTCF ChIP-seq, K562) + `ENCSR000DWA` (input control), John Stamatoyannopoulos lab, UW.

---

## 1. Environment setup

```bash
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
```

**Problem encountered:** `CondaToSNonInteractiveError` — Anaconda requires explicit ToS acceptance for the `main`/`r` channels before creating environments.
**Solution:** run `conda tos accept` for both channels before `conda create`.

**Note:** the official Miniconda URL is `repo.anaconda.com/miniconda/` (not `miniconda3/` — that path returns 404).

---

## 2. Data download

Downloaded manually from the ENCODE portal (encodeproject.org):
- Experiment: `ENCSR000DWE` → file `ENCFF001HTO.fastq.gz` (replicate 1, SE36nt)
- Control: `ENCSR000DWA` → file `ENCFF001HTT.fastq.gz`

```bash
mkdir -p ~/chipseq_exercise/raw_fastq
mv /path/to/ENCFF001HTO.fastq.gz ~/chipseq_exercise/raw_fastq/ctcf.fastq.gz
mv /path/to/ENCFF001HTT.fastq.gz ~/chipseq_exercise/raw_fastq/control.fastq.gz
```

**Note:** SE36nt dataset (single-end, 36nt) — the two rows listed on ENCODE for a given experiment are independent isogenic biological replicates, not paired-end R1/R2. One replicate is enough for a single exercise.

---

## 3. Quality Control (raw data)

```bash
mkdir -p ~/chipseq_exercise/qc_raw
fastqc raw_fastq/ctcf.fastq.gz raw_fastq/control.fastq.gz -o qc_raw/
```

What to check (ignore the rest of the report):
- **Per-base sequence quality** — should stay in the green/yellow zone
- **Sequence Duplication Levels** — being elevated is normal for TF ChIP-seq (real enrichment)
- **Adapter Content** — if clean, trimming can be skipped

**Result obtained:** ~47% duplication, adapters and quality OK → trimming skipped, moved straight to alignment.

---

## 4. Genome indexing (chr21, not the full genome)

```bash
mkdir -p ~/chipseq_exercise/genome
cd ~/chipseq_exercise/genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz
gunzip chr21.fa.gz
bowtie2-build chr21.fa chr21_index
```

**Choice:** used only chr21 (not the full human genome) to stay within the available RAM (~3-4GB free) and time constraints. The full genome index requires ~3.2-4GB RAM and 30-60 min to build; chr21 takes seconds and negligible memory. CTCF binds sites distributed across all chromosomes, so chr21 still produces real, interpretable peaks — a legitimate simplification for a learning exercise, not for production work.

---

## 5. Alignment

```bash
mkdir -p ~/chipseq_exercise/aligned
bowtie2 -x genome/chr21_index -U raw_fastq/ctcf.fastq.gz -S aligned/ctcf.sam --no-unal -p 4
bowtie2 -x genome/chr21_index -U raw_fastq/control.fastq.gz -S aligned/control.sam --no-unal -p 4
```

**Results:**
| File | Total reads | Overall alignment rate |
|-------|-------------|------------------------|
| CTCF | 32,740,518 | 8.17% |
| Control | 18,123,856 | 16.62% |

A low alignment rate in absolute terms is expected: chr21 is only ~1.5% of the human genome, so most reads belong to other chromosomes that weren't indexed.

**Observed issue (non-blocking):** warning `Failed to launch x86-64-v3 version, staying with default` — the CPU doesn't support the optimized instruction set, so bowtie2 automatically falls back to the standard binary. No impact on the result, only a minor speed loss.

---

## 6. Filtering and deduplication

```bash
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
```

**Final results:** 417,152 unique reads (CTCF), 479,629 unique reads (control), 0 duplicates removed (plausible with short 36nt reads and MAPQ≥30: few false duplicates from position collisions).

**Conceptual note:** deduplication must happen *after* alignment, never on the raw FASTQ — a true PCR duplicate is identified by an identical genomic position, not by sequence alone.

---

## 7. Peak calling (MACS2)

```bash
mkdir -p ~/chipseq_exercise/peaks
macs2 callpeak \
  -t aligned/ctcf.dedup.bam \
  -c aligned/control.dedup.bam \
  -f BAM -g hs -n ctcf_chr21 \
  --outdir peaks -q 0.05
```

### Critical issue encountered: `ImportError: undefined symbol: __log_finite`

The `macs2` package installed via bioconda (in the original `chipseq` environment, in a dedicated clean environment `macs2env`, and even pinning Python 3.9) consistently failed with the same binary library error (a conflict between the precompiled Cython build and the system math library/glibc). **Not resolved even in a second, clean conda environment** — a sign the issue was in the bioconda package itself, not the local configuration.

**Solution that worked:** bypass conda entirely for MACS2, installing it via `pip` in a dedicated Python virtualenv (local compilation from source instead of a precompiled binary):

```bash
# Use the Python 3.9 already present in an existing conda env, not the system's Python 3.14
# (MACS2 is not guaranteed compatible with very recent Python versions)
~/miniconda3/envs/<env_with_py39>/bin/python3.9 -m venv ~/macs2_pip_env
source ~/macs2_pip_env/bin/activate
pip install macs2
```

**General lesson:** when a bioconda package throws low-level binary errors (`undefined symbol`), recreating the environment is not guaranteed to fix it — it can be a bug in the build itself. `pip install` with local compilation is a reliable fallback, especially for older C/Cython tools like MACS2.

### Output and result

```bash
wc -l peaks/ctcf_chr21_peaks.narrowPeak
head -5 peaks/ctcf_chr21_peaks.narrowPeak
```

- **572 peaks** called on chr21, q-value cutoff 0.05
- Predicted fragment length: **292 bp** (alternatives: 150, 248, 292, 542, 576 bp — dispersion plausibly due to the low depth of a single chromosome)
- Output files: `ctcf_chr21_peaks.narrowPeak` (standard BED format, loadable in IGV/UCSC), `ctcf_chr21_peaks.xls`, `ctcf_chr21_summits.bed`

---

## Summary of issues encountered and solutions

| # | Problem | Cause | Solution |
|---|----------|-------|-----------|
| 1 | `wget` 404 on Miniconda | Wrong URL path (`miniconda3/` instead of `miniconda/`) | Fixed the path in the URL |
| 2 | `CondaToSNonInteractiveError` | ToS not accepted for `main`/`r` channels | `conda tos accept` before creating the environment |
| 3 | Limited RAM (~3.4GB free) | Full human genome index requires ~4GB | Used a single chromosome (chr21) instead of the full genome |
| 4 | Moved files not found | Assumed path differed from the actual download path | `find ~ -iname "*.fastq.gz"` to locate them |
| 5 | `macs2: undefined symbol: __log_finite` | Binary build bug in the bioconda package, unrelated to local config | Installed via `pip` in a virtualenv with Python 3.9, compiled from source |

---

## Tips for next time (short notes)

- **Always QC before aligning** — unnecessary trimming costs time, but skipping it when it's needed ruins the whole downstream.
- **MAPQ ≥ 30 is the ENCODE standard for TF ChIP-seq** — don't lower it without reason.
- **Always deduplicate after alignment**, never on raw FASTQ.
- **The input control is not optional** — without it, MACS2 can't distinguish real enrichment from technical/mappability bias.
- **A single chromosome is acceptable only for exercises/debugging** — a publishable result needs the full genome.
- **If a bioconda package throws `undefined symbol` errors**, don't keep retrying in the same ecosystem for too long: pip/virtualenv is often the faster fix.
- **For real reproducibility**, this command sequence should be converted into a Snakemake or Nextflow workflow, not left as a linear bash script.
