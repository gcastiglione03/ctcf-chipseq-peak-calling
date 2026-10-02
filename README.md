# CTCF ChIP-seq Peak Calling Pipeline

End-to-end pipeline for ChIP-seq analysis, from raw read alignment to identifying binding sites of the transcription factor CTCF, applied to public ENCODE data.

## Goal

Reproduce a standard ChIP-seq analysis locally — alignment, filtering, peak calling — on a subset of real data (chromosome 21), to validate the full pipeline before applying it genome-wide.

## Dataset

- **ChIP (CTCF)**: [ENCSR000DWE](https://www.encodeproject.org/experiments/ENCSR000DWE/) — K562, John Stamatoyannopoulos lab (UW)
- **Input control**: [ENCSR000DWA](https://www.encodeproject.org/experiments/ENCSR000DWA/)
- Format: single-end, 36nt (Illumina GA2)

## Pipeline

```
Raw FASTQ
   │
   ├─ FastQC (quality control)
   │
   ▼
Bowtie2 alignment (vs chr21, hg38)
   │
   ▼
samtools view -q 30 (MAPQ filter, unique alignments only)
   │
   ▼
samtools markdup -r (PCR duplicate removal)
   │
   ▼
MACS2 callpeak (ChIP vs input control, q < 0.05)
   │
   ▼
572 peaks (narrowPeak, BED, summits)
```

## Tech stack

`bowtie2` `samtools` `MACS2` `FastQC` `conda` — environment managed with Miniconda/bioconda.

## Results

| Metric | CTCF | Control |
|---|---|---|
| Total reads | 32,740,518 | 18,123,856 |
| Unique reads post-filter (MAPQ≥30, dedup) | 417,152 | 479,629 |
| Overall alignment rate (chr21 only) | 8.17% | 16.62% |

- **572 peaks** identified on chr21 (q-value < 0.05)
- Predicted fragment length: 292 bp

## Methodological notes

To keep compute time and resources manageable, alignment was restricted to chromosome 21 (hg38) instead of the full genome. CTCF binds sites distributed across all chromosomes, so this subset still produces real, interpretable peaks — but a publishable analysis or final report would require genome-wide alignment.

The MAPQ ≥ 30 filter follows the ENCODE standard for transcription factor ChIP-seq, excluding multi-mapped or low-confidence alignments.

## Technical issues resolved during development

Worth including, since it's part of real bioinformatics work: installing MACS2 via bioconda consistently failed with a binary library error (`undefined symbol: __log_finite`), reproducible even in clean conda environments. Resolved by installing MACS2 via `pip` in a dedicated Python virtualenv (compiled from source). Full details, with all commands and intermediate attempts, in [`chipseq_pipeline_guide.md`](./chipseq_pipeline_guide.md).

## Next steps

- Extend alignment to the full genome
- Peak annotation (e.g. `ChIPseeker`)
- Convert the pipeline into a Snakemake/Nextflow workflow for reproducibility
- Visualize peaks in IGV/UCSC Genome Browser

---

*Portfolio/learning project — public ENCODE data, non-commercial use.*
