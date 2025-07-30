# RNA-Seq Analysis Pipeline: Phyllotis Pulmonary Adaptation

This repository is for data processing and analyses related to Bautista et al. (in review), "Elevational variation in heart mass and suppression of hypoxia-induced right ventricle hypertrophy in Andean leaf-eared mice (Phyllotis)".

Here, we provide the [phenotypic data](Data/Bautista_etal_heartmassdata_upload.xlsx), [rna-seq QC metrics](Data/rnaseq_P.vaccarum_RV_QC_metrics_040325.csv), the raw [featureCounts data](Data/Pvac_readcounts_RV_CountMM.txt), and associatated [metadata]() used in this study and outline the transcriptomic analyses.

The pipeline is as follows:
1. **Read preprocessing and alignment** (Nextflow, FastP, HiSat2)
2. **Gene expression quantification** (featureCounts)
3. **Differential expression and WGCNA** (R)
4. **GO enrichment** and **visual summaries**

---

## Workflow Overview

### Nextflow Pipeline

This pipeline automates read cleaning, alignment, QC, quantification, and reporting.
- There are two files, the workflow and a config file for job submission via slurm.
- [process_RNA-Seq.nf](process_RNA-Seq.nf)
- [process_RNA-seq_slurm.config](process_RNA-seq_slurm.config)

**Steps:**
1. **Fastp** – Adapter trimming and quality filtering
2. **HISAT2 + Picard** – Paired-end read alignment and read group tagging
3. **Qualimap** – BAM-level QC
4. **MultiQC** – Summary reports (Fastp + Qualimap)
5. **featureCounts** – Gene-level read counting

**Inputs:**
- Paired-end FASTQ files (with `_1.fq.gz`/`_2.fq.gz` naming)
- Reference genome index (HISAT2)
- Gene annotation (GTF)

**Command Example:**
```bash
nextflow run process_RNA-Seq.nf -c process_RNA-seq_slurm.config -with-trace
```

**Key Outputs:**
- `results/fastp_cleaned_reads/`: Trimmed reads
- `results/bams/`: Sorted and indexed BAM files
- `results/qualimap/`: Alignment quality reports
- `results/featurecounts/featureCounts_combined.txt`: Final count matrix

---

### R Script: DE and WGCNA Analysis

Once quantified, gene counts are processed in R for sequencing statistics to identify outliers (see: [rnaseq_qc](pman_rnaseq_QC/pman_rnaseq_QC.md)), differential expression, and co-expression network analysis.

**Main Steps:**
- Load count matrix and metadata
- Remove outliers and filter lowly expressed genes
- DE analysis with `edgeR`
- Volcano plot, PCA
- GO enrichment via `gProfiler2`
- Network analysis with `WGCNA`
- Module-trait correlation with Fulton's index and RV mass
- Visualizations: scatterplots, composite plots

R script: [phyllotis_RV_rnaseq_analysis.R](phyllotis_RV_rnaseq_analyses.R)

**Dependencies:**
```r
# Install with BiocManager or install.packages()
library(edgeR)
library(WGCNA)
library(gprofiler2)
library(ggplot2)
library(plotly)
library(ggpmisc)
```

**Key Outputs:**
- `results/DE_analysis/`: DE tables, volcano plots, PCA
- `results/WGCNA/`: modules, membership, hub genes, eigengene plots
- `results/WGCNA/GO_by_module/`: module-specific GO enrichment
- `results/WGCNA/module_trait_plots/`: module-trait regression plots
- `results/WGCNA/genes_by_module/`: per-module gene lists

---
## 📖 Citation

If using this workflow or results, we would love it if you could please cite our paper!

---

## 📬 Contact

For questions or suggestions, contact Nathanae Herrera (ndh04c at gmail.com) or open an issue on this repository.

---
