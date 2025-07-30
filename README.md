# RNA-Seq Analysis Pipeline: Phyllotis Pulmonary Adaptation

This project processes and analyzes RNA-Seq data from *Phyllotis vaccarum* to explore transcriptomic signatures of pulmonary adaptation at high altitude. It includes:

1. **Read preprocessing and alignment** (Nextflow, FastP, HiSat2)
2. **Gene expression quantification** (featureCounts)
3. **Differential expression and WGCNA** (R)
4. **GO enrichment** and **visual summaries**

---

## 🗂️ Workflow Overview

### 🔁 Nextflow Pipeline

This pipeline automates read cleaning, alignment, QC, quantification, and reporting.

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
nextflow run main.nf -c 
```

**Key Outputs:**
- `results/fastp_cleaned_reads/`: Trimmed reads
- `results/bams/`: Sorted and indexed BAM files
- `results/qualimap/`: Alignment quality reports
- `results/featurecounts/featureCounts_combined.txt`: Final count matrix

---

### 📊 R Script: DE and WGCNA Analysis

Once quantified, gene counts are processed in R for sequencing statistics to identify outliers, differential expression, and co-expression network analysis.

**Main Steps:**
- Load count matrix and metadata
- Remove outliers and filter lowly expressed genes
- DE analysis with `edgeR`
- Volcano plot, PCA
- GO enrichment via `gProfiler2`
- Network analysis with `WGCNA`
- Module-trait correlation with Fulton's index and RV mass
- Visualizations: scatterplots, composite plots

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
