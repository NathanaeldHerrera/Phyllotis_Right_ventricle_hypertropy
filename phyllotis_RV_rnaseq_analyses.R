#==========================================================================
# RNA-Seq Differential Expression and WGCNA Analysis: Phyllotis Pulmonary Adaptation  
#==========================================================================

#=================================
# Load Packages and a Clean Start 
#=================================
rm(list = ls())  # Clear all variables from the environment

library(tidyr)
library(dplyr)
library(ggplot2)
library(gprofiler2)
library(edgeR)
library(WGCNA)
library(plotly)
library(ggpmisc)

#-----------------------------#
# Set the current run date for output file naming
run_date <- format(Sys.Date(), "%d%b%Y")

# ==============================
# Import Data and preprocessing
# ==============================

# Load raw gene count matrix and remove metadata columns (2:6)
counts0 <- read.csv("./data/Pvac_readcounts_RV_CountMM.csv", header = TRUE)
counts <- counts0[, -c(2:6)]
rownames(counts) <- counts$Geneid  # Set row names to gene IDs
counts$Geneid <- NULL              # Drop redundant column

# Load sample metadata
id <- read.csv("./data/vaccarum_RNAseq_RV_metadata.csv", header = TRUE)

# Remove technical sequence outlier(s) from metadata and count matrix
excluded <- "GD2096"
id <- id[!id$code %in% excluded, ]
counts <- counts[, colnames(counts) %in% id$code]

# Keep only shared samples between metadata and count matrix
shared_samples <- intersect(id$code, colnames(counts))
id <- id[id$code %in% shared_samples, ]
counts <- counts[, shared_samples]

# Ensure column order in counts matches sample order in metadata
id <- id[match(colnames(counts), id$code), ]
stopifnot(all(colnames(counts) == id$code))  # Check for consistency

# Filter genes with low average expression (mean < 20 across samples)
counts$mean <- rowMeans(counts)
keep_counts <- subset(counts, mean >= 20)
keep_counts$mean <- NULL

# =========================================
# Differential Expression Analysis (edgeR)
# =========================================

# Use locality (population) as grouping factor (due to collinearity elevation and population)
meta <- id
meta$locality <- factor(meta$locality)

# Build design matrix for GLM
design_pop <- model.matrix(~ locality, data = meta)

# Create DGEList object and run normalization and dispersion estimation
y <- DGEList(counts = keep_counts)
y <- calcNormFactors(y)
y <- estimateDisp(y, design_pop)
fit <- glmQLFit(y, design_pop)

# Test for differences between localities (Ojos vs reference)
print(colnames(design_pop))  # Confirm design
qlf_pop <- glmQLFTest(fit, coef = "localityOjos")

# Extract DE genes (FDR <= 0.1)
top_pop <- topTags(qlf_pop, n = Inf)
top_pop_df <- top_pop$table %>% filter(FDR <= 0.1)

# Save DE results
write.csv(top_pop_df,
          file = paste0("results/DE_analysis/sig_DE_population_fdr_logFC_n10_", run_date, ".csv"))

# ==========================
# GO Enrichment (gProfiler)
# ==========================
sig_genes <- rownames(top_pop_df)  # Significant DE genes
all_genes <- rownames(y)           # Background gene list

# Run GO enrichment on DE genes using g:Profiler
if (length(sig_genes) > 0) {
  go_res <- gost(sig_genes, organism = "mmusculus", custom_bg = all_genes,
                 ordered_query = F, exclude_iea = F,
                 correction_method = "g_SCS", significant = TRUE, sources = c("GO","KEGG","REAC","TF","HP"))
  if (!is.null(go_res$result)) {
    go_result_df <- as.data.frame(lapply(go_res$result, function(col) {
      if (is.list(col)) sapply(col, toString) else col
    }))
    write.csv(go_result_df,
              file = paste0("results/DE_analysis/GO_analysis/population_DE_GO_", run_date, ".csv"),
              row.names = FALSE)
  } else {
    message("No significant GO terms returned by gProfiler.")
  }
} else {
  message("No significant DE genes for GO enrichment.")
}

#===============
# Volcano Plot
#===============
# For a two-group DE contrast, the logFC > 0 implies higher expression in Ojos:
# logFC > 0 → Ojos (↑ in Ojos)
# logFC < 0 → Llullaillaco (↑ in Llullaillaco)

# Create a data frame with volcano plot information
# Classify genes as significantly differentially expressed (sig) and by direction of change
volcano_df <- topTags(qlf_pop, n = Inf)$table %>%
  mutate(
    sig = ifelse(FDR <= 0.1 & abs(logFC) >= 1, "yes", "no"),  # flag DE genes
    direction = case_when(  # categorize direction of expression change
      FDR <= 0.1 & logFC > 0  ~ "up",
      FDR <= 0.1 & logFC < 0 ~ "down",
      TRUE ~ "ns"  # not significant
    )
  )

# Count up- and down-regulated genes for labeling
n_up <- sum(volcano_df$direction == "up")
n_down <- sum(volcano_df$direction == "down")

# Build interactive volcano plot
p <- ggplot(volcano_df, aes(x = logFC, y = -log10(FDR), color = sig, text = rownames(volcano_df))) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("no" = "grey70", "yes" = "black")) +
  theme_minimal() +
  labs(
    title = paste0("Volcano Plot: population Effect\nUp: ", n_up, " | Down: ", n_down),
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  )

ggplotly(p)  # render interactive version

# Save static version of volcano plot to PDF
ggsave(
  filename = paste0("results/DE_analysis/volcano_population_DE_n10_", run_date, ".pdf"),
  plot = p, width = 7, height = 6
)

#==================
# PCA of DE Genes
#==================
# Subset expression matrix to significant DE genes
# Then log2-transform using CPM (counts per million)
de_expr <- keep_counts[rownames(keep_counts) %in% rownames(top_pop_df), ]
de_expr_log <- cpm(DGEList(counts = de_expr), log = TRUE, prior.count = 2)

# Perform PCA on transposed matrix (samples as rows)
pca <- prcomp(t(de_expr_log), scale. = TRUE)
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)  # percent variance explained

# Merge PCA coordinates with sample metadata
pca_df <- as.data.frame(pca$x)
pca_df$code <- rownames(pca_df)
pca_df <- left_join(pca_df, id, by = "code")
pca_df$locality <- as.factor(pca_df$locality)

# Plot PCA
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, shape = locality)) +
  geom_point(size = 3, fill = "black", color = "black") +
  scale_shape_manual(values = 21:(20 + length(unique(pca_df$locality)))) +
  theme_bw() +
  labs(
    title = expression(PCA~of~DE~Genes~(FDR <= 0.1)),
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)")
  )

print(pca_plot)

# Save PCA plot
ggsave(
  filename = paste0("results/DE_analysis/pca_population_DE_n10_", run_date, ".pdf"),
  plot = pca_plot, width = 7, height = 6
)
#===================================================
# WGCNA Network Construction and Trait Correlation
#===================================================

# Normalize raw counts to log2 CPM, then transpose to genes-as-columns format
rv.norm <- cpm(y, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)
Expr <- as.data.frame(t(rv.norm))

# Check for genes and samples with too many missing values or zero variance
gsg <- goodSamplesGenes(Expr, verbose = 0)
if (!gsg$allOK) {
  Expr <- Expr[gsg$goodSamples, gsg$goodGenes]
}

# Sample clustering to detect outliers
sampleTree <- hclust(dist(Expr), method = "average")
treecut <- cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
sampleConnectivity <- rowSums(cor(t(Expr))^2)
zconnect <- scale(sampleConnectivity)
outliers <- names(zconnect)[which(zconnect < -2)]  # flag low-connectivity samples

# Save clustering dendrogram to PDF
pdf(file = paste("./results/WGCNA/sampleClustering_rv_n10_", run_date, ".pdf", sep = ""), width = 12, height = 9)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
dev.off()

# Choose soft thresholding power for network construction
powers <- c(1:12, seq(14, 30, 2))
sft <- pickSoftThreshold(Expr, powerVector = powers, networkType = "signed", verbose = 0)

# Save scale-free topology and connectivity plots
pdf(paste("results/WGCNA/rv_beta_signed_plot_n10_", run_date, ".pdf", sep = ""), h = 7, w = 7)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, col = "red")
abline(h = 0.9, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

# Build co-expression modules using WGCNA
Net <- blockwiseModules(
  Expr, power = 25, maxBlockSize = 12000, TOMType = "signed",
  networkType = "signed", minModuleSize = 30, reassignThreshold = 0,
  mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE, saveTOMFileBase = "ExprTOM", verbose = 0,
  randomSeed = 46693
)

# Extract module information
moduleLabels <- Net$colors
moduleColors <- labels2colors(Net$colors)
MEs <- Net$MEs
geneTree <- Net$dendrograms[[1]]

# Save number of genes per module
write.csv(table(moduleLabels), file = paste("results/WGCNA/rv_modules_signed_n10_", run_date, ".csv", sep = ""))

# Identify hub genes in each module
hub_genes <- chooseTopHubInEachModule(Expr, moduleLabels, power = 25)
write.csv(hub_genes, file = paste("results/WGCNA/rv_signed_hubGenes_n10_", run_date, ".csv", sep = ""))

# Plot dendrogram with module colors
pdf(file = paste("results/WGCNA/Rplot_moduleDendrogram_rv_signed_n10_", run_date, ".pdf", sep = ""), h = 9, w = 9)
plotDendroAndColors(Net$dendrograms[[1]], moduleColors[Net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#=================================================
# Module-Trait Correlation and Module Membership
#=================================================
# Prepare trait data and align it with expression samples
# We are looking at two traits, relative right ventricle mass (RV_standardized) and Fulton's index (fi)
traitData <- id
Samples <- rownames(Expr)
traitRows <- match(Samples, traitData$code)
Traits <- traitData[traitRows, c("fi", "RV_standardized")]

# Calculate module eigengenes and number of genes/samples
nGenes <- ncol(Expr)
nSamples <- nrow(Expr)
MEs <- orderMEs(moduleEigengenes(Expr, moduleLabels)$eigengenes)

# Compute correlation and p-values between module eigengenes and traits
moduleTraitCor <- cor(MEs, Traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Create correlation table with raw and FDR-adjusted p-values
corr.table <- as.data.frame(moduleTraitCor)
colnames(corr.table) <- paste0("cor_", colnames(Traits))
corr.table$p_fi <- moduleTraitPvalue[, "fi"]
corr.table$p_rv <- moduleTraitPvalue[, "RV_standardized"]
corr.table$FDR_fi <- p.adjust(corr.table$p_fi, method = "fdr")
corr.table$FDR_rv <- p.adjust(corr.table$p_rv, method = "fdr")

# Save correlation results
write.csv(corr.table, file = paste0("results/WGCNA/Pvac.rv.corrTable_n10_fi_rvmass_", run_date, ".csv"))

# Compute module membership (kME) and p-values for each gene
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(Expr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM.", modNames)
names(MMPvalue) <- paste0("p.MM.", modNames)

# Build gene info table combining module label, membership, and p-values
genes <- names(Expr)
geneInfo0 <- data.frame(Gene = genes, moduleLabel = moduleLabels)
geneInfo0 <- data.frame(geneInfo0, geneModuleMembership, MMPvalue)
geneInfo <- geneInfo0[order(geneInfo0$moduleLabel), ]  # sort by module

# Write gene-level module membership output
write.csv(geneInfo, file = paste0("results/WGCNA/rv_signed_geneModuleMembership_n10_", run_date, ".csv"), row.names = FALSE)

# Generate and save heatmap of module-trait correlations
pdf(paste0("results/WGCNA/pvac_rv_module_trait_relationships_plot_n10_fi_rvmass_", run_date, ".pdf"),
    width = 8.5, height = 11)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(10, 10, 4, 4))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(Traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.2,
               cex.lab.x = 0.8,
               cex.lab.y = 0.4,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

#==========================================================
# ANOVA-based Trait Association and Module-Trait Plotting
#==========================================================

# Remove grey (ME0) module before downstream analysis
ME <- removeGreyME(MEs, greyMEName = paste(moduleColor.getMEprefix(), "0", sep = ""))
ME$code <- rownames(ME)  # add sample code column
ME_info <- merge(ME, id, by = "code")  # join with metadata

# Subset to only module eigengenes (We'll drop sample metadata and traits as well)
ME_info2 <- ME_info[, !(names(ME_info) %in% c("code", "locality", "elevation", "RV_standardized", "fi"))]
modNames <- names(ME_info2)  # list of module names

# Initialize matrix to store raw ANOVA p-values (both traits)
raw_pvals <- matrix(NA, nrow = length(modNames), ncol = 2,
                    dimnames = list(modNames, c("fi", "RV_standardized")))

# Loop through each module and perform ANOVA with each trait
for (i in seq_along(modNames)) {
  rME <- rank(ME_info2[[i]])  # rank-transform for robustness
  aov_fi <- summary(aov(rME ~ ME_info$fi))
  aov_rv <- summary(aov(rME ~ ME_info$RV_standardized))
  raw_pvals[i, "fi"] <- aov_fi[[1]][["Pr(>F)"]][1]  # extract p-values
  raw_pvals[i, "RV_standardized"] <- aov_rv[[1]][["Pr(>F)"]][1]
}

# Adjust p-values using FDR
adjusted_pvals <- apply(raw_pvals, 2, function(p) p.adjust(p, method = "fdr"))

# Format raw and adjusted p-values into data frames
raw_df <- round(as.data.frame(raw_pvals), 4) %>%
  tibble::rownames_to_column("module")

fdr_df <- round(as.data.frame(adjusted_pvals), 4) %>%
  tibble::rownames_to_column("module")

# Merge raw and FDR-adjusted tables by module
anova_df <- raw_df %>%
  rename_with(~ paste0("raw_", .), -module) %>%
  inner_join(fdr_df %>% rename_with(~ paste0("FDR_", .), -module), by = "module")

# Sort by significance for fi
anova_df <- anova_df %>%
  arrange(FDR_fi, module)

# Write table of ANOVA results to disk
write.csv(anova_df, paste0("results/WGCNA/table_rv_wgcna_signed_anova_merged_sorted_n10_", run_date, ".csv"), row.names = FALSE)

#==============================================================================
# Plotting: Module Eigengene vs. Trait with Linear Fit
#==============================================================================

# Set raw p-value threshold for significance
raw_pval_cutoff <- 0.05

# Create output folder for plots
plot_dir <- "results/WGCNA/module_trait_plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Add module name column to correlation table
corr.table$module <- rownames(corr.table)

# Define user-friendly! axis labels for traits
x_axis_labels <- c(
  fi = "Fulton's Index",
  RV_standardized = "Relative right ventricle mass (g)"
)

# Loop over each trait
for (trait in c("fi", "RV_standardized")) {
  
  # Get relevant column name for significance filtering
  pval_col <- if (trait == "fi") "p_fi" else "p_rv"
  
  # Filter modules significantly associated with trait
  sig_modules <- corr.table %>%
    filter(.data[[pval_col]] < raw_pval_cutoff) %>%
    pull(module)
  
  # Plot each significant module
  for (mod in sig_modules) {
    
    # Subset metadata and module expression values
    df <- ME_info %>%
      dplyr::select(locality, elevation, fi, RV_standardized, !!mod := all_of(mod))
    
    # Fetch p-value for title if needed (ended up dropping but I'll leave this here)
    raw_pval <- corr.table[mod, pval_col]
    
    # Format module name for display (e.g., ME17 → RV17)
    mod_display <- gsub("^ME", "RV", mod)
    
    # Build scatter plot with regression and R²
    p <- ggplot(df, aes_string(x = trait, y = mod)) +
      geom_point(color = "black", size = 2) +
      geom_smooth(method = "lm", color = "steelblue", se = TRUE) +
      stat_poly_eq(
        formula = y ~ x,
        aes(label = paste(..rr.label.., sep = "~`,`~")),
        parse = TRUE,
        size = 5
      ) +
      labs(
        title = paste0(mod_display,' ', x_axis_labels[[trait]]),
        x = x_axis_labels[[trait]],
        y = "Module expression"
      ) +
      theme_minimal(base_size = 14)
    
    # Save plot as PDF
    ggsave(
      filename = file.path(plot_dir, paste0(mod_display, "_vs_", trait, "_", run_date, ".pdf")),
      plot = p, width = 5.5, height = 4.5
    )
  }
}

# ============================================================================
# Plotting: Ordered 5x5 Composite Plot for Fulton's Index
# ============================================================================
library(gridExtra)
library(grid)
library(dplyr)
library(stringr)
library(ggplot2)

# Threshold for significance
raw_pval_cutoff <- 0.05

# Output path
plot_dir <- "results/WGCNA/module_trait_plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Trait column name
trait <- "fi"
pval_col <- "p_fi"
cor_col <- "cor_fi"

# Add module column to correlation table
corr.table$module <- rownames(corr.table)

# Filter significant modules for FI
sig_corr <- corr.table %>%
  filter(.data[[pval_col]] < raw_pval_cutoff)

# Convert ME to RV-style labels and order
sig_corr <- sig_corr %>%
  mutate(
    mod_display = gsub("^ME", "RV", module),
    mod_num = as.numeric(str_extract(mod_display, "\\d+")),
    cor_sign = ifelse(.data[[cor_col]] >= 0, "positive", "negative")
  )

ordered_mods <- sig_corr %>%
  arrange(desc(cor_sign), mod_num)

# Create letter labels
panel_labels <- LETTERS[1:nrow(ordered_mods)]

# Container for plots
plot_list <- list()

# Build plots
for (i in seq_len(nrow(ordered_mods))) {
  mod <- ordered_mods$module[i]
  mod_display <- ordered_mods$mod_display[i]
  panel_label <- panel_labels[i]
  
  df <- ME_info %>%
    dplyr::select(locality, elevation, fi, RV_standardized, !!mod := all_of(mod))
  
  # Compute R² from linear model
  lm_fit <- lm(df[[mod]] ~ df[[trait]])
  r2 <- summary(lm_fit)$r.squared
  r2_label <- paste0("R² = ", format(round(r2, 2), nsmall = 2))
  
  p <- ggplot(df, aes_string(x = trait, y = mod)) +
    geom_point(aes(shape = locality), color = "black", size = 1.8) +
    geom_smooth(method = "lm", color = "steelblue", se = TRUE) +
    annotate("text", x = Inf, y = Inf, label = r2_label,
             hjust = 1.1, vjust = 1.3, size = 3) +
    annotate("text", x = -Inf, y = 0.7, label = paste0(panel_label, ") ", mod_display),
             hjust = -0.1, vjust = 1.1, size = 3.5, fontface = "bold") +
    coord_cartesian(ylim = c(-0.8, 0.8)) +
    scale_shape_manual(values = c("Ojos" = 15, "Llullaillaco" = 16)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.title = element_blank()
    )
  
  plot_list[[i]] <- p
}

# Save composite plot with shared axes
composite_path <- file.path(plot_dir, paste0("FI_composite_ordered_", run_date, ".pdf"))
pdf(composite_path, width = 11, height = 8)

grid.arrange(
  grobs = plot_list,
  ncol = 5,
  top = NULL,
  left = textGrob("Module Expression (Eigengene)", rot = 90, gp = gpar(fontsize = 12)),
  bottom = textGrob("Fulton's Index", gp = gpar(fontsize = 12))
)

dev.off()

 #================================================================
# GO enrichment for significant WGCNA modules (from corr.table)
#================================================================

# Identify significant modules (raw p < 0.05 for either trait)
# we have so many modules we do not have any significant modules after FDR correction...
corr.table$module <- rownames(corr.table)
sig_modules <- corr.table %>%
  filter(p_fi < 0.05 | p_rv < 0.05) %>%
  pull(module)

# Extract genes for each module
module_genes <- lapply(sig_modules, function(mod) {
  names(moduleLabels)[moduleLabels == gsub("ME", "", mod)]
})
names(module_genes) <- sig_modules

# Run GO analysis per module
all_go_results <- list()

for (mod in sig_modules) {
  genes <- module_genes[[mod]]
  if (length(genes) >= 10) {
    go <- gost(genes,
               organism = "mmusculus",
               custom_bg = all_genes,
               correction_method = "fdr",
               significant = TRUE,
               sources = c("GO", "KEGG"))
    if (!is.null(go$result)) {
      filtered <- go$result %>% filter(p_value < 0.5)
      out_df <- as.data.frame(lapply(filtered, function(col) {
        if (is.list(col)) sapply(col, toString) else col
      }))
      out_df$module <- mod
      all_go_results[[mod]] <- out_df
      write.csv(out_df,
                file = paste0("results/WGCNA/GO_by_module/", mod, "_GO.csv"),
                row.names = FALSE)
    } else {
      message(paste("No significant GO terms for", mod))
    }
  } else {
    message(paste("Too few genes in module", mod, "- skipping GO"))
  }
}

# Save combined GO results table
if (length(all_go_results) > 0) {
  combined_go <- bind_rows(all_go_results)
  
  # Replace query values (e.g., "ME17") with "RV17"
  combined_go$query <- gsub("^ME", "RV", combined_go$module)
  
  write.csv(combined_go,
            file = paste0("results/WGCNA/GO_by_module/combined_GO_results_n10_", run_date, ".csv"),
            row.names = FALSE)
} else {
  message("No GO results to combine.")
}

#======================================================
# Save Combined Gene List for All Significant Modules
#======================================================

# Initialize list to store results for each module
all_module_genes <- list()

# Loop through each significant module
for (mod in sig_modules) {
  
  # Convert ME → RV (RV for right ventricle. Just feels right)
  mod_display <- gsub("^ME", "RV", mod)
  
  # Get genes for the current module
  genes_in_mod <- names(moduleLabels)[moduleLabels == gsub("ME", "", mod)]
  
  # Store in list
  all_module_genes[[mod_display]] <- data.frame(
    Module = mod_display,
    Gene = genes_in_mod,
    stringsAsFactors = FALSE
  )
}

# Combine all into a single data frame
combined_gene_df <- do.call(rbind, all_module_genes)

# Save to CSV
write.csv(combined_gene_df,
          file = paste0("results/WGCNA/genes_by_module/combined_RV_module_genes_", run_date, ".csv"),
          row.names = FALSE)

#=======================================================================
# Intersect DE outlier list with Velotta et al. 2018 Evolution DE list
#=======================================================================

# Velotta DE gene list intersection:
DE_gene_list <- read.csv("Velotta_2018_Evolution_DE_GeneList.csv", header = TRUE)

# Intersection of mus Ensembl IDs from Velotta DE list
intersect_genes <- unique(DE_gene_list$mus_ortholog)  # Adjust column name if needed

# Step 2: Add gene ID as column to top_pop_df
top_pop_df$Ensembl_ID <- rownames(top_pop_df)

# Step 3: Filter for shared genes
shared_DE_genes <- top_pop_df %>%
  filter(Ensembl_ID %in% intersect_genes)

# Step 4: Save the intersecting genes
write.csv(shared_DE_genes,
          file = paste0("./results/DE_analysis/shared_DE_genes_with_velotta_2018_", run_date, ".csv"),
          row.names = FALSE)

# =======================================
# GO Enrichment of Velotta gene overlap
# ======================================

sig_genes <- rownames(shared_DE_genes)  # Significant DE genes
all_genes <- rownames(y)           # Background gene list

# Run GO enrichment on DE genes using g:Profiler
if (length(sig_genes) > 0) {
  go_res <- gost(sig_genes, organism = "mmusculus", custom_bg = all_genes,
                 ordered_query = F, exclude_iea = F,
                 correction_method = "fdr", significant = TRUE, sources = c("GO","KEGG"))
  if (!is.null(go_res$result)) {
    go_result_df <- as.data.frame(lapply(go_res$result, function(col) {
      if (is.list(col)) sapply(col, toString) else col
    }))
    write.csv(go_result_df,
              file = paste0("GO_analysis/population_DE_GO_", run_date, ".csv"),
              row.names = FALSE)
  } else {
    message("No significant GO terms returned by gProfiler.")
  }
} else {
  message("No significant DE genes for GO enrichment.")
}

### FIN! ###
