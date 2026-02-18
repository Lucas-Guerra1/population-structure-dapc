# ==============================================================================
# DAPC - COMPLETE PIPELINE FOR POPULATION GENETICS
# Discriminant Analysis of Principal Components (DAPC)
#
# Based on:
#   - Jombart et al. (2010) BMC Genetics. doi:10.1186/1471-2156-11-94
#   - Miller et al. (2020) Molecular Ecology Resources
#
# Dependencies: vcfR, adegenet, viridis
# Input:        .vcf or .vcf.gz file
# Output:       PDFs with plots, CSV with results, TXT with summary
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. PACKAGES
# ------------------------------------------------------------------------------
if (!requireNamespace("vcfR",     quietly = TRUE)) install.packages("vcfR")
if (!requireNamespace("adegenet", quietly = TRUE)) install.packages("adegenet")
if (!requireNamespace("viridis",  quietly = TRUE)) install.packages("viridis")

library(vcfR)
library(adegenet)
library(viridis)

# ------------------------------------------------------------------------------
# 1. LOAD VCF DATA
# ------------------------------------------------------------------------------
cat("Select the VCF file in the window that will open...\n")

vcf <- tryCatch(
  read.vcfR(file.choose()),
  error = function(e) stop("Error reading VCF file: ", e$message)
)

genlight_obj <- vcfR2genlight(vcf)

cat("\n=== ORIGINAL DATA ===\n")
cat("Individuals:", nInd(genlight_obj), "\n")
cat("SNPs:       ", nLoc(genlight_obj), "\n\n")

# ------------------------------------------------------------------------------
# 2. QUALITY CONTROL (QC)
# ------------------------------------------------------------------------------
cat("=== APPLYING QC FILTERS ===\n")

# 2.1 Remove individuals with > 50% missing data
ind_missing <- rowSums(is.na(as.matrix(genlight_obj))) / nLoc(genlight_obj)
cat("Individuals removed (> 50% missing):", sum(ind_missing >= 0.5), "\n")
genlight_obj <- genlight_obj[ind_missing < 0.5, ]

# 2.2 Remove loci with > 10% missing data
loc_missing <- glNA(genlight_obj, alleleAsUnit = FALSE)
cat("SNPs removed (> 10% missing):       ", sum(loc_missing >= 0.1), "\n")
genlight_obj <- genlight_obj[, loc_missing < 0.1]

# 2.3 Remove loci with MAF < 0.05 (Minor Allele Frequency)
maf <- glMean(genlight_obj)
maf[maf > 0.5] <- 1 - maf[maf > 0.5]
cat("SNPs removed (MAF < 0.05):          ", sum(maf <= 0.05), "\n")
genlight_obj <- genlight_obj[, maf > 0.05]

cat("\n=== DATA AFTER QC ===\n")
cat("Individuals:", nInd(genlight_obj), "\n")
cat("SNPs:       ", nLoc(genlight_obj), "\n\n")

# ------------------------------------------------------------------------------
# 3. FIND OPTIMAL K VIA BIC (K-means)
# ------------------------------------------------------------------------------
cat("=== DETERMINING OPTIMAL NUMBER OF CLUSTERS (K) ===\n")

# Maximum K: square root of N or 10 (whichever is smaller)
max_k <- min(10, floor(sqrt(nInd(genlight_obj))))
cat("Testing K from 1 to", max_k, "\n")

# Number of PCs for K-means: conservative rule of N/3
n_pca_initial <- floor(nInd(genlight_obj) / 3)
cat("Using", n_pca_initial, "PCs for K-means\n")

grp <- find.clusters(
  genlight_obj,
  max.n.clust    = max_k,
  n.pca          = n_pca_initial,
  choose.n.clust = FALSE,   # non-interactive
  criterion      = "diffNgroup",
  stat           = "BIC",
  n.iter         = 1e6,
  n.start        = 10
)

k_optimal <- which.min(grp$Kstat)
cat("\nOptimal K (lowest BIC):", k_optimal, "\n")
cat("Clusters identified:", length(unique(grp$grp)), "\n\n")

# Save BIC plot
pdf("BIC_K_selection.pdf", width = 7, height = 5)
plot(
  grp$Kstat, type = "b",
  xlab = "Number of clusters (K)",
  ylab = "BIC",
  main = "BIC for K-means Clustering",
  pch = 19, col = "steelblue", lwd = 2
)
abline(v = k_optimal, col = "red", lty = 2, lwd = 2)
text(k_optimal, min(grp$Kstat),
     labels = paste("K =", k_optimal),
     pos = 4, col = "red", cex = 1.2)
dev.off()
cat("Plot saved: BIC_K_selection.pdf\n\n")

# ------------------------------------------------------------------------------
# 4. CROSS-VALIDATION — DETERMINE OPTIMAL NUMBER OF PCs (xvalDapc)
# ------------------------------------------------------------------------------
cat("=== CROSS-VALIDATION (xvalDapc) ===\n")
cat("Testing up to 300 PCs with 30 replicates...\n")

set.seed(999)  # reproducibility

# Open PDF before the function to capture the automatic plot
pdf("cross_validation_PCs.pdf", width = 8, height = 6)
xval <- xvalDapc(
  tab(genlight_obj, NA.method = "mean"),
  grp$grp,
  n.pca.max    = 300,
  training.set = 0.9,
  result       = "groupMean",
  center       = TRUE,
  scale        = FALSE,
  n.pca        = NULL,
  n.rep        = 30,
  xval.plot    = TRUE
)
dev.off()
cat("Plot saved: cross_validation_PCs.pdf\n\n")

# Safe extraction of optimal number of PCs
optimal_n_pca <- as.numeric(xval$`Number of PCs Achieving Lowest MSE`)
cat("Optimal number of PCs:", optimal_n_pca, "\n")

# Access MSA by name (string) to avoid incorrect positional indexing
msa_value <- xval$`Mean Successful Assignment by Number of PCs`[
  as.character(optimal_n_pca)
]
cat("Mean Successful Assignment (MSA):", round(msa_value, 4), "\n\n")

# ------------------------------------------------------------------------------
# 5. FINAL DAPC WITH OPTIMIZED PARAMETERS
# ------------------------------------------------------------------------------
cat("=== RUNNING FINAL DAPC ===\n")

# Number of discriminants: maximum is K - 1
n_da <- min(length(unique(grp$grp)) - 1, 3)
cat("Number of DAs used:", n_da, "\n")

dapc_result <- dapc(
  genlight_obj,
  grp$grp,
  n.pca = optimal_n_pca,
  n.da  = n_da
)

# ------------------------------------------------------------------------------
# 6. VISUALIZATIONS
# ------------------------------------------------------------------------------
cat("=== GENERATING FIGURES ===\n")

cores <- viridis(length(unique(grp$grp)))

# A. Scatter Plot (discriminant axes)
pdf("DAPC_scatter.pdf", width = 10, height = 8)
scatter(
  dapc_result,
  col      = cores,
  scree.da = TRUE,
  bg       = "white",
  pch      = 20,
  cell     = 0,
  cstar    = 1,
  solid    = 0.5,
  cex      = 3,
  clab     = 0,
  leg      = TRUE
)
dev.off()
cat("Plot saved: DAPC_scatter.pdf\n")

# B. Compoplot (assignment probabilities per individual)
pdf("DAPC_compoplot.pdf", width = 12, height = 6)
compoplot(dapc_result, col = cores, lab = "", main = "Cluster Assignment Probabilities per Individual")
dev.off()
cat("Plot saved: DAPC_compoplot.pdf\n\n")

# ------------------------------------------------------------------------------
# 7. EXPORT RESULTS
# ------------------------------------------------------------------------------
cat("=== EXPORTING RESULTS ===\n")

# C. Table: assigned cluster + posterior probabilities
results_df <- data.frame(
  Sample           = indNames(genlight_obj),
  Assigned_Cluster = dapc_result$assign,
  dapc_result$posterior
)
write.csv(results_df, "DAPC_results.csv", row.names = FALSE)
cat("Table saved: DAPC_results.csv\n")

# D. Text summary
sink("DAPC_summary.txt")
cat("==================== DAPC ANALYSIS SUMMARY ====================\n")
cat("Date:                    ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total individuals:       ", nInd(genlight_obj), "\n")
cat("Total SNPs (post-QC):    ", nLoc(genlight_obj), "\n")
cat("Optimal K (BIC):         ", k_optimal, "\n")
cat("Optimal PCs (xvalDapc):  ", optimal_n_pca, "\n")
cat("Discriminant axes (DAs): ", n_da, "\n")
cat("Mean Successful Assignment:", round(msa_value, 4), "\n")
cat("===============================================================\n")
sink()
cat("Summary saved: DAPC_summary.txt\n")

cat("\n✓ Pipeline complete! Check the files in your working directory.\n")

cat("\n✓ Pipeline concluído! Verifique os arquivos no diretório de trabalho.\n")
