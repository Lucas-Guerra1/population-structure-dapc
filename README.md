# DAPC Pipeline — Population Genetics

Complete R pipeline for **Discriminant Analysis of Principal Components (DAPC)**, applied to population genetic structure from SNP data in VCF format.

---

## Overview

DAPC is a multivariate method free of population-level assumptions (such as Hardy-Weinberg equilibrium) that identifies and visualizes genetic clusters by maximizing between-group variance while minimizing within-group variance. This pipeline automates the QC, parameter selection, and visualization steps.

```
VCF → QC → Optimal K (BIC) → Optimal PCs (xvalDapc) → DAPC → Figures + CSV
```

---

## Dependencies

| Package    | Minimum version | Usage                        |
|------------|----------------|------------------------------|
| `vcfR`     | ≥ 1.14         | Reading VCF files            |
| `adegenet` | ≥ 2.1          | DAPC, find.clusters, xvalDapc|
| `viridis`  | ≥ 0.6          | Perceptual color palette     |

Installation:
```r
install.packages(c("vcfR", "adegenet", "viridis"))
```

---

## How to Use

1. Clone the repository and open `dapc_pipeline.R` in RStudio
2. Set the working directory (where outputs will be saved):
   ```r
   setwd("path/to/your/directory")
   ```
3. Run the script — a file selection window will open for you to choose the `.vcf` file

---

## Pipeline Steps

### 1. Quality Control (QC)
| Filter | Criterion |
|--------|-----------|
| Missing per individual | Removes samples with > 50% missing data |
| Missing per locus | Removes SNPs with > 10% missing data |
| MAF | Removes SNPs with minor allele frequency < 0.05 |

### 2. K Selection (find.clusters + BIC)
- Tests K from 1 to `min(10, √N)`
- Uses the **BIC** criterion with `diffNgroup` to select the lowest-cost K
- Initial number of PCs: `⌊N/3⌋` (conservative)

### 3. Cross-Validation (xvalDapc)
- Determines the optimal number of PCs for the final DAPC
- 30 replicates with 90% training data
- Criterion: lowest mean MSE

### 4. Final DAPC
- Runs with optimized K and n.pca
- Number of discriminant axes: `min(K − 1, 3)`

---

## Generated Outputs

| File | Description |
|------|-------------|
| `BIC_K_selection.pdf` | BIC curve for K selection |
| `cross_validation_PCs.pdf` | MSE by number of PCs (cross-validation) |
| `DAPC_scatter.pdf` | Scatter plot of discriminant axes |
| `DAPC_compoplot.pdf` | Assignment probabilities per individual |
| `DAPC_results.csv` | Assigned cluster + posterior probabilities per sample |
| `DAPC_summary.txt` | Summary of analysis parameters and metrics |

---

## References

- Jombart T, Devillard S, Balloux F (2010). *Discriminant analysis of principal components: a new method for the analysis of genetically structured populations*. BMC Genetics, 11:94. [doi:10.1186/1471-2156-11-94](https://doi.org/10.1186/1471-2156-11-94)
- Jombart T, Ahmed I (2011). *adegenet 1.3-1: new tools for the analysis of genome-wide SNP data*. Bioinformatics, 27:3070–3071.
- Miller JM et al. (2020). Molecular Ecology Resources.

---

## Notes

- The script is **non-interactive**: all parameters are determined automatically
- `set.seed(999)` ensures reproducibility of the cross-validation
- For large datasets (> 5,000 SNPs), the `xvalDapc` step may be slow — reduce `n.rep` if needed
