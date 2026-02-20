# DAPC Pipeline — Population Genetics

A complete R pipeline for Discriminant Analysis of Principal Components (DAPC), applied to population genetic structure inference from SNP data in VCF format.

---

## Table of Contents

- [Overview](#overview)
- [Workflow](#workflow)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Output Files](#output-files)
- [Notes](#notes)
- [Limitations](#limitations)
- [References](#references)
- [License](#license)

---

## Overview

DAPC is a multivariate statistical method free of population-level assumptions (such as Hardy-Weinberg equilibrium) that identifies and visualizes genetic clusters by maximizing between-group variance while minimizing within-group variance.

This pipeline automates quality control, optimal parameter selection (K and number of PCs), and figure generation, requiring only a `.vcf` or `.vcf.gz` file as input.

---

## Workflow

```
VCF → QC → Optimal K (BIC) → Optimal PCs (xvalDapc) → DAPC → Figures + CSV
```

---

## Dependencies

| Package | Minimum version | Purpose |
|---|---|---|
| `vcfR` | ≥ 1.14 | Reading VCF files |
| `adegenet` | ≥ 2.1 | DAPC, `find.clusters`, `xvalDapc` |
| `viridis` | ≥ 0.6 | Perceptual color palettes |

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/your-username/dapc-pipeline.git
cd dapc-pipeline
```

### 2. Install R dependencies

```r
install.packages(c("vcfR", "adegenet", "viridis"))
```

---

## Usage

1. Open `dapc_pipeline.R` in RStudio or run it from the terminal
2. Set the working directory (where output files will be saved):

```r
setwd("path/to/your/output/directory")
```

3. Run the script — a file selection dialog will open for you to choose the `.vcf` or `.vcf.gz` input file:

```bash
Rscript dapc_pipeline.R
```

---

## Pipeline Steps

### 1. Quality Control (QC)

| Filter | Criterion |
|---|---|
| Missing per individual | Removes samples with > 50% missing data |
| Missing per locus | Removes SNPs with > 10% missing data |
| MAF | Removes SNPs with minor allele frequency < 0.05 |

### 2. K Selection (`find.clusters` + BIC)

- Tests K from 1 to min(10, √N)
- Uses the BIC criterion with `diffNgroup` to select the lowest-cost K
- Initial number of PCs: ⌊N/3⌋ (conservative rule)

### 3. Cross-Validation (`xvalDapc`)

- Determines the optimal number of PCs for the final DAPC
- 30 replicates with 90% training data
- Selection criterion: lowest mean MSE

### 4. Final DAPC

- Runs with optimized K and n.pca
- Number of discriminant axes: min(K − 1, 3)

---

## Output Files

| File | Description |
|---|---|
| `BIC_K_selection.pdf` | BIC curve for K selection |
| `cross_validation_PCs.pdf` | MSE by number of PCs (cross-validation) |
| `DAPC_scatter.pdf` | Scatter plot of discriminant axes |
| `DAPC_compoplot.pdf` | Assignment probabilities per individual |
| `DAPC_results.csv` | Assigned cluster and posterior probabilities per sample |
| `DAPC_summary.txt` | Summary of analysis parameters and metrics |

---

## Notes

- The script is **non-interactive**: all parameters are determined automatically
- `set.seed(999)` ensures reproducibility of cross-validation results
- For large datasets (> 5,000 SNPs), the `xvalDapc` step may be slow — reduce `n.rep` if needed
- The script handles diploid data only; polyploid VCF files are not supported

---

## Limitations

- Designed for **diploid organisms** only; polyploid genotypes are not handled
- Minimum recommended sample size: ~10 individuals per expected cluster
- Performance degrades with very low sample sizes relative to K
- The `xvalDapc` step can be computationally intensive for large SNP datasets

---

## References

Jombart T, Devillard S, Balloux F (2010). Discriminant analysis of principal components: a new method for the analysis of genetically structured populations. *BMC Genetics*, 11:94. [doi:10.1186/1471-2156-11-94](https://doi.org/10.1186/1471-2156-11-94)

Jombart T, Ahmed I (2011). adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. *Bioinformatics*, 27:3070–3071. [doi:10.1093/bioinformatics/btr521](https://doi.org/10.1093/bioinformatics/btr521)

Miller JM, Picco C, Ciarniello M, de Haan A (2020). Population structure and demographic history of an alpine ungulate. *Molecular Ecology Resources*, 20(5):1390–1406. [doi:10.1111/1755-0998.13195](https://doi.org/10.1111/1755-0998.13195)

---

## License

MIT License — feel free to use, modify, and distribute with attribution.
