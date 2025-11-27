# Alternative Splicing Analysis with Vast-tools Data

This folder contains R scripts for analyzing alternative splicing events using outputs from [Vast-tools](https://github.com/vastgroup/vast-tools), focusing on the FET protein family (EWSR1, FUS, TAF15) knockout vs. wildtype (WT) samples.

## Folder Structure

```
project-root/
├── Splicing_Script.R
├── Splicing_plots.R
├── data/                # [Place small/sample data here, DO NOT add raw or large data]
├── results/             # Output Excel and image files (auto-ignored by .gitignore)
├── .gitignore
└── README.md
```

## Contents

* **Splicing\_Script.R**
  Performs filtering of splicing events from Vast-tools output, splits by event type (EX, INT, ALTA, ALTD), extracts associated genes, and runs Gene Ontology (GO) enrichment analysis for each gene set.
  Results are saved as Excel files for each event type and sample.

* **Splicing\_plots.R**
  Loads filtered splicing events, creates intersections (UpSet plots) of events shared among EWSR1, FUS, and TAF15 knockouts, and exports detailed event lists and summary Excel files for downstream analysis or visualization.

## How to Use

### 1. Prerequisites

* R (>= 4.0.0)
* [Vast-tools](https://github.com/vastgroup/vast-tools) tab output files for your splicing samples
* Recommended: [RStudio](https://www.rstudio.com/)

**Required R packages:**

* `dplyr`
* `openxlsx`
* `UpSetR`
* `clusterProfiler`
* `org.Hs.eg.db`
* `writexl`

Install any missing packages in R:

```r
install.packages(c("dplyr", "openxlsx", "UpSetR", "writexl"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
```

### 2. Directory Setup

* Edit the `base_dir` variable in both scripts to match the location of your Vast-tools output files, e.g.:

  ```r
  base_dir <- "D:/FETseq/vast-tools_results"
  ```

### 3. Running the Scripts

#### a) **Filtering and GO Analysis**

* Open `Splicing_Script.R` in RStudio.
* Run the script.
  It will:

  * Filter events by minimum `MV.dPsi._at_0.95`
  * Split by event type (EX, INT, ALTA, ALTD)
  * Extract gene symbols for each set
  * Run GO enrichment and save results as Excel files

#### b) **Intersection Analysis and Plotting**

* Open `Splicing_plots.R` in RStudio.
* Run the script.
  It will:

  * Generate UpSet plots for each splicing event type
  * Export summary Excel files listing intersecting events, associated genes, MV values, and genomic coordinates