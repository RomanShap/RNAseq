
# Differential Expression Pipeline  
*Kallisto → DESeq2 → PCA/SVD → GO/FGSEA*

This repository contains an R-based pipeline for differential expression (DE) analysis starting from Kallisto quantification. The pipeline is designed to be:

- **Config-driven** (via a single `config.yml`).
- **Modular** (separate numbered scripts for each step).
- **Reusable on a cluster** (headless, Rscript-compatible, no interactive steps).

At a high level, the pipeline:

1. Reads Kallisto `abundance.tsv` files.
2. Aggregates transcript-level estimates to **gene-level** counts.
3. Applies a **bimodal expression filter** and runs DESeq2 for all conditions vs a reference.
4. Performs **PCA, SVD and hierarchical clustering** on filtered counts.
5. Runs **GO over-representation analysis** and **FGSEA**.
6. Computes **overlaps of gene sets and GO terms** across conditions and visualises them.
7. Runs **GO analysis on shared/non-shared gene groups** derived from overlaps.

Everything is orchestrated by a single main script (`FET_main.R` / `Quantification_main.R`).

---

## 1. Prerequisites

### 1.1 Software

- R (4.x recommended)
- Kallisto (quantifications already run; the pipeline assumes existing `abundance.tsv` files)

### 1.2 R packages

Install these (manually or via `renv`):

- Core and tidy: `yaml`, `dplyr`, `readr`, `tibble`, `purrr`, `tidyr`
- Statistics / DE: `DESeq2`
- Annotation: `org.Hs.eg.db`, `AnnotationDbi`, `biomaRt`
- GO and GSEA: `clusterProfiler`, `fgsea`, `msigdbr`
- Plotting: `ggplot2`, `ComplexHeatmap`, `circlize`, `VennDiagram`, `UpSetR`, `ComplexUpset`, `scales`, `grid`, `grDevices`, `forcats`
- I/O: `openxlsx`, `writexl` (or just `openxlsx` as used in the scripts)
- Optional: `limma` (for quantile normalisation in SVD helper)

The main script assumes these are available; it does **not** install them automatically.

---

## 2. Directory layout

### 2.1 Project root (dataset-specific)

This is where you run the pipeline for a given dataset. The main script expects:

project_root/
├── cluster/        # job scripts for the cluster (optional)
├── config/         # config/config.yml (your main YAML config)
├── logs/           # R logs written here
└── output/         # all analysis results, organised by step

If these subdirectories are missing, the main script will create them and then stop, so you can check everything before rerunning.

### 2.2 Pipeline root (this repository)

This is where the pipeline code lives:

pipeline_root/
├── resources/
│   └── rds/                       # cached tx2gene and MSigDB RDS files
└── scripts/
    └── Differential_expression_scripts/
        ├── 01_kallisto_to_counts.R
        ├── 02_deseq_bimodal_filter.R
        ├── 03_pca_filtered_counts.R        # in file 02.1_PCA_SVD.txt
        ├── 03_heatmap_lfc.R                # in file 02.2_heatmap_lfc.txt
        ├── 04_go_overrep.R                 # commented skeleton currently
        ├── 05_fgsea.R
        ├── 06_gene_overlap_plots.R
        ├── 07_go_overlap_plots.R
        ├── 08_go_shared_genes.R
        └── 99_utils.R

`FET_main.R` (or `Quantification_main.R`) lives in the pipeline root and:

- sets up paths,
- reads the config,
- defines global objects,
- and sources all helper scripts in this folder (ensuring `99_utils.R` is sourced first).

---

## 3. Configuration (`config/config.yml`)

The pipeline is driven by a single YAML file, typically:

project_root/config/config.yml

Example structure:

kallisto_dir: "/path/to/KallistoBootstrap"

rds_dir: "resources/rds"

sets:
  EWSR1: [EWSR1_1, EWSR1_2, EWSR1_3]
  FUS:   [FUS_1,   FUS_2,   FUS_3]
  TAF15: [TAF15_1, TAF15_2, TAF15_3]
  WT:    [WT_1,    WT_2,    WT_3]

de:
  lfc: 1
  padj: 0.05
  q_value: 0.2
  reference: "WT"

go:
  ontology: "ALL"
  simplify_cutoff: 0.6
  count_min: 10
  count_max: 500

colors:
  EWSR1: "#d53031"
  FUS:   "#009E73"
  TAF15: "#0072B2"
  WT:    "#999999"

cache:
  tx2gene: "tx2gene_hsapiens.rds"

msigdb_c5: "msigdb_C5_hs_symbol.rds"

pca:
  enabled: true
  counts_source: "kallisto_bimodal_filtered"
  n_top_var: NULL

  svd:
    enabled: true
    subtract_pc1: true
    max_components: 6

  hclust:
    enabled: true
    n_top_var: NULL
    dist_type: "pearson"
    hclust_method: "complete"

fgsea:
  nperm: 10000
  top_n: 10
  padj_max: 0.05
  abs_nes_min: 0.0
  max_pathways: 10

Key points:

- `kallisto_dir` must contain one subdirectory per sample, each holding an `abundance.tsv` file. The folder name must begin with the sample ID (e.g. `EWSR1_1...`).
- `sets` defines conditions and their samples. These names drive DESeq2 contrasts, group colours, and most downstream plots.
- `de$reference` is the condition used as reference level in DESeq2.
- `colors` must cover every condition in `sets`.
- `cache.tx2gene` and `msigdb_c5` are RDS files that will be created on first run and reused later.

---

## 4. Running the pipeline

From your project root:

cd /path/to/project_root

Rscript /path/to/pipeline_root/FET_main.R config/config.yml

Rscript /path/to/pipeline_root/FET_main.R

What happens:

1. The script determines `project_root` and `pipeline_root`.
2. It checks for required subdirectories and creates them if missing (`cluster`, `config`, `logs`, `output`, plus pipeline-side `resources`, `scripts`).
3. It reads `config/config.yml` into a global `cfg` list.
4. It builds `SAMPLE_INFO` and a named vector of Kallisto paths (`kallisto_paths`) from `cfg$sets` and `kallisto_dir`, with strict checks on folder names.
5. It opens a log file under `logs/` and redirects both standard output and messages there (via `sink()`).
6. It sets global thresholds (e.g., `DE_LFC`, `DE_PADJ`, `DE_QVAL`, `GO_ONTOLOGY`, colours).
7. It discovers all scripts in `scripts/Differential_expression_scripts/`, makes sure `99_utils.R` is first, and then sources all of them.

On a cluster, you typically wrap the `Rscript` call in a batch script under `cluster/` and submit via `sbatch` or equivalent.

---

## 5. Step-by-step description

Below is a conceptual description of what each numbered script does. Exact details are in the scripts themselves.

### 5.1 `01_kallisto_to_counts.R`

Goal: Convert Kallisto transcript-level quantification to gene-level counts.

Inputs:

- `kallisto_paths`
- `SAMPLE_INFO`
- `tx2gene_rds`
- `output_dir`

Main operations:

1. Read all `abundance.tsv` files.
2. Build or load a `tx2gene` mapping:
   - If `tx2gene_rds` exists, load it.
   - Otherwise, query Ensembl via `biomaRt` to map transcripts to genes, then save to RDS.
3. Aggregate transcript counts to gene counts for each sample.
4. Construct a gene × sample matrix: `COUNT_DATA_GENE`.
5. Save counts (RDS and XLSX) to `output/01_counts/`.

### 5.2 `02_deseq_bimodal_filter.R`

Goal: Filter lowly expressed genes with a bimodal model and run DESeq2.

Inputs:

- `COUNT_DATA_GENE`
- `SAMPLE_INFO`
- `DE_LFC`, `DE_PADJ`
- `output_dir`
- `cfg$de$reference`

Main operations:

1. Choose reference condition:
   - Use `cfg$de$reference` if provided.
   - Otherwise, fall back to `"WT"`, `"NT"`, or the first condition.
2. For each sample, apply `bimodal_filter_counts()`:
   - Fit a 2-component Gaussian mixture to log2(counts + 1).
   - Identify the “background” component.
   - Use a high percentile of that component as an expression threshold.
   - Mark genes as “expressed” if they exceed this threshold.
3. Keep genes that are expressed in at least one condition, forming `FILTERED_COUNTS`.
4. Create a `DESeqDataSet` from `FILTERED_COUNTS`, with `design = ~ condition`.
5. Run DESeq2 and extract results for each non-reference condition vs reference.
6. Store results in a list `DE_RESULTS` and save each contrast as:
   - `output/02_deseq/DE_<COND>_vs_<REF>.xlsx`
   - `output/02_deseq/DE_<COND>_vs_<REF>.rds`
7. Save `FILTERED_COUNTS` to RDS and XLSX in `output/02_deseq/`.

### 5.3 `03_heatmap_lfc.R`

Goal: Visualise log2 fold changes of shared DE genes across contrasts.

Inputs:

- `DE_RESULTS`
- `DE_LFC`, `DE_PADJ`
- `COLOR_KO`
- `output_dir`

Main operations:

1. Combine DE results across contrasts.
2. Select genes passing:
   - `padj < DE_PADJ`
   - `|log2FC| > DE_LFC`
3. Build a log2FC matrix (`lfc_mat`: genes × contrasts).
4. Generate a heatmap of LFC values, optionally:
   - clustering rows,
   - printing numeric LFC values inside cells only if the matrix is small.
5. Save the heatmap to `output/04_heatmap/`.

### 5.4 `03_pca_filtered_counts.R`

Goal: PCA, SVD and hierarchical clustering on count data.

Inputs:

- `cfg`
- `output_dir`
- `FILTERED_COUNTS` and/or `COUNT_DATA_GENE`
- `SAMPLE_INFO`
- Colour mapping

Main operations:

1. Decide which counts to use:
   - If `cfg$pca$counts_source == "kallisto_bimodal_filtered"`, load `filtered_counts_matrix.rds`.
   - Otherwise, use `COUNT_DATA_GENE`.
2. Convert counts to CPM and (optionally) select the `n_top_var` most variable genes.
3. Align `sample_info` with columns of the count matrix.
4. Run PCA via `pca_from_counts()`:
   - log2-transform,
   - centre and scale,
   - compute principal components.
5. Plot PC1 vs PC2 using `plot_pca_scores()` and save to `output/03_pca/`.
6. If `cfg$pca$hclust$enabled`:
   - Compute sample–sample distances,
   - Perform hierarchical clustering,
   - Save a dendrogram as `Sample_hclust_from_expression.png`.
7. If `cfg$pca$svd$enabled`:
   - Log2-transform the counts for SVD.
   - Optionally apply quantile normalisation.
   - Run SVD and compute explained variance per component.
   - Optionally subtract PC1 and recompute SVD on the residuals.
   - Save SVD results and plots in `output/03_pca/`.

### 5.5 `04_go_overrep.R`

Goal: GO over-representation (enrichGO) for DEGs.

Typical operations:

1. For each contrast:
   - Select all / upregulated / downregulated genes based on `DE_PADJ`, `DE_LFC`, `DE_QVAL`.
2. Convert gene IDs to Entrez using `AnnotationDbi` and `org.Hs.eg.db`.
3. Run `clusterProfiler::enrichGO()` for each set.
4. Simplify terms with `clusterProfiler::simplify()` using `GO_SIMPLIFY_CUTOFF`.
5. Filter by `GO_COUNT_MIN` and `GO_COUNT_MAX`.
6. Save results as XLSX tables and barplots (top N terms) to `output/04_go_overrep/`.

### 5.6 `05_fgsea.R`

Goal: GSEA using FGSEA for each contrast.

Inputs:

- `DE_RESULTS`
- `msigdb_c5_rds`
- `output_dir`
- FGSEA parameters

Main operations:

1. Load or download MSigDB C5 gene sets using `load_or_build_msigdb_c5()`.
2. For each contrast:
   - Build a named numeric vector of gene-level statistics (`ranks`), using only genes with known symbols.
   - Remove duplicates and NAs.
3. Run `fgsea::fgsea()` with:
   - `pathways = pathways_list`
   - `stats = ranks`
   - `nperm = FGSEA_NPERM`
4. Store FGSEA results in `FGSEA_RESULTS` and save to `output/06_fgsea/fgsea_results_list.rds`.
5. Generate:
   - Dotplots of the top N pathways per contrast.
   - An NES heatmap across contrasts, using `build_fgsea_nes_matrix()` and `make_fgsea_nes_heatmap()`.

### 5.7 `06_gene_overlap_plots.R`

Goal: Overlaps of gene sets across contrasts (UpSet + Venn).

Inputs:

- `DE_RESULTS`
- `DE_LFC`, `DE_PADJ`
- `output_dir`

Main operations:

1. For each contrast, extract sets of genes:
   - All DEGs,
   - Upregulated,
   - Downregulated.
2. Label genes by symbol when available, falling back to Ensembl IDs.
3. Build binary membership matrices from these sets using `build_binary_matrix_from_sets()`.
4. Produce UpSet plots for various combinations and export them to `output/07_gene_overlap/`.
5. Write the gene sets themselves as RDS:
   - `gene_sets_all_by_condition.rds`
   - `gene_sets_pos_by_condition.rds`
   - `gene_sets_neg_by_condition.rds`

### 5.8 `07_go_overlap_plots.R`

Goal: Overlaps of GO terms across contrasts.

Inputs:

- GO results from `04_go_overrep.R`
- `GO_COUNT_MIN`, `GO_COUNT_MAX`
- `output_dir`

Main operations:

1. From GO results, build sets of GO term descriptions for each condition and direction.
2. Use `build_binary_matrix_from_sets()` to generate membership matrices for GO terms.
3. Create UpSet plots and (optionally) z-score heatmaps of GO enrichment across contrasts.
4. Export overlap matrices and plots to `output/08_go_overlap/`.

### 5.9 `08_go_shared_genes.R`

Goal: GO analysis on shared/non-shared DEG groups across conditions.

Inputs:

- `output_dir`
- `DE_PADJ`, `DE_QVAL`
- `GO_ONTOLOGY`, `GO_SIMPLIFY_CUTOFF`, `GO_COUNT_MIN`, `GO_COUNT_MAX`
- RDS files from `07_gene_overlap`:
  - `gene_sets_all_by_condition.rds`
  - `gene_sets_pos_by_condition.rds`
  - `gene_sets_neg_by_condition.rds`

Main operations:

1. Load gene sets by condition.
2. Build exact membership groups using `build_exact_membership_groups()`:
   - Group names like `"EWSR1&FUS"`, `"FUS&TAF15"`, `"EWSR1"`, etc.
3. For each group (all/pos/neg), run GO over-representation using `run_enrich_go_for_result()`.
4. Export:
   - GO tables to XLSX,
   - Barplots of top N GO terms per group,
   under `output/09_go_shared_genes/`.

---

## 6. Output overview

All outputs are written under `project_root/output/`:

- `01_counts/`
  - Gene-level counts from Kallisto.

- `02_deseq/`
  - `DE_<COND>_vs_<REF>.xlsx` and `.rds` for each contrast.
  - `filtered_counts_matrix.rds` and `.xlsx`.
  - Optional bimodal histograms per sample.

- `03_pca/`
  - PCA objects (`pca_filtered_counts.rds`).
  - PCA plots (PC1 vs PC2).
  - Sample dendrograms.
  - SVD plots and summary tables (if enabled).

- `04_heatmap/`
  - LFC heatmaps of shared DEGs.

- `04_go_overrep/`
  - GO enrichment tables and barplots per contrast/direction.

- `06_fgsea/`
  - `fgsea_results_list.rds`.
  - Dotplots of enriched pathways.
  - NES heatmaps.

- `07_gene_overlap/`
  - Gene set RDS (all/pos/neg).
  - UpSet plots and Venn diagrams (if used).

- `08_go_overlap/`
  - GO overlap matrices and plots.

- `09_go_shared_genes/`
  - GO results for shared/non-shared DEG groups, with tables and barplots.

---

## 7. Common failure modes and checks

1. Sample–folder mismatches

   - Every sample in `cfg$sets` must have exactly one matching subfolder in `kallisto_dir` whose name starts with the sample ID.
   - If zero or multiple matches are found, the main script stops.

2. Config or directory missing

   - If `config/config.yml` is missing or unreadable, the script stops early.
   - If required `project_root` or `pipeline_root` subdirectories are missing, they are created and the script exits with a message.

3. Missing or inconsistent metadata

   - `SAMPLE_INFO` must contain all sample IDs present in the count matrix for DESeq2 and PCA.
   - If they do not match, DE or PCA scripts stop with an explicit error.

4. Too few genes or samples

   - PCA/SVD require at least 2 genes and 2 samples after filtering.
   - If filtering leaves too few, the PCA step aborts.

5. Annotation problems in GO/FGSEA

   - If genes cannot be mapped to Entrez IDs or if no gene sets pass thresholds, GO or FGSEA functions may return `NULL`.
   - Plotting helpers check for empty inputs and skip plotting rather than failing.

---

This README describes how the pipeline is structured and how each part fits together. For implementation details and edge cases, consult the comments in the individual scripts, especially `99_utils.R` and the numbered step scripts in `scripts/Differential_expression_scripts/`.
