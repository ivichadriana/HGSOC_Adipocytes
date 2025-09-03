# HGSOC Adipocytes Deconvolution

**Deconvolution of bulk RNA-seq of High Grade Serous Ovarian Carcinoma (HGSOC) incorporating adipocytes for survival analysis**

High Grade Serous Ovarian Carcinoma (HGSOC) is characterized by heterogeneity at the cellular level. Deconvolution methods can estimate cell-type composition from bulk RNA-seq data using single-cell references. However, adipocytes are often underrepresented in single-cell RNA-seq. This project integrates single-nucleus RNA-seq adipocyte data with HGSOC single-cell references and applies [InstaPrism](https://github.com/humengying0907/InstaPrism) for deconvolution of bulk RNA-seq and microarray datasets from Black and White patient cohorts: “SchildkrautB” and “SchildkrautW”.

## Features

- Preprocessing of bulk RNA-seq.
- TRanscriptomic subtypes: consensusOV.
- Construction of integrated single-cell/single-nucleus reference (adipocytes).
- Cox Proportional Hazard Survival analysis.
- Clinical metadata comparisons (race, age, BMI, FIGO stage, etc.)

## Running the full pipeline

1. Data Preparation

  Place all raw data files into `input_data/` following the structure described in the [Data Requirements](#data-requirements) section. Use `scripts/unzip_input_data.py` to decompress archives if necessary (not necessary if running through run_pipeline.sh)

2. **Clone the repository:**

   ```bash
   git clone https://github.com/ivichadriana/HGSOC_Adipocytes.git
   cd HGSOC_Adipocytes
   ```

3. **Run pipeline**

   ```bash
   bash run_pipeline.sh
   ```

  This script will create the environment (currently supports Linux or Mac) and run the full analysis.
  You can also navigate through scripts and run each step (they are numbered).

## Workflow Overview
0. **Create appropiate environments**: `0_create_environment.sh` creates conda R and Python environments.
1. **Decompress raw data**: `unzip_input_data.py` 
2. **Process bulk datasets**: `1_process_data_and_subtypes.R` filters, transforms, and clusters data.
3. **Prepare reference matrix**: `2_prepare_deconvolution.R` generates combined single-cell/nucleus reference InstaPrism.
4. **Run deconvolution**: `3_run_deconvolution.R` applies InstaPrism on bulk datasets for proportion estimates.
5. **Notebooks with analysis**: `notebooks/analysis_X.ipynb` Each notebook contains a distinct analysis. See file names for details.

## Notebooks

Interactive analyses are provided in `notebooks/`:

- `analysis_X.ipynb`: run analysis X where X is comparison (e.g., proportions vs. survival)

## Results

All generated files from deconvolution are saved under `output_data/`, 
and the visualizations and analysis are saved in the notebooks/ folder.

## Data Requirements

### Bulk HGSOC (Schildkraut) datasets and clinical data.

- **SchildkrautB** (Black patients; clinical metadata & raw counts)
- **SchildkrautW** (White patients; clinical metadata & raw counts)

> Available upon request. Please email [casey.greene@cuanschutz.edu](mailto\:casey.greene@cuanschutz.edu) for access.

### Reference gene lists & clustering metadata

#### greenelab/hgsc\_characterization

- **Reference mappings**\
  `reference_data/ensembl_hgnc_entrez.tsv`
- **Gene list & cluster assignments**\
  `data/way_pipeline_results_10removed_NeoRemoved_inclWhites/1.DataInclusion-Data-Genes/GlobalMAD_genelist.csv`\
  `data/way_pipeline_results_10removed_NeoRemoved_inclWhites/2.Clustering_DiffExprs-Tables-ClusterMembership/FullClusterMembership.csv`

### Single‑cell RNA‑seq (HGSOC reference for deconvolution)

**GEO accession:** [GSE217517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517)\
Replicates 1–8 (GSM6720925–GSM6720932) each include:

```
GSMxxxxxxx_single_cell_barcodes_<n>.tsv.gz  
GSMxxxxxxx_single_cell_features_<n>.tsv.gz  
GSMxxxxxxx_single_cell_matrix_<n>.mtx.gz
```

### Cell‑type labels

From **greenelab/deconvolution\_pilot**:

```
data/cell_labels/
  ├── 2251_labels.txt
  ├── 2267_labels.txt
  ├── 2283_labels.txt
  ├── 2293_labels.txt
  ├── 2380_labels.txt
  ├── 2428_labels.txt
  ├── 2467_labels.txt
  └── 2497_labels.txt
```

### Single‑nucleus adipocyte RNA‑seq

**GEO accession:** [GSE176171](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176171)\
Samples GSM5359325–GSM5820686 (e.g.):

```
GSM5359325_Hs_OAT_01-1.dge.tsv.gz  
…  
GSM5820686_Hs_SAT_11-1.dge.tsv.gz
```

---

**Directory layout:**

```text
project-root/
├── renv/
├── scripts/
├── input_data/    ← place all downloaded files here  
└── output_data/   ← script-generated results
```

## Contributing and Issues.

Please open issues or pull requests for improvements. We aim to be responsive in the Issues section.

## License

This project is licensed under the [BSD-3-Clause License](LICENSE).
