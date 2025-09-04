# ------------------------------------------------------------------------------
# Description:
#   Builds InstaPrism reference objects from the unified sc/sn reference
#   (produced by 1_get_data.R), then deconvolves bulk expression datasets
#   twice—(a) including adipocytes and (b) excluding adipocytes—to estimate
#   cell-type proportions. Reproducibility is enforced via fixed RNG seeds and
#   parallel-safe RNG settings.
#
# Inputs:
#   - output_data/sc_sn_reference_data/reference_expr_final.csv
#   - output_data/sc_sn_reference_data/reference_expr_final_genes.csv
#   - output_data/sc_sn_reference_data/reference_cell_types_final.csv
#   - input_data/intersect_3ds.csv              # gene symbols to remove (adipocyte run only)
#   - output_data/bulk_datasets/*_filtered_asImported.csv   # bulk RNA-seq (non-logged)
#
# Outputs (written under output_data/instaprism/):
#   - instaprism_reference_objects/
#       • instaprism_reference_with_adipos           # phi.cs matrix (CSV)
#       • instaprism_reference_no_adipos             # phi.cs matrix (CSV)
#   - instaprism_outputs/
#       • instaprism_output_<DATASET>_with_adipocytes.csv
#       • instaprism_output_<DATASET>_no_adipocytes.csv
#     (each file: samples × cell types; estimated proportions θ)
#
# Main Steps:
#   0) Load libraries; set parallel-safe RNG (L’Ecuyer-CMRG) and register workers
#   1) Read sc/sn reference matrices & labels; drop intersect_3ds genes (adipocyte run only)
#   2) Downsample reference to ≤500 cells per cell type (for speed/memory)
#   3) Build two InstaPrism reference objects: with and without adipocytes
#   4) Read bulk datasets (e.g., SchildkrautB, SchildkrautW)
#   5) Run InstaPrism deconvolution for each dataset with both references
#   6) Save theta (cell-type proportion) matrices to disk
# ------------------------------------------------------------------------------

# Clear R's memory to make room
rm(list = ls())
gc()

## ── Activate correct library path ────────────────────────────────────────
if (Sys.getenv("CONDA_DEFAULT_ENV") != "env_deconv_R") {
  stop("Please run `conda activate env_deconv_R` (or add it to your job script) ",
       "before launching 1_get_data.R")
}
## Show current library location for sanity
message("Using library path: ", .libPaths()[1])

# Load required libraries
library(InstaPrism)
library(data.table)
library(Matrix)
library(here)

################################################################################
### reproducibility
################################################################################
SEED <- 88

set.seed(SEED)       # base‑R RNG                       :contentReference[oaicite:0]{index=0}
RNGkind(kind = "L'Ecuyer-CMRG")   # parallel‑safe RNG   :contentReference[oaicite:1]{index=1}

library(BiocParallel)
register(
  MulticoreParam(workers = 4, RNGseed = SEED, progressbar = TRUE)
)
################################################################################

# Be sure to open this script along with a corresponding project
# that is located in the same directory as the scripts and the 
# input_data folder and enclosed files.
output_data <- file.path(here(), "output_data")
instaprism_folder <- file.path(output_data, "instaprism")

# Ensure instaprism_folder directory exists
dir.create(instaprism_folder, recursive = TRUE, showWarnings = FALSE)

##########################################################
# 0) Genes to remove (for adipocyte reference only)
##########################################################
genes_to_remove <- fread(
  file.path(here(), "input_data", "intersect_3ds.csv"),
  header = FALSE,        # your file has no header row
  data.table = FALSE
)[, 2]                   # second column holds the symbols
genes_to_remove <- unique(na.omit(genes_to_remove))

##########################################################
# 1) Read sn+sc reference files
##########################################################

# Read expression matrix
scExpr_all <- Matrix::readMM(file.path(output_data,"sc_sn_reference_data/reference_expr_final.csv"))
rownames(scExpr_all) <- (fread(file.path(output_data,"sc_sn_reference_data/reference_expr_final_genes.csv"),
                               header = TRUE, data.table = FALSE))$"x"
colnames(scExpr_all) <- (1:ncol(scExpr_all))

# Read cell type label (key)
cell_type_labels <- fread(file.path(output_data,"sc_sn_reference_data/reference_cell_types_final.csv"),
                          data.table = FALSE)
cell_type_labels <- cell_type_labels[,"cellType"]

#  Remove the unwanted genes *only for the adipocyte run* intersection genes based on previous work
scExpr_all <- scExpr_all[
  !rownames(scExpr_all) %in% genes_to_remove, ]

##########################################################
# 2) Select 500 cells of each cell type from the large 
#    reference dataset
##########################################################

# Function to select a random subset of each cell type to decrease the size of the single cell/nucleus reference matrix
select_cell_subsets <- function(cell_type_labels, n){
  cell_types <- unique(cell_type_labels) # Getting a list of cell types present
  all_selected_cells <- c()
  for(cell in cell_types){
    all <- grep(cell, cell_type_labels) # Get vector indices for cell type
    if(length(all) < n){
      all_selected_cells <- c(all_selected_cells, all)
      print(paste0("Fewer than ",n," cells of type ",cell," so will not subset. Total cell count ",length(all)))
      next
    }
    selected_cells <- sample(all, n) # Select a random subset of these vector inidices
    all_selected_cells <- c(all_selected_cells, selected_cells)
  }
  return(all_selected_cells) # Indices for selected cells
}

# Using the function defined above to subset the genes x cells matrix and the cell type labels vector
selected_cells <- select_cell_subsets(cell_type_labels, n = 500)
scExpr_subset_all <- scExpr_all[, selected_cells]
cell_type_labels_subset_all <- cell_type_labels[selected_cells]

# Removing large objects from R's memory to save room 
rm(scExpr_all, cell_type_labels)
gc()

# Now create a version with no adipocytes
adipos <- grep("Adipocytes", cell_type_labels_subset_all)
scExpr_subset_no_adipos <- scExpr_subset_all[,-adipos]
cell_type_labels_subset_no_adipos <- cell_type_labels_subset_all[-adipos]

##########################################################
# 3) Create and save reference objects for input
#    into InstaPrism
##########################################################
# Create reference objects for input into InstaPrism
refPhi_obj_all = refPrepare(sc_Expr = scExpr_subset_all,
                        cell.type.labels = cell_type_labels_subset_all,
                        cell.state.labels = cell_type_labels_subset_all)
refPhi_obj_no_adipos = refPrepare(sc_Expr = scExpr_subset_no_adipos,
                            cell.type.labels = cell_type_labels_subset_no_adipos,
                            cell.state.labels = cell_type_labels_subset_no_adipos)

# Write these reference objects as files
dir.create(file.path(instaprism_folder,"instaprism_reference_objects"), recursive = TRUE, showWarnings = FALSE)
write.csv(slot(refPhi_obj_all,"phi.cs"),
        file.path(instaprism_folder, "instaprism_reference_objects/instaprism_reference_with_adipos"))
write.csv(slot(refPhi_obj_no_adipos,"phi.cs"),
        file.path(instaprism_folder, "instaprism_reference_objects/instaprism_reference_no_adipos"))

##########################################################
# 4) Read bulk dataset files
##########################################################

# List of datasets
dataset_list <- c("SchildkrautB", "SchildkrautW")

# Reading the bulk datasets
# Instaprism requires non-log-transformed bulk data, so it is performed on the original, non-transformed data
# for the bulk RNA sequencing datasets and on the 2^(…) transformed data for the microarray datasets.
for (ds in dataset_list){
  if (ds == "SchildkrautB" | ds == "SchildkrautW"){
    bulk_expr <- fread(file.path(output_data, paste0("bulk_datasets/",ds,"_filtered_asImported.csv")),
          header = TRUE, data.table = FALSE)
  } else {
    bulk_expr <- fread(file.path(output_data, paste0("bulk_datasets/",ds,"_filtered_transformed.csv")),
                       header = TRUE, data.table = FALSE)
  }
  # Setting column 1 (sample IDs) as the rownames
  bulk_expr <- data.frame(bulk_expr[,-1], row.names=bulk_expr[,1]) 
  # Assigning object and transposing so that rows = genes and columns = samples
  assign(paste0("bulk_expr_",ds), t(bulk_expr)) 
}

##########################################################
# 5) Run InstaPrism twice on each bulk dataset, both with
#    and without adipocytes
##########################################################

# Run Instaprism with adipocytes
for (ds in dataset_list){
  print(paste0("Running InstaPrism (with adipocytes) on ", ds))
  instaprism_output <- InstaPrism(bulk_Expr = get(paste0("bulk_expr_", ds)), verbose = T, input_type="refPhi_cs",
                                  refPhi_cs = refPhi_obj_all, n.iter = 2000)
  assign(paste0("instaprism_output_", ds, "_with_adipocytes"),
         instaprism_output)
}

# Run Instaprism without adipocytes
for (ds in dataset_list){
  print(paste0("Running InstaPrism (no adipocytes) on ", ds))
  instaprism_output <- InstaPrism(bulk_Expr = get(paste0("bulk_expr_", ds)),  verbose = T, input_type="refPhi_cs",
                                  refPhi_cs = refPhi_obj_no_adipos, n.iter = 2000)
  assign(paste0("instaprism_output_", ds, "_no_adipocytes"),
         instaprism_output)
}

# Creating directory to store Instaprism outputs
dir.create(file.path(instaprism_folder,"instaprism_outputs"), recursive = TRUE, showWarnings = FALSE)

# Writing final files (with adipocytes)
for (ds in dataset_list) {
  write.csv(t((get(paste0("instaprism_output_", ds, "_with_adipocytes"))@Post.ini.ct@theta)),
            file.path(instaprism_folder, paste0("instaprism_outputs/instaprism_output_", ds,"_with_adipocytes.csv")),
            row.names=TRUE)
}

# Writing final files (no adipocytes)
for (ds in dataset_list) {
  write.csv(t((get(paste0("instaprism_output_", ds, "_no_adipocytes"))@Post.ini.ct@theta)),
            file.path(instaprism_folder, paste0("instaprism_outputs/instaprism_output_", ds,"_no_adipocytes.csv")),
            row.names=TRUE)
}
