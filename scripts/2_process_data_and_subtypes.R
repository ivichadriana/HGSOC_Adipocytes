##########################################################################################
### 1_process_data_and_subtypes.R
### 
### This script reads the bulk RNA sequencing datasets and filters them
### to only include genes present in one common gene mapping list.
### – this is used for InstaPrism deconvolution. It also transforms the bulk RNA sequencing
### data using log10(***+1) - this is used for clustering. All of these matrices are saved in 
### a uniform format containing on sample ID (rows) and genes (columns);
### file names are appended with either “asImported” or
### “transformed.” It also saves any metadata information (ex. prior clustering) about
### the samples in a separate file for reference.
##########################################################################################

# These are large files that take up a lot of memory - clear R's memory to make room
rm(list = ls())
gc()

## ── Activate correct library path ────────────────────────────────────────
if (Sys.getenv("CONDA_DEFAULT_ENV") != "env_deconv_R") {
  stop("Please run `conda activate env_deconv_R` (or add it to your job script) ", 
    "before launching 1_process_data_and_subtypes.R")

}
## Show current library location for sanity
message("Using library path: ", .libPaths()[1])

# Load Required Libraries
library(data.table)
library(dplyr)
library(here)
library(consensusOV)

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

# Set base directories
# Be sure to open this script along with a corresponding project
# that is located in the same directory as the scripts and the 
# input_data folder and enclosed files.
input_data <- file.path(here(), "input_data")
output_data <- file.path(here(), "output_data/bulk_datasets")

# Ensure output_data directory exists
dir.create(output_data, recursive = TRUE, showWarnings = FALSE)

# Load Gene List
MAD_genes <- data.frame(fread(file.path(input_data,"GlobalMAD_genelist.csv")))
colnames(MAD_genes) <- "hgnc_symbol"

##########################################################
# 1) Helper reading/merging functions 
##########################################################

# read Schildkraut data from file
read_format_expr <- function(in_file, metadata_table){
  # 1) Read file with fill=TRUE to handle irregular columns
  rnaseq_expr_df <- fread(in_file,
                          header = TRUE,
                          data.table = FALSE,
                          fill   = TRUE)
  
  # 2) Remove any trailing empty columns
  empty_cols <- sapply(rnaseq_expr_df, function(x) all(is.na(x)))
  if (any(empty_cols)) {
    rnaseq_expr_df <- rnaseq_expr_df[, !empty_cols, drop = FALSE]
  }
  
  # 3) The first column is gene IDs
  gene_ids <- rnaseq_expr_df[, 1]
  rnaseq_expr_df <- rnaseq_expr_df[, -1, drop = FALSE]
  
  # 4) Transpose
  rnaseq_expr_df <- data.frame(t(rnaseq_expr_df))
  colnames(rnaseq_expr_df) <- gene_ids
  
  # 5) Create sample ID column
  rnaseq_expr_df$ID <- gsub("^Sample_", "", rownames(rnaseq_expr_df))
  
  # 6) Merge with metadata by ID
  full_df <- merge(metadata_table, rnaseq_expr_df, by = "ID")
  
  # Return the combined data.frame as [[1]] and gene IDs as [[2]]
  return(list(full_df, gene_ids))
}

##########################################################
# 2) Filter to MAD_genes & write out as CSV plus metadata
#    We'll define a function that saves both the 'as
#    imported' version and the 'transformed' version.
##########################################################

save_dual_versions <- function(expr_merged, dataset_name, transform_type = c("rnaseq", "microarray")) {
  expr_merged <- as.data.frame(expr_merged)
  
  # 1) Save as-imported (raw/log2) version
  save_filtered(expr_merged, dataset_name, suffix = "asImported")
  
  # 2) Create a “transformed” version
  transform_type <- match.arg(transform_type)
  if (transform_type == "rnaseq") {
    expr_transformed <- transform_rnaseq(expr_merged)
    
    # 3) Save the transformed version
    save_filtered(expr_transformed, dataset_name, suffix = "transformed")
  } 
  if (transform_type == "microarray") {
    expr_transformed <- transform_microarray(expr_merged)
  
  # 3) Save the transformed version
  save_filtered(expr_transformed, dataset_name, suffix = "transformed")
}}

# This sub-function saves ONLY ID + gene columns in the expression CSV,
# and everything else in a separate metadata file.
save_filtered <- function(expr_merged, dataset_name, suffix) {
  # 1) Gene columns = those in the MAD_genes list
  valid_genes <- intersect(colnames(expr_merged), MAD_genes$hgnc_symbol)
  
  # 2) The final expression CSV has only: ID + (valid_genes)
  final_order <- c("ID", valid_genes)
  final_order <- intersect(final_order, colnames(expr_merged))
  expr_filtered <- expr_merged[, final_order, drop = FALSE]

  # Write expression CSV
  expr_output_file <- file.path(output_data,
    paste0(dataset_name, "_filtered_", suffix, ".csv")
  )
  write.csv(expr_filtered, expr_output_file, row.names = FALSE)
  
  # 3) Metadata file: everything else (excluding ID + valid_genes)
  meta_cols <- setdiff(colnames(expr_merged), final_order)
  metadata_only <- expr_merged[, meta_cols, drop = FALSE]
  
  meta_output_file <- file.path(output_data,
    paste0(dataset_name, "_metadata_", suffix, ".csv")
  )
  write.csv(metadata_only, meta_output_file, row.names = FALSE)
}

##########################################
# Helper transformations for each data type
##########################################

# 1) For RNA-seq => log10(... + 1)
transform_rnaseq <- function(expr_df) {
  df <- expr_df
  # We'll define 'ID' as non-gene. 
  # Every other column that is also in 'MAD_genes' is numeric. 
  # So we can just do:
  numeric_cols <- intersect(colnames(df), MAD_genes$hgnc_symbol)
  
  for (cc in numeric_cols) {
    df[[cc]] <- as.numeric(df[[cc]])
    df[[cc]] <- log10(df[[cc]] + 1)
  }
  return(df)
}

# 2) For microarray => 2^(...) to revert from log2
transform_microarray <- function(expr_df) {
  df <- expr_df
  # Again, only exponentiate gene columns
  numeric_cols <- intersect(colnames(df), MAD_genes$hgnc_symbol)
  
  for (cc in numeric_cols) {
    df[[cc]] <- as.numeric(df[[cc]])
    df[[cc]] <- 2^(df[[cc]])
  }
  return(df)
}

##########################################################
# 3) Load microarray metadata
##########################################################

clust_file <- file.path(input_data, "FullClusterMembership.csv")

if (file.exists(clust_file)) {
  clust_df <- data.frame(fread(clust_file))
  colnames(clust_df)[1] <- "ID"
} else {
  stop("Could not find FullClusterMembership.csv for microarray metadata.")
}

##########################################################
# 4) Process, transform, and save datasets
##########################################################

## A) SchildkrautB (RNA-seq)
message("Processing SchildkrautB (RNA-seq)")
schildB_metadata <- fread(file.path(input_data, "main_AA_metadata_table.tsv"))
schildB_res <- read_format_expr(
  in_file = file.path(input_data, "salmon_raw_counts_for_way_pipeline.tsv"),
  metadata_table = schildB_metadata
)
save_dual_versions(
  expr_merged = schildB_res[[1]],
  dataset_name = "SchildkrautB",
  transform_type = "rnaseq"
)

## B) SchildkrautW (RNA-seq)
message("Processing SchildkrautW (RNA-seq)")
schildW_metadata <- fread(file.path(input_data, "main_white_metadata_table.tsv"))
schildW_res <- read_format_expr(
  in_file = file.path(input_data, "salmon_raw_counts_for_way_pipeline_whites.tsv"),
  metadata_table = schildW_metadata
)
save_dual_versions(
  expr_merged = schildW_res[[1]],
  dataset_name = "SchildkrautW",
  transform_type = "rnaseq"
)

message("All datasets transformed and saved as both asImported & transformed versions!")

######### Now Subtyping #########

# Set base directories
# Be sure to open this script along with a corresponding project
# that is located in the same directory as the scripts and the 
# input_data folder and enclosed files.
input_data <- file.path(here(), "input_data")
output_data <- file.path(here(), "output_data")
clustering_results_dir <- file.path(output_data, "bulk_data_clustering_and_subtypes")

# Ensure clustering_results_dir directory exists
dir.create(clustering_results_dir, recursive = TRUE, showWarnings = FALSE)

# Set seed for reproducibility
set.seed(5)

# Load gene_map, needed for real Entrez IDs for consensusOV
gene_map_file <- file.path(input_data, "ensembl_hgnc_entrez.tsv")
if (!file.exists(gene_map_file)) {
  stop("gene_map_file not found!")
}
gene_map <- data.frame(fread(gene_map_file)) 
# We expect columns like: ensembl_gene_id, hgnc_symbol, entrezgene_id, etc.

##########################################################
# 1) Helper reading, clustering, and subtyping functions 
##########################################################

# read_dataset_counts():
# Reads a CSV (with ID + gene columns) that has *rows=samples* after we set rownames=ID.
# This function reads raw counts (for RNAseq data - asImported), and
# 2^(...) scaled "pseudocounts" data (for microarray - transformed).
read_dataset_counts <- function(dataset_name) {
  if (dataset_name == "SchildkrautB" | dataset_name == "SchildkrautW"){
    fpath <- file.path(output_data, paste0("bulk_datasets/", dataset_name, "_filtered_asImported.csv"))
    df <- data.frame(fread(fpath))
    rownames(df) <- df$ID
    df$ID <- NULL
    return(df)  # row=sample, col=gene
  } else {
    fpath <- file.path(output_data, paste0("bulk_datasets/", dataset_name, "_filtered_transformed.csv"))
    df <- data.frame(fread(fpath))
    rownames(df) <- df$ID
    df$ID <- NULL
    return(df)
  }
}

# read_dataset_log():
# Reads a CSV (with ID + gene columns) that has *rows=samples* after we set rownames=ID.
# This function reads log10(...+1) transformed data (for RNAseq data - transformed), and
# raw log2 data (for microarray - asImported).
read_dataset_log <- function(dataset_name) {
  if (dataset_name == "SchildkrautB" | dataset_name == "SchildkrautW"){
    fpath <- file.path(output_data, paste0("bulk_datasets/", dataset_name, "_filtered_transformed.csv"))
    df <- data.frame(fread(fpath))
    rownames(df) <- df$ID
    df$ID <- NULL
    return(df)  # row=sample, col=gene
  } else {
    fpath <- file.path(output_data, paste0("bulk_datasets/", dataset_name, "_filtered_asImported.csv"))
    df <- data.frame(fread(fpath))
    rownames(df) <- df$ID
    df$ID <- NULL
    return(df)
  }
}

# run_consensusOV():
# We will run cOV subtyping on raw counts/2^(...) scaled data
# cOV expects row=genes, col=samples, plus *actual Entrez IDs* for geneIDs.
# We'll do the same transpose, then map gene symbols -> entrez IDs
run_consensusOV <- function(samp_x_gene) {

  # put everything on ~log₂ scale 
  samp_x_gene <- log2(samp_x_gene + 1)
  # transpose => row=genes, col=samples
  gene_x_samp <- t(samp_x_gene)
  
  # rownames(gene_x_samp) are gene symbols. We must find their entrez IDs in gene_map
  # We'll keep only genes that appear in gene_map. 
  # Because cOV needs a numeric vector of geneIDs in the same order as the rows.
  all_genes <- rownames(gene_x_samp)
  # subset gene_map to these symbols
  map_sub <- subset(gene_map, hgnc_symbol %in% all_genes & !is.na(entrezgene_id))
  # remove duplicates if any
  map_sub <- unique(map_sub[, c("hgnc_symbol","entrezgene_id")])
  # reorder map_sub to match rownames(gene_x_samp)
  map_sub <- map_sub[match(all_genes, map_sub$hgnc_symbol), ]
  
  # Some genes might be missing from map_sub => NA in entrezgene_id => exclude them
  keep_idx <- which(!is.na(map_sub$entrezgene_id))
  # filter gene_x_samp to only those rows
  gene_x_samp_filt <- gene_x_samp[keep_idx, , drop=FALSE]
  gene_ids_W <- map_sub$entrezgene_id[keep_idx]
  
  if (nrow(gene_x_samp_filt) < 2) {
    # too few genes => can't run cOV meaningfully
    # return all NA
    return(data.frame(
      ID=rownames(samp_x_gene),
      consensusOV=rep(NA_character_, nrow(samp_x_gene))
    ))
  }
  
  # Now run get.subtypes
  subtypes_res <- get.subtypes(gene_x_samp_filt, gene_ids_W, method='consensus')
  
  # Generate a dataframe that aligns sample ID with consensus subtype
  out_df <- data.frame(ID = rownames(subtypes_res$rf.probs),
                       consensusOVsubtype = subtypes_res$consensusOV.subtypes)
  
  # Ensure that no samples were accidentally filtered out
  if (nrow(out_df) != nrow(samp_x_gene)) {
    warning("consensusOV did not generate a result for all samples.")
  }
  
  return(out_df)
}

# combine_clusterings():
# merges multiple data frames with the same 'ID' column.
combine_clusterings <- function(...) {
  dfs <- list(...)
  out <- Reduce(function(x, y) merge(x, y, by='ID', all=TRUE), dfs)
  return(out)
}

##################################################
# 3) Main Script
##################################################

dataset_list <- c("SchildkrautB", "SchildkrautW")

all_clustering_res <- list()

for (ds in dataset_list) {
  message("Clustering dataset: ", ds)
  
  # 1) Load counts and log transformed (row=sample, col=gene)
  expr_df_counts <- read_dataset_counts(ds)
  expr_df_log <- read_dataset_log(ds)
  
  # 2) consensusOV
  cov_out <- run_consensusOV(expr_df_counts)
  
  # 3) Combine
  combined <- cov_out
  combined$Dataset <- ds
  
  # 4) Write per-dataset
  outf <- file.path(clustering_results_dir, paste0(ds, "_clustering_labels.csv"))
  write.csv(combined, outf, row.names=FALSE)
  
  all_clustering_res[[ds]] <- combined
}

final_df <- dplyr::bind_rows(all_clustering_res)
final_out <- file.path(clustering_results_dir, "all_datasets_clustering_labels.csv")
write.csv(final_df, final_out, row.names=FALSE)

message("Done! Clustering results written to ", final_out)
