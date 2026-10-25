# ---- Install/load required packages ----
required_packages <- c("readxl", "readr", "dplyr", "purrr", "stringr", "fs")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(readxl)
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(fs)

# ---- Set target folder (contains the 167 folders) ----
root_folder <- "Path to File"

if (!dir.exists(root_folder)) {
  stop("Folder not found: ", root_folder)
}

# ---- Set output location (final combined file goes directly here) ----
output_folder <- "Path to File"

# ---- Get list of all subfolders ----
subfolders <- dir_ls(root_folder, type = "directory")
cat("Number of subfolders found:", length(subfolders), "\n")

# ---- Function to find the datasheet inside a folder ----
find_datasheet <- function(folder) {
  files <- dir_ls(folder, regexp = "\\.(xlsx|xls|csv)$")
  if (length(files) == 0) return(NA_character_)
  files[1]
}

# ---- Function to read a file regardless of format ----
read_datasheet <- function(filepath) {
  ext <- tolower(path_ext(filepath))
  tryCatch({
    if (ext %in% c("xlsx", "xls")) {
      read_excel(filepath)
    } else if (ext == "csv") {
      read_csv(filepath, show_col_types = FALSE)
    } else {
      NULL
    }
  }, error = function(e) {
    message("Error reading file: ", filepath, " -> ", e$message)
    NULL
  })
}

# ---- Function to extract the Gene column from one datasheet ----
extract_genes <- function(df, folder_name) {
  if (is.null(df)) return(NULL)
  
  gene_col <- names(df)[str_detect(tolower(names(df)), "gene")]
  
  if (length(gene_col) == 0) {
    message("No 'Gene' column found in folder: ", folder_name)
    return(NULL)
  }
  
  gene_col <- gene_col[1]
  
  df %>%
    select(Gene = all_of(gene_col)) %>%
    mutate(Gene = as.character(Gene)) %>%
    filter(!is.na(Gene), Gene != "")
}

# ---- Loop through all 167 folders, collecting genes into one combined list ----
all_genes <- map_dfr(subfolders, function(folder) {
  folder_name <- path_file(folder)
  file_path <- find_datasheet(folder)
  
  if (is.na(file_path)) {
    message("No datasheet found in: ", folder_name)
    return(NULL)
  }
  
  df <- read_datasheet(file_path)
  extract_genes(df, folder_name)
})

# ---- Remove duplicates across ALL folders combined ----
unique_genes <- all_genes %>%
  distinct(Gene, .keep_all = TRUE)

cat("\nTotal gene entries collected (before dedup):", nrow(all_genes), "\n")
cat("Total unique genes (after dedup):", nrow(unique_genes), "\n")

# ---- Save the single final combined file ----
out_file <- file.path(output_folder, "Combined_Unique_Genes.csv")
write_csv(unique_genes, out_file)

cat("\nSaved final combined gene list to:", out_file, "\n")

