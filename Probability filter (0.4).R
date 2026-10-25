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

# ---- Set source folder (contains the 167 folders) ----
root_folder <- "Path to File"

if (!dir.exists(root_folder)) {
  stop("Folder not found: ", root_folder)
}

# ---- Set output parent folder (where 167 output subfolders will be created) ----
output_root <- "Path to File"
dir_create(output_root)

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

# ---- Function to filter rows where Probability >= 0.4 ----
filter_probability <- function(df, folder_name, cutoff = 0.4) {
  if (is.null(df)) return(NULL)
  
  # Find the column name that matches "probability" (case-insensitive, flexible)
  prob_col <- names(df)[str_detect(tolower(names(df)), "probab")]
  
  if (length(prob_col) == 0) {
    message("No 'Probability' column found in folder: ", folder_name)
    return(NULL)
  }
  
  prob_col <- prob_col[1]  # use first match if multiple
  
  df_filtered <- df %>%
    filter(as.numeric(.data[[prob_col]]) >= cutoff)
  
  df_filtered
}

# ---- Main loop: process all 167 folders and save each result separately ----
summary_log <- map_dfr(subfolders, function(folder) {
  folder_name <- path_file(folder)
  file_path <- find_datasheet(folder)
  
  if (is.na(file_path)) {
    message("No datasheet found in: ", folder_name)
    return(tibble(Folder = folder_name, Status = "No datasheet found", Rows_Extracted = 0))
  }
  
  df <- read_datasheet(file_path)
  df_filtered <- filter_probability(df, folder_name, cutoff = 0.4)
  
  if (is.null(df_filtered)) {
    return(tibble(Folder = folder_name, Status = "No Probability column found", Rows_Extracted = 0))
  }
  
  # Create corresponding output subfolder
  out_subfolder <- file.path(output_root, folder_name)
  dir_create(out_subfolder)
  
  if (nrow(df_filtered) == 0) {
    message("No rows with Probability >= 0.4 in: ", folder_name)
    out_file <- file.path(out_subfolder, paste0(folder_name, "_filtered.csv"))
    write_csv(df_filtered, out_file)  # saves empty (headers only) for consistency
    return(tibble(Folder = folder_name, Status = "0 rows passed cutoff", Rows_Extracted = 0))
  }
  
  out_file <- file.path(out_subfolder, paste0(folder_name, "_filtered.csv"))
  write_csv(df_filtered, out_file)
  
  tibble(Folder = folder_name, Status = "Success", Rows_Extracted = nrow(df_filtered))
})

# ---- Save a summary log ----
write_csv(summary_log, file.path(output_root, "processing_summary.csv"))

cat("\nDone. Output folders created in:", output_root, "\n")
print(summary_log)

