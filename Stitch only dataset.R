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

# ---- Set source folder (contains the 202 folders) ----
root_folder <- "Path to File"

if (!dir.exists(root_folder)) {
  stop("Folder not found: ", root_folder)
}

# ---- Set destination folder (will contain only the STITCH folders) ----
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

# ---- Function to check if a datasheet has ANY row with STITCH in the source column ----
has_stitch <- function(df) {
  if (is.null(df)) return(FALSE)
  
  source_col <- names(df)[str_detect(tolower(names(df)), "source")]
  if (length(source_col) == 0) return(FALSE)
  
  source_col <- source_col[1]
  
  any(str_detect(str_to_upper(as.character(df[[source_col]])), "STITCH"), na.rm = TRUE)
}

# ---- Scan all folders and copy the STITCH ones directly ----
copy_log <- map_dfr(subfolders, function(folder) {
  folder_name <- path_file(folder)
  file_path <- find_datasheet(folder)
  
  if (is.na(file_path)) {
    return(tibble(Folder = folder_name, Has_STITCH = FALSE, Status = "No datasheet found"))
  }
  
  df <- read_datasheet(file_path)
  stitch_present <- has_stitch(df)
  
  if (stitch_present) {
    dest <- file.path(output_root, folder_name)
    dir_copy(folder, dest, overwrite = TRUE)
    return(tibble(Folder = folder_name, Has_STITCH = TRUE, Status = "Copied"))
  } else {
    return(tibble(Folder = folder_name, Has_STITCH = FALSE, Status = "Skipped (no STITCH)"))
  }
})

# ---- Save a log of what happened ----
write_csv(copy_log, file.path(root_folder, "..", "stitch_copy_log.csv"))

cat("\nDone.\n")
cat("Folders copied:", sum(copy_log$Status == "Copied"), "\n")
cat("Folders skipped:", sum(copy_log$Status != "Copied"), "\n")
cat("New folder location:", normalizePath(output_root), "\n")
print(copy_log)

