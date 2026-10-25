library(readr)
library(dplyr)

input_folder <- "Path to File"

all_files <- list.files(input_folder, pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
cat("Files found:", length(all_files), "\n")

if (length(all_files) == 0) {
  stop("No .csv files found — check input_folder path.")
}

# Helper: scan the first N rows (no header assumption) to find which row
# actually contains "gene" as a cell value, and return that row number
find_header_row <- function(f, max_scan_rows = 10) {
  preview <- tryCatch(
    read_csv(f, col_names = FALSE, n_max = max_scan_rows, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(preview)) return(NA_integer_)
  
  for (r in seq_len(nrow(preview))) {
    row_vals <- tolower(trimws(as.character(preview[r, ])))
    row_vals_clean <- gsub("[^a-z]", "", row_vals)
    if (any(row_vals_clean == "gene", na.rm = TRUE)) {
      return(r)
    }
  }
  return(NA_integer_)
}

find_gene_col <- function(col_names) {
  cleaned <- gsub("[^a-z]", "", tolower(col_names))
  match_idx <- which(cleaned == "gene")
  if (length(match_idx) == 0) return(NA_character_)
  col_names[match_idx[1]]
}

raw_data <- bind_rows(lapply(all_files, function(f) {
  
  hdr_row <- find_header_row(f)
  
  if (is.na(hdr_row)) {
    cat("Could not locate a 'gene' header in first 10 rows of:", basename(f), "\n")
    row1_preview <- tryCatch(read_csv(f, col_names = FALSE, n_max = 1, show_col_types = FALSE),
                             error = function(e) NULL)
    if (!is.null(row1_preview)) {
      cat("   Row 1 contents:", paste(row1_preview[1, ], collapse = " | "), "\n")
    }
    return(NULL)
  }
  
  df <- tryCatch(
    read_csv(f, skip = hdr_row - 1, show_col_types = FALSE),
    error = function(e) {
      cat("Could not read:", basename(f), "-", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(df) || nrow(df) == 0) {
    cat("Empty or unreadable file (after header fix):", basename(f), "\n")
    return(NULL)
  }
  
  gene_col <- find_gene_col(names(df))
  if (is.na(gene_col)) {
    cat("Header row found but column still not matched in:", basename(f),
        "| columns were:", paste(names(df), collapse = ", "), "\n")
    return(NULL)
  }
  
  # Use the parent subfolder name as the metabolite identifier
  metabolite_val <- basename(dirname(f))
  
  data.frame(
    metabolite = metabolite_val,
    gene       = df[[gene_col]],
    stringsAsFactors = FALSE
  )
}))

cat("Class of raw_data:", class(raw_data), "| Rows:", nrow(raw_data), "\n")

if (is.null(raw_data) || nrow(raw_data) == 0) {
  stop("raw_data is empty — check diagnostic messages above for header issues.")
}

# ---- Clean NA/blank gene names ----
g <- trimws(raw_data$gene)
na_like <- is.na(g) | g == "" | grepl("^n[\\s\\.\\/]*a$", g, ignore.case = TRUE, perl = TRUE)
data_clean <- raw_data[!na_like, , drop = FALSE]
data_clean$gene <- trimws(data_clean$gene)
cat("Rows after cleaning:", nrow(data_clean), "\n")

if (nrow(data_clean) == 0) stop("data_clean is empty after cleaning.")

# ---- Count how many distinct metabolites each gene is related to ----
metabolite_counts <- data_clean %>%
  distinct(metabolite, gene) %>%
  count(gene, name = "n_metabolites") %>%
  arrange(desc(n_metabolites))

cat("Unique genes found:", nrow(metabolite_counts), "\n")
print(head(metabolite_counts, 20))

# ---- Save results ----
out_dir <- "Path to File"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write.csv(metabolite_counts,
          file.path(out_dir, "metabolite_counts_per_gene.csv"),
          row.names = FALSE)

cat("Saved to:", file.path(out_dir, "metabolite_counts_per_gene.csv"), "\n")
