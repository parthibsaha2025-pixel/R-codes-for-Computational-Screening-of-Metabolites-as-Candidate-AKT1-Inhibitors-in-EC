# ---- Install/load required packages ----
required_packages <- c("readxl", "readr", "dplyr", "stringr", "VennDiagram", "grid")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(readxl)
library(readr)
library(dplyr)
library(stringr)
library(VennDiagram)
library(grid)

# ---- Set file paths ----
genecards_path <- "Path to File"
combined_genes_path <- "Path to File" 

if (!file.exists(genecards_path)) stop("File not found: ", genecards_path)
if (!file.exists(combined_genes_path)) stop("File not found: ", combined_genes_path)

# ---- Set output folder ----
output_folder <- "Path to File"

# ---- Read both files ----
genecards_df <- read_excel(genecards_path)
combined_genes_df <- read_csv(combined_genes_path, show_col_types = FALSE)

# ---- Locate the relevant columns (flexible matching) ----
symbol_col <- names(genecards_df)[str_detect(tolower(names(genecards_df)), "symbol")]
if (length(symbol_col) == 0) {
  stop("No 'Symbol' column found in GeneCards EC.xlsx. Columns present: ",
       paste(names(genecards_df), collapse = ", "))
}
symbol_col <- symbol_col[1]

gene_col <- names(combined_genes_df)[str_detect(tolower(names(combined_genes_df)), "gene")]
if (length(gene_col) == 0) {
  stop("No 'Gene' column found in Combined_Unique_Genes.csv. Columns present: ",
       paste(names(combined_genes_df), collapse = ", "))
}
gene_col <- gene_col[1]

# ---- Clean up gene lists (trim whitespace, remove NA/blank, keep unique) ----
genecards_symbols <- genecards_df[[symbol_col]] %>%
  as.character() %>%
  str_trim() %>%
  na.omit() %>%
  unique()
genecards_symbols <- genecards_symbols[genecards_symbols != ""]

combined_genes <- combined_genes_df[[gene_col]] %>%
  as.character() %>%
  str_trim() %>%
  na.omit() %>%
  unique()
combined_genes <- combined_genes[combined_genes != ""]

cat("Number of unique symbols in GeneCards EC:", length(genecards_symbols), "\n")
cat("Number of unique genes in Combined_Unique_Genes:", length(combined_genes), "\n")

# ---- Perform intersection ----
intersected_genes <- intersect(genecards_symbols, combined_genes)

cat("Number of intersecting genes:", length(intersected_genes), "\n")

# ---- Store intersection as a dataframe ----
intersection_df <- tibble(Gene = sort(intersected_genes))

# ---- Save intersection dataframe ----
out_file <- file.path(output_folder, "Intersected_Genes.csv")
write_csv(intersection_df, out_file)
cat("Saved intersected gene list to:", out_file, "\n")

# ---- Generate Venn diagram ----
venn_file <- file.path(output_folder, "Venn_Diagram_GeneIntersection.png")

venn.plot <- venn.diagram(
  x = list(
    "GeneCards EC" = genecards_symbols,
    "Combined Unique Genes" = combined_genes
  ),
  filename = NULL,
  fill = c("#4E79A7", "#F28E2B"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = c(-20, 20),
  main = "Intersection of GeneCards EC and Combined Unique Genes"
)

png(venn_file, width = 1600, height = 1600, res = 300)
grid.draw(venn.plot)
dev.off()

cat("Saved Venn diagram to:", venn_file, "\n")

# ---- Clean up log files that VennDiagram package generates automatically ----
log_files <- list.files(output_folder, pattern = "^VennDiagram.*\\.log$", full.names = TRUE)
if (length(log_files) > 0) file.remove(log_files)

