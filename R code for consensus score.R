library(readr)
library(dplyr)

# ---- Config ----
NODE_TABLE_PATH <- "Path to File"
CYTOHUBBA_PATH  <- "Path to File"
OUTPUT_PATH     <- "Path to File"
TOP_N           <- 30

if (!file.exists(NODE_TABLE_PATH)) stop("File not found: ", NODE_TABLE_PATH)
if (!file.exists(CYTOHUBBA_PATH)) stop("File not found: ", CYTOHUBBA_PATH)

# ---- Read and map gene symbols ----
node_table <- read_csv(NODE_TABLE_PATH, show_col_types = FALSE)
hubba      <- read_csv(CYTOHUBBA_PATH, show_col_types = FALSE)

names(node_table) <- trimws(names(node_table))
names(hubba) <- trimws(names(hubba))

cat("Columns in STRING node table:\n")
print(names(node_table))
cat("\nColumns in CytoHubba score file:\n")
print(names(hubba))

# ---- Check expected columns exist before mapping ----
if (!all(c("name", "display name") %in% names(node_table))) {
  stop("Expected columns 'name' and 'display name' not found in node table. ",
       "Available columns: ", paste(names(node_table), collapse = ", "))
}

if (!"node_name" %in% names(hubba)) {
  stop("Expected column 'node_name' not found in CytoHubba file. ",
       "Available columns: ", paste(names(hubba), collapse = ", "),
       "\nUpdate the join key below to match the actual identifier column name.")
}

gene_map <- node_table %>%
  dplyr::select(name, `display name`) %>%
  dplyr::rename(node_name = name, gene = `display name`)

hubba_named <- hubba %>% left_join(gene_map, by = "node_name")

cat("\nRows in cytoHubba export:", nrow(hubba), "\n")
cat("Rows mapped to a gene symbol:", sum(!is.na(hubba_named$gene)), "\n")

# ---- Flag top-N membership per method, then count consensus ----
hubba_methods <- c("MCC","DMNC","MNC","Degree","EPC","BottleNeck","EcCentricity",
                   "Closeness","Radiality","Betweenness","Stress","ClusteringCoefficient")

# Only keep methods actually present in the file
hubba_methods <- intersect(hubba_methods, names(hubba_named))
cat("\nScoring methods found and used:", paste(hubba_methods, collapse = ", "), "\n")

if (length(hubba_methods) == 0) {
  stop("None of the expected CytoHubba method columns were found. ",
       "Available columns: ", paste(names(hubba_named), collapse = ", "))
}

for (m in hubba_methods) {
  flag_col <- paste0(m, "_top", TOP_N)
  hubba_named[[flag_col]] <- rank(-as.numeric(hubba_named[[m]]), ties.method = "min") <= TOP_N
}

flag_cols <- paste0(hubba_methods, "_top", TOP_N)
hubba_named$cytohubba_consensus <- rowSums(hubba_named[, flag_cols], na.rm = TRUE)

cytohubba_summary <- hubba_named %>%
  dplyr::select(gene, node_name, all_of(hubba_methods), cytohubba_consensus) %>%
  arrange(desc(cytohubba_consensus))

# Write to a NEW file (not overwriting your original cytoHubba export)
write_csv(cytohubba_summary, OUTPUT_PATH)

cat("\nSummary saved to:", OUTPUT_PATH, "\n")

cat("\nTop 20 by consensus (out of", length(hubba_methods), "methods):\n")
print(head(cytohubba_summary[, c("gene","cytohubba_consensus")], 20))

