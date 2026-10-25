library(readr)
library(dplyr)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)

# =========================================================
# SECTION — GO/KEGG enrichment + EC-pathway plausibility
# (Intersected_Genes.csv version)
# =========================================================

## ---- 0. Paths -------------------------------------------------------
# This is a CSV file, so use read_csv(), not read_excel().
input_path  <- "Path to File"
output_dir  <- "Path to File"

## ---- 1. Read input ---------------------------------------------------
core_targets <- read_csv(input_path, show_col_types = FALSE)
names(core_targets) <- trimws(names(core_targets))

# Your file's gene column is named "Gene"
gene_col <- "Gene"

if (!gene_col %in% names(core_targets)) {
  stop("Column '", gene_col, "' not found. Available columns are: ",
       paste(names(core_targets), collapse = ", "))
}

# Clean gene symbols: trim whitespace, drop blanks/NAs, dedupe
gene_symbols <- core_targets[[gene_col]] %>%
  as.character() %>%
  trimws() %>%
  unique()
gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]

cat("Number of unique input gene symbols:", length(gene_symbols), "\n")
if (length(gene_symbols) == 0) {
  stop("No valid gene symbols found in the '", gene_col, "' column — check the input file.")
}

## ---- 2. Map SYMBOL -> ENTREZID --------------------------------------
id_map <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

cat("Mapped", nrow(id_map), "out of", length(gene_symbols), "gene symbols to Entrez IDs.\n")

# Show which symbols failed to map (helps catch typos/aliases/case issues)
unmapped <- setdiff(gene_symbols, id_map$SYMBOL)
if (length(unmapped) > 0) {
  cat("Unmapped gene symbols (", length(unmapped), "):\n", sep = "")
  print(unmapped)
}

entrez_ids <- unique(id_map$ENTREZID)

if (length(entrez_ids) == 0) {
  stop("bitr() mapped 0 genes to Entrez IDs. Check that:\n",
       "  1) gene_symbols contains real HGNC gene symbols (not e.g. Ensembl IDs),\n",
       "  2) there's no hidden whitespace/encoding issue in the '", gene_col, "' column,\n",
       "  3) org.Hs.eg.db is the correct organism database for your species.")
}

## ---- 3. GO (BP) enrichment ------------------------------------------
go_bp <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

if (is.null(go_bp) || nrow(as.data.frame(go_bp)) == 0) {
  warning("No significantly enriched GO BP terms found at the given cutoffs.")
  go_bp_df <- data.frame()
} else {
  go_bp_df <- as.data.frame(go_bp)
}

write.csv(go_bp_df, file.path(output_dir, "GO_BP_enrichment.csv"), row.names = FALSE)
saveRDS(go_bp_df, file.path(output_dir, "GO_BP_enrichment.rds"))

## ---- 4. KEGG enrichment ----------------------------------------------
kegg <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

if (is.null(kegg) || nrow(as.data.frame(kegg)) == 0) {
  warning("No significantly enriched KEGG pathways found at the given cutoffs.")
  kegg_df <- data.frame()
} else {
  kegg_readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  kegg_df <- as.data.frame(kegg_readable)
}

write.csv(kegg_df, file.path(output_dir, "KEGG_enrichment.csv"), row.names = FALSE)
saveRDS(kegg_df, file.path(output_dir, "KEGG_enrichment.rds"))

## ---- 5. EC-pathway plausibility ---------------------------------------
# Your specified EC pathways of interest (now 5 total):
EC_PATHWAYS <- c(
  "PI3K-Akt signaling pathway",
  "MAPK signaling pathway",
  "Estrogen signaling pathway",
  "Endometrial cancer",
  "Apoptosis"
)

cat("\nKEGG pathways available in your enrichment results:\n")
print(kegg_df$Description)

# Check which of your EC_PATHWAYS actually appear in the enrichment results
missing_pathways <- setdiff(EC_PATHWAYS, kegg_df$Description)
if (length(missing_pathways) > 0) {
  warning("The following EC_PATHWAYS were not found among the enriched KEGG results ",
          "(they will contribute 0 hits): ", paste(missing_pathways, collapse = ", "))
}

if (nrow(kegg_df) == 0) {
  warning("kegg_df is empty — EC-pathway plausibility scores will all be zero.")
}

pathway_membership <- lapply(EC_PATHWAYS, function(p) {
  row <- kegg_df[kegg_df$Description == p, ]
  if (nrow(row) == 0) return(character(0))
  strsplit(row$geneID, "/")[[1]]
})
names(pathway_membership) <- EC_PATHWAYS

ec_pathway_count <- sapply(gene_symbols, function(g) {
  sum(sapply(pathway_membership, function(genes) g %in% genes))
})

ec_pathway_which <- sapply(gene_symbols, function(g) {
  hits <- names(pathway_membership)[sapply(pathway_membership, function(genes) g %in% genes)]
  paste(hits, collapse = "; ")
})

pathway_summary <- data.frame(
  gene             = gene_symbols,
  ec_pathway_count = ec_pathway_count,
  ec_pathways      = ec_pathway_which,
  stringsAsFactors = FALSE
) %>% arrange(desc(ec_pathway_count))

write.csv(pathway_summary, file.path(output_dir, "EC_Pathway_Plausibility.csv"), row.names = FALSE)
saveRDS(pathway_summary, file.path(output_dir, "EC_Pathway_Plausibility.rds"))

save(go_bp_df, kegg_df, pathway_summary,
     file = file.path(output_dir, "enrichment_results.RData"))

cat("\nTop pathway-plausibility genes:\n")
print(head(pathway_summary, 15))

# =========================================================
# SECTION 6 - STORE RESULTS AS SEPARATE DATAFRAMES
# =========================================================

GO_BP_df <- go_bp_df
KEGG_df <- kegg_df
EC_Pathway_df <- pathway_summary

cat("\nDataframe Summary:\n")
cat("GO_BP_df rows:", nrow(GO_BP_df), "\n")
cat("KEGG_df rows:", nrow(KEGG_df), "\n")
cat("EC_Pathway_df rows:", nrow(EC_Pathway_df), "\n")

# =========================================================
# STORE EVERYTHING IN A MASTER LIST
# =========================================================

Enrichment_Results <- list(
  GO_BP = GO_BP_df,
  KEGG = KEGG_df,
  EC_Pathway = EC_Pathway_df
)

save(Enrichment_Results, file = file.path(output_dir, "Enrichment_Results.RData"))
save(GO_BP_df, file = file.path(output_dir, "GO_BP_df.RData"))
save(KEGG_df, file = file.path(output_dir, "KEGG_df.RData"))
save(EC_Pathway_df, file = file.path(output_dir, "EC_Pathway_df.RData"))

# =========================================================
# ACCESS DATAFRAMES LATER
# =========================================================
# load(file.path(output_dir, "Enrichment_Results.RData"))
# Enrichment_Results$GO_BP
# Enrichment_Results$KEGG
# Enrichment_Results$EC_Pathway

