library(readr)
library(dplyr)

# ---- Config ----
CYTOHUBBA_PATH   <- "Path to File"
PATHWAY_PATH     <- "Path to File"
CONVERGENCE_PATH <- "Path to File"
OUT_DIR          <- "Path to File"
OUT_PATH         <- file.path(OUT_DIR, "Hub_Gene_Composite_Score.csv")

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
  cat("Created output folder:", OUT_DIR, "\n")
}

CANDIDATE_GENES <- c("AKT1", "JUN", "TP53", "MAPK1", "RELA")

CYTOHUBBA_MAX <- 12
PATHWAY_MAX   <- 6

# ---- Load the three scored tables ----
if (!file.exists(CYTOHUBBA_PATH)) stop("File not found: ", CYTOHUBBA_PATH)
if (!file.exists(PATHWAY_PATH)) stop("File not found: ", PATHWAY_PATH)
if (!file.exists(CONVERGENCE_PATH)) stop("File not found: ", CONVERGENCE_PATH)

cytohubba_summary <- read_csv(CYTOHUBBA_PATH, show_col_types = FALSE)
pathway_summary   <- read_csv(PATHWAY_PATH, show_col_types = FALSE)
convergence       <- read_csv(CONVERGENCE_PATH, show_col_types = FALSE)

# Print raw column names BEFORE anything else — this tells us the truth
cat("\n=== RAW COLUMN NAMES ===\n")
cat("cytohubba_summary:", paste(names(cytohubba_summary), collapse = " | "), "\n")
cat("pathway_summary:  ", paste(names(pathway_summary), collapse = " | "), "\n")
cat("convergence:      ", paste(names(convergence), collapse = " | "), "\n\n")

# ---- Auto-detect and rename the identifier column to "gene" in each table ----
rename_id_to_gene <- function(df, df_name) {
  cleaned_names <- gsub("[^a-z]", "", tolower(names(df)))
  match_idx <- which(cleaned_names %in% c("gene", "commonname", "genesymbol", "genename"))
  if (length(match_idx) == 0) {
    stop("[", df_name, "] No gene/common-name-like column found. Columns are: ",
         paste(names(df), collapse = ", "))
  }
  actual_col <- names(df)[match_idx[1]]
  cat("[", df_name, "] using column '", actual_col, "' as gene identifier\n", sep = "")
  df %>% rename(gene = !!actual_col)
}

cytohubba_summary <- rename_id_to_gene(cytohubba_summary, "cytohubba_summary")
pathway_summary   <- rename_id_to_gene(pathway_summary, "pathway_summary")
convergence       <- rename_id_to_gene(convergence, "convergence")

# Normalize metabolite convergence against the TRUE dataset-wide max
metabolite_max <- max(convergence$n_metabolites, na.rm = TRUE)
cat("\nMetabolite convergence normalized against dataset max:", metabolite_max, "\n")

# ---- Build the composite table ----
composite_table <- tibble(gene = CANDIDATE_GENES) %>%
  left_join(cytohubba_summary %>% dplyr::select(gene, cytohubba_consensus), by = "gene") %>%
  left_join(pathway_summary   %>% dplyr::select(gene, ec_pathway_count, ec_pathways), by = "gene") %>%
  left_join(convergence       %>% dplyr::select(gene, n_metabolites), by = "gene") %>%
  mutate(
    cytohubba_norm   = cytohubba_consensus / CYTOHUBBA_MAX,
    ec_pathway_norm  = ec_pathway_count / PATHWAY_MAX,
    metabolite_norm  = n_metabolites / metabolite_max,
    composite_score  = cytohubba_norm + ec_pathway_norm + metabolite_norm
  ) %>%
  arrange(desc(composite_score))

print(composite_table %>%
        dplyr::select(gene, cytohubba_consensus, ec_pathway_count, n_metabolites, composite_score))

write_csv(composite_table, OUT_PATH)
cat("\nSaved to:", OUT_PATH, "\n")

