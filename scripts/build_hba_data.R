# =============================================================================
# Build pre-indexed HbaDb annotation dataset combining:
#   1. database.xlsx (11,466 H. bacteriophora genes x 24 eggNOG/KEGG/PFAM/CAZy/MetaCyc columns)
#   2. org.Hbacteriophora.eg.sqlite (Bioconductor OrgDb: GO BP/CC/MF, KEGG, SYMBOL, GENENAME)
#   3. Optional GO.db term names (if installed) for human-readable GO enrichment & gene cards
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(DBI)
  library(RSQLite)
  library(dplyr)
})

app_dir <- "E:/Helpings/IARI, NEW DELHI/Wesly_work/Database_Hba"
xlsx_path <- file.path(app_dir, "database.xlsx")
sqlite_path <- "E:/Helpings/IARI, NEW DELHI/Wesly_work/R_package/org.Hbacteriophora.eg.sqlite"
out_dir <- file.path(app_dir, "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Reading database.xlsx...\n")
ann <- as.data.frame(readxl::read_excel(xlsx_path), stringsAsFactors = FALSE)
cat("Loaded database.xlsx:", nrow(ann), "rows x", ncol(ann), "cols\n")

# Clean "-" / NA strings consistently while keeping display tidy
clean_dash <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == "" | trimws(x) == "NA" | trimws(x) == "Not Available"] <- "-"
  trimws(x)
}

for (col in setdiff(colnames(ann), c("evalue", "score"))) {
  ann[[col]] <- clean_dash(ann[[col]])
}
ann$evalue <- suppressWarnings(as.numeric(ann$evalue))
ann$score  <- suppressWarnings(as.numeric(ann$score))

# Load org.Hbacteriophora.eg.sqlite to enrich with ontology-specific GO & OrgDb fields
cat("Reading org.Hbacteriophora.eg.sqlite...\n")
con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
org_meta <- DBI::dbGetQuery(con, "SELECT * FROM metadata")
org_genes <- DBI::dbGetQuery(con, "SELECT _id, GID FROM genes")
org_info <- DBI::dbGetQuery(con, "SELECT _id, SYMBOL, GENENAME FROM gene_info")
org_go <- DBI::dbGetQuery(con, "SELECT _id, GO, ONTOLOGY FROM go")
org_kegg <- DBI::dbGetQuery(con, "SELECT _id, KEGGPATHWAY FROM kegg")
DBI::dbDisconnect(con)

org_go <- merge(org_go, org_genes, by = "_id")
org_kegg <- merge(org_kegg, org_genes, by = "_id")
org_info <- merge(org_info, org_genes, by = "_id")

# Build per-ontology GO collapsed strings
collapse_by_gid <- function(df, val_col) {
  if (nrow(df) == 0) return(setNames(character(0), character(0)))
  sp <- split(df[[val_col]], df$GID)
  vapply(sp, function(v) paste(sort(unique(v)), collapse = ","), character(1))
}

go_bp_map <- collapse_by_gid(org_go[org_go$ONTOLOGY == "BP", ], "GO")
go_cc_map <- collapse_by_gid(org_go[org_go$ONTOLOGY == "CC", ], "GO")
go_mf_map <- collapse_by_gid(org_go[org_go$ONTOLOGY == "MF", ], "GO")
go_cnt_map <- tapply(org_go$GO, org_go$GID, function(v) length(unique(v)))
kegg_org_map <- collapse_by_gid(org_kegg, "KEGGPATHWAY")

ann$GO_BP <- ifelse(ann$GID %in% names(go_bp_map), unname(go_bp_map[ann$GID]), "-")
ann$GO_CC <- ifelse(ann$GID %in% names(go_cc_map), unname(go_cc_map[ann$GID]), "-")
ann$GO_MF <- ifelse(ann$GID %in% names(go_mf_map), unname(go_mf_map[ann$GID]), "-")
ann$GO_Count <- ifelse(ann$GID %in% names(go_cnt_map), as.integer(unname(go_cnt_map[ann$GID])), 0L)
ann$OrgDb_KEGG <- ifelse(ann$GID %in% names(kegg_org_map), unname(kegg_org_map[ann$GID]), "-")

# Reorder columns cleanly with GID, Preferred_name (SYMBOL), Description first
front_cols <- c("GID", "Preferred_name", "Description", "COG_category",
                "GO_Count", "GO_BP", "GO_CC", "GO_MF", "GOs",
                "KEGG_Pathway", "OrgDb_KEGG", "KEGG_ko", "KEGG_Module",
                "PFAMs", "PFam", "PFAM Accession", "EC", "CAZy",
                "MetaCyc/Reactome-Pathway", "evalue", "score",
                "seed_ortholog", "max_annot_lvl", "eggNOG_OGs",
                "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "BiGG_Reaction")
col_order <- unique(c(intersect(front_cols, colnames(ann)), colnames(ann)))
ann <- ann[, col_order, drop = FALSE]

# Build compact term -> GID lookup lists for instant ORA enrichment (<0.2s)
make_term2gene <- function(gids, terms) {
  df <- data.frame(term = terms, gid = gids, stringsAsFactors = FALSE)
  df <- unique(df[!is.na(df$term) & df$term != "" & df$term != "-", ])
  split(df$gid, df$term)
}

t2g_bp <- make_term2gene(org_go$GID[org_go$ONTOLOGY == "BP"], org_go$GO[org_go$ONTOLOGY == "BP"])
t2g_cc <- make_term2gene(org_go$GID[org_go$ONTOLOGY == "CC"], org_go$GO[org_go$ONTOLOGY == "CC"])
t2g_mf <- make_term2gene(org_go$GID[org_go$ONTOLOGY == "MF"], org_go$GO[org_go$ONTOLOGY == "MF"])
t2g_kegg <- make_term2gene(org_kegg$GID, paste0("map", org_kegg$KEGGPATHWAY))

# PFAM term2gene from PFAMs column
pfam_rows <- do.call(rbind, lapply(seq_len(nrow(ann)), function(i) {
  v <- ann$PFAMs[i]
  if (is.na(v) || v == "-") return(NULL)
  parts <- unique(trimws(unlist(strsplit(v, ","))))
  parts <- parts[nzchar(parts) & parts != "-"]
  if (!length(parts)) return(NULL)
  data.frame(term = parts, gid = ann$GID[i], stringsAsFactors = FALSE)
}))
t2g_pfam <- split(pfam_rows$gid, pfam_rows$term)

# COG term2gene from COG_category column
cog_rows <- do.call(rbind, lapply(seq_len(nrow(ann)), function(i) {
  v <- ann$COG_category[i]
  if (is.na(v) || v == "-") return(NULL)
  letters_v <- unique(unlist(strsplit(gsub("[^A-Z]", "", v), "")))
  if (!length(letters_v)) return(NULL)
  data.frame(term = letters_v, gid = ann$GID[i], stringsAsFactors = FALSE)
}))
t2g_cog <- split(cog_rows$gid, cog_rows$term)

# Optional human-readable GO term dictionary from GO.db if installed
go_term_names <- setNames(character(0), character(0))
if (requireNamespace("GO.db", quietly = TRUE)) {
  cat("Extracting GO term names from GO.db...\n")
  all_go_ids <- unique(org_go$GO)
  go_res <- tryCatch(
    AnnotationDbi::select(GO.db::GO.db, keys = all_go_ids, columns = c("TERM", "ONTOLOGY"), keytype = "GOID"),
    error = function(e) NULL
  )
  if (!is.null(go_res) && nrow(go_res) > 0) {
    go_res <- go_res[!duplicated(go_res$GOID) & !is.na(go_res$TERM), ]
    go_term_names <- setNames(go_res$TERM, go_res$GOID)
    cat("Resolved", length(go_term_names), "of", length(all_go_ids), "GO term names.\n")
  }
}

# Build PFAM short-name -> full description lookup from ann
pfam_desc_map <- setNames(character(0), character(0))
has_both_pfam <- ann$PFAMs != "-" & ann$PFam != "-" & !grepl(",", ann$PFAMs) & !grepl(";", ann$PFam)
if (any(has_both_pfam)) {
  sub_p <- unique(ann[has_both_pfam, c("PFAMs", "PFam")])
  sub_p <- sub_p[!duplicated(sub_p$PFAMs), ]
  pfam_desc_map <- setNames(sub_p$PFam, sub_p$PFAMs)
}

hba_bundle <- list(
  annotation    = ann,
  org_metadata  = org_meta,
  t2g           = list(
    GO_BP = t2g_bp,
    GO_CC = t2g_cc,
    GO_MF = t2g_mf,
    KEGG  = t2g_kegg,
    PFAM  = t2g_pfam,
    COG   = t2g_cog
  ),
  go_term_names = go_term_names,
  pfam_desc_map = pfam_desc_map,
  stats         = list(
    n_genes          = nrow(ann),
    n_cols           = ncol(ann),
    n_symbols        = sum(ann$Preferred_name != "-"),
    n_descriptions   = sum(ann$Description != "-"),
    n_go_genes       = sum(ann$GO_Count > 0L),
    n_go_terms       = length(unique(org_go$GO)),
    n_go_bp_terms    = length(t2g_bp),
    n_go_cc_terms    = length(t2g_cc),
    n_go_mf_terms    = length(t2g_mf),
    n_kegg_genes     = sum(ann$OrgDb_KEGG != "-"),
    n_kegg_pathways  = length(t2g_kegg),
    n_pfam_genes     = sum(ann$PFAMs != "-" | ann$PFam != "-"),
    n_pfam_domains   = length(t2g_pfam),
    n_metacyc_genes  = sum(ann[["MetaCyc/Reactome-Pathway"]] != "-"),
    n_ec_genes       = sum(ann$EC != "-"),
    n_cazy_genes     = sum(ann$CAZy != "-"),
    n_cog_genes      = sum(ann$COG_category != "-"),
    build_date       = as.character(Sys.Date())
  )
)

out_file <- file.path(out_dir, "hba_annotation.rds")
saveRDS(hba_bundle, out_file, compress = "xz")
cat("Saved", out_file, "| size:", round(file.info(out_file)$size / 1024^2, 2), "MB\n")
cat("Summary stats:\n")
print(unlist(hba_bundle$stats))
