# =============================================================================
# HbaDb: Heterorhabditis bacteriophora Functional Genome Database
# Companion Interactive Shiny Portal for Bioconductor's org.Hbacteriophora.eg.db
# Bioconductor DOI: https://doi.org/10.18129/B9.bioc.org.Hbacteriophora.eg.db
# =============================================================================

library(shiny)
library(bslib)
library(bsicons)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)
library(writexl)

# -----------------------------------------------------------------------------
# 1. Load Pre-Indexed Annotation Bundle (with graceful database.xlsx fallback)
# -----------------------------------------------------------------------------
load_hba_bundle <- function() {
  rds_path <- file.path("data", "hba_annotation.rds")
  if (file.exists(rds_path)) {
    return(readRDS(rds_path))
  }
  if (file.exists("database.xlsx") && requireNamespace("readxl", quietly = TRUE)) {
    ann <- as.data.frame(readxl::read_excel("database.xlsx"), stringsAsFactors = FALSE)
    if (!"GO_BP" %in% colnames(ann)) ann$GO_BP <- "-"
    if (!"GO_CC" %in% colnames(ann)) ann$GO_CC <- "-"
    if (!"GO_MF" %in% colnames(ann)) ann$GO_MF <- "-"
    if (!"GO_Count" %in% colnames(ann)) ann$GO_Count <- 0L
    if (!"OrgDb_KEGG" %in% colnames(ann)) ann$OrgDb_KEGG <- ann$KEGG_Pathway
    return(list(
      annotation    = ann,
      org_metadata  = data.frame(
        name  = c("ORGANISM", "Taxonomy ID", "Db type", "CENTRALID"),
        value = c("Heterorhabditis bacteriophora", "37862", "OrgDb", "GID"),
        stringsAsFactors = FALSE
      ),
      t2g           = list(),
      go_term_names = setNames(character(0), character(0)),
      pfam_desc_map = setNames(character(0), character(0)),
      stats         = list(
        n_genes         = nrow(ann),
        n_cols          = ncol(ann),
        n_symbols       = sum(ann$Preferred_name != "-" & !is.na(ann$Preferred_name)),
        n_descriptions  = sum(ann$Description != "-" & !is.na(ann$Description)),
        n_go_genes      = sum(ann$GOs != "-" & !is.na(ann$GOs)),
        n_go_terms      = NA_integer_,
        n_go_bp_terms   = NA_integer_,
        n_go_cc_terms   = NA_integer_,
        n_go_mf_terms   = NA_integer_,
        n_kegg_genes    = sum(ann$KEGG_Pathway != "-" & !is.na(ann$KEGG_Pathway)),
        n_kegg_pathways = NA_integer_,
        n_pfam_genes    = sum(ann$PFAMs != "-" & !is.na(ann$PFAMs)),
        n_pfam_domains  = NA_integer_,
        n_metacyc_genes = sum(ann[["MetaCyc/Reactome-Pathway"]] != "-" & !is.na(ann[["MetaCyc/Reactome-Pathway"]])),
        n_ec_genes      = sum(ann$EC != "-" & !is.na(ann$EC)),
        n_cazy_genes    = sum(ann$CAZy != "-" & !is.na(ann$CAZy)),
        n_cog_genes     = sum(ann$COG_category != "-" & !is.na(ann$COG_category)),
        build_date      = as.character(Sys.Date())
      )
    ))
  }
  stop("Neither data/hba_annotation.rds nor database.xlsx was found in the app directory.")
}

hba_bundle      <- load_hba_bundle()
annotation_data <- hba_bundle$annotation
hba_stats       <- hba_bundle$stats
t2g_list        <- hba_bundle$t2g
go_term_names   <- hba_bundle$go_term_names
pfam_desc_map   <- hba_bundle$pfam_desc_map

# Standard NCBI COG functional category descriptions
COG_LABELS <- c(
  A = "A: RNA processing and modification",
  B = "B: Chromatin structure and dynamics",
  C = "C: Energy production and conversion",
  D = "D: Cell cycle control, cell division, chromosome partitioning",
  E = "E: Amino acid transport and metabolism",
  F = "F: Nucleotide transport and metabolism",
  G = "G: Carbohydrate transport and metabolism",
  H = "H: Coenzyme transport and metabolism",
  I = "I: Lipid transport and metabolism",
  J = "J: Translation, ribosomal structure and biogenesis",
  K = "K: Transcription",
  L = "L: Replication, recombination and repair",
  M = "M: Cell wall/membrane/envelope biogenesis",
  N = "N: Cell motility",
  O = "O: Post-translational modification, protein turnover, chaperones",
  P = "P: Inorganic ion transport and metabolism",
  Q = "Q: Secondary metabolites biosynthesis, transport and catabolism",
  S = "S: Function unknown",
  T = "T: Signal transduction mechanisms",
  U = "U: Intracellular trafficking, secretion, and vesicular transport",
  V = "V: Defense mechanisms",
  W = "W: Extracellular structures",
  Y = "Y: Nuclear structure",
  Z = "Z: Cytoskeleton"
)

# Column presets for clean browsing
COLUMN_PRESETS <- list(
  core = c("GID", "Preferred_name", "Description", "COG_category", "GO_Count",
           "KEGG_Pathway", "PFAMs", "PFAM Accession", "evalue", "score"),
  go_func = c("GID", "Preferred_name", "Description", "COG_category",
              "GO_Count", "GO_BP", "GO_CC", "GO_MF", "GOs"),
  kegg_met = c("GID", "Preferred_name", "Description", "EC", "KEGG_ko",
               "KEGG_Pathway", "OrgDb_KEGG", "KEGG_Module", "KEGG_Reaction",
               "BRITE", "KEGG_TC", "CAZy", "MetaCyc/Reactome-Pathway"),
  pfam_dom = c("GID", "Preferred_name", "Description", "PFAMs", "PFam",
               "PFAM Accession", "seed_ortholog", "evalue", "score",
               "max_annot_lvl", "eggNOG_OGs"),
  all = colnames(annotation_data)
)

# Example gene list (40 kinase / signal transduction & transporter genes for instant ORA demo)
example_gids <- head(
  annotation_data$GID[
    grepl("T", annotation_data$COG_category) &
      annotation_data$GO_Count >= 10L &
      annotation_data$OrgDb_KEGG != "-"
  ],
  40L
)
if (length(example_gids) < 10L) example_gids <- head(annotation_data$GID, 30L)

# Helper: fast hypergeometric Over-Representation Analysis (ORA)
run_ora_enrichment <- function(query_gids, term2gene, universe_gids,
                               term_names = NULL, min_gs = 3L,
                               p_cutoff = 0.05) {
  query_gids <- intersect(unique(query_gids), universe_gids)
  k <- length(query_gids)
  N <- length(universe_gids)
  if (k == 0L || length(term2gene) == 0L) return(data.frame())

  rows <- lapply(names(term2gene), function(term) {
    term_genes <- intersect(term2gene[[term]], universe_gids)
    M <- length(term_genes)
    if (M < min_gs) return(NULL)
    overlap <- intersect(query_gids, term_genes)
    x <- length(overlap)
    if (x == 0L) return(NULL)
    pval <- stats::phyper(x - 1L, M, N - M, k, lower.tail = FALSE)
    data.frame(
      TermID     = term,
      Count      = x,
      QueryTotal = k,
      BgCount    = M,
      BgTotal    = N,
      GeneRatio  = x / k,
      BgRatio    = M / N,
      FoldEnrich = (x / k) / (M / N),
      PValue     = pval,
      Genes      = paste(overlap, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, rows)
  if (is.null(res) || nrow(res) == 0L) return(data.frame())

  res$Padj <- stats::p.adjust(res$PValue, method = "BH")
  res$Description <- if (!is.null(term_names) && length(term_names) > 0L) {
    ifelse(res$TermID %in% names(term_names), unname(term_names[res$TermID]), res$TermID)
  } else {
    res$TermID
  }
  res <- res[order(res$Padj, res$PValue, -res$Count), , drop = FALSE]
  res <- res[res$Padj <= p_cutoff, , drop = FALSE]
  rownames(res) <- NULL
  res[, c("TermID", "Description", "Count", "GeneRatio", "BgRatio",
          "FoldEnrich", "PValue", "Padj", "Genes"), drop = FALSE]
}

# Helper: render clickable HTML badge pills
make_link_badges <- function(raw_str, split_pat = ",", url_fn = NULL,
                             label_fn = NULL, badge_cls = "bg-primary-subtle text-primary border border-primary-subtle",
                             max_items = 25L) {
  if (is.null(raw_str) || is.na(raw_str) || raw_str == "-" || !nzchar(trimws(raw_str))) {
    return(tags$span(class = "text-muted fst-italic", "Not annotated"))
  }
  items <- unique(trimws(unlist(strsplit(as.character(raw_str), split_pat))))
  items <- items[nzchar(items) & items != "-" & items != "NA"]
  if (length(items) == 0L) {
    return(tags$span(class = "text-muted fst-italic", "Not annotated"))
  }
  overflow <- max(0L, length(items) - max_items)
  show_items <- head(items, max_items)
  pills <- lapply(show_items, function(it) {
    lbl <- if (!is.null(label_fn)) label_fn(it) else it
    if (!is.null(url_fn)) {
      href <- url_fn(it)
      tags$a(
        href = href, target = "_blank", rel = "noopener noreferrer",
        class = paste("badge text-decoration-none me-1 mb-1", badge_cls),
        style = "font-weight: 500; font-size: 0.78rem;",
        lbl
      )
    } else {
      tags$span(
        class = paste("badge me-1 mb-1", badge_cls),
        style = "font-weight: 500; font-size: 0.78rem;",
        lbl
      )
    }
  })
  if (overflow > 0L) {
    pills <- c(pills, list(tags$span(
      class = "badge bg-secondary-subtle text-secondary me-1 mb-1",
      sprintf("+%d more", overflow)
    )))
  }
  tagList(pills)
}

decode_cog <- function(cog_str) {
  if (is.na(cog_str) || cog_str == "-") return("Unassigned")
  chars <- unique(unlist(strsplit(gsub("[^A-Z]", "", cog_str), "")))
  if (length(chars) == 0L) return(cog_str)
  paste(ifelse(chars %in% names(COG_LABELS), unname(COG_LABELS[chars]), chars), collapse = " | ")
}

pct_fmt <- function(n, tot) sprintf("%.1f%%", 100 * n / max(1, tot))

# Precompute genome-wide analytics tables once at startup for instant plotting
build_analytics_cache <- function(ann, t2g, stats, pfam_map) {
  # 1. Coverage summary
  cov_df <- data.frame(
    Database = c(
      "Functional Description", "COG Category", "GO Terms (OrgDb)",
      "PFAM Domains", "KEGG Orthology (KO)", "MetaCyc / Reactome",
      "Gene Symbol (SYMBOL)", "KEGG Pathways", "EC Numbers", "CAZy Enzymes"
    ),
    Genes = c(
      stats$n_descriptions, stats$n_cog_genes, stats$n_go_genes,
      stats$n_pfam_genes, sum(ann$KEGG_ko != "-"), stats$n_metacyc_genes,
      stats$n_symbols, stats$n_kegg_genes, stats$n_ec_genes, stats$n_cazy_genes
    ),
    stringsAsFactors = FALSE
  )
  cov_df$Percent <- round(100 * cov_df$Genes / stats$n_genes, 1)
  cov_df <- cov_df[order(cov_df$Genes), ]

  # 2. COG distribution
  cog_counts <- if (!is.null(t2g$COG)) {
    vapply(t2g$COG, length, integer(1))
  } else {
    integer(0)
  }
  cog_df <- if (length(cog_counts) > 0L) {
    d <- data.frame(
      Code  = names(cog_counts),
      Label = ifelse(names(cog_counts) %in% names(COG_LABELS),
                     unname(COG_LABELS[names(cog_counts)]), names(cog_counts)),
      Genes = as.integer(cog_counts),
      stringsAsFactors = FALSE
    )
    d[order(d$Genes), ]
  } else {
    data.frame(Code = character(0), Label = character(0), Genes = integer(0))
  }

  # 3. Top 20 PFAM domains
  pfam_counts <- if (!is.null(t2g$PFAM)) vapply(t2g$PFAM, length, integer(1)) else integer(0)
  pfam_df <- if (length(pfam_counts) > 0L) {
    top_p <- sort(pfam_counts, decreasing = TRUE)[seq_len(min(20L, length(pfam_counts)))]
    d <- data.frame(
      Domain = names(top_p),
      Desc   = ifelse(names(top_p) %in% names(pfam_map), unname(pfam_map[names(top_p)]), names(top_p)),
      Genes  = as.integer(top_p),
      stringsAsFactors = FALSE
    )
    d[order(d$Genes), ]
  } else {
    data.frame(Domain = character(0), Desc = character(0), Genes = integer(0))
  }

  # 4. Top 20 KEGG pathways
  kegg_counts <- if (!is.null(t2g$KEGG)) vapply(t2g$KEGG, length, integer(1)) else integer(0)
  kegg_df <- if (length(kegg_counts) > 0L) {
    top_k <- sort(kegg_counts, decreasing = TRUE)[seq_len(min(20L, length(kegg_counts)))]
    d <- data.frame(
      Pathway = names(top_k),
      Genes   = as.integer(top_k),
      stringsAsFactors = FALSE
    )
    d[order(d$Genes), ]
  } else {
    data.frame(Pathway = character(0), Genes = integer(0))
  }

  # 5. GO ontology breakdown
  go_ont_df <- data.frame(
    Ontology = c("Biological Process (BP)", "Cellular Component (CC)", "Molecular Function (MF)"),
    Genes    = c(sum(ann$GO_BP != "-"), sum(ann$GO_CC != "-"), sum(ann$GO_MF != "-")),
    Terms    = c(stats$n_go_bp_terms %||% 0L, stats$n_go_cc_terms %||% 0L, stats$n_go_mf_terms %||% 0L),
    stringsAsFactors = FALSE
  )

  list(cov = cov_df, cog = cog_df, pfam = pfam_df, kegg = kegg_df, go_ont = go_ont_df)
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0L && !all(is.na(a))) a else b
analytics_cache <- build_analytics_cache(annotation_data, t2g_list, hba_stats, pfam_desc_map)

# -----------------------------------------------------------------------------
# 2. UI Definition (Modern Bootstrap 5 via bslib)
# -----------------------------------------------------------------------------
app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#0F3D5E",
  secondary = "#475569",
  success = "#0F766E",
  info = "#0284C7",
  "font-size-base" = "0.93rem",
  "border-radius" = "0.5rem",
  "navbar-bg" = "#0F3D5E"
)

ui <- page_navbar(
  id = "main_nav",
  title = tags$span(
    class = "d-flex align-items-center gap-2 fw-bold",
    bs_icon("database-fill-gear", size = "1.25em"),
    "HbaDb",
    tags$span(
      class = "badge bg-info text-white fw-normal d-none d-md-inline",
      style = "font-size: 0.7rem;",
      "org.Hbacteriophora.eg.db"
    )
  ),
  window_title = "HbaDb | Heterorhabditis bacteriophora Genome Annotation Database",
  theme = app_theme,
  fillable = FALSE,
  header = tags$head(
    tags$style(HTML("
      .hero-banner {
        background: linear-gradient(135deg, #0F3D5E 0%, #155E75 60%, #0E7490 100%);
        color: #FFFFFF;
        border-radius: 12px;
        padding: 28px 32px;
        margin-bottom: 20px;
        box-shadow: 0 6px 20px rgba(15, 61, 94, 0.16);
      }
      .hero-banner a.badge-link {
        color: #FFFFFF;
        text-decoration: none;
        background: rgba(255, 255, 255, 0.16);
        border: 1px solid rgba(255, 255, 255, 0.35);
        padding: 4px 12px;
        border-radius: 999px;
        font-size: 0.8rem;
        display: inline-flex;
        align-items: center;
        gap: 6px;
        transition: background 0.15s ease;
      }
      .hero-banner a.badge-link:hover {
        background: rgba(255, 255, 255, 0.28);
      }
      .gene-detail-card {
        border-left: 4px solid #0284C7 !important;
      }
      .code-block {
        background: #0F172A;
        color: #E2E8F0;
        padding: 14px 16px;
        border-radius: 8px;
        font-family: 'Cascadia Code', 'Fira Code', Consolas, monospace;
        font-size: 0.84rem;
        overflow-x: auto;
        white-space: pre;
      }
      .author-card {
        border-top: 3px solid #0F3D5E;
      }
    "))
  ),

  # ===========================================================================
  # TAB 1: HOME / OVERVIEW DASHBOARD
  # ===========================================================================
  nav_panel(
    title = "Home",
    value = "home",
    icon = bs_icon("house-door-fill"),

    div(
      class = "container-fluid py-3",

      # Hero Banner
      div(
        class = "hero-banner",
        div(
          class = "d-flex flex-wrap justify-content-between align-items-start gap-3",
          div(
            style = "max-width: 780px;",
            tags$h2(
              class = "fw-bold mb-2",
              tags$em("Heterorhabditis bacteriophora"), " Functional Genome Database (HbaDb)"
            ),
            tags$p(
              class = "mb-3 opacity-90",
              style = "font-size: 1.0rem;",
              "Interactive web portal and enrichment workbench powered by Bioconductor's ",
              tags$strong("org.Hbacteriophora.eg.db"),
              " annotation package — integrating eggNOG orthology, ontology-resolved Gene Ontology (BP/CC/MF), ",
              "KEGG pathways, PFAM protein domains, EC numbers, CAZy enzymes, and MetaCyc/Reactome pathways across ",
              tags$strong(format(hba_stats$n_genes, big.mark = ",")), " genes."
            ),
            div(
              class = "d-flex flex-wrap gap-2",
              tags$a(
                class = "badge-link",
                href = "https://doi.org/10.18129/B9.bioc.org.Hbacteriophora.eg.db",
                target = "_blank",
                bs_icon("box-seam"), "Bioconductor DOI: 10.18129/B9.bioc.org.Hbacteriophora.eg.db"
              ),
              tags$a(
                class = "badge-link",
                href = "https://github.com/kabilanbio/org.Hbacteriophora.eg.db",
                target = "_blank",
                bs_icon("github"), "org.Hbacteriophora.eg.db"
              ),
              tags$span(
                class = "badge-link",
                bs_icon("tag"), "NCBI TaxID: 37862 | AnnotationHub: AH121717"
              )
            )
          ),
          div(
            class = "d-flex flex-column gap-2",
            actionButton("go_explorer", "Explore Genes & Annotations", icon = icon("search"), class = "btn btn-light fw-semibold text-primary"),
            actionButton("go_enrichment", "Batch Lookup & Enrichment", icon = icon("layer-group"), class = "btn btn-outline-light fw-semibold"),
            actionButton("go_analytics", "Genome-Wide Analytics", icon = icon("chart-bar"), class = "btn btn-outline-light fw-semibold")
          )
        )
      ),

      # Live KPI Value Boxes
      layout_columns(
        col_widths = c(2, 2, 2, 2, 2, 2),
        fill = FALSE,
        value_box(
          title = "Annotated Genes",
          value = format(hba_stats$n_genes, big.mark = ","),
          showcase = bs_icon("fingerprint"),
          theme = "primary",
          tags$small("Custom GIDs (Hba_*)")
        ),
        value_box(
          title = "Gene Symbols",
          value = format(hba_stats$n_symbols, big.mark = ","),
          showcase = bs_icon("tags"),
          theme = "info",
          tags$small(paste0(pct_fmt(hba_stats$n_symbols, hba_stats$n_genes), " named symbols"))
        ),
        value_box(
          title = "GO-Annotated",
          value = format(hba_stats$n_go_genes, big.mark = ","),
          showcase = bs_icon("diagram-3"),
          theme = "success",
          tags$small(paste0(format(hba_stats$n_go_terms, big.mark = ","), " unique GO terms"))
        ),
        value_box(
          title = "KEGG Pathways",
          value = format(hba_stats$n_kegg_genes, big.mark = ","),
          showcase = bs_icon("bezier2"),
          theme = "secondary",
          tags$small(paste0(hba_stats$n_kegg_pathways, " pathways mapped"))
        ),
        value_box(
          title = "PFAM Domains",
          value = format(hba_stats$n_pfam_genes, big.mark = ","),
          showcase = bs_icon("bricks"),
          theme = "primary",
          tags$small(paste0(format(hba_stats$n_pfam_domains, big.mark = ","), " unique domains"))
        ),
        value_box(
          title = "MetaCyc / Reactome",
          value = format(hba_stats$n_metacyc_genes, big.mark = ","),
          showcase = bs_icon("share"),
          theme = "info",
          tags$small(paste0(pct_fmt(hba_stats$n_metacyc_genes, hba_stats$n_genes), " pathway coverage"))
        )
      ),

      # About + Quick Guide + Citation
      layout_columns(
        col_widths = c(7, 5),
        class = "mt-3",

        card(
          card_header(class = "fw-bold", bs_icon("info-circle", class = "me-2"), "About HbaDb & Bioconductor Integration"),
          card_body(
            tags$p(
              tags$em("Heterorhabditis bacteriophora"),
              " is an entomopathogenic nematode widely studied as a biological control agent against insect pests ",
              "and as a model for symbiosis with ", tags$em("Photorhabdus"), " bacteria. ",
              "HbaDb unifies genome-wide structural and functional annotations from ",
              tags$code("org.Hbacteriophora.eg.db"), " (Bioconductor OrgDb) with rich eggNOG-mapper, ",
              "WormBase ParaSite, InterPro/PFAM, KEGG, EC, CAZy, and MetaCyc/Reactome annotations."
            ),
            tags$h6(class = "fw-bold mt-3", "Key Capabilities"),
            tags$ul(
              tags$li(tags$strong("Interactive Gene Explorer:"), " Search across all 29 annotation fields, filter by COG category, e-value, or database presence, and click any row to inspect clickable QuickGO, KEGG, PFAM/InterPro, EC, and WormBase ParaSite links."),
              tags$li(tags$strong("Batch Gene Lookup & ORA Enrichment:"), " Paste or upload a list of differentially expressed ", tags$code("Hba_*"), " genes or symbols to extract full annotations and run instant hypergeometric enrichment (BH-FDR) across GO (BP/CC/MF), KEGG, PFAM, and COG."),
              tags$li(tags$strong("Genome-Wide Analytics:"), " Interactive visualizations summarizing annotation coverage, COG functional distributions, GO ontology breakdown, and top PFAM/KEGG families."),
              tags$li(tags$strong("Reproducible Bioconductor R Code:"), " Auto-generated R scripts using ", tags$code("org.Hbacteriophora.eg.db"), ", ", tags$code("AnnotationHub"), ", and ", tags$code("clusterProfiler"), ".")
            )
          )
        ),

        card(
          card_header(class = "fw-bold", bs_icon("journal-bookmark", class = "me-2"), "Citation & Provenance"),
          card_body(
            tags$p("If you use HbaDb or ", tags$code("org.Hbacteriophora.eg.db"), " in your research, please cite:"),
            tags$blockquote(
              class = "p-3 rounded border-start border-4 border-primary bg-light-subtle",
              style = "font-size: 0.88rem;",
              "Johnson LW, Bhat CG, Dillman A, Somvanshi VS, Kabilan S (2025). ",
              tags$em("org.Hbacteriophora.eg.db: Genome-wide Annotation for Heterorhabditis bacteriophora."),
              " R/Bioconductor package. DOI: ",
              tags$a(
                href = "https://doi.org/10.18129/B9.bioc.org.Hbacteriophora.eg.db",
                target = "_blank",
                "10.18129/B9.bioc.org.Hbacteriophora.eg.db"
              )
            ),
            tags$hr(),
            div(
              class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
              div(
                tags$small(class = "text-muted d-block", paste("Annotation build:", hba_stats$build_date)),
                tags$small(class = "text-muted d-block", paste("Total fields:", hba_stats$n_cols, "| OrgDb schema: 2.1"))
              ),
              tags$div(
                class = "visitor-map",
                tags$a(
                  href = "https://mapmyvisitors.com/web/1bxbm",
                  target = "_blank",
                  title = "Visitor tracker",
                  tags$img(
                    src = "https://mapmyvisitors.com/map.png?d=JjosJ6JAVG9bsS0dycAv8Wp7IEc9rPmDuZmbxnj2dyA&cl=ffffff",
                    style = "border: none; width: 135px; border-radius: 6px;"
                  )
                )
              )
            )
          )
        )
      )
    )
  ),

  # ===========================================================================
  # TAB 2: GENE & ANNOTATION EXPLORER
  # ===========================================================================
  nav_panel(
    title = "Gene Explorer",
    value = "explorer",
    icon = bs_icon("table"),

    layout_sidebar(
      sidebar = sidebar(
        title = "Search & Filter",
        width = 330,

        textInput("search", "Keyword Search (All Columns):", placeholder = "e.g., Hba_00004, kinase, GO:0035556, PF00211"),
        checkboxInput("exact_match", "Exact match (case-insensitive)", FALSE),

        selectInput(
          "column_select", "Column-Specific Filter:",
          choices = c("All Columns", colnames(annotation_data)),
          selected = "All Columns"
        ),
        conditionalPanel(
          condition = "input.column_select != 'All Columns'",
          textInput("column_value", "Filter Value:", placeholder = "Enter value or substring...")
        ),

        selectizeInput(
          "cog_filter", "COG Functional Category:",
          choices = c("All Categories" = "ALL", setNames(names(COG_LABELS), unname(COG_LABELS))),
          selected = "ALL",
          multiple = TRUE
        ),

        selectInput(
          "evalue_filter", "eggNOG E-value Threshold:",
          choices = c(
            "All e-values" = "all",
            "<= 1e-5"      = "1e-5",
            "<= 1e-10"     = "1e-10",
            "<= 1e-20"     = "1e-20",
            "<= 1e-50"     = "1e-50"
          ),
          selected = "all"
        ),

        checkboxGroupInput(
          "presence_filter", "Require Annotation In:",
          choices = c(
            "Gene Symbol"        = "symbol",
            "GO Terms"           = "go",
            "KEGG Pathway"       = "kegg",
            "PFAM Domain"        = "pfam",
            "EC Number"          = "ec",
            "CAZy Enzyme"        = "cazy",
            "MetaCyc / Reactome" = "metacyc"
          )
        ),

        div(
          class = "d-grid gap-2 mt-2",
          actionButton("search_btn", "Apply Filters", icon = icon("filter"), class = "btn-primary fw-semibold"),
          actionButton("reset_btn", "Reset All Filters", icon = icon("rotate-left"), class = "btn-outline-secondary btn-sm")
        ),

        tags$hr(),
        div(
          class = "d-grid gap-2",
          downloadButton("downloadDataCSV", "Download Filtered (.csv)", class = "btn-success btn-sm"),
          downloadButton("downloadDataXLSX", "Download Filtered (.xlsx)", class = "btn-outline-success btn-sm")
        )
      ),

      # Main content: Gene Inspector + DataTable
      uiOutput("gene_detail_card"),

      card(
        full_screen = TRUE,
        card_header(
          div(
            class = "d-flex flex-wrap justify-content-between align-items-center gap-2",
            div(
              class = "d-flex align-items-center gap-2",
              tags$span(class = "fw-bold", bs_icon("grid-3x3-gap", class = "me-1"), "Annotation Table"),
              uiOutput("result_summary_badge", inline = TRUE)
            ),
            div(
              class = "d-flex align-items-center gap-2",
              selectInput(
                "col_preset", NULL,
                choices = c(
                  "Core Overview (10 cols)"       = "core",
                  "GO & Functional (9 cols)"      = "go_func",
                  "KEGG & Metabolic (13 cols)"    = "kegg_met",
                  "PFAM & Orthology (11 cols)"    = "pfam_dom",
                  "All 29 Columns"                = "all"
                ),
                selected = "core",
                width = "230px"
              )
            )
          )
        ),
        card_body(
          tags$small(class = "text-muted d-block mb-2", bs_icon("cursor"), " Tip: Click any row to open the interactive Gene Detail Card with clickable QuickGO, KEGG, PFAM/InterPro, EC, and WormBase ParaSite links."),
          DTOutput("annotationTable")
        )
      )
    )
  ),

  # ===========================================================================
  # TAB 3: BATCH GENE LOOKUP & FUNCTIONAL ENRICHMENT (ORA)
  # ===========================================================================
  nav_panel(
    title = "Batch Lookup & Enrichment",
    value = "enrichment",
    icon = bs_icon("stars"),

    layout_sidebar(
      sidebar = sidebar(
        title = "Gene List & ORA Settings",
        width = 340,

        textAreaInput(
          "batch_genes_text",
          "Paste Gene IDs (GID) or Gene Symbols:",
          rows = 7,
          placeholder = "One per line, comma, or space separated:\nHba_00001\nHba_00004\nHba_00009\nwht-7"
        ),
        fileInput("batch_genes_file", "Or Upload .txt / .csv Gene List:", accept = c(".txt", ".csv", ".tsv")),

        div(
          class = "d-flex gap-2 mb-3",
          actionButton("btn_load_example_genes", "Load Example (40 Genes)", class = "btn-outline-info btn-sm flex-fill"),
          actionButton("btn_clear_batch", "Clear", class = "btn-outline-secondary btn-sm")
        ),

        tags$hr(),
        selectInput(
          "ora_db", "Enrichment Database:",
          choices = c(
            "GO: Biological Process (BP)" = "GO_BP",
            "GO: Molecular Function (MF)" = "GO_MF",
            "GO: Cellular Component (CC)" = "GO_CC",
            "KEGG Pathways"               = "KEGG",
            "PFAM Protein Domains"        = "PFAM",
            "COG Functional Categories"   = "COG"
          ),
          selected = "GO_BP"
        ),
        numericInput("ora_fdr", "BH-Adjusted FDR Cutoff:", value = 0.05, min = 0.0001, max = 1, step = 0.01),
        numericInput("ora_min_gs", "Min Background Term Size:", value = 3, min = 1, max = 100, step = 1),
        sliderInput("ora_top_n", "Top Enriched Terms to Plot:", min = 5, max = 30, value = 15, step = 1),

        div(
          class = "d-grid",
          actionButton("btn_run_batch", "Run Batch Lookup & ORA", icon = icon("play"), class = "btn-primary fw-semibold")
        )
      ),

      uiOutput("batch_match_summary"),

      navset_card_tab(
        full_screen = TRUE,
        nav_panel(
          "Enrichment Dotplot",
          card_body(
            plotlyOutput("ora_dotplot", height = "520px")
          )
        ),
        nav_panel(
          "Enrichment Table",
          card_body(
            div(
              class = "mb-2 d-flex gap-2",
              downloadButton("download_ora_csv", "Download Enrichment (.csv)", class = "btn-success btn-sm"),
              downloadButton("download_ora_xlsx", "Download Enrichment (.xlsx)", class = "btn-outline-success btn-sm")
            ),
            DTOutput("ora_table")
          )
        ),
        nav_panel(
          "Matched Gene Annotations",
          card_body(
            div(
              class = "mb-2 d-flex gap-2",
              downloadButton("download_batch_csv", "Download Matched Genes (.csv)", class = "btn-primary btn-sm"),
              downloadButton("download_batch_xlsx", "Download Matched Genes (.xlsx)", class = "btn-outline-primary btn-sm")
            ),
            DTOutput("batch_genes_table")
          )
        )
      )
    )
  ),

  # ===========================================================================
  # TAB 4: GENOME-WIDE ANALYTICS
  # ===========================================================================
  nav_panel(
    title = "Genome Analytics",
    value = "analytics",
    icon = bs_icon("bar-chart-line-fill"),

    div(
      class = "container-fluid py-3",
      layout_columns(
        col_widths = c(6, 6),
        card(
          full_screen = TRUE,
          card_header(class = "fw-bold", "Annotation Coverage Across Functional Databases"),
          card_body(plotlyOutput("plot_coverage", height = "390px"))
        ),
        card(
          full_screen = TRUE,
          card_header(class = "fw-bold", "Gene Ontology (GO) Breakdown (org.Hbacteriophora.eg.db)"),
          card_body(plotlyOutput("plot_go_ont", height = "390px"))
        )
      ),
      layout_columns(
        col_widths = c(12),
        class = "mt-3",
        card(
          full_screen = TRUE,
          card_header(class = "fw-bold", "COG Functional Category Distribution (All 11,466 Genes)"),
          card_body(plotlyOutput("plot_cog_dist", height = "520px"))
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        class = "mt-3",
        card(
          full_screen = TRUE,
          card_header(class = "fw-bold", "Top 20 PFAM Protein Domains"),
          card_body(plotlyOutput("plot_top_pfam", height = "460px"))
        ),
        card(
          full_screen = TRUE,
          card_header(class = "fw-bold", "Top 20 KEGG Pathways"),
          card_body(plotlyOutput("plot_top_kegg", height = "460px"))
        )
      )
    )
  ),

  # ===========================================================================
  # TAB 5: BIOCONDUCTOR PACKAGE & REPRODUCIBLE R CODE
  # ===========================================================================
  nav_panel(
    title = "Bioconductor & R Code",
    value = "bioc_r",
    icon = bs_icon("code-slash"),

    div(
      class = "container-fluid py-3",
      layout_columns(
        col_widths = c(7, 5),

        card(
          card_header(class = "fw-bold", bs_icon("terminal", class = "me-2"), "Using org.Hbacteriophora.eg.db in R / Bioconductor"),
          card_body(
            tags$p(
              "All core gene identifiers (", tags$code("GID"), "), gene symbols (", tags$code("SYMBOL"), "), ",
              "gene names (", tags$code("GENENAME"), "), Gene Ontology terms (", tags$code("GO"), ", ", tags$code("ONTOLOGY"), "), ",
              "and KEGG pathways (", tags$code("KEGGPATHWAY"), ") are distributed programmatically through the Bioconductor ",
              tags$strong("org.Hbacteriophora.eg.db"), " package."
            ),
            tags$h6(class = "fw-bold mt-3", "1. Install from Bioconductor & Retrieve via AnnotationHub"),
            tags$div(
              class = "code-block",
'if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("org.Hbacteriophora.eg.db", "AnnotationHub"))

library(AnnotationHub)
library(AnnotationDbi)

ah <- AnnotationHub()
q  <- query(ah, c("org.Hbacteriophora.eg.db", "Heterorhabditis bacteriophora"))
orgdb <- q[[1]] # AH121717

keytypes(orgdb)
columns(orgdb)'
            ),
            tags$h6(class = "fw-bold mt-3", "2. Dynamic R Script (Synced With Your Active Genes)"),
            tags$p(class = "text-muted small", "Select a row in Gene Explorer or run a Batch Lookup to auto-populate your GIDs below:"),
            uiOutput("dynamic_r_code_ui"),
            downloadButton("download_r_script", "Download Reproducible R Script (.R)", class = "btn-primary btn-sm mt-2")
          )
        ),

        card(
          card_header(class = "fw-bold", bs_icon("box2-heart", class = "me-2"), "OrgDb Metadata & Links"),
          card_body(
            tags$ul(
              class = "list-unstyled mb-3",
              tags$li(class = "mb-2", tags$strong("Bioconductor DOI: "),
                      tags$a(href = "https://doi.org/10.18129/B9.bioc.org.Hbacteriophora.eg.db", target = "_blank", "10.18129/B9.bioc.org.Hbacteriophora.eg.db")),
              tags$li(class = "mb-2", tags$strong("AnnotationHub ID: "), tags$code("AH121717")),
              tags$li(class = "mb-2", tags$strong("Package GitHub: "),
                      tags$a(href = "https://github.com/kabilanbio/org.Hbacteriophora.eg.db", target = "_blank", "github.com/kabilanbio/org.Hbacteriophora.eg.db")),
              tags$li(class = "mb-2", tags$strong("HbaDb Web App GitHub: "),
                      tags$a(href = "https://github.com/kabilanbio/HbaDb", target = "_blank", "github.com/kabilanbio/HbaDb")),
              tags$li(class = "mb-2", tags$strong("License: "), tags$code("Artistic-2.0"))
            ),
            tags$h6(class = "fw-bold", "SQLite OrgDb Internal Metadata"),
            tableOutput("orgdb_meta_tbl")
          )
        )
      )
    )
  ),

  # ===========================================================================
  # TAB 6: TEAM, CITATION & CONTACT
  # ===========================================================================
  nav_panel(
    title = "Team & Contact",
    value = "contact",
    icon = bs_icon("people-fill"),

    div(
      class = "container-fluid py-3",
      card(
        card_header(class = "fw-bold bg-primary text-white", bs_icon("people", class = "me-2"), "Development Team & Contributors"),
        card_body(
          layout_columns(
            col_widths = c(4, 4, 4),
            card(
              class = "author-card",
              card_body(
                tags$h5(class = "fw-bold mb-1", "Lalson Wesly Johnson"),
                tags$span(class = "badge bg-primary-subtle text-primary mb-2", "Author"),
                tags$p(class = "small mb-1", bs_icon("envelope", class = "me-1"), "WESLYJOHNZEN@gmail.com"),
                tags$p(class = "small mb-0", bs_icon("person-badge", class = "me-1"),
                       tags$a(href = "https://orcid.org/0000-0001-9121-1927", target = "_blank", "ORCID: 0000-0001-9121-1927"))
              )
            ),
            card(
              class = "author-card",
              card_body(
                tags$h5(class = "fw-bold mb-1", "Kabilan S"),
                tags$span(class = "badge bg-success-subtle text-success mb-2", "Author & Maintainer"),
                tags$p(class = "small mb-1", bs_icon("envelope", class = "me-1"), "kabilan151414@gmail.com"),
                tags$p(class = "small mb-0", bs_icon("person-badge", class = "me-1"),
                       tags$a(href = "https://orcid.org/0009-0002-9422-0938", target = "_blank", "ORCID: 0009-0002-9422-0938"))
              )
            ),
            card(
              class = "author-card",
              card_body(
                tags$h5(class = "fw-bold mb-1", "Chaitra Ganapati Bhat"),
                tags$span(class = "badge bg-primary-subtle text-primary mb-2", "Author"),
                tags$p(class = "small mb-1", bs_icon("envelope", class = "me-1"), "bhatchaitragm@gmail.com"),
                tags$p(class = "small mb-0", bs_icon("person-badge", class = "me-1"),
                       tags$a(href = "https://orcid.org/0000-0003-1497-890X", target = "_blank", "ORCID: 0000-0003-1497-890X"))
              )
            ),
            card(
              class = "author-card",
              card_body(
                tags$h5(class = "fw-bold mb-1", "Adler Dillman"),
                tags$span(class = "badge bg-primary-subtle text-primary mb-2", "Author"),
                tags$p(class = "small mb-1", bs_icon("envelope", class = "me-1"), "adlerd@ucr.edu"),
                tags$p(class = "small mb-0", bs_icon("person-badge", class = "me-1"),
                       tags$a(href = "https://orcid.org/0000-0001-7171-4332", target = "_blank", "ORCID: 0000-0001-7171-4332"))
              )
            ),
            card(
              class = "author-card",
              card_body(
                tags$h5(class = "fw-bold mb-1", "Vishal Singh Somvanshi"),
                tags$span(class = "badge bg-primary-subtle text-primary mb-2", "Author"),
                tags$p(class = "small mb-0", bs_icon("envelope", class = "me-1"), "vishalsomvanshi@gmail.com")
              )
            ),
            card(
              class = "author-card",
              card_body(
                tags$h5(class = "fw-bold mb-1", "Affiliations"),
                tags$p(class = "small mb-1", tags$strong("ICAR–Indian Agricultural Research Institute (IARI)"), ", New Delhi, India"),
                tags$p(class = "small mb-0", tags$strong("University of California, Riverside"), ", California, USA")
              )
            )
          ),

          tags$hr(),
          tags$h5(class = "fw-bold", "BibTeX Citation"),
          tags$div(
            class = "code-block",
'@Manual{org.Hbacteriophora.eg.db,
  title  = {org.Hbacteriophora.eg.db: Genome-wide Annotation for Heterorhabditis bacteriophora},
  author = {Lalson Wesly Johnson and Chaitra Ganapati Bhat and Adler Dillman and Vishal Singh Somvanshi and Kabilan S},
  year   = {2025},
  note   = {Bioconductor OrgDb annotation package},
  doi    = {10.18129/B9.bioc.org.Hbacteriophora.eg.db},
  url    = {https://bioconductor.org/packages/org.Hbacteriophora.eg.db}
}'
          )
        )
      )
    )
  ),

  nav_spacer(),
  nav_item(input_dark_mode(id = "dark_mode", mode = "light")),

  footer = tags$footer(
    class = "py-3 mt-4 border-top text-center text-muted small",
    tags$span(
      "HbaDb | Powered by Bioconductor ",
      tags$a(href = "https://doi.org/10.18129/B9.bioc.org.Hbacteriophora.eg.db", target = "_blank", "org.Hbacteriophora.eg.db"),
      " | ICAR-IARI New Delhi & UC Riverside"
    )
  )
)

# -----------------------------------------------------------------------------
# 3. Server Logic
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  # Quick navigation buttons on Home tab
  observeEvent(input$go_explorer,   { nav_select("main_nav", "explorer") })
  observeEvent(input$go_enrichment, { nav_select("main_nav", "enrichment") })
  observeEvent(input$go_analytics,  { nav_select("main_nav", "analytics") })

  # Reactive filtered data for Tab 2
  filtered_data <- reactiveVal(annotation_data)

  apply_explorer_filters <- function() {
    df <- annotation_data

    # 1. Global search (case-insensitive, properly supporting exact vs substring)
    kw <- trimws(input$search %||% "")
    if (nzchar(kw)) {
      kw_low <- tolower(kw)
      exact <- isTRUE(input$exact_match)
      hit <- Reduce(`|`, lapply(df, function(col) {
        vals <- tolower(as.character(col))
        vals[is.na(vals)] <- ""
        if (exact) vals == kw_low else grepl(kw_low, vals, fixed = TRUE)
      }))
      df <- df[hit, , drop = FALSE]
    }

    # 2. Column-specific filter
    col_sel <- input$column_select %||% "All Columns"
    col_val <- trimws(input$column_value %||% "")
    if (!identical(col_sel, "All Columns") && nzchar(col_val) && col_sel %in% colnames(df)) {
      vals <- tolower(as.character(df[[col_sel]]))
      vals[is.na(vals)] <- ""
      hit <- if (isTRUE(input$exact_match)) {
        vals == tolower(col_val)
      } else {
        grepl(tolower(col_val), vals, fixed = TRUE)
      }
      df <- df[hit, , drop = FALSE]
    }

    # 3. COG functional category filter
    cogs <- input$cog_filter
    if (!is.null(cogs) && length(cogs) > 0L && !("ALL" %in% cogs)) {
      pat <- paste(cogs, collapse = "|")
      df <- df[grepl(pat, df$COG_category), , drop = FALSE]
    }

    # 4. E-value threshold
    ev_thr <- input$evalue_filter %||% "all"
    if (!identical(ev_thr, "all")) {
      thr_num <- suppressWarnings(as.numeric(ev_thr))
      if (is.finite(thr_num)) {
        df <- df[!is.na(df$evalue) & df$evalue <= thr_num, , drop = FALSE]
      }
    }

    # 5. Annotation presence checkboxes
    pres <- input$presence_filter
    if (length(pres) > 0L) {
      if ("symbol"  %in% pres) df <- df[df$Preferred_name != "-", , drop = FALSE]
      if ("go"      %in% pres) df <- df[df$GO_Count > 0L | df$GOs != "-", , drop = FALSE]
      if ("kegg"    %in% pres) df <- df[df$OrgDb_KEGG != "-" | df$KEGG_Pathway != "-", , drop = FALSE]
      if ("pfam"    %in% pres) df <- df[df$PFAMs != "-" | df$PFam != "-", , drop = FALSE]
      if ("ec"      %in% pres) df <- df[df$EC != "-", , drop = FALSE]
      if ("cazy"    %in% pres) df <- df[df$CAZy != "-", , drop = FALSE]
      if ("metacyc" %in% pres) df <- df[df[["MetaCyc/Reactome-Pathway"]] != "-", , drop = FALSE]
    }

    filtered_data(df)
  }

  observeEvent(input$search_btn, { apply_explorer_filters() })

  observeEvent(input$reset_btn, {
    updateTextInput(session, "search", value = "")
    updateCheckboxInput(session, "exact_match", value = FALSE)
    updateSelectInput(session, "column_select", selected = "All Columns")
    updateTextInput(session, "column_value", value = "")
    updateSelectizeInput(session, "cog_filter", selected = "ALL")
    updateSelectInput(session, "evalue_filter", selected = "all")
    updateCheckboxGroupInput(session, "presence_filter", selected = character(0))
    filtered_data(annotation_data)
  })

  output$result_summary_badge <- renderUI({
    n_cur <- nrow(filtered_data())
    n_tot <- nrow(annotation_data)
    cls <- if (n_cur == 0L) "bg-danger-subtle text-danger" else "bg-primary-subtle text-primary"
    tags$span(
      class = paste("badge border", cls),
      sprintf("%s / %s genes", format(n_cur, big.mark = ","), format(n_tot, big.mark = ","))
    )
  })

  output$annotationTable <- renderDT({
    df <- filtered_data()
    preset_key <- input$col_preset %||% "core"
    vis_cols <- intersect(COLUMN_PRESETS[[preset_key]] %||% colnames(df), colnames(df))
    hide_idx <- which(!(colnames(df) %in% vis_cols)) - 1L

    datatable(
      df,
      selection = "single",
      rownames = FALSE,
      filter = "top",
      extensions = c("Buttons"),
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = "Blfrtip",
        buttons = c("copy", "csv", "excel", "colvis"),
        lengthMenu = list(c(15, 50, 100, 250), c("15", "50", "100", "250")),
        columnDefs = if (length(hide_idx) > 0L) {
          list(list(visible = FALSE, targets = as.list(hide_idx)))
        } else {
          list()
        }
      ),
      class = "compact stripe hover"
    )
  }, server = TRUE)

  # Interactive Gene Detail Card on Row Click
  selected_gene_row <- reactive({
    sel <- input$annotationTable_rows_selected
    df <- filtered_data()
    if (is.null(sel) || length(sel) == 0L || sel[1] > nrow(df)) return(NULL)
    df[sel[1], , drop = FALSE]
  })

  output$gene_detail_card <- renderUI({
    row <- selected_gene_row()
    if (is.null(row)) return(NULL)

    gid <- row$GID[1]
    sym <- row$Preferred_name[1]
    desc <- row$Description[1]

    go_label_fn <- function(go_id) {
      nm <- go_term_names[go_id]
      if (!is.na(nm) && nzchar(nm)) sprintf("%s (%s)", go_id, nm) else go_id
    }

    card(
      class = "gene-detail-card mb-3",
      card_header(
        div(
          class = "d-flex flex-wrap justify-content-between align-items-center gap-2",
          div(
            tags$span(class = "fw-bold fs-5 text-primary me-2", gid),
            if (sym != "-") tags$span(class = "badge bg-info text-white me-2", paste("SYMBOL:", sym)),
            tags$span(class = "text-muted", if (desc != "-") desc else "Uncharacterized protein")
          ),
          tags$a(
            href = paste0("https://parasite.wormbase.org/Multi/Search/Results?q=", utils::URLencode(gid)),
            target = "_blank", rel = "noopener noreferrer",
            class = "btn btn-outline-primary btn-sm",
            bs_icon("box-arrow-up-right", class = "me-1"), "Search WormBase ParaSite"
          )
        )
      ),
      card_body(
        layout_columns(
          col_widths = c(4, 4, 4),
          div(
            tags$p(class = "mb-1 small", tags$strong("COG Category: "), decode_cog(row$COG_category[1])),
            tags$p(class = "mb-1 small", tags$strong("Seed Ortholog: "), tags$code(row$seed_ortholog[1])),
            tags$p(class = "mb-0 small", tags$strong("E-value / Score: "),
                   sprintf("%s / %s", format(row$evalue[1], scientific = TRUE, digits = 3), row$score[1]))
          ),
          div(
            tags$p(class = "mb-1 small fw-bold", "KEGG Pathways & Orthology"),
            make_link_badges(
              row$OrgDb_KEGG[1],
              url_fn = function(x) paste0("https://www.genome.jp/pathway/map", sub("^map|^ko", "", x)),
              label_fn = function(x) paste0("map", sub("^map|^ko", "", x)),
              badge_cls = "bg-success-subtle text-success border border-success-subtle"
            ),
            tags$div(class = "mt-1", make_link_badges(
              row$KEGG_ko[1],
              url_fn = function(x) paste0("https://www.genome.jp/entry/", sub("^ko:", "", x)),
              badge_cls = "bg-secondary-subtle text-secondary border"
            ))
          ),
          div(
            tags$p(class = "mb-1 small fw-bold", "PFAM Domains, EC & CAZy"),
            make_link_badges(
              row[["PFAM Accession"]][1],
              split_pat = "[,;]",
              url_fn = function(x) paste0("https://www.ebi.ac.uk/interpro/entry/pfam/", x),
              badge_cls = "bg-warning-subtle text-dark border border-warning-subtle"
            ),
            tags$div(class = "mt-1", make_link_badges(
              row$EC[1],
              url_fn = function(x) paste0("https://enzyme.expasy.org/EC/", x),
              label_fn = function(x) paste0("EC:", x),
              badge_cls = "bg-info-subtle text-info border border-info-subtle"
            )),
            if (row$CAZy[1] != "-") tags$div(class = "mt-1", make_link_badges(row$CAZy[1], badge_cls = "bg-dark-subtle text-dark"))
          )
        ),
        tags$hr(class = "my-2"),
        layout_columns(
          col_widths = c(4, 4, 4),
          div(
            tags$p(class = "mb-1 small fw-bold text-primary", "GO Biological Process (BP)"),
            make_link_badges(
              row$GO_BP[1],
              url_fn = function(x) paste0("https://www.ebi.ac.uk/QuickGO/term/", x),
              label_fn = go_label_fn,
              max_items = 10L
            )
          ),
          div(
            tags$p(class = "mb-1 small fw-bold text-success", "GO Molecular Function (MF)"),
            make_link_badges(
              row$GO_MF[1],
              url_fn = function(x) paste0("https://www.ebi.ac.uk/QuickGO/term/", x),
              label_fn = go_label_fn,
              badge_cls = "bg-success-subtle text-success border border-success-subtle",
              max_items = 10L
            )
          ),
          div(
            tags$p(class = "mb-1 small fw-bold text-info", "GO Cellular Component (CC)"),
            make_link_badges(
              row$GO_CC[1],
              url_fn = function(x) paste0("https://www.ebi.ac.uk/QuickGO/term/", x),
              label_fn = go_label_fn,
              badge_cls = "bg-info-subtle text-info border border-info-subtle",
              max_items = 10L
            )
          )
        )
      )
    )
  })

  # Download handlers (Tab 2)
  output$downloadDataCSV <- downloadHandler(
    filename = function() paste0("HbaDb_filtered_", Sys.Date(), ".csv"),
    content  = function(file) utils::write.csv(filtered_data(), file, row.names = FALSE)
  )
  output$downloadDataXLSX <- downloadHandler(
    filename = function() paste0("HbaDb_filtered_", Sys.Date(), ".xlsx"),
    content  = function(file) writexl::write_xlsx(list(HbaDb_Filtered = filtered_data()), path = file)
  )

  # ---------------------------------------------------------------------------
  # Tab 3: Batch Gene Lookup & Hypergeometric ORA Enrichment
  # ---------------------------------------------------------------------------
  observeEvent(input$btn_load_example_genes, {
    updateTextAreaInput(session, "batch_genes_text", value = paste(example_gids, collapse = "\n"))
  })
  observeEvent(input$btn_clear_batch, {
    updateTextAreaInput(session, "batch_genes_text", value = "")
  })

  parse_batch_input <- reactive({
    tokens <- character(0)
    txt <- input$batch_genes_text %||% ""
    if (nzchar(trimws(txt))) {
      tokens <- c(tokens, unlist(strsplit(txt, "[,;\\s]+", perl = TRUE)))
    }
    if (!is.null(input$batch_genes_file) && file.exists(input$batch_genes_file$datapath)) {
      lines <- readLines(input$batch_genes_file$datapath, warn = FALSE)
      tokens <- c(tokens, unlist(strsplit(paste(lines, collapse = "\n"), "[,;\\s]+", perl = TRUE)))
    }
    tokens <- unique(trimws(tokens))
    tokens[nzchar(tokens)]
  })

  batch_results <- eventReactive(list(input$btn_run_batch, input$ora_db, input$ora_fdr, input$ora_min_gs), {
    tokens <- parse_batch_input()
    if (length(tokens) == 0L) return(NULL)

    tok_low <- tolower(tokens)
    gid_low <- tolower(annotation_data$GID)
    sym_low <- tolower(annotation_data$Preferred_name)

    matched_idx <- which((gid_low %in% tok_low) | (sym_low != "-" & sym_low %in% tok_low))
    matched_df <- annotation_data[matched_idx, , drop = FALSE]

    matched_tokens <- unique(c(
      tokens[tok_low %in% gid_low],
      tokens[tok_low %in% sym_low[sym_low != "-"]]
    ))
    unmatched_tokens <- setdiff(tokens, matched_tokens)

    ora_db_key <- input$ora_db %||% "GO_BP"
    t2g <- t2g_list[[ora_db_key]]
    term_dict <- switch(
      ora_db_key,
      GO_BP = go_term_names,
      GO_CC = go_term_names,
      GO_MF = go_term_names,
      PFAM  = pfam_desc_map,
      COG   = COG_LABELS,
      NULL
    )

    ora_df <- run_ora_enrichment(
      query_gids    = matched_df$GID,
      term2gene     = t2g,
      universe_gids = annotation_data$GID,
      term_names    = term_dict,
      min_gs        = as.integer(input$ora_min_gs %||% 3L),
      p_cutoff      = as.numeric(input$ora_fdr %||% 0.05)
    )

    list(
      input_n    = length(tokens),
      matched_df = matched_df,
      unmatched  = unmatched_tokens,
      ora_df     = ora_df,
      db_key     = ora_db_key
    )
  }, ignoreNULL = FALSE)

  output$batch_match_summary <- renderUI({
    res <- batch_results()
    if (is.null(res)) {
      return(div(
        class = "alert alert-info mb-3",
        bs_icon("lightbulb", class = "me-2"),
        "Paste gene IDs (e.g., ", tags$code("Hba_00004"), ") or click ",
        tags$strong("Load Example (40 Genes)"), " then ", tags$strong("Run Batch Lookup & ORA"), "."
      ))
    }
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      class = "mb-3",
      value_box(title = "Input Queries", value = res$input_n, showcase = bs_icon("list-ul"), theme = "secondary"),
      value_box(title = "Matched Genes", value = nrow(res$matched_df), showcase = bs_icon("check-circle"), theme = "success"),
      value_box(title = "Unmatched", value = length(res$unmatched), showcase = bs_icon("question-circle"), theme = if (length(res$unmatched) > 0L) "warning" else "light"),
      value_box(title = "Enriched Terms", value = nrow(res$ora_df), showcase = bs_icon("stars"), theme = "primary")
    )
  })

  output$ora_dotplot <- renderPlotly({
    res <- batch_results()
    validate(need(!is.null(res), "Enter genes on the left and click 'Run Batch Lookup & ORA'."))
    df <- res$ora_df
    validate(need(nrow(df) > 0L, "No significantly enriched terms at the current FDR cutoff. Try relaxing the FDR cutoff or selecting another database."))

    top_n <- min(as.integer(input$ora_top_n %||% 15L), nrow(df))
    plot_df <- head(df, top_n)
    plot_df$ShortLabel <- ifelse(
      nchar(plot_df$Description) > 48,
      paste0(substr(plot_df$Description, 1, 45), "..."),
      plot_df$Description
    )
    plot_df$ShortLabel <- make.unique(plot_df$ShortLabel)
    plot_df$ShortLabel <- factor(plot_df$ShortLabel, levels = rev(plot_df$ShortLabel))
    plot_df$NegLog10FDR <- -log10(pmax(plot_df$Padj, 1e-300))
    plot_df$Hover <- sprintf(
      "<b>%s</b><br>%s<br>Count: %d | GeneRatio: %.3f<br>Fold Enrichment: %.2fx<br>P-value: %.2e | FDR: %.2e",
      plot_df$TermID, plot_df$Description, plot_df$Count,
      plot_df$GeneRatio, plot_df$FoldEnrich, plot_df$PValue, plot_df$Padj
    )

    p <- ggplot(plot_df, aes(x = GeneRatio, y = ShortLabel, size = Count, color = NegLog10FDR, text = Hover)) +
      geom_point(alpha = 0.9) +
      scale_color_gradient(low = "#38BDF8", high = "#0F3D5E", name = "-log10(FDR)") +
      scale_size_continuous(range = c(4, 12), name = "Gene Count") +
      labs(
        title = paste("Top Enriched Terms —", res$db_key),
        x = "Gene Ratio (Overlap / Query Genes)",
        y = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))

    ggplotly(p, tooltip = "text") %>% layout(margin = list(l = 220))
  })

  output$ora_table <- renderDT({
    res <- batch_results()
    req(res)
    df <- res$ora_df
    if (nrow(df) == 0L) return(datatable(data.frame(Message = "No enriched terms found at current FDR cutoff."), options = list(dom = "t"), rownames = FALSE))
    df$GeneRatio  <- round(df$GeneRatio, 4)
    df$BgRatio    <- round(df$BgRatio, 4)
    df$FoldEnrich <- round(df$FoldEnrich, 2)
    df$PValue     <- signif(df$PValue, 3)
    df$Padj       <- signif(df$Padj, 3)
    datatable(df, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE), class = "compact stripe hover")
  })

  output$batch_genes_table <- renderDT({
    res <- batch_results()
    req(res)
    vis <- intersect(COLUMN_PRESETS$core, colnames(res$matched_df))
    datatable(res$matched_df[, vis, drop = FALSE], rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE), class = "compact stripe hover")
  })

  output$download_ora_csv <- downloadHandler(
    filename = function() paste0("HbaDb_ORA_", input$ora_db, "_", Sys.Date(), ".csv"),
    content  = function(file) utils::write.csv(batch_results()$ora_df, file, row.names = FALSE)
  )
  output$download_ora_xlsx <- downloadHandler(
    filename = function() paste0("HbaDb_ORA_", input$ora_db, "_", Sys.Date(), ".xlsx"),
    content  = function(file) writexl::write_xlsx(list(Enrichment = batch_results()$ora_df), path = file)
  )
  output$download_batch_csv <- downloadHandler(
    filename = function() paste0("HbaDb_matched_genes_", Sys.Date(), ".csv"),
    content  = function(file) utils::write.csv(batch_results()$matched_df, file, row.names = FALSE)
  )
  output$download_batch_xlsx <- downloadHandler(
    filename = function() paste0("HbaDb_matched_genes_", Sys.Date(), ".xlsx"),
    content  = function(file) writexl::write_xlsx(list(Matched_Genes = batch_results()$matched_df), path = file)
  )

  # ---------------------------------------------------------------------------
  # Tab 4: Genome-Wide Analytics Plots
  # ---------------------------------------------------------------------------
  output$plot_coverage <- renderPlotly({
    df <- analytics_cache$cov
    df$Database <- factor(df$Database, levels = df$Database)
    plot_ly(
      df, x = ~Percent, y = ~Database, type = "bar", orientation = "h",
      marker = list(color = "#0E7490"),
      text = ~sprintf("%s genes (%.1f%%)", format(Genes, big.mark = ","), Percent),
      hoverinfo = "text+y"
    ) %>% layout(xaxis = list(title = "Genome Coverage (%)", range = c(0, 100)),
                 yaxis = list(title = ""), margin = list(l = 165))
  })

  output$plot_go_ont <- renderPlotly({
    df <- analytics_cache$go_ont
    plot_ly(df, x = ~Ontology, y = ~Genes, type = "bar", name = "Annotated Genes",
            marker = list(color = "#0F3D5E")) %>%
      add_trace(y = ~Terms, name = "Unique GO Terms", marker = list(color = "#0EA5E9")) %>%
      layout(barmode = "group", yaxis = list(title = "Count"), xaxis = list(title = ""),
             legend = list(orientation = "h", x = 0.2, y = 1.12))
  })

  output$plot_cog_dist <- renderPlotly({
    df <- analytics_cache$cog
    df$Label <- factor(df$Label, levels = df$Label)
    plot_ly(
      df, x = ~Genes, y = ~Label, type = "bar", orientation = "h",
      marker = list(color = "#155E75"),
      text = ~sprintf("%s genes", format(Genes, big.mark = ",")),
      hoverinfo = "text+y"
    ) %>% layout(xaxis = list(title = "Number of Genes"), yaxis = list(title = ""), margin = list(l = 330))
  })

  output$plot_top_pfam <- renderPlotly({
    df <- analytics_cache$pfam
    df$Domain <- factor(df$Domain, levels = df$Domain)
    plot_ly(
      df, x = ~Genes, y = ~Domain, type = "bar", orientation = "h",
      marker = list(color = "#0F766E"),
      text = ~sprintf("%s (%s): %d genes", Domain, Desc, Genes),
      hoverinfo = "text"
    ) %>% layout(xaxis = list(title = "Number of Genes"), yaxis = list(title = ""), margin = list(l = 140))
  })

  output$plot_top_kegg <- renderPlotly({
    df <- analytics_cache$kegg
    df$Pathway <- factor(df$Pathway, levels = df$Pathway)
    plot_ly(
      df, x = ~Genes, y = ~Pathway, type = "bar", orientation = "h",
      marker = list(color = "#0284C7"),
      text = ~sprintf("KEGG %s: %d genes", Pathway, Genes),
      hoverinfo = "text"
    ) %>% layout(xaxis = list(title = "Number of Genes"), yaxis = list(title = ""), margin = list(l = 100))
  })

  # ---------------------------------------------------------------------------
  # Tab 5: Dynamic Reproducible R Script & OrgDb Metadata
  # ---------------------------------------------------------------------------
  output$orgdb_meta_tbl <- renderTable({
    hba_bundle$org_metadata
  }, striped = TRUE, bordered = TRUE, width = "100%")

  dynamic_r_code <- reactive({
    gids <- NULL
    if (!is.null(batch_results()) && nrow(batch_results()$matched_df) > 0L) {
      gids <- head(batch_results()$matched_df$GID, 20L)
    } else if (!is.null(selected_gene_row())) {
      gids <- selected_gene_row()$GID[1]
    } else {
      gids <- head(annotation_data$GID, 5L)
    }
    gids_vec <- paste(sprintf('"%s"', gids), collapse = ", ")
    sprintf(
'# Reproducible R workflow using Bioconductor org.Hbacteriophora.eg.db
# DOI: https://doi.org/10.18129/B9.bioc.org.Hbacteriophora.eg.db

library(AnnotationHub)
library(AnnotationDbi)

ah    <- AnnotationHub()
orgdb <- ah[["AH121717"]]

my_genes <- c(%s)

# 1. Retrieve Gene Symbols & Descriptions
gene_info <- select(
  orgdb,
  keys    = my_genes,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "GID"
)
print(gene_info)

# 2. Retrieve GO Annotations (Biological Process, Molecular Function, Cellular Component)
go_annot <- select(
  orgdb,
  keys    = my_genes,
  columns = c("GO", "ONTOLOGY"),
  keytype = "GID"
)
head(go_annot)

# 3. Retrieve KEGG Pathways
kegg_annot <- select(
  orgdb,
  keys    = my_genes,
  columns = "KEGGPATHWAY",
  keytype = "GID"
)
head(kegg_annot)

# 4. Optional: GO Over-Representation Analysis with clusterProfiler
# library(clusterProfiler)
# ego <- enrichGO(gene = my_genes, OrgDb = orgdb, keyType = "GID", ont = "BP")
# dotplot(ego)
',
      gids_vec
    )
  })

  output$dynamic_r_code_ui <- renderUI({
    tags$div(class = "code-block", dynamic_r_code())
  })

  output$download_r_script <- downloadHandler(
    filename = function() paste0("HbaDb_orgHbacteriophora_workflow_", Sys.Date(), ".R"),
    content  = function(file) writeLines(dynamic_r_code(), con = file)
  )
}

# Run the application
shinyApp(ui = ui, server = server)
