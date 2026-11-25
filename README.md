# HbaDb: *Heterorhabditis bacteriophora* Functional Genome Database

[![Bioconductor OrgDb](https://img.shields.io/badge/Bioconductor-org.Hbacteriophora.eg.db-1A81C2.svg)](https://doi.org/10.18129/B9.bioc.org.Hbacteriophora.eg.db)
[![AnnotationHub](https://img.shields.io/badge/AnnotationHub-AH121717-0F766E.svg)](https://bioconductor.org/packages/org.Hbacteriophora.eg.db)
[![Live Shiny App](https://img.shields.io/badge/ShinyApps-HbaDb-0284C7.svg)](https://dabiniasri.shinyapps.io/HbaDb/)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)

**HbaDb** is an interactive R Shiny web portal and functional enrichment workbench for the entomopathogenic nematode ***Heterorhabditis bacteriophora*** (NCBI TaxID: `37862`), built as a companion web interface to Bioconductor's [**`org.Hbacteriophora.eg.db`**](https://doi.org/10.18129/B9.bioc.org.Hbacteriophora.eg.db) `OrgDb` annotation package.

---

## Database Coverage (`11,466` Genes × `29` Annotation Fields)

| Annotation Layer | Annotated Genes | Coverage | Unique Terms / Families |
|:---|---:|---:|---:|
| Custom Gene Identifiers (`GID`, `Hba_*`) | 11,466 | 100.0% | 11,466 |
| Functional Descriptions (`GENENAME` / `Description`) | 9,600 | 83.7% | — |
| COG Functional Categories | 9,600 | 83.7% | 24 categories |
| Gene Ontology (`GO_BP`, `GO_CC`, `GO_MF` via `OrgDb`) | 9,470 | 82.6% | 17,439 (`12,679` BP / `1,626` CC / `3,134` MF) |
| PFAM Protein Domains | 7,199 | 62.8% | 2,875 domains |
| KEGG Orthology (`KEGG_ko`) | 6,327 | 55.2% | — |
| Named Gene Symbols (`SYMBOL` / `Preferred_name`) | 5,865 | 51.2% | — |
| MetaCyc / Reactome Pathways | 5,452 | 47.5% | — |
| KEGG Pathways (`KEGG_Pathway` / `OrgDb_KEGG`) | 4,049 | 35.3% | 379 pathways |
| EC Numbers | 2,461 | 21.5% | — |
| CAZy Carbohydrate-Active Enzymes | 154 | 1.3% | — |

---

## Key Features

1. **Overview Dashboard (`Home`)**: Live KPI summary cards, Bioconductor DOI / AnnotationHub links, and citation reference.
2. **Gene & Annotation Explorer (`Gene Explorer`)**:
   - Case-insensitive global and column-specific search (substring or exact match).
   - Multi-facet filters for decoded **COG functional categories**, **eggNOG e-value cutoffs**, and **annotation presence** (Symbol, GO, KEGG, PFAM, EC, CAZy, MetaCyc/Reactome).
   - Column presets (`Core Overview`, `GO & Functional`, `KEGG & Metabolic`, `PFAM & Domains`, `All 29 Columns`) and `.csv` / `.xlsx` export.
   - Click any row to open the **Gene Detail Card** with clickable external links to **QuickGO** (with human-readable GO term names), **KEGG Pathways/KO**, **InterPro/PFAM**, **ExPASy EC**, and **WormBase ParaSite**.
3. **Batch Gene Lookup & Over-Representation Analysis (`Batch Lookup & Enrichment`)**:
   - Paste or upload lists of `Hba_*` gene IDs or gene symbols (includes a 1-click 40-gene example).
   - Instant hypergeometric **ORA enrichment** (BH-adjusted FDR) across **GO BP**, **GO MF**, **GO CC**, **KEGG Pathways**, **PFAM Domains**, and **COG Categories**, with interactive `plotly` dotplots and `.csv` / `.xlsx` downloads.
4. **Genome-Wide Analytics (`Genome Analytics`)**: Interactive `plotly` charts for database coverage, GO ontology breakdown, COG distribution, top 20 PFAM domains, and top 20 KEGG pathways.
5. **Bioconductor Integration (`Bioconductor & R Code`)**: Installation commands, `AnnotationHub` (`AH121717`) usage, `OrgDb` schema metadata, and a dynamic reproducible `.R` script generator synced with your active genes.

---

## Running Locally

### Prerequisites
```r
install.packages(c("shiny", "bslib", "bsicons", "DT", "dplyr", "ggplot2", "plotly", "writexl"))
```

### Launch the App
```r
git clone https://github.com/kabilanbio/HbaDb.git
cd HbaDb
R -e "shiny::runApp('.')"
```

### Rebuilding `data/hba_annotation.rds` (Optional)
The pre-indexed `data/hba_annotation.rds` (`3.06 MB`) is included in the repository so the app launches immediately. To regenerate it from raw `database.xlsx`, `org.Hbacteriophora.eg.sqlite`, and `GO.db`, run:
```r
source("scripts/build_hba_data.R")
```

---

## Citation

If you use **HbaDb** or **`org.Hbacteriophora.eg.db`** in your research, please cite:

> Johnson LW, Bhat CG, Dillman A, Somvanshi VS, Kabilan S (2025). *org.Hbacteriophora.eg.db: Genome-wide Annotation for Heterorhabditis bacteriophora*. R/Bioconductor package. DOI: [10.18129/B9.bioc.org.Hbacteriophora.eg.db](https://doi.org/10.18129/B9.bioc.org.Hbacteriophora.eg.db)

## Authors & Maintainer

- **Lalson Wesly Johnson** (Author) — [ORCID: 0000-0001-9121-1927](https://orcid.org/0000-0001-9121-1927)
- **Chaitra Ganapati Bhat** (Author) — [ORCID: 0000-0003-1497-890X](https://orcid.org/0000-0003-1497-890X)
- **Adler Dillman** (Author) — [ORCID: 0000-0001-7171-4332](https://orcid.org/0000-0001-7171-4332)
- **Vishal Singh Somvanshi** (Author)
- **Kabilan S** (Author & Maintainer, `kabilan151414@gmail.com`) — [ORCID: 0009-0002-9422-0938](https://orcid.org/0009-0002-9422-0938)
