library(dplyr)
library(tidyr)
library(stringr)

go <- final_GO_enrichment_file[,c(1,2)]
colnames(go) <- c("query", "GOs")

go1 <- go %>%
  group_by(query) %>%
  summarise(GOs = paste(GOs, collapse = ","), .groups = "drop")

file <- Hbac_annotation_new[,c(1,10)]

# Merge the two dataframes
merged_go <- full_join(go1, file, by = "query", suffix = c("_1", "_2"))

# Combine both GOs columns, remove NAs, split by comma, remove duplicates, and paste back
merged_go <- merged_go %>%
  rowwise() %>%
  mutate(
    GOs_combined = paste(c(GOs_1, GOs_2), collapse = ","),
    GOs_final = paste(unique(na.omit(str_split(GOs_combined, ",")[[1]])), collapse = ",")
  ) %>%
  ungroup() %>%
  select(query, GOs = GOs_final)

cleaned_go <- merged_go %>%
  mutate(
    GOs = GOs %>%
      str_replace_all(",\\-", "") %>%   # Remove ',-'
      str_replace_all("\\-,", "") %>%   # Remove '-,'
      str_replace_all(",{2,}", ",") %>% # Replace multiple commas with single comma
      str_replace_all("^,|,$", "")      # Remove leading or trailing commas
  )

data <- Hbac_annotation_new
data_updated <- data %>%
  left_join(cleaned_go, by = "query", suffix = c("", ".new")) %>%
  rowwise() %>%
  mutate(
    # Combine both columns, handle NAs safely
    combined = paste(na.omit(c(GOs, GOs.new)), collapse = ","),

    # Split by comma, remove unwanted values, trim whitespace
    cleaned = combined %>%
      str_split(",") %>%
      unlist() %>%
      str_trim() %>%
      .[. != "" & . != "NA" & . != "-"] %>%
      unique() %>%
      paste(collapse = ","),

    # Final assignment with fallback to "-"
    GOs = ifelse(cleaned == "", "-", cleaned)
  ) %>%
  ungroup() %>%
  select(-GOs.new, -combined, -cleaned)

colnames(data_updated)[1] <- "GID"
writexl::write_xlsx(data_updated, "./database.xlsx")

div(style = "text-align: center; margin: 30px auto; max-width: 250px;",  # Container for centering
    tags$div(
      class = "visitor-map",
      style = "display: flex; justify-content: center;",  # Additional centering styles
      tags$a(
        href = "https://mapmyvisitors.com/web/1bxbm",
        target = "_blank",
        title = "Visit tracker",
        tags$img(
          src = "https://mapmyvisitors.com/map.png?d=JjosJ6JAVG9bsS0dycAv8Wp7IEc9rPmDuZmbxnj2dyA&cl=ffffff",
          style = "border: none; width: 150px;"
        )
      )
    )
)
