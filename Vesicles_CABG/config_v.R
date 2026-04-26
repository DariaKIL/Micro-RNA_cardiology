# Transplantation/config_c.R

PHENOTYPE_FILE <- here("Vesicles_CABG", "data", "phenotableV.tsv")

OUTPUT_DIR <- "results"
FIGURES_DIR <- "figures"

# Create directories if they do not exist
dir.create(OUTPUT_DIR, showWarnings = FALSE)
dir.create(FIGURES_DIR, showWarnings = FALSE)


# Color palette for groups
GROUP_COLORS <- c(
  "before_16"  = "#4DAF4A",
  "before_150" = "#FF7F00",
  "after_16"   = "#6E4A99",
  "after_150"  = "#377EB8"
)

CONTRASTS <- list(
  time = c("time", "before", "after"),
  type = c("type", "16", "150"),
  patient = c("patient", "29", "15", "17"),
  condition = c("condition", "after_150", "after_16", "before_150", "before_16")
)

set.seed(42)

cat("✓ Configuration loaded successfully\n")