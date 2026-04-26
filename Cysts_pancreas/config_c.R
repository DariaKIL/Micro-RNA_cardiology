# Transplantation/config_c.R

PHENOTYPE_FILE <- here("Cysts_pancreas", "data", "phenotable.tsv")

OUTPUT_DIR <- "results"
FIGURES_DIR <- "figures"

# Create directories if they do not exist
dir.create(OUTPUT_DIR, showWarnings = FALSE)
dir.create(FIGURES_DIR, showWarnings = FALSE)


# Color palette for groups
GROUP_COLORS <- c(
  "mucinous" = "#E41A1C",
  "serous" = "#377EB8"
)


CONTRASTS <- c("group", "mucinous", "serous")

set.seed(42)

cat("✓ Configuration loaded successfully\n")