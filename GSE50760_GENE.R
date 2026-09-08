#GSE50760 genes
counts <- read.table(
  "GSE50760_GENE_counts_matrix.final.txt",
  header      = TRUE,
  row.names   = 1,
  sep         = "\t",
  check.names = FALSE
)

counts <- counts[-1, ]
counts <- as.matrix(counts)
mode(counts) <- "numeric"
dim(counts)

#Metadata
library(GEOquery)
library(dplyr)

gse <- getGEO(
  "GSE50760",
  GSEMatrix = TRUE
)

gse  <- getGEO("GSE50760", GSEMatrix = TRUE)

meta <- pData(gse[[1]])

metadata <- data.frame(
  GSM        = meta$geo_accession,
  tissue_raw = meta$`tissue:ch1`,
  ajcc_stage = sub("ajcc stage: ", "", meta$`ajcc stage:ch1`),
  patient    = factor(sub(".*(AMC_[0-9]+).*", "\\1", meta$title)),
  stringsAsFactors = FALSE
)

metadata$tissue <- recode(
  metadata$tissue_raw,
  "primary colorectal cancer"                     = "primary_colon",
  "normal-looking surrounding colonic epithelium" = "normal_colon",
  "metastatic colorectal cancer to the liver"     = "liver_metastasis"
)

metadata$cond <- factor(
  metadata$tissue,
  levels = c("normal_colon", "primary_colon", "liver_metastasis")
)

rownames(metadata) <- colnames(counts)
all(rownames(metadata) == colnames(counts))
table(metadata$cond, useNA = "ifany")


#QC: Tumor purity

library(ggplot2)
library(ggrepel)

purity <- read.table("GSE50760_purity_results.tsv", header = TRUE)

metadata$purity <- purity$purity[
  match(rownames(metadata), rownames(purity))
]

metadata$cond <- factor(
  metadata$tissue,
  levels = c("normal_colon", "primary_colon", "liver_metastasis")
)

purity_table <- metadata %>%
  group_by(cond) %>%
  summarise(
    Samples       = n(),
    Min_Purity    = round(min(purity),    3),
    Median_Purity = round(median(purity), 3),
    Max_Purity    = round(max(purity),    3)
  )

metadata %>%
  dplyr::select(GSM, cond, patient, purity) %>%
  arrange(purity)

# Remove low-purity sample
sample_to_remove <- "SRR975557"
keep             <- colnames(counts) != sample_to_remove
counts_qc        <- counts[, keep]
metadata_qc      <- metadata[keep, ]

purity_table_qc <- metadata_qc %>%
  group_by(cond) %>%
  summarise(
    Samples       = n(),
    Min_Purity    = round(min(purity),    3),
    Median_Purity = round(median(purity), 3),
    Max_Purity    = round(max(purity),    3)
  )

# TPM threshold — STEP 1: Compute TPM from counts_qc
library(biomaRt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(biomaRt)

# Connect directly to human genes
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

ensg_clean <- sub(
  "\\..*$",
  "",
  rownames(counts_qc)
)

gene_coords <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "start_position",
    "end_position"
  ),
  filters = "ensembl_gene_id",
  values = unique(ensg_clean),
  mart = mart
) %>%
  mutate(
    length_bp = end_position - start_position + 1
  ) %>%
  group_by(
    ensembl_gene_id
  ) %>%
  summarise(
    length_bp = max(length_bp),
    .groups = "drop"
  )

matched_lengths <- gene_coords$length_bp[
  match(ensg_clean, gene_coords$ensembl_gene_id)
]

valid         <- !is.na(matched_lengths)
counts_valid  <- counts_qc[valid, ]
lengths_valid <- matched_lengths[valid]

rpk <- sweep(counts_valid, 1, lengths_valid / 1000, FUN = "/")
tpm <- sweep(rpk, 2, colSums(rpk) / 1e6, FUN = "/")


# Retain genes with >= 5 counts in >= 20% of samples
min_fraction <- 0.20
min_samples  <- ceiling(min_fraction * ncol(counts_qc))

keep_genes <- rowSums(counts_qc >= 5) >= min_samples

counts_filtered <- counts_qc[keep_genes, ]

# Remove PCA outlier 
samples_to_remove <- rownames(metadata)[
  metadata$patient == "AMC_10" &
    metadata$cond == "liver_metastasis"
]

counts_qc2 <- counts_filtered[
  ,
  !(colnames(counts_filtered) %in% samples_to_remove)
]

metadata_qc2 <- metadata[
  colnames(counts_qc2),
]

stopifnot(
  all(colnames(counts_qc2) == rownames(metadata_qc2))
)

metadata_qc2$cond <- factor(
  metadata_qc2$cond,
  levels = c(
    "normal_colon",
    "primary_colon",
    "liver_metastasis"
  )
)



##Check the expression of the tigger7_dup2343 nearest gene --> LINC02418 
library(tibble)

# Gene of interest
genes_of_interest <- "LINC02418"


# Get Ensembl ID
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)

gene_info <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol"
  ),
  filters = "hgnc_symbol",
  values = genes_of_interest,
  mart = mart
)

print(gene_info)


# Match gene to TPM matrix
ensg <- gene_info$ensembl_gene_id[1]

tpm_qc <- tpm[, rownames(metadata_qc2)]

stopifnot(
  all(colnames(tpm_qc) == rownames(metadata_qc2))
)

gene_row <- grep(
  paste0("^", ensg),
  rownames(tpm_qc),
  value = TRUE
)[1]

print(gene_row)

# ------------------------------------------------------------
# Long TPM dataframe
# ------------------------------------------------------------

df_gene <- as.data.frame(
  tpm_qc[
    gene_row,
    ,
    drop = FALSE
  ]
) %>%
  rownames_to_column("ensembl") %>%
  mutate(
    gene = "LINC02418"
  ) %>%
  pivot_longer(
    cols = -c(ensembl, gene),
    names_to = "sample",
    values_to = "TPM"
  ) %>%
  left_join(
    metadata_qc2 %>%
      rownames_to_column("sample") %>%
      dplyr::select(
        sample,
        cond
      ),
    by = "sample"
  )

# ------------------------------------------------------------
# Summary table
# ------------------------------------------------------------

gene_summary <- df_gene %>%
  group_by(cond) %>%
  summarise(
    Mean_TPM = round(mean(TPM), 2),
    Median_TPM = round(median(TPM), 2),
    Max_TPM = round(max(TPM), 2),
    .groups = "drop"
  )

print(gene_summary)

# ------------------------------------------------------------
# Boxplot
# ------------------------------------------------------------

ggplot(
  df_gene,
  aes(
    x = cond,
    y = TPM + 0.01,
    fill = cond
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  scale_y_log10() +
  scale_fill_manual(
    values = c(
      normal_colon     = "#2CA02C",
      primary_colon    = "#1F77B4",
      liver_metastasis = "#D62728"
    )
  ) +
  scale_x_discrete(
    labels = c(
      normal_colon     = "Normal Colon",
      primary_colon    = "Primary Colon",
      liver_metastasis = "Liver Metastasis"
    )
  ) +
  labs(
    title = "LINC02418 Expression Across Tumor Progression",
    x = NULL,
    y = "Expression (TPM, log10 scale)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      face = "bold"
    )
  )
