#GSE207194

#Load the TE count matrix

#Read the TE count matrix from file (rows = TE subfamilies, columns = samples)
counts <- read.table(
  "GSE207194_TE_counts_matrix.final.txt",
  header      = TRUE,
  row.names   = 1,
  sep         = "\t",
  check.names = FALSE
)

#Check the dimensions of the matrix
dim(counts)  #1180 TEs x 61 samples

#2) Load and prepare the metadata

#Load the GEOquery library to download sample metadata from GEO
library(GEOquery)

#Download the GEO ExpressionSet for dataset GSE50760
gse  <- getGEO("GSE207194", GSEMatrix = TRUE)

#Extract the sample metadata table (rows = samples, columns = annotations)
meta <- pData(gse[[1]])

#Build a clean metadata dataframe with relevant columns
metadata <- data.frame(
  GSM      = meta$geo_accession,
  title    = meta$title,
  tissue   = sub("^tissue: ", "", meta$`tissue:ch1`),
  treatment = sub("^treatment: ", "", meta$`treatment:ch1`),
  paired_sample = meta$`paired_sample:ch1`,
  stringsAsFactors = FALSE
)


#Set rownames of metadata to match column names of the count matrix
rownames(metadata) <- colnames(counts)

#Verify that metadata rows and count matrix columns are aligned
all(rownames(metadata) == colnames(counts))


#Rename the tissues
tissue_map <- c(
  "COLON" = "primary_colon",
  "LIVER" = "liver_metastasis"
)

metadata$tissue <- tissue_map[metadata$tissue]


# Create a patient column 
metadata$patient <- sub("^([A-Z]{2}).*", "\\1", metadata$title)


all(colnames(counts) == rownames(metadata))


# 3) Quality Control
# 3.1 Tumor purity QC

library(dplyr)
library(ggplot2)
library(ggrepel)

# Load tumor purity estimates
purity <- read.table(
  "GSE207194_purity_results.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Convert metadata to standard dataframe
metadata <- as.data.frame(metadata)

# Add purity values to metadata
metadata$purity <- as.numeric(
  purity$purity[
    match(rownames(metadata), purity$X)
  ]
)

# Create ordered factor for tissue condition
metadata$cond <- factor(
  metadata$tissue,
  levels = c("primary_colon", "liver_metastasis")
)

# Summarize purity statistics
purity_table <- metadata %>%
  group_by(cond) %>%
  summarise(
    Samples       = n(),
    Min_Purity    = round(min(purity, na.rm = TRUE), 3),
    Median_Purity = round(median(purity, na.rm = TRUE), 3),
    Max_Purity    = round(max(purity, na.rm = TRUE), 3)
  )

# Add SRR column for plotting
metadata$SRR <- rownames(metadata)

# Identify lowest purity sample
low_purity <- metadata[which.min(metadata$purity), ]

# Plot tumor purity distribution
ggplot(metadata, aes(x = cond, y = purity)) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "grey90",
    colour = "black"
  ) +
  geom_jitter(
    aes(color = cond),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  geom_text_repel(
    data = low_purity,
    aes(label = SRR),
    color = "black",
    size = 3,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = c(
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    )
  ) +
  scale_x_discrete(
    labels = c(
      primary_colon    = "Primary colon",
      liver_metastasis = "Liver metastasis"
    )
  ) +
  labs(
    title = "Tumor purity distribution by condition",
    x     = "Condition",
    y     = "Tumor purity (PUREE)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# 3.2 Library size QC

# Calculate total read counts per sample
lib_size <- colSums(counts)

# Create dataframe
lib_df <- data.frame(
  sample   = names(lib_size),
  lib_size = lib_size,
  cond     = metadata$cond,
  stringsAsFactors = FALSE
)

# Basic statistics
summary(lib_size)[c("Min.", "Median", "Max.")]

# Order samples for plotting
lib_df <- lib_df %>%
  mutate(
    cond = factor(
      cond,
      levels = c("primary_colon", "liver_metastasis")
    )
  ) %>%
  arrange(cond, lib_size)

# Plot library sizes
ggplot(
  lib_df,
  aes(
    x = factor(sample, levels = sample),
    y = lib_size,
    fill = cond
  )
) +
  geom_bar(stat = "identity") +
  
  scale_fill_manual(
    values = c(
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    ),
    labels = c(
      primary_colon    = "Primary colon",
      liver_metastasis = "Liver metastasis"
    )
  ) +
  
  scale_y_continuous(labels = scales::comma) +
  
  labs(
    title = "GSE207194 — Library size per sample",
    x     = "Sample",
    y     = "Total TE counts",
    fill  = "Condition"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 8
    ),
    legend.position = "right"
  )


# Correlation plots: tumor purity vs sequencing depth

library(ggplot2)
library(ggrepel)
library(scales)

# Calculate sequencing depth per sample
depth <- colSums(counts)

# Build correlation dataframe
corr_df <- data.frame(
  sample  = names(depth),
  depth   = depth,
  purity  = metadata$purity,
  cond    = metadata$cond,
  patient = metadata$patient,
  stringsAsFactors = FALSE
)

# Compute Pearson correlation
r    <- cor(corr_df$depth,
            corr_df$purity,
            use = "complete.obs")

r2   <- r^2

pval <- cor.test(
  corr_df$depth,
  corr_df$purity
)$p.value

pval_label <- ifelse(
  pval < 0.001,
  "p < 0.001",
  paste0("p = ", round(pval, 3))
)

#Plot
ggplot(corr_df,
       aes(x = depth,
           y = purity,
           color = cond)) +
  
  geom_point(
    size = 3.5,
    alpha = 0.85
  ) +
  
  scale_x_continuous(
    labels = comma
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "grey40",
    fill = "grey85",
    linewidth = 0.8
  ) +
  
  geom_text_repel(
    aes(label = patient),
    size = 3,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    ),
    labels = c(
      primary_colon    = "Primary colon",
      liver_metastasis = "Liver metastasis"
    )
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.5,
    label = paste0(
      "r = ", round(r, 3), "\n",
      "R² = ", round(r2, 3), "\n",
      pval_label
    ),
    size = 3.5,
    color = "grey20"
  ) +
  
  labs(
    title = "Tumor purity vs sequencing depth",
    x     = "Sequencing depth (total TE counts)",
    y     = "Tumor purity (PUREE)",
    color = "Condition"
  ) +
  
  theme_bw()

# Plot 2: Correlation faceted by condition

ggplot(corr_df,
       aes(x = depth,
           y = purity,
           color = cond)) +
  
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  
  scale_x_continuous(
    labels = comma
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "grey40",
    fill = "grey85",
    linewidth = 0.8
  ) +
  
  geom_text_repel(
    aes(label = patient),
    size = 2.5,
    show.legend = FALSE
  ) +
  
  facet_wrap(
    ~ cond,
    labeller = labeller(
      cond = c(
        primary_colon    = "Primary colon",
        liver_metastasis = "Liver metastasis"
      )
    )
  ) +
  
  scale_color_manual(
    values = c(
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    ),
    labels = c(
      primary_colon    = "Primary colon",
      liver_metastasis = "Liver metastasis"
    )
  ) +
  
  labs(
    title = "Tumor purity vs sequencing depth by condition",
    x     = "Sequencing depth (total TE counts)",
    y     = "Tumor purity (PUREE)",
    color = "Condition"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "none"
  )


# 4) DESeq2 — PCA and outlier detection

library(DESeq2)
library(ggplot2)
library(ggrepel)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ patient + cond
)

# Run DESeq2 pipeline
dds <- DESeq(dds)

# Inspect size factors
sizeFactors(dds)

# Check normalization range
range(sizeFactors(dds))

# Store size factors
sf_df <- data.frame(
  sample      = names(sizeFactors(dds)),
  size_factor = sizeFactors(dds),
  cond        = metadata$cond
)

# Plot size factors
ggplot(sf_df,
       aes(x = sample,
           y = size_factor,
           fill = cond)) +
  
  geom_bar(stat = "identity") +
  
  scale_fill_manual(
    values = c(
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    )
  ) +
  
  labs(
    title = "DESeq2 size factors per sample",
    x     = "Sample",
    y     = "Size factor",
    fill  = "Condition"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  )

# Dispersion plot
plotDispEsts(
  dds,
  main = "Dispersion estimates — TEs"
)


# VST for PCA/QC
vsd <- vst(dds, blind = TRUE)

# Perform PCA
pca_res <- prcomp(t(assay(vsd)))

# Variance explained
pct_var <- round(
  100 * pca_res$sdev^2 / sum(pca_res$sdev^2),
  1
)

# Build PCA dataframe
df_pca <- data.frame(
  sample  = colnames(vsd),
  cond    = metadata$cond,
  patient = metadata$patient,
  pca_res$x
)

# PC1 vs PC2
ggplot(
  df_pca,
  aes(
    PC1,
    PC2,
    colour = cond,
    label = patient
  )
) +
  
  geom_point(size = 4) +
  
  geom_text_repel(
    size = 3,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    ),
    labels = c(
      primary_colon    = "Primary colon",
      liver_metastasis = "Liver metastasis"
    )
  ) +
  
  labs(
    title = "PC1 vs PC2 of VST-normalized TE counts",
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)"),
    colour = "Condition"
  ) +
  
  theme_bw()


## 5. Remove PCA outlier patient (HP) and rerun DESeq2 ####

# Identify all samples from patient HP
outlier_samples <- rownames(metadata)[
  metadata$patient == "HP"
]

outlier_samples

# Remove the whole patient from the count matrix
counts_qc <- counts[, !colnames(counts) %in% outlier_samples]

# Update metadata
metadata_qc <- metadata[
  !rownames(metadata) %in% outlier_samples,
]

# Reorder metadata to perfectly match counts
metadata_qc <- metadata_qc[colnames(counts_qc), ]

# Verify alignment
all(colnames(counts_qc) == rownames(metadata_qc))

# Create DESeq2 object
dds_qc <- DESeqDataSetFromMatrix(
  countData = counts_qc,
  colData   = metadata_qc,
  design    = ~ patient + cond
)

# Run DESeq2
dds_qc <- DESeq(dds_qc)

# VST for PCA
vsd_qc <- vst(dds_qc, blind = TRUE)

# PCA
pca_qc <- prcomp(t(assay(vsd_qc)))

# Variance explained
pct_var_qc <- round(
  100 * pca_qc$sdev^2 / sum(pca_qc$sdev^2),
  1
)

# PCA dataframe
df_pca_qc <- data.frame(
  sample  = colnames(vsd_qc),
  cond    = metadata_qc$cond,
  patient = metadata_qc$patient,
  pca_qc$x
)

# PC1 vs PC2 after removal
ggplot(
  df_pca_qc,
  aes(PC1, PC2, colour = cond, label = patient)
) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  
  scale_color_manual(
    values = c(
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    ),
    labels = c(
      primary_colon    = "Primary colon",
      liver_metastasis = "Liver metastasis"
    )
  ) +
  
  labs(
    title = "PC1 vs PC2 after removal of HP patient",
    x = paste0("PC1 (", pct_var_qc[1], "%)"),
    y = paste0("PC2 (", pct_var_qc[2], "%)"),
    colour = "Condition"
  ) +
  
  theme_bw()



#Compute the signature scores
library(scuttle)

# TPM matrix (after outlier removal)
tpm_mat <- calculateTPM(
  counts_qc,
  lengths = rep(1, nrow(counts_qc))
)

# Load TE signatures 
TE_68  <- scan("tumor_associated_68_TEs.txt", what = character())
TE_141 <- scan("metastasis_associated_141_TEs.txt", what = character())

# Add treatment group to QC metadata
metadata_qc$treatment_group <- ifelse(
  metadata_qc$treatment == "TRG4-5",
  "Treated",
  "Untreated"
)

metadata_qc$treatment_group <- factor(
  metadata_qc$treatment_group,
  levels = c("Untreated", "Treated")
)

# 1. 68‑TE signature 
TE_common_68 <- intersect(TE_68, rownames(tpm_mat))

expr_68 <- tpm_mat[TE_common_68, ]
score_68 <- colMeans(expr_68, na.rm = TRUE)

df_68 <- data.frame(
  sample          = names(score_68),
  score           = score_68,
  tissue          = metadata_qc$tissue,
  treatment_group = metadata_qc$treatment_group,
  stringsAsFactors = FALSE
)

write.csv(
  df_68,
  "TE_signature68_scores_GSE207194_QC_TREATED_UNTREATED.csv",
  row.names = FALSE
)

# 2. 141‑TE signature
TE_common_141 <- intersect(TE_141, rownames(tpm_mat))

expr_141 <- tpm_mat[TE_common_141, ]
score_141 <- colMeans(expr_141, na.rm = TRUE)

df_141 <- data.frame(
  sample          = names(score_141),
  score           = score_141,
  tissue          = metadata_qc$tissue,
  treatment_group = metadata_qc$treatment_group,
  stringsAsFactors = FALSE
)

write.csv(
  df_141,
  "TE_signature141_scores_GSE207194_QC_TREATED_UNTREATED.csv",
  row.names = FALSE
)

