#GSE124535_liver
#Load the TE count matrix

#Read the TE count matrix from file (rows = TE subfamilies, columns = samples)
counts <- read.table(
  "GSE124535_liver_TE_counts_matrix.final.txt",
  header      = TRUE,
  row.names   = 1,
  sep         = "\t",
  check.names = FALSE
)

#Check the dimensions of the matrix
dim(counts)  #1180 TEs x 70 samples

#2) Load and prepare the metadata
#Load the GEOquery library to download sample metadata from GEO
library(GEOquery)

#Download the GEO ExpressionSet for dataset GSE50760
gse  <- getGEO("GSE124535", GSEMatrix = TRUE)

#Extract the sample metadata table (rows = samples, columns = annotations)
meta <- pData(gse[[1]])


# Build metadata dataframe
metadata <- data.frame(
  GSM          = meta$geo_accession,
  tissue       = meta$`tissue:ch1`,
  tumor_stage  = meta$`tumor stage:ch1`,
  sample       = meta$title,
  stringsAsFactors = FALSE
)

# Clean tissue labels
library(dplyr)
metadata$tissue <- recode(
  metadata$tissue,
  "non-tumor" = "normal_liver",
  "tumor"     = "HCC"
)

#Set rownames of metadata to match column names of the count matrix
rownames(metadata) <- colnames(counts)

#Verify that metadata rows and count matrix columns are aligned
all(rownames(metadata) == colnames(counts))


#3) Quality Control

# 3.1 Tumor purity QC

#Load libraries for data manipulation and visualization
library(dplyr)
library(ggplot2)
library(ggrepel)

#Load tumor purity estimates
purity <- read.table(
  "GSE124535_liver_purity_results.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

#Add purity values to metadata
metadata$purity <- purity$purity[
  match(rownames(metadata), rownames(purity))
]

#Create ordered factor for tissue condition
metadata$cond <- factor(
  metadata$tissue,
  levels = c("normal_liver", "HCC")
)

#Summarize purity statistics per tissue group
purity_table <- metadata %>%
  group_by(cond) %>%
  summarise(
    Samples       = n(),
    Min_Purity    = round(min(purity), 3),
    Median_Purity = round(median(purity), 3),
    Max_Purity    = round(max(purity), 3)
  )

purity_table


#Identify lowest purity samples
metadata %>%
  dplyr::select(GSM, sample, cond, purity) %>%
  arrange(purity)


#Plot tumor purity distribution
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
  
  scale_color_manual(
    values = c(
      normal_liver = "#4C72B0",
      HCC           = "#C44E52"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      normal_liver = "Normal liver",
      HCC          = "HCC"
    )
  ) +
  
  labs(
    title = "Tumor purity distribution by condition",
    x = "Condition",
    y = "Tumor purity (PUREE)"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "none"
  )

##Change the variables name
counts_qc <- counts
metadata_qc <- metadata


#3.2 Library size QC

# Calculate total read counts per sample (library size)
lib_size <- colSums(counts_qc)

# Create a dataframe with sample IDs, library sizes, and condition
lib_df <- data.frame(
  sample   = names(lib_size),
  lib_size = lib_size,
  cond     = metadata_qc$cond,
  stringsAsFactors = FALSE
)

# Print basic library size statistics
summary(lib_size)[c("Min.", "Median", "Max.")]

# Order samples by condition and library size for plotting
lib_df <- lib_df %>%
  mutate(
    cond = factor(cond, levels = c("normal_liver", "HCC"))
  ) %>%
  arrange(cond, lib_size)

# Plot library size per sample, colored by tissue condition
ggplot(
  lib_df,
  aes(x = factor(sample, levels = sample),
      y = lib_size,
      fill = cond)
) +
  
  geom_bar(stat = "identity") +
  
  scale_fill_manual(
    values = c(
      normal_liver = "#55A868",
      HCC          = "#4C72B0"
    ),
    labels = c(
      normal_liver = "Normal liver",
      HCC          = "HCC"
    )
  ) +
  
  labs(
    title = "GSE124535 — Library size per sample",
    x     = "Sample",
    y     = "Total read counts",
    fill  = "Condition"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    legend.position = "right"
  )


#Correlation plots(purity of samples x depth of sequencing)
metadata_qc$patient <- metadata_qc$sample

library(ggplot2)
library(ggrepel)

# Calculate library size
depth <- colSums(counts_qc)

# Build correlation dataframe
corr_df <- data.frame(
  sample  = names(depth),
  depth   = depth,
  purity  = metadata_qc$purity,
  cond    = metadata_qc$cond,
  patient = metadata_qc$sample,
  stringsAsFactors = FALSE
)

# Pearson correlation
r <- cor(
  corr_df$depth,
  corr_df$purity,
  use = "complete.obs"
)

r2 <- r^2

pval <- cor.test(
  corr_df$depth,
  corr_df$purity
)$p.value

pval_label <- ifelse(
  pval < 0.001,
  "p < 0.001",
  paste0("p = ", round(pval, 3))
)


# Combined correlation plot
ggplot(corr_df,
       aes(x = depth, y = purity, color = cond)) +
  
  geom_point(size = 3.5, alpha = 0.85) +
  
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
  
  scale_color_manual(
    values = c(
      normal_liver = "#55A868",
      HCC          = "#4C72B0"
    ),
    labels = c(
      normal_liver = "Normal liver",
      HCC          = "HCC"
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
    x = "Sequencing depth",
    y = "Tumor purity (PUREE)",
    color = "Condition"
  ) +
  
  theme_bw()


# Faceted correlation plot
ggplot(corr_df,
       aes(x = depth, y = purity, color = cond)) +
  
  geom_point(size = 3, alpha = 0.85) +
  
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
        normal_liver = "Normal liver",
        HCC          = "HCC"
      )
    )
  ) +
  
  scale_color_manual(
    values = c(
      normal_liver = "#55A868",
      HCC          = "#4C72B0"
    )
  ) +
  
  labs(
    title = "Tumor purity vs sequencing depth by condition",
    x = "Sequencing depth",
    y = "Tumor purity (PUREE)",
    color = "Condition"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "none"
  )


# Tumor purity vs sequencing depth (HCC only)

library(dplyr)
library(ggplot2)
library(ggrepel)

# Calculate sequencing depth
depth <- colSums(counts_qc)

# Build dataframe
corr_df <- data.frame(
  sample  = names(depth),
  depth   = depth,
  purity  = metadata_qc$purity,
  cond    = metadata_qc$cond,
  patient = metadata_qc$sample,
  stringsAsFactors = FALSE
)

# Keep only HCC samples
corr_df_tumor <- corr_df %>%
  filter(cond == "HCC")

# Compute correlation
r <- cor(
  corr_df_tumor$depth,
  corr_df_tumor$purity,
  use = "complete.obs"
)

r2 <- r^2

pval <- cor.test(
  corr_df_tumor$depth,
  corr_df_tumor$purity
)$p.value

pval_label <- ifelse(
  pval < 0.001,
  "p < 0.001",
  paste0("p = ", round(pval, 3))
)

# Plot
ggplot(
  corr_df_tumor,
  aes(x = depth, y = purity)
) +
  
  geom_point(
    color = "#4C72B0",
    size = 3.5,
    alpha = 0.85
  ) +
  
  scale_x_continuous(
    labels = scales::comma
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
    size = 2.8
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
    title = "Tumor purity vs sequencing depth (HCC only)",
    x = "Sequencing depth (total read counts)",
    y = "Tumor purity (PUREE)"
  ) +
  
  theme_bw()


### Fix the pairing 
metadata_qc$patient <- sub("[PT]$", "", metadata_qc$sample)

#4) DESeq2 - PCA and QC

# Load libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Create paired patient variable
metadata_qc$patient <- sub("[PT]$", "", metadata_qc$sample)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_qc,
  colData   = metadata_qc,
  design    = ~ patient + cond
)

# Run DESeq2
dds <- DESeq(dds)

# Inspect size factors
sizeFactors(dds)

# Check size factor range
range(sizeFactors(dds))

# Store size factors for plotting
sf_df <- data.frame(
  sample      = names(sizeFactors(dds)),
  size_factor = sizeFactors(dds),
  cond        = metadata_qc$cond
)

# Plot size factors
ggplot(
  sf_df,
  aes(x = sample, y = size_factor, fill = cond)
) +
  
  geom_bar(stat = "identity") +
  
  scale_fill_manual(
    values = c(
      normal_liver = "#55A868",
      HCC          = "#4C72B0"
    ),
    labels = c(
      normal_liver = "Normal liver",
      HCC          = "HCC"
    )
  ) +
  
  theme_bw() +
  
  labs(
    title = "DESeq2 size factors per sample",
    x     = "Sample",
    y     = "Size factor",
    fill  = "Condition"
  ) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )


# Plot dispersion estimates
plotDispEsts(
  dds,
  main = "Dispersion estimates — TEs"
)


# Variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# PCA
pca_res <- prcomp(t(assay(vsd)))

# Percent variance explained
pct_var <- round(
  100 * pca_res$sdev^2 / sum(pca_res$sdev^2),
  1
)

# PCA dataframe
df_pca <- data.frame(
  sample  = colnames(vsd),
  cond    = metadata_qc$cond,
  patient = metadata_qc$patient,
  pca_res$x
)

# PC1 vs PC2
ggplot(
  df_pca,
  aes(PC1, PC2, colour = cond, label = patient)
) +
  
  geom_point(size = 3) +
  
  geom_text_repel(
    size = 2.5,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      normal_liver = "#55A868",
      HCC          = "#4C72B0"
    ),
    labels = c(
      normal_liver = "Normal liver",
      HCC          = "HCC"
    )
  ) +
  
  labs(
    title = "PC1 vs PC2 of VST-normalized TE counts",
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)"),
    colour = "Condition"
  ) +
  
  theme_bw()

##Change the variables
counts_qc2 <- counts_qc
metadata_qc2 <- metadata_qc
dds2 <- dds

# 6. Extract differential expression results ####

#Extract DE results using the cleaned DESeq2 object
res_HCC <- results(
  dds2,
  contrast = c("cond", "HCC", "normal_liver")
)

#Print summary of significant results
summary(res_HCC)

#Show top hits ranked by adjusted p-value
head(res_HCC[order(res_HCC$padj), ])

#Build a summary dataframe of up/downregulated TE counts
df_summary <- data.frame(
  contrast  = rep("HCC vs Normal Liver", each = 2),
  direction = c("Up", "Down"),
  count = c(
    sum(
      res_HCC$padj < 0.05 &
        res_HCC$log2FoldChange > 0,
      na.rm = TRUE
    ),
    sum(
      res_HCC$padj < 0.05 &
        res_HCC$log2FoldChange < 0,
      na.rm = TRUE
    )
  )
)

#Set contrast order for plotting
df_summary$contrast <- factor(
  df_summary$contrast,
  levels = c("HCC vs Normal Liver")
)

#Plot number of up- and downregulated TE subfamilies
ggplot(
  df_summary,
  aes(x = contrast, y = count, fill = direction)
) +
  
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.7),
    width = 0.6
  ) +
  
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.7),
    vjust = -0.4,
    size = 4
  ) +
  
  scale_fill_manual(
    values = c(
      "Up"   = "firebrick",
      "Down" = "steelblue"
    )
  ) +
  
  theme_bw() +
  
  labs(
    title    = "Differentially Expressed TE Subfamilies",
    subtitle = "HCC vs Normal Liver (padj < 0.05)",
    x        = "",
    y        = "Number of TE subfamilies",
    fill     = "Regulation"
  )

#### 8. Filter significant TE subfamilies (padj < 0.05) ####

#Subset results keeping only TEs with adjusted p-value < 0.05
sig_HCC <- subset(res_HCC, padj < 0.05)

#Count total significant TE subfamilies
nrow(sig_HCC)

#Split significant TEs by direction of change

#HCC vs Normal Liver
up_HCC   <- sig_HCC[sig_HCC$log2FoldChange > 0, ]
down_HCC <- sig_HCC[sig_HCC$log2FoldChange < 0, ]

#Print counts of up and downregulated TEs
nrow(up_HCC)
nrow(down_HCC)

#Build summary dataframe for plotting
df_sig <- data.frame(
  comparison = rep("HCC vs Normal Liver", each = 2),
  direction  = c("Up", "Down"),
  count = c(
    nrow(up_HCC),
    nrow(down_HCC)
  )
)

#Set comparison order for plotting
df_sig$comparison <- factor(
  df_sig$comparison,
  levels = c("HCC vs Normal Liver")
)




##Save the signature scores
library(scuttle)

#TPM matrix
tpm_mat <- calculateTPM(
  counts_qc2,
  lengths = rep(1, nrow(counts_qc2))
)

#Load the TE signatures from the discovery cohort
TE_68  <- scan("tumor_associated_68_TEs.txt", what = character())
TE_141 <- scan("metastasis_associated_141_TEs.txt", what = character())

# 68 tumor-associated signature
TE_common_68 <- intersect(TE_68, rownames(tpm_mat))

expr_68 <- tpm_mat[TE_common_68, ]
score_68 <- colMeans(expr_68, na.rm = TRUE)

df_68 <- data.frame(
  sample = names(score_68),
  score  = score_68,
  tissue = metadata_qc2$tissue,   # normal_liver / HCC
  stringsAsFactors = FALSE
)

write.csv(
  df_68,
  "TE_signature68_scores_GSE124535.csv",
  row.names = FALSE
)

# 141 metastasis-associated signature
TE_common_141 <- intersect(TE_141, rownames(tpm_mat))

expr_141 <- tpm_mat[TE_common_141, ]
score_141 <- colMeans(expr_141, na.rm = TRUE)

df_141 <- data.frame(
  sample = names(score_141),
  score  = score_141,
  tissue = metadata_qc2$tissue,
  stringsAsFactors = FALSE
)

write.csv(
  df_141,
  "TE_signature141_scores_GSE124535.csv",
  row.names = FALSE
)
