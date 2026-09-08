#GSE50760

#1) Data loading and metadata preparation

#Read the TE count matrix from file (rows = TE subfamilies, columns = samples)
counts <- read.table(
  "GSE50760_TE_counts_matrix.final.txt",
  header      = TRUE,
  row.names   = 1,
  sep         = "\t",
  check.names = FALSE
)

#Check the dimensions of the matrix
dim(counts)  #1180 TEs x 54 samples

#2) Load and prepare the metadata

#Download the GEO ExpressionSet for dataset GSE50760
library(GEOquery)

gse <- getGEO(
  "GSE50760",
  GSEMatrix = TRUE,
  destdir = getwd()
)

#Extract the sample metadata table (rows = samples, columns = annotations)
meta <- pData(gse[[1]])

#Build a clean metadata dataframe 
metadata <- data.frame(
  GSM        = meta$geo_accession,
  tissue     = sub("tissue: ",     "", meta$`tissue:ch1`),
  ajcc_stage = sub("ajcc stage: ", "", meta$`ajcc stage:ch1`),
  patient    = factor(sub(".*(AMC_[0-9]+).*", "\\1", meta$title)),
  stringsAsFactors = FALSE
)

#Rename tissue labels to shorter, cleaner names
metadata$tissue <- c(
  "primary colorectal cancer"                     = "primary_colon",
  "normal-looking surrounding colonic epithelium" = "normal_colon",
  "metastatic colorectal cancer to the liver"     = "liver_metastasis"
)[metadata$tissue]


#Set rownames of metadata to match column names of the count matrix
rownames(metadata) <- colnames(counts)

#Verify that metadata rows and count matrix columns are aligned
all(rownames(metadata) == colnames(counts))





#2) Quality control (QC)

# a) Tumor purity QC

#Load libraries for data manipulation and visualization
library(dplyr)
library(ggplot2)
library(ggrepel)

#Load tumor purity estimates derived from PUREE (Python/cluster Gustave Roussy)
purity <- read.table("GSE50760_purity_results.tsv", header = TRUE)

#Add purity values to metadata, matched by sample name
metadata$purity <- purity$purity[
  match(rownames(metadata), rownames(purity))
]

#Create ordered factor for tissue condition (used in plots and models)
metadata$cond <- factor(
  metadata$tissue,
  levels = c("normal_colon", "primary_colon", "liver_metastasis")
)

#Summarize purity statistics per tissue group (min, median, max)
purity_table <- metadata %>%
  group_by(cond) %>%
  summarise(
    Samples       = n(),
    Min_Purity    = round(min(purity),    3),
    Median_Purity = round(median(purity), 3),
    Max_Purity    = round(max(purity),    3)
  )

#Identify the sample with the lowest purity value
metadata %>%
  dplyr::select(GSM, cond, patient, purity) %>%
  arrange(purity)

# Create tumor-only metadata for plotting
metadata_tumor <- metadata %>%
  dplyr::filter(cond %in% c("primary_colon", "liver_metastasis"))


# Plot tumor purity distribution per condition (tumor samples only)
ggplot(metadata_tumor, aes(x = cond, y = purity)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", colour = "black") +
  geom_jitter(
    aes(color = cond),
    width = 0.15, size = 3, alpha = 0.8
  ) +
  geom_text_repel(
    data = metadata_tumor[rownames(metadata_tumor) == "SRR975557", ],
    aes(label = rownames(metadata_tumor[rownames(metadata_tumor) == "SRR975557", ])),
    color = "black", size = 3, fontface = "bold"
  ) +
  scale_color_manual(
    values = c(
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    )
  ) +
  scale_x_discrete(
    labels = c(
      primary_colon    = "primary colon",
      liver_metastasis = "liver metastasis"
    )
  ) +
  labs(
    title = "Tumor purity distribution in tumor samples",
    x     = "Condition",
    y     = "Tumor purity (PUREE)"
  ) +
  theme_bw() +
  theme(legend.position = "none")



#Define the low-purity sample to remove (metastatic)
sample_to_remove <- "SRR975557"

#Create logical vector: TRUE = keep, FALSE = remove
keep <- colnames(counts) != sample_to_remove

#Remove the low-purity sample from the count matrix
counts_qc <- counts[, keep]

#Remove the low-purity sample from the metadata
metadata_qc <- metadata[keep, ]

#Recalculate purity summary after QC
purity_table_qc <- metadata_qc %>%
  group_by(cond) %>%
  summarise(
    Samples       = n(),
    Min_Purity    = round(min(purity),    3),
    Median_Purity = round(median(purity), 3),
    Max_Purity    = round(max(purity),    3)
  )


#b) Library size QC

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
    cond = factor(cond, levels = c("normal_colon", "primary_colon", "liver_metastasis"))
  ) %>%
  arrange(cond, lib_size)

# Plot library size per sample, colored by tissue condition
ggplot(
  lib_df,
  aes(x = factor(sample, levels = sample), y = lib_size, fill = cond)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c(
      normal_colon     = "#55A868",
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    ),
    labels = c(
      normal_colon     = "normal colon",
      primary_colon    = "primary colon",
      liver_metastasis = "liver metastasis"
    )
  ) +
  labs(
    title = "GSE50760 — Library size per sample (grouped by tissue)",
    x     = "Sample",
    y     = "Total read counts",
    fill  = "Tissue"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

#Correlation plots(purity of samples x depth of sequencing)

# Tumor purity vs sequencing depth (tumor samples only)

# Calculate sequencing depth
depth <- colSums(counts_qc)

# Build full dataframe
corr_df <- data.frame(
  sample  = names(depth),
  depth   = depth,
  purity  = metadata_qc$purity,
  cond    = metadata_qc$cond,
  patient = metadata_qc$patient,
  stringsAsFactors = FALSE
)

# Keep only tumor samples
corr_df_tumor <- corr_df %>%
  filter(cond %in% c("primary_colon", "liver_metastasis"))

# Compute correlation (tumors only)
r    <- cor(corr_df_tumor$depth, corr_df_tumor$purity, use = "complete.obs")
r2   <- r^2
pval <- cor.test(corr_df_tumor$depth, corr_df_tumor$purity)$p.value
pval_label <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", round(pval, 3)))

# Plot 1: Combined
ggplot(corr_df_tumor, aes(x = depth, y = purity, color = cond)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_x_continuous(labels = scales::comma)+
  geom_smooth(method = "lm", se = TRUE,
              color = "grey40", fill = "grey85", linewidth = 0.8) +
  geom_text_repel(aes(label = patient),
                  size = 2.8, show.legend = FALSE) +
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
    x = Inf, y = Inf,
    hjust = 1.1, vjust = 1.5,
    label = paste0(
      "r = ", round(r, 3), "\n",
      "R² = ", round(r2, 3), "\n",
      pval_label
    ),
    size = 3.5, color = "grey20"
  ) +
  labs(
    title = "Tumor purity vs sequencing depth (tumor samples only)",
    x = "Sequencing depth (total read counts)",
    y = "Tumor purity (PUREE)",
    color = "Condition"
  ) +
  theme_bw()



# Plot 2: Faceted by condition
ggplot(corr_df_tumor, aes(x = depth, y = purity, color = cond)) +
  
  geom_point(size = 3, alpha = 0.85) +
  scale_x_continuous(labels = scales::comma)+
  geom_smooth(method = "lm", se = TRUE,
              color = "grey40", fill = "grey85", linewidth = 0.8) +
  geom_text_repel(aes(label = patient),
                  size = 2.5, show.legend = FALSE) +
  facet_wrap(
    ~ cond,
    labeller = labeller(cond = c(
      primary_colon    = "Primary colon",
      liver_metastasis = "Liver metastasis"
    ))
  ) +
  scale_color_manual(
    values = c(
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    )
  ) +
  labs(
    title = "Tumor purity vs sequencing depth by condition (tumor samples only)",
    x = "Sequencing depth (total read counts)",
    y = "Tumor purity (PUREE)"
  ) +
  theme_bw() +
  theme(legend.position = "none")




#3) DESeq2 Normalization, variance stabilization, PCA and outlier removal

#a) DESeq2 for PCA and outlier detection

#Load the DESeq2 library for differential expression analysis
library(DESeq2)

#Create DESeq2 object: design controls for patient (paired) and tests condition
dds <- DESeqDataSetFromMatrix(
  countData = counts_qc,
  colData   = metadata_qc,
  design    = ~ patient + cond
)

#Run the full DESeq2 pipeline: normalization, dispersion estimation, model fitting
dds <- DESeq(dds)

#Extract and inspect size factors (normalization factors per sample)
sizeFactors(dds)

#Check the range of size factors (should be close to 1 for well-normalized data)
range(sizeFactors(dds))

#Store size factors in a dataframe for plotting
sf_df <- data.frame(
  sample      = names(sizeFactors(dds)),
  size_factor = sizeFactors(dds)
)

#Plot size factors per sample to check normalization consistency
ggplot(sf_df, aes(x = sample, y = size_factor)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_bw() +
  labs(
    title = "DESeq2 size factors per sample",
    x     = "Sample",
    y     = "Size factor"
  ) +
  theme(axis.text.x = element_text(angle = 90))

# Plot dispersion estimates to assess variance-mean relationship
plotDispEsts(dds, main = "Dispersion estimates — TEs")

#Apply variance stabilizing transformation (VST) for PCA and visualization
vsd <- vst(dds, blind = TRUE)

#Perform PCA on the VST-normalized data
pca_res <- prcomp(t(assay(vsd)))

#Compute percentage of variance explained by each principal component
pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

#Build a dataframe for PCA plotting
df_pca <- data.frame(
  sample  = colnames(vsd),
  cond    = metadata_qc$cond,
  patient = metadata_qc$patient,
  pca_res$x
)

#Plot PC1 vs PC2, colored by condition and labeled by patient
ggplot(df_pca, aes(PC1, PC2, colour = cond, label = patient)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5, show.legend = FALSE) +
  scale_color_manual(
    values = c(
      normal_colon     = "#55A868",
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    ),
    breaks = c(
      "normal_colon",
      "primary_colon",
      "liver_metastasis"
    ),
    labels = c(
      "Normal colon",
      "Primary colon",
      "Liver metastases"
    )
  ) +
  labs(
    title = "PCA of TE subfamily expression (GSE50760)",
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)"),
    colour = "Condition"
  ) +
  theme_bw()
##From the PCA, these is an outlier (metastatic AMC_10)


##Remove PCA outlier (AMC_10 metastasis) and rerun DESeq2

#Identify the metastasis sample from patient AMC_10 (PCA outlier)
outlier_sample <- rownames(metadata_qc)[
  metadata_qc$patient == "AMC_10" &
    metadata_qc$cond    == "liver_metastasis"
]

#Remove the outlier sample from the count matrix
counts_qc2 <- counts_qc[, colnames(counts_qc) != outlier_sample]

#Update metadata to match the filtered count matrix
metadata_qc2 <- metadata_qc[colnames(counts_qc2), ]

#Create a new DESeq2 object with the cleaned dataset
dds2 <- DESeqDataSetFromMatrix(
  countData = counts_qc2,
  colData   = metadata_qc2,
  design    = ~ patient + cond
)

#Run DESeq2 pipeline again on the cleaned dataset
dds2 <- DESeq(dds2)

#Apply VST on the cleaned dataset for PCA visualization
vsd2 <- vst(dds2, blind = TRUE)

#Perform PCA on the cleaned VST data
pca2 <- prcomp(t(assay(vsd2)))

#Compute % of variance explained by each PC
pct_var2 <- round(100 * pca2$sdev^2 / sum(pca2$sdev^2), 1)

#Build PCA dataframe for the cleaned dataset
df_pca2 <- data.frame(
  sample  = colnames(vsd2),
  cond    = metadata_qc2$cond,
  patient = metadata_qc2$patient,
  pca2$x
)

#Plot PC1 vs PC2 after outlier removal
# Ensure condition is a factor (important for legend order)
df_pca2$cond <- factor(
  df_pca2$cond,
  levels = c("normal_colon", "primary_colon", "liver_metastasis")
)

ggplot(df_pca2, aes(PC1, PC2, colour = cond, label = patient)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  
  scale_color_manual(
    values = c(
      normal_colon     = "#55A868",
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    ),
    breaks = c(
      "normal_colon",
      "primary_colon",
      "liver_metastasis"
    ),
    labels = c(
      "Normal colon",
      "Primary colon",
      "Liver metastases"
    ),
    name = "Condition"
  ) +
  
  labs(
    title = "PCA of TE Subfamily Expression (GSE50760)",
    x = paste0("PC1 (", pct_var2[1], "%)"),
    y = paste0("PC2 (", pct_var2[2], "%)")
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    axis.text = element_text(
      size = 13
    ),
    legend.title = element_text(
      size = 15,
      face = "bold"
    ),
    legend.text = element_text(
      size = 13
    )
  )


#b) Tumor purity analysis to check its role on the PC1 variation, and specificallly for the primary distribution 

# Add purity into your PCA dataframe (using df_pca2, after outlier removal)
df_pca2$purity <- metadata_qc2$purity[match(df_pca2$sample, rownames(metadata_qc2))]

# Correlation between PC1 and purity, separately per condition -> gives the r, R2, and p-value
cor_by_cond <- df_pca2 %>%
  group_by(cond) %>%
  summarise(
    n     = n(),
    r_PC1 = cor(PC1, purity),
    r2_PC1 = r_PC1^2,
    p_PC1 = cor.test(PC1, purity)$p.value
  )
print(cor_by_cond)

# Check if the purity-PC1 relationship differ significantly by condition -> this is the purity x condition interaction test
model_interaction <- lm(PC1 ~ purity * cond, data = df_pca2)
summary(model_interaction)
anova(model_interaction)

# Does purity explain expression independent of condition, at the TE level? -> likelihood ratio test comparing model with vs without purity
metadata_qc2$purity_scaled <- scale(metadata_qc2$purity)

dds_purity <- DESeqDataSetFromMatrix(
  countData = counts_qc2,
  colData   = metadata_qc2,
  design    = ~ patient + purity_scaled + cond
)
dds_purity <- DESeq(dds_purity)

dds_lrt <- DESeq(dds_purity, test = "LRT", reduced = ~ patient + cond)
res_lrt <- results(dds_lrt)
summary(res_lrt)

# Plot the relationship (PC1 vs purity, colored by condition)
ggplot(df_pca2, aes(x = purity, y = PC1, colour = cond)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_text_repel(aes(label = patient), size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = c(
    normal_colon = "#55A868",
    primary_colon = "#4C72B0",
    liver_metastasis = "#C44E52"
  )) +
  labs(
    title = "PC1 vs tumor purity",
    x = "Tumor purity (PUREE)",
    y = "PC1"
  ) +
  theme_bw()



#4) Extract differential expression results
#Extract DE results for each pairwise comparison using the cleaned DESeq2 model
res_primary <- results(dds2, contrast = c("cond", "primary_colon",    "normal_colon"))
res_meta    <- results(dds2, contrast = c("cond", "liver_metastasis",  "normal_colon"))
res_mvp     <- results(dds2, contrast = c("cond", "liver_metastasis",  "primary_colon"))

#Print summary of significant results for each comparison
summary(res_primary)
summary(res_meta)
summary(res_mvp)

#Show top hits ranked by adjusted p-value for each comparison
head(res_primary[order(res_primary$padj), ])
head(res_meta[order(res_meta$padj), ])
head(res_mvp[order(res_mvp$padj), ])

#Build a summary dataframe of up/downregulated TE counts per comparison
df_summary <- data.frame(
  contrast  = rep(c("Primary Colon vs Normal Colon", "Liver Metastasis vs Normal Colon", "Liver Metastasis vs Primary Colon"), each = 2),
  direction = rep(c("Up", "Down"), 3),
  count = c(
    sum(res_primary$padj < 0.05 & res_primary$log2FoldChange > 0, na.rm = TRUE),
    sum(res_primary$padj < 0.05 & res_primary$log2FoldChange < 0, na.rm = TRUE),
    sum(res_meta$padj    < 0.05 & res_meta$log2FoldChange    > 0, na.rm = TRUE),
    sum(res_meta$padj    < 0.05 & res_meta$log2FoldChange    < 0, na.rm = TRUE),
    sum(res_mvp$padj     < 0.05 & res_mvp$log2FoldChange     > 0, na.rm = TRUE),
    sum(res_mvp$padj     < 0.05 & res_mvp$log2FoldChange     < 0, na.rm = TRUE)
  )
)

#Set contrast order for plotting
df_summary$contrast <- factor(
  df_summary$contrast,
  levels = c("Primary Colon vs Normal Colon", "Liver Metastasis vs Normal Colon", "Liver Metastasis vs Primary Colon")
)

#Plot number of up- and downregulated TE subfamilies per comparison
ggplot(df_summary, aes(x = contrast, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.7),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = c("Up" = "firebrick", "Down" = "steelblue")) +
  theme_bw() +
  labs(
    title    = "Differentially Expressed TE Subfamilies Across Tumor Progression (GSE50760)",
    subtitle = "padj < 0.05",
    x        = "",
    y        = "Number of TE subfamilies",
    fill     = "Regulation"
  )



#5) Filter significant TE subfamilies (padj < 0.05) ####

#Subset results keeping only TEs with adjusted p-value < 0.05
sig_primary <- subset(res_primary, padj < 0.05)
sig_meta    <- subset(res_meta,    padj < 0.05)
sig_mvp     <- subset(res_mvp,     padj < 0.05)

#Count total significant TE subfamilies per comparison
nrow(sig_primary)
nrow(sig_meta)
nrow(sig_mvp)

#Split significant TEs by direction of change (upregulated vs downregulated)
#Primary colon vs Normal colon
up_primary   <- sig_primary[sig_primary$log2FoldChange > 0, ]
down_primary <- sig_primary[sig_primary$log2FoldChange < 0, ]

#Liver metastasis vs Normal colon
up_meta      <- sig_meta[sig_meta$log2FoldChange > 0, ]
down_meta    <- sig_meta[sig_meta$log2FoldChange < 0, ]

#Liver metastasis vs Primary colon
up_mvp       <- sig_mvp[sig_mvp$log2FoldChange > 0, ]
down_mvp     <- sig_mvp[sig_mvp$log2FoldChange < 0, ]

#Print counts of up and downregulated TEs per comparison
nrow(up_primary);   nrow(down_primary)
nrow(up_meta);      nrow(down_meta)
nrow(up_mvp);       nrow(down_mvp)

#Build summary dataframe of up/down counts per comparison for plotting
df_sig <- data.frame(
  comparison = rep(c("Primary Colon vs Normal Colon", "Liver Metastasis vs Normal Colon", "Liver Metastasis vs Primary Colon"), each = 2),
  direction  = rep(c("Up", "Down"), 3),
  count = c(
    nrow(up_primary), nrow(down_primary),
    nrow(up_meta),    nrow(down_meta),
    nrow(up_mvp),     nrow(down_mvp)
  )
)

#Set comparison order for plotting
df_sig$comparison <- factor(
  df_sig$comparison,
  levels = c("Primary Colon vs Normal Colon", "Liver Metastasis vs Normal Colon", "Liver Metastasis vs Primary Colon")
)



#6) Volcano Plots for ALL TEs (up and down regulated)
library(ggplot2)

make_volcano_TE <- function(res, title) {
  
  df <- as.data.frame(res)
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  
  df$status <- factor(
    ifelse(df$padj < 0.05 & df$log2FoldChange > 0, "Up",
           ifelse(df$padj < 0.05 & df$log2FoldChange < 0, "Down", "NS")),
    levels = c("Up", "Down", "NS")
  )
  
  # Count TE subfamilies per direction
  n_up   <- sum(df$status == "Up")
  n_down <- sum(df$status == "Down")
  
  # Text to display on plot
  label_text <- paste0(
    "Up: ", n_up, "\n",
    "Down: ", n_down
  )
  
  print(
    ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
      geom_point(size = 2.2, alpha = 0.85) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_manual(
        values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey75")
      ) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = label_text,
        hjust = 1.05, vjust = 1.2,
        size = 4,
        color = "black"
      ) +
      labs(
        title = title,
        x = "Log2 fold change",
        y = "-log10(padj)",
        color = "Regulation"
      ) +
      theme_bw(base_size = 11)
  )
}

###Generate the 3 volcano plots
make_volcano_TE(
  res_primary,
  "Primary Colon vs Normal Colon — All differentially expressed TE subfamilies (padj < 0.05)"
)

make_volcano_TE(
  res_meta,
  "Liver Metastasis vs Normal Colon — All differentially expressed TE subfamilies (padj < 0.05)"
)

make_volcano_TE(
  res_mvp,
  "Liver Metastasis vs Primary Colon — All differentially expressed TE subfamilies (padj < 0.05)"
)



#7) Venn diagram — overlap of significant TE subfamilies ####

#Load the library
library(ggvenn)

#Make a VEnn with all these
TE_meta_vs_primary_up <- rownames(
  res_mvp[
    !is.na(res_mvp$padj) &
      res_mvp$padj < 0.05 &
      res_mvp$log2FoldChange > 0,
  ]
)

meta_vs_primary_up <- TE_meta_vs_primary_up   

primary_vs_normal_up <- rownames(
  res_primary[
    !is.na(res_primary$padj) &
      res_primary$padj < 0.05 &
      res_primary$log2FoldChange > 0,
  ]
)

meta_vs_normal_up <- rownames(
  res_meta[
    !is.na(res_meta$padj) &
      res_meta$padj < 0.05 &
      res_meta$log2FoldChange > 0,
  ]
)


#Make the plot
ggvenn(
  list(
    "Primary vs normal"      = primary_vs_normal_up,
    "Metastasis vs normal"  = meta_vs_normal_up,
    "Metastasis vs primary" = meta_vs_primary_up
  ),
  fill_color = c(
    "#6BAED6",  # muted steel blue
    "#FDAE6B",  # soft burnt orange
    "#9E9AC8"   # dusty lavender
  ),
  stroke_size = 0.6,
  set_name_size = 4,
  text_size = 4,
  show_percentage = FALSE
) +
  labs(
    title = "Relationship between overexpressed TE subfamilies across tumor progression comparisons ",
    subtitle = "Primary vs normal, metastasis vs normal, and metastasis vs primary (padj < 0.05, log2FC > 0)"
  )


#8) Define the two signatures of interest: tumor-associated (65+3) and metastasis-associated (106+35)

#Tumor-associated
tumor_associated <- intersect(primary_vs_normal_up, meta_vs_normal_up)
length(tumor_associated)   # should be 68

#Metastasis-associated
# 35 metastasis‑specific (met vs primary AND met vs normal AND NOT primary)
meta_specific <- intersect(meta_vs_normal_up, meta_vs_primary_up)
meta_specific <- setdiff(meta_specific, primary_vs_normal_up)

# 106 metastasis‑only vs normal
meta_only_normal <- setdiff(meta_vs_normal_up, primary_vs_normal_up)

# Combine both
metastasis_associated <- union(meta_only_normal, meta_specific)
length(metastasis_associated)   # should be 141




#9) TE class composition of TE subfamilies with RepeatMasker (1180) as background. For the tumor-associated (68)
#and metastasis-associated (141) signatures

# Create the function
get_class_composition <- function(te_set, signature_name){
  
  te_class <- sapply(
    strsplit(te_set, ":"),
    function(x) if(length(x) >= 3) x[3] else "Other"
  )
  
  te_class <- toupper(te_class)
  
  te_class[!te_class %in% c(
    "DNA",
    "LTR",
    "LINE",
    "SINE"
  )] <- "Other"
  
  df <- as.data.frame(
    table(te_class)
  )
  
  colnames(df) <- c(
    "TE_class",
    "Count"
  )
  
  df$Percent <- 100 * df$Count / sum(df$Count)
  
  df$Signature <- signature_name
  
  df
}

#RepeatMasker background
library(rtracklayer)

rmsk <- import("hg38_rmsk_TE.gtf")

rmsk_df <- as.data.frame(rmsk)

rm_control <- rmsk_df %>%
  
  distinct(
    gene_name,
    class_id
  ) %>%
  
  mutate(
    TE_class = toupper(class_id)
  ) %>%
  
  mutate(
    TE_class = ifelse(
      TE_class %in% c(
        "DNA",
        "LTR",
        "LINE",
        "SINE"
      ),
      TE_class,
      "Other"
    )
  ) %>%
  
  dplyr::count(
    TE_class,
    name = "Count"
  ) %>%
  
  mutate(
    Percent = 100 * Count / sum(Count),
    Signature = "RepeatMasker (1180)"
  )


#Tumor-associated signature
TE_tumor_all <- intersect(meta_vs_normal_up, primary_vs_normal_up)

comp_68 <- get_class_composition(
  TE_tumor_all,
  "Tumor-associated signature (68)"
)

#Metastasis-associated signature
comp_141 <- get_class_composition(
  metastasis_associated,
  "Metastasis-associated signature (141)"
)

#Combine
plot_df <- bind_rows(
  rm_control,
  comp_68,
  comp_141
)

#Order
plot_df$TE_class <- factor(
  plot_df$TE_class,
  levels = c(
    "DNA",
    "LINE",
    "LTR",
    "SINE",
    "Other"
  )
)

plot_df$Signature <- factor(
  plot_df$Signature,
  levels = c(
    "RepeatMasker (1180)",
    "Tumor-associated signature (68)",
    "Metastasis-associated signature (141)"
  )
)

#Colors
class_colors <- c(
  DNA   = "#8172B2",
  LINE  = "#4C72B0",
  LTR   = "#55A868",
  SINE  = "#C44E52",
  Other = "#8C8C8C"
)

#Plot
ggplot(
  plot_df,
  aes(
    x = Signature,
    y = Percent,
    fill = TE_class
  )
) +
  
  geom_col(
    width = 0.85,
    colour = "white"
  ) +
  
  geom_text(
    aes(
      label = ifelse(
        Percent >= 3,
        paste0(round(Percent, 1), "%"),
        ""
      )
    ),
    position = position_stack(vjust = 0.5),
    colour = "white",
    fontface = "bold",
    size = 3.5
  ) +
  
  scale_fill_manual(
    values = class_colors,
    breaks = c(
      "DNA",
      "LINE",
      "LTR",
      "SINE",
      "Other"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  
  labs(
    title = "TE Class Composition Across Signatures",
    subtitle = "Percentage of TE subfamilies belonging to each TE class",
    x = NULL,
    y = "Percentage of TE subfamilies",
    fill = "TE Class"
  ) +
  
  theme_bw(base_size = 13) +
  
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      face = "bold",
      size = 11
    ),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.position = "right"
  )


#10) Enrichment analysis of the TE class composition to see it they are statistically enriched/depleted
#compared to the genomic background (RepeatMasker)

#Make a function
extract_class <- function(te_vector){
  
  cls <- sapply(
    strsplit(te_vector, ":"),
    function(x) if(length(x) >= 3) x[3] else "Other"
  )
  
  cls <- toupper(cls)
  
  cls[!cls %in% c(
    "LTR",
    "LINE",
    "SINE",
    "DNA"
  )] <- "Other"
  
  cls
}

#Background

background_tes <- rownames(counts_qc2)

bg_class <- extract_class(
  background_tes
)

#Function
run_enrichment <- function(
    signature_tes,
    signature_name
){
  
  sig_class <- extract_class(
    signature_tes
  )
  
  classes <- c(
    "LTR",
    "LINE",
    "SINE",
    "DNA",
    "Other"
  )
  
  out <- lapply(
    classes,
    function(cls){
      
      a <- sum(sig_class == cls)
      b <- sum(sig_class != cls)
      
      c <- sum(bg_class == cls) - a
      d <- sum(bg_class != cls) - b
      
      ft <- fisher.test(
        matrix(
          c(a,b,c,d),
          nrow = 2
        )
      )
      
      data.frame(
        Signature = signature_name,
        Class = cls,
        Count = a,
        OddsRatio = unname(ft$estimate),
        Pvalue = ft$p.value
      )
    }
  )
  
  do.call(
    rbind,
    out
  )
}

#Run tests
enrichment_results <- rbind(
  
  run_enrichment(
    TE_tumor_all,
    "Tumor-associated (68)"
  ),
  
  run_enrichment(
    metastasis_associated,
    "Metastasis-associated (141)"
  )
  
)

#FDR correction
enrichment_results$FDR <- p.adjust(
  enrichment_results$Pvalue,
  method = "BH"
)

#Sort
enrichment_results <- enrichment_results[
  order(
    enrichment_results$Signature,
    enrichment_results$FDR
  ),
]

print(enrichment_results)


#11) External validation with the other datasets.
#For the external validation we created the mean TPM of the 2 signatures and we compared it with the
#mean TPM of these signatures in the validation cohorts

#Save the tumor-associated signature TPM

# TPM matrix
tpm_mat <- calculateTPM(
  counts_qc2,
  lengths = rep(1, nrow(counts_qc2))
)

expr_68 <- tpm_mat[TE_tumor_all, ]

score_68 <- colMeans(
  expr_68,
  na.rm = TRUE
)

df_68 <- data.frame(
  sample = names(score_68),
  score  = score_68,
  cond   = metadata_qc2$cond
)

write.csv(
  df_68,
  "TE_signature68_scores_GSE50760.csv",
  row.names = FALSE
)


#Save the metastasis-associated signature TPM
expr_141 <- tpm_mat[metastasis_associated, ]

score_141 <- colMeans(
  expr_141,
  na.rm = TRUE
)

df_141 <- data.frame(
  sample = names(score_141),
  score  = score_141,
  cond   = metadata_qc2$cond
)

write.csv(
  df_141,
  "TE_signature141_scores_GSE50760.csv",
  row.names = FALSE
)



#11) External validation with the other datasets.
#For the external validation we created the mean TPM of the 2 signatures and we compared it with the
#mean TPM of these signatures in the validation cohorts

#Save the tumor-associated signature TPM

expr_68 <- tpm_mat[TE_tumor_all, ]

score_68 <- colMeans(
  expr_68,
  na.rm = TRUE
)

df_68 <- data.frame(
  sample = names(score_68),
  score  = score_68,
  cond   = metadata_qc2$cond
)

write.csv(
  df_68,
  "TE_signature68_scores_GSE50760.csv",
  row.names = FALSE
)


#Save the metastasis-associated signature TPM
expr_141 <- tpm_mat[metastasis_associated, ]

score_141 <- colMeans(
  expr_141,
  na.rm = TRUE
)

df_141 <- data.frame(
  sample = names(score_141),
  score  = score_141,
  cond   = metadata_qc2$cond
)

write.csv(
  df_141,
  "TE_signature141_scores_GSE50760.csv",
  row.names = FALSE
)


#11) External validation with the other datasets.
#For the external validation we created the mean TPM of the 2 signatures and we compared it with the
#mean TPM of these signatures in the validation cohorts

#Save the tumor-associated signature TPM

expr_68 <- tpm_mat[TE_tumor_all, ]

score_68 <- colMeans(
  expr_68,
  na.rm = TRUE
)

df_68 <- data.frame(
  sample = names(score_68),
  score  = score_68,
  cond   = metadata_qc2$cond
)

write.csv(
  df_68,
  "TE_signature68_scores_GSE50760.csv",
  row.names = FALSE
)


#Save the metastasis-associated signature TPM
expr_141 <- tpm_mat[metastasis_associated, ]

score_141 <- colMeans(
  expr_141,
  na.rm = TRUE
)

df_141 <- data.frame(
  sample = names(score_141),
  score  = score_141,
  cond   = metadata_qc2$cond
)

write.csv(
  df_141,
  "TE_signature141_scores_GSE50760.csv",
  row.names = FALSE
)

#And as txt
write.table(
  tumor_associated,
  "tumor_associated_68_TEs.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  metastasis_associated,
  "metastasis_associated_141_TEs.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


##Different scripts were created for the signatures scores of the validation cohorts. 
#We will create combined boxplots to see the expression of the signature scores

# Load datasets

df50760_141 <- read.csv("TE_signature141_scores_GSE50760.csv")
df144259_141 <- read.csv("TE_signature141_scores_GSE144259.csv")
df207194_141 <- read.csv("TE_signature141_scores_GSE207194_QC_TREATED_UNTREATED.csv")
df124535_141 <- read.csv("TE_signature141_scores_GSE124535.csv")

df50760_68 <- read.csv("TE_signature68_scores_GSE50760.csv")
df144259_68 <- read.csv("TE_signature68_scores_GSE144259.csv")
df207194_68 <- read.csv("TE_signature68_scores_GSE207194_QC_TREATED_UNTREATED.csv")
df124535_68 <- read.csv("TE_signature68_scores_GSE124535.csv")

#Add cohort labels
add_cohort <- function(df, cohort) {
  df$cohort <- cohort
  df
}

df50760_141 <- add_cohort(df50760_141, "GSE50760")
df144259_141 <- add_cohort(df144259_141, "GSE144259")
df207194_141 <- add_cohort(df207194_141, "GSE207194")
df124535_141 <- add_cohort(df124535_141, "GSE124535")

df50760_68 <- add_cohort(df50760_68, "GSE50760")
df144259_68 <- add_cohort(df144259_68, "GSE144259")
df207194_68 <- add_cohort(df207194_68, "GSE207194")
df124535_68 <- add_cohort(df124535_68, "GSE124535")

# Tissue recoding 
recode_tissue <- function(df) {
  
  
  if ("cond" %in% names(df)) {
    recode_map <- c(
      normal_colon     = "Normal",
      primary_colon    = "Primary",
      liver_metastasis = "Liver metastasis",
      normal_liver     = "Normal liver",
      HCC              = "HCC"
    )
    df$tissue <- recode_map[df$cond]
    return(df)
  }
  
  if ("tissue" %in% names(df)) {
    df$tissue <- recode(df$tissue,
                        primary_colon    = "Primary",
                        liver_metastasis = "Liver metastasis",
                        normal_liver     = "Normal liver",
                        HCC              = "HCC")
    return(df)
  }
  
  stop("ERROR: dataframe has neither 'cond' nor 'tissue'.")
}

df50760_141 <- recode_tissue(df50760_141)
df144259_141 <- recode_tissue(df144259_141)
df207194_141 <- recode_tissue(df207194_141)
df124535_141 <- recode_tissue(df124535_141)

df50760_68 <- recode_tissue(df50760_68)
df144259_68 <- recode_tissue(df144259_68)
df207194_68 <- recode_tissue(df207194_68)
df124535_68 <- recode_tissue(df124535_68)


# Recode treatment for GSE207194
df207194_141$treatment <- recode(df207194_141$treatment_group,
                                 Untreated = "Untreated",
                                 Treated   = "Treated")

df207194_68$treatment <- recode(df207194_68$treatment_group,
                                Untreated = "Untreated",
                                Treated   = "Treated")

df50760_141$treatment <- NA
df144259_141$treatment <- NA
df124535_141$treatment <- NA

df50760_68$treatment <- NA
df144259_68$treatment <- NA
df124535_68$treatment <- NA

# Merge datasets
all141 <- bind_rows(df50760_141, df144259_141, df207194_141, df124535_141)
all68  <- bind_rows(df50760_68,  df144259_68,  df207194_68,  df124535_68)


# Build x-axis groups
build_groups <- function(df) {
  df$x_group <- NA_character_
  
  idx_other <- df$cohort != "GSE207194"
  df$x_group[idx_other] <- paste(df$cohort[idx_other], df$tissue[idx_other], sep = "_")
  
  idx_u <- df$cohort == "GSE207194" & df$treatment == "Untreated"
  df$x_group[idx_u] <- paste("GSE207194_Untreated", df$tissue[idx_u], sep = "_")
  
  idx_t <- df$cohort == "GSE207194" & df$treatment == "Treated"
  df$x_group[idx_t] <- paste("GSE207194_Treated", df$tissue[idx_t], sep = "_")
  
  df$x_group <- factor(df$x_group, levels = c(
    "GSE50760_Normal",
    "GSE50760_Primary",
    "GSE50760_Liver metastasis",
    "GSE144259_Normal",
    "GSE144259_Primary",
    "GSE144259_Liver metastasis",
    "GSE207194_Untreated_Primary",
    "GSE207194_Untreated_Liver metastasis",
    "GSE207194_Treated_Primary",
    "GSE207194_Treated_Liver metastasis",
    "GSE124535_Normal liver",
    "GSE124535_HCC"
  ))
  
  df %>% filter(!is.na(x_group))
}

all141 <- build_groups(all141)
all68  <- build_groups(all68)


# Colors 
color_map <- c(
  "Normal"           = "#2E8B57",  # green
  "Primary"          = "#1F78B4",  # blue
  "Liver metastasis" = "#D73027",  # red
  "Normal liver"     = "#E6AB02",  # yellow
  "HCC"              = "#984EA3"   # purple
)

# Plot function (LOG10 scale, log ticks, custom y label)
make_plot <- function(df, title, y_label) {
  df$tissue <- factor(df$tissue,
                      levels = c("Normal", "Primary", "Liver metastasis", "Normal liver", "HCC"))
  
  top_labels <- data.frame(
    label = c("GSE50760", "GSE144259", "GSE207194\nUntreated", "GSE207194\nTreated", "GSE124535"),
    x_mid = c(2, 5, 7.5, 9.5, 11.5)
  )
  
  ggplot(df, aes(x = x_group, y = score, fill = tissue)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.75, color = "black") +
    scale_y_log10() +
    annotation_logticks(sides = "l") +
    geom_vline(xintercept = c(3.5, 6.5, 10.5), linetype = "dashed", color = "black") +
    geom_vline(xintercept = 8.5, linetype = "dashed", color = "#2E8B57", linewidth = 1) +
    geom_text(data = top_labels,
              aes(x = x_mid, y = max(df$score, na.rm = TRUE) * 1.1, label = label),
              inherit.aes = FALSE, fontface = "bold", size = 4) +
    scale_fill_manual(values = color_map,
                      breaks = c("Normal", "Primary", "Liver metastasis", "Normal liver", "HCC"),
                      labels = c("Normal colon", "Primary colon", "Liver metastasis",
                                 "Normal liver", "HCC"),
                      name = "Tissue") +
    scale_x_discrete(labels = function(x) sub(".*_", "", x)) +
    theme_bw() +
    labs(title = title, x = "", y = y_label) +
    theme(axis.text.x = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "horizontal",
          plot.title = element_text(face = "bold", size = 14))
}


# Final plots 
p141 <- make_plot(
  all141,
  title   = "141-TE Metastasis-Associated Signature Across Cohorts",
  y_label = "Mean TPM across 141 TE subfamilies (log10 scale)"
)

p68 <- make_plot(
  all68,
  title   = "68-TE Tumor-Associated Signature Across Cohorts",
  y_label = "Mean TPM across 68 TE subfamilies (log10 scale)"
)

print(p141)
print(p68)
