#GSE50760_Individual TEs

#1) Load the data
counts <- read.delim(
  "TE_individual_counts_matrix.txt",
  check.names = FALSE
)

rownames(counts) <- counts$TE
counts$TE <- NULL


#2) Metadata
library(GEOquery)
library(dplyr)

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
  "primary colorectal cancer"                    = "primary_colon",
  "normal-looking surrounding colonic epithelium" = "normal_colon",
  "metastatic colorectal cancer to the liver"    = "liver_metastasis"
)

metadata$cond <- factor(
  metadata$tissue,
  levels = c("normal_colon", "primary_colon", "liver_metastasis")
)

rownames(metadata) <- colnames(counts)


#3) Remove the low-purity sample (SRR975557)

counts_qc   <- counts[, colnames(counts) != "SRR975557"]
metadata_qc <- metadata[rownames(metadata) != "SRR975557", ]

stopifnot(all(colnames(counts_qc) == rownames(metadata_qc)))

#4) Check the library size
library(ggplot2)

lib_size <- colSums(counts_qc)

lib_df <- data.frame(
  sample   = names(lib_size),
  lib_size = lib_size,
  cond     = metadata_qc$cond,
  stringsAsFactors = FALSE
) %>%
  mutate(cond = factor(
    cond,
    levels = c("normal_colon", "primary_colon", "liver_metastasis")
  )) %>%
  arrange(cond, lib_size)

summary(lib_size)

#Plot
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
      normal_colon     = "Normal colon",
      primary_colon    = "Primary colon",
      liver_metastasis = "Liver metastasis"
    )
  ) +
  labs(
    title = "Individual TE loci - library size per sample",
    x     = "Sample",
    y     = "Total TE counts",
    fill  = "Tissue"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#5) Expression filtering assessment

library(ggpubr)

thresholds  <- c(1, 2, 5, 10)
min_samples <- ceiling(0.20 * ncol(counts_qc))
n_total     <- nrow(counts_qc)

filter_summary <- do.call(
  rbind,
  lapply(thresholds, function(th) {
    
    keep <- rowSums(counts_qc >= th) >= min_samples
    
    data.frame(
      `Filtering threshold` = paste0(
        "\u2265", th,
        ifelse(th == 1, " count", " counts"),
        " in \u226520% of samples"
      ),
      `TE loci retained` = format(sum(keep), big.mark = ","),
      `TE loci removed`  = format(n_total - sum(keep), big.mark = ","),
      `Retention (%)`    = round(100 * sum(keep) / n_total, 1),
      check.names = FALSE
    )
  })
)

tbl <- ggtexttable(
  filter_summary,
  rows = NULL,
  theme = ttheme(
    "light",
    base_size = 13
  )
)

tbl <- annotate_figure(
  tbl,
  top = text_grob(
    "Impact of expression filtering on individual TE loci",
    face = "bold",
    size = 15
  )
)

tbl


#6) Filter with >= 5 counts in >= 20% of samples (for the expression)
min_samples <- ceiling(0.20 * ncol(counts_qc))
keep        <- rowSums(counts_qc >= 5) >= min_samples
counts_filt <- counts_qc[keep, ]

#Print how many they remain after the filtering 
cat(
  "TE loci retained after filtering:",
  nrow(counts_filt),
  "of",
  nrow(counts_qc),
  "\n"
)

#7) PCA for outlier detection

library(DESeq2)
library(ggrepel)

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filt),
  colData   = metadata_qc,
  design    = ~ patient + cond
)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

vst_mat    <- assay(vsd)
pca        <- prcomp(t(vst_mat), scale. = FALSE)
percentVar <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

pca_df <- data.frame(
  PC1     = pca$x[, 1],
  PC2     = pca$x[, 2],
  cond    = metadata_qc$cond,
  patient = metadata_qc$patient
)

ggplot(pca_df, aes(PC1, PC2, color = cond)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = patient), size = 3, max.overlaps = Inf) +
  theme_bw() +
  scale_color_manual(
    values = c(
      normal_colon     = "#55A868",
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    ),
    labels = c("Normal colon", "Primary colon", "Liver metastasis")
  ) +
  labs(
    title    = "PCA of individual TE loci (VST-normalized)",
    subtitle = paste0(
      "PC1: ", round(percentVar[1], 1),
      "% | PC2: ", round(percentVar[2], 1), "%"
    ),
    x     = paste0("PC1 (", round(percentVar[1], 1), "%)"),
    y     = paste0("PC2 (", round(percentVar[2], 1), "%)"),
    color = "Tissue"
  )


#9) Remove the AMC_10 outlier
sample_to_remove <- rownames(metadata_qc)[
  metadata_qc$patient == "AMC_10" &
    metadata_qc$cond == "liver_metastasis"
]

counts_qc <- counts_qc[, !colnames(counts_qc) %in% sample_to_remove]

metadata_qc <- metadata_qc[
  !rownames(metadata_qc) %in% sample_to_remove,
]

counts_filt <- counts_filt[, rownames(metadata_qc)]

stopifnot(all(colnames(counts_filt) == rownames(metadata_qc)))


#10) PCA after the outleir removal
stopifnot(all(colnames(counts_filt) == rownames(metadata_qc)))

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filt),
  colData   = metadata_qc,
  design    = ~ patient + cond
)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

vst_mat <- assay(vsd)

pca <- prcomp(
  t(vst_mat),
  center = TRUE,
  scale. = FALSE
)

percentVar <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

pca_df <- data.frame(
  PC1     = pca$x[, 1],
  PC2     = pca$x[, 2],
  cond    = metadata_qc$cond,
  patient = metadata_qc$patient
)

ggplot(pca_df, aes(PC1, PC2, color = cond)) +
  geom_point(size = 4) +
  geom_text_repel(
    aes(label = patient),
    size = 4
  ) +
  scale_color_manual(
    values = c(
      normal_colon     = "#55A868",
      primary_colon    = "#4C72B0",
      liver_metastasis = "#C44E52"
    ),
    labels = c(
      "Normal colon",
      "Primary colon",
      "Liver metastasis"
    )
  ) +
  theme_bw() +
  labs(
    title = "PCA of Individual TE Loci (GSE50760)",
    x = paste0("PC1 (", round(percentVar[1], 1), "%)"),
    y = paste0("PC2 (", round(percentVar[2], 1), "%)"),
    color = "Tissue"
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      face = "bold"
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



#11) Differential expression analysis (DESeq2)
library(DESeq2)

stopifnot(all(colnames(counts_filt) == rownames(metadata_qc)))

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filt),
  colData   = metadata_qc,
  design    = ~ patient + cond
)

dds <- DESeq(dds)

# Differential expression contrasts
res_primary <- results(
  dds,
  contrast = c("cond", "primary_colon", "normal_colon")
)

res_meta <- results(
  dds,
  contrast = c("cond", "liver_metastasis", "normal_colon")
)

res_mvp <- results(
  dds,
  contrast = c("cond", "liver_metastasis", "primary_colon")
)

summary(res_primary)
summary(res_meta)
summary(res_mvp)


# Summary of differentially expressed TE loci
# (padj < 0.05 & |log2FC| > 2)
df_summary <- data.frame(
  contrast = rep(
    c(
      "Primary vs Normal",
      "Metastasis vs Normal",
      "Metastasis vs Primary"
    ),
    each = 2
  ),
  direction = rep(c("Up", "Down"), 3),
  count = c(
    sum(res_primary$padj < 0.05 & res_primary$log2FoldChange > 2, na.rm = TRUE),
    sum(res_primary$padj < 0.05 & res_primary$log2FoldChange < -2, na.rm = TRUE),
    sum(res_meta$padj < 0.05 & res_meta$log2FoldChange > 2, na.rm = TRUE),
    sum(res_meta$padj < 0.05 & res_meta$log2FoldChange < -2, na.rm = TRUE),
    sum(res_mvp$padj < 0.05 & res_mvp$log2FoldChange > 2, na.rm = TRUE),
    sum(res_mvp$padj < 0.05 & res_mvp$log2FoldChange < -2, na.rm = TRUE)
  )
)

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
    vjust = -0.3,
    size = 4
  ) +
  scale_fill_manual(values = c(Up = "#D62728", Down = "#1F77B4")) +
  theme_bw() +
  labs(
    title = "Differentially expressed TE loci",
    subtitle = paste0(
      "Paired analysis (patient-adjusted)\n",
      "Filtered loci: \u22655 counts in \u226520% of samples | ",
      "padj < 0.05 & |log2FC| > 2"
    ),
    x = "",
    y = "Number of TE loci",
    fill = "Regulation"
  )


# Volcano plots with |log2FC| > 2 threshold
make_volcano <- function(res, title) {
  
  df <- as.data.frame(res)
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  
  df$status <- factor(
    ifelse(df$padj < 0.05 & df$log2FoldChange > 2, "Up",
           ifelse(df$padj < 0.05 & df$log2FoldChange < -2, "Down", "NS")),
    levels = c("Up", "Down", "NS")
  )
  
  n_up   <- sum(df$status == "Up")
  n_down <- sum(df$status == "Down")
  
  ggplot(df, aes(log2FoldChange, -log10(padj), color = status)) +
    geom_point(size = 2, alpha = 0.75) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
    annotate(
      "text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
      size = 4, fontface = "bold",
      label = paste0(
        "Up: ", format(n_up, big.mark = ","), "\n",
        "Down: ", format(n_down, big.mark = ",")
      )
    ) +
    scale_color_manual(values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey75")) +
    theme_bw(base_size = 11) +
    labs(
      title = title,
      subtitle = "padj < 0.05 & |log2FC| > 2",
      x = "Log2 fold change",
      y = "-log10(adjusted p-value)",
      color = "Regulation"
    )
}

make_volcano(res_primary, "Volcano plot - Primary Colon vs Normal Colon")
make_volcano(res_meta,    "Volcano plot - Liver Metastasis vs Normal Colon")
make_volcano(res_mvp,     "Volcano plot - Liver Metastasis vs Primary Colon")



#From the volcano plots we can see that maybe there are two populations in the upregulated side
#We will use the shrinkage method to check the confidence of the log2fc

# LOG2FC shrinkage analysis
library(apeglm)
library(ashr)
library(ggplot2)


# Check available coefficients
resultsNames(dds)

# Shrinkage
#Liver vs Normal
res_meta_shrink <- lfcShrink(
  dds,
  coef = "cond_liver_metastasis_vs_normal_colon",
  type = "apeglm"
)

#Liver vs Primary
res_mvp_shrink <- lfcShrink(
  dds,
  contrast = c(
    "cond",
    "liver_metastasis",
    "primary_colon"
  ),
  type = "ashr"
)

# Primary vs Normal
res_primary_shrink <- lfcShrink(
  dds,
  coef = "cond_primary_colon_vs_normal_colon",
  type = "apeglm"
)


# Volcano function, to see the volcanos after the shrinkage
make_volcano <- function(res, title){
  
  df <- as.data.frame(res)
  
  df <- df[
    !is.na(df$padj) &
      !is.na(df$log2FoldChange),
  ]
  
  df$status <- factor(
    ifelse(
      df$padj < 0.05 & df$log2FoldChange > 2,
      "Up",
      ifelse(
        df$padj < 0.05 & df$log2FoldChange < -2,
        "Down",
        "NS"
      )
    ),
    levels = c("Up","Down","NS")
  )
  
  n_up <- sum(df$status == "Up")
  n_down <- sum(df$status == "Down")
  
  ggplot(
    df,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      colour = status
    )
  ) +
    geom_point(
      size = 1.8,
      alpha = 0.7
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = c(-2,2),
      linetype = "dashed",
      colour = "grey50"
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.1,
      vjust = 1.5,
      label = paste0(
        "Up: ", format(n_up, big.mark = ","), "\n",
        "Down: ", format(n_down, big.mark = ",")
      ),
      fontface = "bold"
    ) +
    scale_colour_manual(
      values = c(
        Up = "#D62728",
        Down = "#1F77B4",
        NS = "grey75"
      )
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = title,
      x = "log2 fold change",
      y = "-log10(adjusted p-value)",
      colour = "Regulation"
    )
}


#Original volcanos
make_volcano(
  res_meta,
  "Metastasis vs Normal (Original)"
)

make_volcano(
  res_mvp,
  "Metastasis vs Primary (Original)"
)

# Shrunk volcanos 
make_volcano(
  res_meta_shrink,
  "Metastasis vs Normal (Shrunk)"
)

make_volcano(
  res_mvp_shrink,
  "Metastasis vs Primary (Shrunk)"
)

##We keep the TEs that remain after the shrinkage



#12) Venn diagram: overlap of overexpressed TE loci (shrunk)
library(ggvenn)
library(ggplot2)

primary_vs_normal_up <- rownames(
  subset(
    as.data.frame(res_primary_shrink),
    !is.na(padj) & padj < 0.05 & log2FoldChange > 2
  )
)

meta_vs_normal_up <- rownames(
  subset(
    as.data.frame(res_meta_shrink),
    !is.na(padj) & padj < 0.05 & log2FoldChange > 2
  )
)

meta_vs_primary_up <- rownames(
  subset(
    as.data.frame(res_mvp_shrink),
    !is.na(padj) & padj < 0.05 & log2FoldChange > 2
  )
)

#Plot
ggvenn(
  list(
    "Primary vs Normal"     = primary_vs_normal_up,
    "Metastasis vs Normal"  = meta_vs_normal_up,
    "Metastasis vs Primary" = meta_vs_primary_up
  ),
  fill_color = c(
    "#1F77B4",
    "#D62728",
    "#9467BD"
  ),
  stroke_size = 0.5,
  set_name_size = 5,
  text_size = 5,
  show_percentage = FALSE
) +
  labs(
    title = "Overlap of Overexpressed TE Loci",
    subtitle = "Shrunk results: padj < 0.05 and log2FC > 2"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )


# Define TE signatures using SHRUNK results

# Tumor-associated
# Upregulated in both Primary vs Normal and Metastasis vs Normal
TE_tumor_associated <- intersect(
  primary_vs_normal_up,
  meta_vs_normal_up
)

# Shared by all three contrasts
TE_shared_all <- Reduce(
  intersect,
  list(
    primary_vs_normal_up,
    meta_vs_normal_up,
    meta_vs_primary_up
  )
)

# Metastasis vs Normal only
TE_metastasis_vs_normal_only <- setdiff(
  meta_vs_normal_up,
  union(
    primary_vs_normal_up,
    meta_vs_primary_up
  )
)

# Metastasis only
# Upregulated in Metastasis vs Normal and Metastasis vs Primary
# but NOT in Primary vs Normal
TE_metastasis_only <- setdiff(
  intersect(
    meta_vs_normal_up,
    meta_vs_primary_up
  ),
  primary_vs_normal_up
)

# Combined metastasis-associated signature
TE_metastasis_associated <- union(
  TE_metastasis_vs_normal_only,
  TE_metastasis_only
)

# Primary-only
TE_primary_only <- setdiff(
  primary_vs_normal_up,
  union(
    meta_vs_normal_up,
    meta_vs_primary_up
  )
)



#13) Signature scores (log10 scale)

#We need to calculate the TPM and we need the exact length of the individual TEs, so we will use one 
#of the sample files
te_length_file <- "C:/Users/styli/Downloads/SRR975551_counts.txt"
te_lengths_raw <- read.delim(
  te_length_file,
  skip = 1,
  check.names = FALSE
)

te_lengths_df <- data.frame(
  TE     = te_lengths_raw$Geneid,
  Length = te_lengths_raw$Length,
  stringsAsFactors = FALSE
)

#Check if there are any missing lengths
n_missing_length <- sum(!(rownames(counts_filt) %in% te_lengths_df$TE))
cat(
  "TE loci in counts_filt with no matching length entry:",
  n_missing_length, "out of", nrow(counts_filt), "\n"
)

#Compute the TPM of the TE loci
te_lengths <- te_lengths_df$Length[
  match(rownames(counts_filt), te_lengths_df$TE)
]

names(te_lengths) <- rownames(counts_filt)
valid <- !is.na(te_lengths)

counts_tpm <- counts_filt[valid, ]
te_lengths <- te_lengths[valid]

rpk <- sweep(
  counts_tpm,
  1,
  te_lengths / 1000,
  "/"
)

tpm <- sweep(
  rpk,
  2,
  colSums(rpk) / 1e6,
  "/"
)


#Function to calculate the TPM score of each signature
plot_signature_score_log <- function(te_set, title_text) {
  
  # Keep only loci present in TPM matrix
  te_set <- intersect(
    te_set,
    rownames(tpm)
  )
  
  expr <- tpm[
    te_set,
    ,
    drop = FALSE
  ]
  
  # Mean TPM across signature loci
  score <- colMeans(
    expr,
    na.rm = TRUE
  )
  
  df <- data.frame(
    sample = names(score),
    score = score,
    stringsAsFactors = FALSE
  )
  
  df$cond <- metadata_qc[
    df$sample,
    "cond"
  ]
  
  ggplot(
    df,
    aes(
      x = cond,
      y = score,
      fill = cond
    )
  ) +
    
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA
    ) +
    
    geom_jitter(
      width = 0.1,
      size = 2
    ) +
    
    scale_fill_manual(
      values = c(
        normal_colon     = "#55A868",
        primary_colon    = "#4C72B0",
        liver_metastasis = "#C44E52"
      )
    ) +
    
    scale_x_discrete(
      labels = c(
        normal_colon     = "Normal Colon",
        primary_colon    = "Primary Colon",
        liver_metastasis = "Liver Metastasis"
      )
    ) +
    
    scale_y_log10(
      expand = expansion(mult = c(0.05, 0.15))
    ) +
    
    labs(
      title = title_text,
      x = NULL,
      y = "Mean TPM across signature loci (log10 scale)"
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      legend.position = "none",
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      
      # KEEP MAJOR GRID LINES
      panel.grid.major = element_line(
        colour = "grey85",
        linewidth = 0.6
      ),
      
      # KEEP MINOR GRID LINES
      panel.grid.minor = element_line(
        colour = "grey92",
        linewidth = 0.4
      )
    )
}


# Tumor-associated signature
p_tumor <- plot_signature_score_log(
  TE_tumor_associated,
  paste0(
    "Tumor-Associated TE Signature (Shrunk)\n(n = ",
    length(TE_tumor_associated),
    " loci)"
  )
)

print(p_tumor)


#Metastasis-associated
p_meta <- plot_signature_score_log(
  TE_metastasis_associated,
  paste0(
    "Metastasis-Associated TE Signature\n(n = ",
    length(TE_metastasis_associated),
    " loci)"
  )
)

print(p_meta)



#14) TE class composition of TE signatures compared to the RepeatMasker
library(tidyr)
library(gridExtra)
library(grid)
library(stringr)

#First we need to create the annotation table
gtf <- read.delim(
  "hg38_rmsk_TE.gtf",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

te_annotation <- gtf %>%
  transmute(
    chrom = V1,
    start = V4,
    end   = V5,
    
    TE = str_match(
      V9,
      "transcript_id ([^;]+)"
    )[,2],
    
    TE_family = str_match(
      V9,
      "family_id ([^;]+)"
    )[,2],
    
    TE_class = str_match(
      V9,
      "class_id ([^;]+)"
    )[,2]
  ) %>%
  mutate(
    across(
      c(TE, TE_family, TE_class),
      trimws
    )
  ) %>%
  distinct(TE, .keep_all = TRUE)


# Build the annotation table
signature_annotation <- bind_rows(
  
  te_annotation %>%
    mutate(Signature = "RepeatMasker"),
  
  te_annotation %>%
    filter(TE %in% TE_tumor_associated) %>%
    mutate(Signature = "Tumor-Associated"),
  
  te_annotation %>%
    filter(TE %in% TE_metastasis_associated) %>%
    mutate(Signature = "Metastasis-Associated")
)

# Class composition
class_composition <- signature_annotation %>%
  dplyr::count(Signature, TE_class) %>%
  group_by(Signature) %>%
  mutate(
    Percent = 100 * n / sum(n)
  ) %>%
  ungroup()


# Keep only major classes
class_composition_plot <- class_composition %>%
  filter(
    TE_class %in% c(
      "DNA",
      "LTR",
      "LINE",
      "SINE"
    )
  )

# Order signatures
class_composition_plot$Signature <- factor(
  class_composition_plot$Signature,
  levels = c(
    "RepeatMasker",
    "Tumor-Associated",
    "Metastasis-Associated"
  )
)


# Order TE classes
class_composition_plot$TE_class <- factor(
  class_composition_plot$TE_class,
  levels = c(
    "DNA",
    "LTR",
    "LINE",
    "SINE"
  )
)

# Colors
class_colors <- c(
  DNA  = "#8172B2",
  LTR  = "#55A868",
  LINE = "#4C72B0",
  SINE = "#C44E52"
)


# Labels with number of loci
sig_labels <- c(
  "RepeatMasker" = paste0(
    "RepeatMasker\n(n=",
    format(nrow(te_annotation), big.mark = ","),
    ")"
  ),
  
  "Tumor-Associated" = paste0(
    "Tumor-Associated\n(n=",
    length(TE_tumor_associated),
    ")"
  ),
  
  "Metastasis-Associated" = paste0(
    "Metastasis-Associated\n(n=",
    length(TE_metastasis_associated),
    ")"
  )
)

# Plot
ggplot(
  class_composition_plot,
  aes(
    x = Signature,
    y = Percent,
    fill = TE_class
  )
) +
  geom_col(
    width = 0.8,
    colour = "white",
    linewidth = 0.4
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
    size = 4
  ) +
  scale_fill_manual(
    values = class_colors
  ) +
  scale_x_discrete(
    labels = sig_labels
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    title = "TE Class Composition of Tumor- and Metastasis-Associated Signatures",
    subtitle = "Compared with the RepeatMasker background",
    x = NULL,
    y = "Percentage of TE loci",
    fill = "TE Class"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(
      face = "bold",
      size = 11
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )


#15) Statistical test for the class composition vs repeatmasker

# Fisher's Exact Test
# Tumor-Associated Signature vs RepeatMasker Background
classes <- c(
  "SINE",
  "LINE",
  "LTR",
  "DNA"
)

tumor_class_enrichment <- lapply(
  classes,
  function(cl) {
    
    # In signature + class
    a <- sum(
      te_annotation$TE %in% TE_tumor_associated &
        te_annotation$TE_class == cl
    )
    
    # In signature + not class
    b <- sum(
      te_annotation$TE %in% TE_tumor_associated &
        te_annotation$TE_class != cl
    )
    
    # Not in signature + class
    c <- sum(
      !(te_annotation$TE %in% TE_tumor_associated) &
        te_annotation$TE_class == cl
    )
    
    # Not in signature + not class
    d <- sum(
      !(te_annotation$TE %in% TE_tumor_associated) &
        te_annotation$TE_class != cl
    )
    
    mat <- matrix(
      c(a, b,
        c, d),
      nrow = 2,
      byrow = TRUE
    )
    
    ft <- fisher.test(mat)
    
    data.frame(
      TE_Class = cl,
      Tumor_Count = a,
      Tumor_Percent = round(
        100 * a / length(TE_tumor_associated),
        1
      ),
      RepeatMasker_Percent = round(
        100 * sum(te_annotation$TE_class == cl) /
          nrow(te_annotation),
        1
      ),
      Odds_Ratio = round(
        unname(ft$estimate),
        3
      ),
      Pvalue = ft$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(
      Pvalue,
      method = "BH"
    )
  ) %>%
  arrange(FDR)

tumor_class_enrichment


# Fisher's Exact Test
# Metastasis-Associated Signature vs RepeatMasker Background

classes <- c(
  "SINE",
  "LINE",
  "LTR",
  "DNA"
)

metastasis_class_enrichment <- lapply(
  classes,
  function(cl) {
    
    a <- sum(
      te_annotation$TE %in% TE_metastasis_associated &
        te_annotation$TE_class == cl
    )
    
    b <- sum(
      te_annotation$TE %in% TE_metastasis_associated &
        te_annotation$TE_class != cl
    )
    
    c <- sum(
      !(te_annotation$TE %in% TE_metastasis_associated) &
        te_annotation$TE_class == cl
    )
    
    d <- sum(
      !(te_annotation$TE %in% TE_metastasis_associated) &
        te_annotation$TE_class != cl
    )
    
    ft <- fisher.test(
      matrix(
        c(a, b,
          c, d),
        nrow = 2,
        byrow = TRUE
      )
    )
    
    data.frame(
      TE_Class = cl,
      Metastasis_Count = a,
      Metastasis_Percent = round(
        100 * a / length(TE_metastasis_associated),
        1
      ),
      RepeatMasker_Percent = round(
        100 * sum(te_annotation$TE_class == cl) /
          nrow(te_annotation),
        1
      ),
      Odds_Ratio = round(
        unname(ft$estimate),
        3
      ),
      Pvalue = ft$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(
      Pvalue,
      method = "BH"
    )
  ) %>%
  arrange(FDR)

metastasis_class_enrichment



#16) Genomic Distribution

#Homer: create bed files from the final DESeq2 Signatures

# All annotated TE loci 
write.table(
  te_annotation[, c("chrom", "start", "end", "TE")],
  file = "TEs.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


# Tumor-Associated
te_annotation %>%
  filter(TE %in% TE_tumor_associated) %>%
  select(chrom, start, end, TE) %>%
  write.table(
    file = "Tumor_Associated.bed",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

# Metastasis-Associated
te_annotation %>%
  filter(TE %in% TE_metastasis_associated) %>%
  select(chrom, start, end, TE) %>%
  write.table(
    file = "Metastasis_Associated.bed",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

# RepeatMasker
write.table(
  te_annotation[, c("chrom", "start", "end", "TE")],
  file = "TEs.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

##17) The HOMER was performed in the HPC and the resulting files where imported back to R fow downstream
#analysis

repeatmasker <- read.delim("RepeatMasker_homer.txt", check.names = FALSE) %>%
  mutate(Signature = "RepeatMasker")

tumor <- read.delim("Tumor_Associated_homer.txt", check.names = FALSE) %>%
  mutate(Signature = "Tumor-Associated")

met_assoc <- read.delim("Metastasis_Associated_homer_FIXED.txt", check.names = FALSE) %>%
  mutate(Signature = "Metastasis-Associated")

#Combine them
homer_df <- bind_rows(
  repeatmasker,
  tumor,
  met_assoc)

homer_df <- homer_df %>%
  mutate(
    Location = case_when(
      str_detect(Annotation, "^promoter-TSS") ~ "Promoter",
      str_detect(Annotation, "^intron") ~ "Intron",
      str_detect(Annotation, "^exon") ~ "Exon",
      str_detect(Annotation, "^3' UTR") ~ "3' UTR",
      str_detect(Annotation, "^5' UTR") ~ "5' UTR",
      str_detect(Annotation, "^TTS") ~ "TTS",
      str_detect(Annotation, "^Intergenic") ~ "Intergenic",
      str_detect(Annotation, "^non-coding") ~ "Non-coding",
      TRUE ~ "Other"
    )
  )

location_summary <- homer_df %>%
  dplyr::count(Signature, Location) %>%
  group_by(Signature) %>%
  mutate(
    Percent = 100 * n / sum(n)
  ) %>%
  ungroup()

print(location_summary, n = Inf)


#Genomic distribution of TE signatures
location_summary_small <- location_summary %>%
  filter(
    Signature %in% c(
      "RepeatMasker",
      "Tumor-Associated",
      "Metastasis-Associated"
    ),
    !is.na(Location),
    Location != "Other"
  )

# Add total N per signature
totals <- location_summary_small %>%
  group_by(Signature) %>%
  summarise(
    N = sum(n),
    .groups = "drop"
  )

location_summary_small <- location_summary_small %>%
  select(-any_of(c("N", "Signature_lab"))) %>%
  left_join(
    totals,
    by = "Signature"
  ) %>%
  mutate(
    Signature_lab = paste0(
      Signature,
      "\n(n=",
      format(N, big.mark = ","),
      ")"
    )
  )

# Signature order
sig_levels <- c(
  "RepeatMasker",
  "Tumor-Associated",
  "Metastasis-Associated"
)

location_summary_small$Signature <- factor(
  location_summary_small$Signature,
  levels = sig_levels
)

label_lookup <- location_summary_small %>%
  distinct(
    Signature,
    Signature_lab
  ) %>%
  arrange(
    match(
      Signature,
      sig_levels
    )
  )

location_summary_small$Signature_lab <- factor(
  location_summary_small$Signature_lab,
  levels = label_lookup$Signature_lab
)


# Location order
loc_levels <- c(
  "Promoter",
  "5' UTR",
  "Exon",
  "Intron",
  "3' UTR",
  "TTS",
  "Intergenic",
  "Non-coding"
)

location_summary_small$Location <- factor(
  location_summary_small$Location,
  levels = loc_levels
)

location_summary_small <- location_summary_small %>%
  filter(!is.na(Location))


# Colors
loc_colors <- c(
  "Promoter"   = "#E41A1C",
  "5' UTR"     = "#FF7F00",
  "Exon"       = "#4DAF4A",
  "Intron"     = "#377EB8",
  "3' UTR"     = "#984EA3",
  "TTS"        = "#F781BF",
  "Intergenic" = "#999999",
  "Non-coding" = "#A65628"
)


# Barplot
p1 <- ggplot(
  location_summary_small,
  aes(
    x = Signature_lab,
    y = Percent,
    fill = Location
  )
) +
  geom_col(
    width = 0.75,
    color = "white",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = loc_colors,
    drop = FALSE
  ) +
  labs(
    title = "Genomic Distribution of TE Signatures",
    subtitle = "RepeatMasker background compared with Tumor- and Metastasis-Associated signatures",
    x = NULL,
    y = "Percentage of loci",
    fill = "Genomic Location"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      colour = "grey30"
    ),
    axis.text.x = element_text(
      face = "bold"
    ),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )

p1


##18) Statistical test for the genomic locations of each signature vs RepeatMasker
# Fisher's exact test

# Add Location category to each HOMER table
add_location <- function(df){
  
  df %>%
    mutate(
      Location = case_when(
        str_detect(Annotation, "^promoter-TSS") ~ "Promoter",
        str_detect(Annotation, "^5' UTR") ~ "5' UTR",
        str_detect(Annotation, "^exon") ~ "Exon",
        str_detect(Annotation, "^intron") ~ "Intron",
        str_detect(Annotation, "^3' UTR") ~ "3' UTR",
        str_detect(Annotation, "^TTS") ~ "TTS",
        str_detect(Annotation, "^Intergenic") ~ "Intergenic",
        str_detect(Annotation, "^non-coding") ~ "Non-coding",
        TRUE ~ NA_character_
      )
    )
}

tumor <- add_location(tumor)
repeatmasker <- add_location(repeatmasker)

# Genomic compartments
locations <- c(
  "Promoter",
  "5' UTR",
  "Exon",
  "Intron",
  "3' UTR",
  "TTS",
  "Intergenic",
  "Non-coding"
)

# Fisher enrichment: Tumor-associated signature vs RepeatMasker
tumor_location_enrichment <- lapply(
  locations,
  function(loc){
    
    a <- sum(tumor$Location == loc, na.rm = TRUE)
    
    b <- sum(
      !is.na(tumor$Location)
      & tumor$Location != loc
    )
    
    c <- sum(
      repeatmasker$Location == loc,
      na.rm = TRUE
    )
    
    d <- sum(
      !is.na(repeatmasker$Location)
      & repeatmasker$Location != loc
    )
    
    ft <- fisher.test(
      matrix(
        c(a, b,
          c, d),
        nrow = 2,
        byrow = TRUE
      )
    )
    
    data.frame(
      Location = loc,
      
      Tumor_Count = a,
      
      Tumor_Percent =
        100 * a / sum(!is.na(tumor$Location)),
      
      RepeatMasker_Count = c,
      
      RepeatMasker_Percent =
        100 * c / sum(!is.na(repeatmasker$Location)),
      
      Odds_Ratio =
        unname(ft$estimate),
      
      Log2_OR =
        log2(unname(ft$estimate)),
      
      Pvalue =
        ft$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(
      Pvalue,
      method = "BH"
    ),
    Direction = case_when(
      Odds_Ratio > 1 ~ "Enriched",
      Odds_Ratio < 1 ~ "Depleted",
      TRUE ~ "No change"
    )
  ) %>%
  arrange(FDR)

tumor_location_enrichment


# Fisher's exact test
# Metastasis-Associated vs RepeatMasker

# Add Location category
add_location <- function(df){
  
  df %>%
    mutate(
      Location = case_when(
        str_detect(Annotation, "^promoter-TSS") ~ "Promoter",
        str_detect(Annotation, "^5' UTR") ~ "5' UTR",
        str_detect(Annotation, "^exon") ~ "Exon",
        str_detect(Annotation, "^intron") ~ "Intron",
        str_detect(Annotation, "^3' UTR") ~ "3' UTR",
        str_detect(Annotation, "^TTS") ~ "TTS",
        str_detect(Annotation, "^Intergenic") ~ "Intergenic",
        str_detect(Annotation, "^non-coding") ~ "Non-coding",
        TRUE ~ NA_character_
      )
    )
}

met_assoc <- add_location(met_assoc)
repeatmasker <- add_location(repeatmasker)


# Genomic compartments
locations <- c(
  "Promoter",
  "5' UTR",
  "Exon",
  "Intron",
  "3' UTR",
  "TTS",
  "Intergenic",
  "Non-coding"
)

# Fisher enrichment
met_location_enrichment <- lapply(
  locations,
  function(loc) {
    
    a <- sum(
      met_assoc$Location == loc,
      na.rm = TRUE
    )
    
    b <- sum(
      !is.na(met_assoc$Location) &
        met_assoc$Location != loc
    )
    
    c <- sum(
      repeatmasker$Location == loc,
      na.rm = TRUE
    )
    
    d <- sum(
      !is.na(repeatmasker$Location) &
        repeatmasker$Location != loc
    )
    
    ft <- fisher.test(
      matrix(
        c(a, b,
          c, d),
        nrow = 2,
        byrow = TRUE
      )
    )
    
    data.frame(
      Location = loc,
      
      Metastasis_Count = a,
      
      Metastasis_Percent =
        round(
          100 * a /
            sum(!is.na(met_assoc$Location)),
          2
        ),
      
      RepeatMasker_Count = c,
      
      RepeatMasker_Percent =
        round(
          100 * c /
            sum(!is.na(repeatmasker$Location)),
          2
        ),
      
      Odds_Ratio =
        round(
          unname(ft$estimate),
          3
        ),
      
      Log2_OR =
        round(
          log2(unname(ft$estimate)),
          3
        ),
      
      Pvalue =
        ft$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(
      Pvalue,
      method = "BH"
    ),
    Direction = case_when(
      Odds_Ratio > 1 ~ "Enriched",
      Odds_Ratio < 1 ~ "Depleted",
      TRUE ~ "No change"
    )
  ) %>%
  arrange(FDR)

met_location_enrichment

# Significant enrichments/depletions only

met_location_enrichment %>%
  filter(FDR < 0.05) %>%
  select(
    Location,
    Odds_Ratio,
    Log2_OR,
    FDR,
    Direction
  )




#19) Evolutionary age of TE loci. RepeatMasker milliDiv-based age estimates

# Load RepeatMasker table
rmsk <- read.delim(
  "rmsk.txt",
  header = FALSE
)

# Create age map
# Age (Myr) = milliDiv / 2.2
# milliDiv is column V3 in standard UCSC rmsk.txt (17-column schema)
#V1 = bin
#V2 = swScore
#V3 = milliDiv --> This is what we need to use
#V4 = milliDel
#V5 = milliIns

rmsk_map <- rmsk %>%
  transmute(
    subfamily = V11,
    age_myr = V3 / 2.2   # milliDiv, not swScore
  ) %>%
  group_by(subfamily) %>%
  mutate(
    dup_id = row_number(),
    TE = paste0(
      subfamily,
      "_dup",
      dup_id
    )
  ) %>%
  ungroup() %>%
  select(
    TE,
    subfamily,
    age_myr
  )


# Attach age to annotated TE loci
te_age <- te_annotation %>%
  select(
    TE,
    TE_family,
    TE_class
  ) %>%
  left_join(
    rmsk_map,
    by = "TE"
  )


# Check mapping success
cat(
  "Mapped loci:",
  sum(!is.na(te_age$age_myr)),
  "/",
  nrow(te_age),
  "\n"
)


# Tumor-Associated
tumor_age <- te_age %>%
  filter(
    TE %in% TE_tumor_associated
  ) %>%
  mutate(
    Signature = "Tumor-Associated"
  )

# Metastasis-Associated
meta_assoc_age <- te_age %>%
  filter(
    TE %in% TE_metastasis_associated
  ) %>%
  mutate(
    Signature = "Metastasis-Associated"
  )


# RepeatMasker background
rm_age <- te_age %>%
  filter(
    !is.na(age_myr)
  ) %>%
  mutate(
    Signature = "RepeatMasker"
  )

# Combine all age data

age_df <- bind_rows(
  rm_age,
  tumor_age,
  meta_assoc_age
)


# Summary statistics
age_summary <- age_df %>%
  group_by(Signature) %>%
  summarise(
    n = n(),
    Median_Age = median(age_myr, na.rm = TRUE),
    Mean_Age = mean(age_myr, na.rm = TRUE),
    Min_Age = min(age_myr, na.rm = TRUE),
    Max_Age = max(age_myr, na.rm = TRUE),
    .groups = "drop"
  )

print(age_summary)



# Violin plot with median agesa nd significance stars

# Data
plot_df <- age_df %>%
  filter(
    Signature %in% c(
      "RepeatMasker",
      "Tumor-Associated",
      "Metastasis-Associated"
    ),
    TE_class %in% c(
      "SINE",
      "LINE",
      "LTR",
      "DNA"
    ),
    !is.na(age_myr)
  )

plot_df$Signature <- factor(
  plot_df$Signature,
  levels = c(
    "RepeatMasker",
    "Tumor-Associated",
    "Metastasis-Associated"
  )
)

plot_df$TE_class <- factor(
  plot_df$TE_class,
  levels = c(
    "SINE",
    "LINE",
    "LTR",
    "DNA"
  )
)


# Wilcoxon tests
tumor_stats <- plot_df %>%
  filter(
    Signature %in% c(
      "RepeatMasker",
      "Tumor-Associated"
    )
  ) %>%
  group_by(TE_class) %>%
  summarise(
    p = wilcox.test(age_myr ~ Signature)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    Comparison = "Tumor-Associated"
  )

meta_stats <- plot_df %>%
  filter(
    Signature %in% c(
      "RepeatMasker",
      "Metastasis-Associated"
    )
  ) %>%
  group_by(TE_class) %>%
  summarise(
    p = wilcox.test(age_myr ~ Signature)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    Comparison = "Metastasis-Associated"
  )

stats_df <- bind_rows(
  tumor_stats,
  meta_stats
)

stats_df$FDR <- p.adjust(
  stats_df$p,
  method = "BH"
)

stats_df$label <- case_when(
  stats_df$FDR < 0.001 ~ "***",
  stats_df$FDR < 0.01  ~ "**",
  stats_df$FDR < 0.05  ~ "*",
  TRUE                 ~ "ns"
)

# Label positions
label_pos <- plot_df %>%
  group_by(TE_class) %>%
  summarise(
    ymax = max(age_myr, na.rm = TRUE),
    .groups = "drop"
  )

#Significance stars
stats_df <- left_join(
  stats_df,
  label_pos,
  by = "TE_class"
)

stats_df$x <- ifelse(
  stats_df$Comparison == "Tumor-Associated",
  2,
  3
)


# Median age labels
median_labels <- plot_df %>%
  group_by(
    TE_class,
    Signature
  ) %>%
  summarise(
    median_age = median(
      age_myr,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  left_join(
    label_pos,
    by = "TE_class"
  )

median_labels$display <- sprintf(
  "%.1f",
  median_labels$median_age
)

#Plot
plot_age <- ggplot(
  plot_df,
  aes(
    x = Signature,
    y = age_myr,
    fill = Signature
  )
) +
  
  geom_violin(
    trim = TRUE,
    scale = "area",
    adjust = 1.5,
    colour = "black",
    linewidth = 0.4,
    alpha = 0.85
  ) +
  
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ) +
  
  # Median ages (all three groups)
  geom_text(
    data = median_labels,
    aes(
      x = Signature,
      y = ymax * 1.02,
      label = display
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4.5
  ) +
  
  # Wilcoxon significance stars
  geom_text(
    data = stats_df,
    aes(
      x = x,
      y = ymax * 1.10,
      label = label
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 6
  ) +
  
  facet_wrap(
    ~ TE_class,
    nrow = 1,
    scales = "free_y"
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.02, 0.22)
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "RepeatMasker" = "#BDBDBD",
      "Tumor-Associated" = "#F8766D",
      "Metastasis-Associated" = "#00BFC4"
    )
  ) +
  
  labs(
    title = "Evolutionary Age of TE Loci",
    subtitle = "Median ages and Wilcoxon significance versus RepeatMasker",
    x = NULL,
    y = "Evolutionary Age (Myr)",
    fill = NULL
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 20
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 14
    ),
    strip.background = element_rect(
      fill = "grey90",
      colour = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 16
    ),
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 12
    ),
    axis.title.y = element_text(
      face = "bold"
    )
  )

print(plot_age)




#20) Unsupervised Heatmap for the 2 signatures in z-score (including te class annotation and age)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

# Sample order
sample_order <- rownames(
  metadata_qc[
    order(metadata_qc$cond),
  ]
)

# Tissue annotation
tissue_labels <- metadata_qc[
  sample_order,
  "cond"
]

tissue_labels <- recode(
  as.character(tissue_labels),
  "normal_colon"     = "Normal Colon",
  "primary_colon"    = "Primary Colon",
  "liver_metastasis" = "Liver Metastasis"
)

tissue_labels <- factor(
  tissue_labels,
  levels = c(
    "Normal Colon",
    "Primary Colon",
    "Liver Metastasis"
  )
)

ha <- HeatmapAnnotation(
  
  Tissue = tissue_labels,
  
  col = list(
    Tissue = c(
      "Normal Colon"     = "#55A868",
      "Primary Colon"    = "#4C72B0",
      "Liver Metastasis" = "#C44E52"
    )
  ),
  
  annotation_name_gp = gpar(
    fontsize = 12,
    fontface = "bold"
  ),
  
  annotation_legend_param = list(
    Tissue = list(
      title = "Tissue",
      at = c(
        "Normal Colon",
        "Primary Colon",
        "Liver Metastasis"
      )
    )
  )
  
)

#Color scale
heat_colors <- colorRamp2(
  c(-2, 0, 2),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)

#Tumor-associated TE-loci
tumor_mat <- tpm[
  rownames(tpm) %in% TE_tumor_associated,
]

tumor_mat <- log2(
  tumor_mat + 1
)

tumor_mat_z <- t(
  scale(
    t(tumor_mat)
  )
)

tumor_mat_z[
  is.na(tumor_mat_z)
] <- 0

tumor_mat_z <- tumor_mat_z[
  ,
  sample_order
]

# TE CLASS + AGE INFORMATION
row_info <- te_annotation %>%
  
  filter(
    TE %in% rownames(tumor_mat_z)
  ) %>%
  
  distinct(
    TE,
    .keep_all = TRUE
  ) %>%
  
  select(
    TE,
    TE_class
  ) %>%
  
  left_join(
    
    te_age %>%
      
      select(
        TE,
        age_myr
      ),
    
    by = "TE"
    
  )

row_info <- row_info[
  match(
    rownames(tumor_mat_z),
    row_info$TE
  ),
]

# CLASS COLOURS
classes <- sort(
  unique(row_info$TE_class)
)

class_cols <- c(
  
  "DNA"        = "#E41A1C",
  "LINE"       = "#377EB8",
  "LTR"        = "#4DAF4A",
  "SINE"       = "#FF7F00",
  "Satellite"  = "#A65628",
  "LTR?"       = "#F781BF",
  "Retroposon" = "#984EA3",
  "Unknown"    = "#999999"
  
)


# Row annotations
ha_row <- rowAnnotation(
  
  Class = row_info$TE_class,
  
  `Age (Myr)` = row_info$age_myr,
  
  col = list(
    
    Class = class_cols,
    
    `Age (Myr)` = colorRamp2(
      
      c(
        min(row_info$age_myr, na.rm = TRUE),
        median(row_info$age_myr, na.rm = TRUE),
        max(row_info$age_myr, na.rm = TRUE)
      ),
      
      c(
        "#440154",
        "#21908C",
        "#F97306"
      )
      
    )
    
  ),
  
  annotation_width = unit(
    c(6, 6),
    "mm"
  )
  
)

#Heatmap
ht_tumor <- Heatmap(
  
  tumor_mat_z,
  
  name = "Z-score on TPM",
  
  top_annotation = ha,
  
  left_annotation = ha_row,
  
  cluster_rows = TRUE,
  
  cluster_columns = TRUE,
  
  show_row_names = FALSE,
  
  show_column_names = FALSE,
  
  col = heat_colors,
  
  column_title =
    "Tumor-Associated Individual TE Loci",
  
  column_title_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),
  
  row_title = paste(
    nrow(tumor_mat_z),
    "TE loci"
  ),
  
  row_title_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  )
  
)


#Print
draw(
  ht_tumor,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)



##Metastasis-associated heatmap

meta_mat <- tpm[
  rownames(tpm) %in% TE_metastasis_associated,
]

meta_mat <- log2(
  meta_mat + 1
)

meta_mat_z <- t(
  scale(
    t(meta_mat)
  )
)

meta_mat_z[
  is.na(meta_mat_z)
] <- 0

meta_mat_z <- meta_mat_z[
  ,
  sample_order
]

#TE class and age information
row_info_meta <- te_annotation %>%
  
  filter(
    TE %in% rownames(meta_mat_z)
  ) %>%
  
  distinct(
    TE,
    .keep_all = TRUE
  ) %>%
  
  select(
    TE,
    TE_class
  ) %>%
  
  left_join(
    
    te_age %>%
      
      select(
        TE,
        age_myr
      ),
    
    by = "TE"
    
  )

row_info_meta <- row_info_meta[
  match(
    rownames(meta_mat_z),
    row_info_meta$TE
  ),
]

#Row annotations 
ha_row_meta <- rowAnnotation(
  
  Class = row_info_meta$TE_class,
  
  `Age (Myr)` = row_info_meta$age_myr,
  
  col = list(
    
    Class = class_cols,
    
    `Age (Myr)` = colorRamp2(
      
      c(
        min(row_info_meta$age_myr, na.rm = TRUE),
        median(row_info_meta$age_myr, na.rm = TRUE),
        max(row_info_meta$age_myr, na.rm = TRUE)
      ),
      
      c(
        "#440154",
        "#21908C",
        "#F97306"
      )
      
    )
    
  ),
  
  annotation_width = unit(
    c(6, 6),
    "mm"
  )
  
)



#Heatmap
ht_meta <- Heatmap(
  
  meta_mat_z,
  
  name = "Z-score on TPM",
  
  top_annotation = ha,
  
  left_annotation = ha_row_meta,
  
  cluster_rows = TRUE,
  
  cluster_columns = TRUE,
  
  show_row_names = FALSE,
  
  show_column_names = FALSE,
  
  col = heat_colors,
  
  column_title =
    "Metastasis-Associated Individual TE Loci",
  
  column_title_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),
  
  row_title = paste(
    nrow(meta_mat_z),
    "TE loci"
  ),
  
  row_title_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  )
  
)

#Print
draw(
  ht_meta,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)






##21) Immunotherapy candidates.
#We will continue only with the tumor-associated signature


#There are ALT-TEs in the signature and they first need to be removed before we search for candidates

#ALT-TEs IDs
alt_tes <- te_annotation %>%
  filter(
    grepl("_alt$", chrom)
  ) %>%
  pull(TE)



# Remove ALT loci from  homer files
repeatmasker <- repeatmasker %>%
  filter(!(.[[1]] %in% alt_tes))

tumor <- tumor %>%
  filter(!(.[[1]] %in% alt_tes))

met_assoc <- met_assoc %>%
  filter(!(.[[1]] %in% alt_tes))

# Combine HOMER tables
homer_df <- bind_rows(
  repeatmasker,
  tumor,
  met_assoc)

# Assign genomic locations
homer_df <- homer_df %>%
  mutate(
    Location = case_when(
      str_detect(Annotation, "^promoter-TSS") ~ "Promoter",
      str_detect(Annotation, "^intron") ~ "Intron",
      str_detect(Annotation, "^exon") ~ "Exon",
      str_detect(Annotation, "^3' UTR") ~ "3' UTR",
      str_detect(Annotation, "^5' UTR") ~ "5' UTR",
      str_detect(Annotation, "^TTS") ~ "TTS",
      str_detect(Annotation, "^Intergenic") ~ "Intergenic",
      str_detect(Annotation, "^non-coding") ~ "Non-coding",
      TRUE ~ "Other"
    )
  )


# Summary table
location_summary <- homer_df %>%
  dplyr::count(Signature, Location) %>%
  group_by(Signature) %>%
  mutate(
    Percent = 100 * n / sum(n)
  ) %>%
  ungroup()

print(location_summary, n = Inf)

# Combine homer tables
homer_df <- bind_rows(
  repeatmasker,
  tumor,
  met_assoc)

#  Assign genomic location
homer_df <- homer_df %>%
  mutate(
    Location = case_when(
      str_detect(Annotation, "^promoter-TSS") ~ "Promoter",
      str_detect(Annotation, "^intron") ~ "Intron",
      str_detect(Annotation, "^exon") ~ "Exon",
      str_detect(Annotation, "^3' UTR") ~ "3' UTR",
      str_detect(Annotation, "^5' UTR") ~ "5' UTR",
      str_detect(Annotation, "^TTS") ~ "TTS",
      str_detect(Annotation, "^Intergenic") ~ "Intergenic",
      str_detect(Annotation, "^non-coding") ~ "Non-coding",
      TRUE ~ "Other"
    )
  )

#Location summary
location_summary <- homer_df %>%
  filter(
    !is.na(Location),
    Location != "Other"
  ) %>%
  dplyr::count(
    Signature,
    Location
  ) %>%
  group_by(Signature) %>%
  mutate(
    Percent = 100 * n / sum(n)
  ) %>%
  ungroup()


#GTEx file
gtex <- readRDS("gtex.rds")


# GTEx TPM matrix
gtex_mat <- gtex %>%
  select(
    TE,
    Sample,
    Count
  ) %>%
  pivot_wider(
    names_from = Sample,
    values_from = Count,
    values_fill = 0
  )

te_lengths <- te_lengths_df$Length[
  match(
    gtex_mat$TE,
    te_lengths_df$TE
  )
]

valid <- !is.na(te_lengths)

gtex_mat <- gtex_mat[valid, ]
te_lengths <- te_lengths[valid]

gtex_counts <- as.matrix(
  gtex_mat[, -1]
)

rownames(gtex_counts) <- gtex_mat$TE

rpk <- sweep(
  gtex_counts,
  1,
  te_lengths / 1000,
  "/"
)

gtex_tpm <- sweep(
  rpk,
  2,
  colSums(rpk) / 1e6,
  "/"
)

# GTEx liver
gtex_liver_samples <- gtex %>%
  filter(Tissue == "Liver") %>%
  pull(Sample) %>%
  unique()

gtex_liver_tpm <- gtex_tpm[
  ,
  intersect(
    gtex_liver_samples,
    colnames(gtex_tpm)
  ),
  drop = FALSE
]

gtex_liver_mean <- data.frame(
  TE = rownames(gtex_liver_tpm),
  Mean_Liver_TPM = rowMeans(
    gtex_liver_tpm,
    na.rm = TRUE
  )
)

#GTEx colon
gtex_colon_samples <- gtex %>%
  filter(
    Tissue %in% c(
      "Colon",
      "Colon_Sigmoid",
      "Colon_Transverse"
    )
  ) %>%
  pull(Sample) %>%
  unique()

gtex_colon_tpm <- gtex_tpm[
  ,
  intersect(
    gtex_colon_samples,
    colnames(gtex_tpm)
  ),
  drop = FALSE
]

gtex_colon_mean <- data.frame(
  TE = rownames(gtex_colon_tpm),
  Mean_Colon_TPM = rowMeans(
    gtex_colon_tpm,
    na.rm = TRUE
  )
)



#Sample groups 
normal_samples <- rownames(
  metadata_qc[metadata_qc$cond == "normal_colon", ]
)

primary_samples <- rownames(
  metadata_qc[metadata_qc$cond == "primary_colon", ]
)

met_samples <- rownames(
  metadata_qc[metadata_qc$cond == "liver_metastasis", ]
)


##Remove the ALT from the signatures
TE_tumor_associated <- setdiff(TE_tumor_associated, alt_tes)
TE_metastasis_associated <- setdiff(TE_metastasis_associated, alt_tes)


#Filter 1: >=70% (TPM >1) in the metastatic samples of the discovery cohort (GSE50760)
tumor_candidates <- data.frame(
  TE = setdiff(TE_tumor_associated, alt_tes),
  stringsAsFactors = FALSE
)

# Metastasis prevalence (TPM > 1)
tumor_candidates$Pct_Metastasis_TPM1 <- apply(
  tpm[
    tumor_candidates$TE,
    met_samples,
    drop = FALSE
  ],
  1,
  function(x) {
    mean(x > 1, na.rm = TRUE) * 100
  }
)

# Keep only TEs present in >70% metastases
tumor_f1 <- tumor_candidates %>%
  filter(
    Pct_Metastasis_TPM1 > 70
  )

# Number of normal colon samples expressing the TE (TPM > 1)

tumor_f1$Normal_N_TPM1 <- apply(
  tpm[
    tumor_f1$TE,
    normal_samples,
    drop = FALSE
  ],
  1,
  function(x) {
    sum(x > 1, na.rm = TRUE)
  }
)

#Candidate sets
tumor_norm0 <- tumor_f1 %>%
  filter(Normal_N_TPM1 == 0)

tumor_norm1 <- tumor_f1 %>%
  filter(Normal_N_TPM1 <= 1)

tumor_norm2 <- tumor_f1 %>%
  filter(Normal_N_TPM1 <= 2)

tumor_norm3 <- tumor_f1 %>%
  filter(Normal_N_TPM1 <= 3)

#Summary
cat(
  "\nCandidates after >70% metastasis filter:",
  nrow(tumor_f1), "\n",
  
  "\nRemaining candidates:\n",
  "0/18 normals  :", nrow(tumor_norm0), "\n",
  "<=1/18 normals :", nrow(tumor_norm1), "\n",
  "<=2/18 normals :", nrow(tumor_norm2), "\n",
  "<=3/18 normals :", nrow(tumor_norm3), "\n"
)


#22) From this we selected the 1/18 filtering which gave us 131 candidates and we visually inspected
#them in the gtex to try and identify the best
#TPM values on log2scale and the output is 131 tiff files 
library(ggtext)
library(scales)

#Candidates
candidates <- tumor_norm1$TE

# Output directory
outdir <- "TumorNorm1_GTEx_Plots"

if (!dir.exists(outdir)) {
  dir.create(outdir)
}

# GTEx annotation
gtex_info <- gtex %>%
  distinct(
    Sample,
    Tissue
  )

# Loop through candidates
for(te_of_interest in candidates){
  
  cohort_df <- data.frame(
    TPM = c(
      as.numeric(tpm[te_of_interest, normal_samples]),
      as.numeric(tpm[te_of_interest, primary_samples]),
      as.numeric(tpm[te_of_interest, met_samples])
    ),
    
    Tissue = c(
      rep("Normal colon", length(normal_samples)),
      rep("Primary colon", length(primary_samples)),
      rep("Liver metastasis", length(met_samples))
    )
  )
  
  gtex_df <- data.frame(
    Sample = colnames(gtex_tpm),
    TPM = as.numeric(
      gtex_tpm[
        te_of_interest,
        colnames(gtex_tpm)
      ]
    )
  ) %>%
    left_join(
      gtex_info,
      by = "Sample"
    ) %>%
    select(
      TPM,
      Tissue
    )
  
  plot_df <- bind_rows(
    cohort_df,
    gtex_df
  )
  
  gtex_tissues <- c(
    "Colon",
    "Liver",
    sort(
      setdiff(
        unique(gtex_df$Tissue),
        c("Colon", "Liver")
      )
    )
  )
  
  tissue_order <- c(
    "Normal colon",
    "Primary colon",
    "Liver metastasis",
    gtex_tissues
  )
  
  plot_df$Tissue <- factor(
    plot_df$Tissue,
    levels = tissue_order
  )
  
  tissue_stats <- plot_df %>%
    group_by(Tissue) %>%
    summarise(
      N = n(),
      Pct_Expressed = round(
        100 * mean(TPM > 1),
        1
      ),
      .groups = "drop"
    )
  
  bold_tissues <- c(
    "Normal colon",
    "Primary colon",
    "Liver metastasis",
    "Colon",
    "Liver"
  )
  
  tissue_labels <- sapply(
    tissue_stats$Tissue,
    function(x){
      
      n_text <- paste0(
        "n=",
        tissue_stats$N[
          tissue_stats$Tissue == x
        ]
      )
      
      pct_text <- paste0(
        tissue_stats$Pct_Expressed[
          tissue_stats$Tissue == x
        ],
        "%"
      )
      
      if(x %in% bold_tissues){
        
        paste0(
          "<b>", x,
          "</b><br><b>(",
          n_text,
          " | ",
          pct_text,
          ")</b>"
        )
        
      } else {
        
        paste0(
          x,
          "<br>(",
          n_text,
          " | ",
          pct_text,
          ")"
        )
        
      }
    }
  )
  
  names(tissue_labels) <- tissue_stats$Tissue
  
  plot_df$TPM_plot <- plot_df$TPM
  plot_df$TPM_plot[plot_df$TPM_plot == 0] <- 0.1
  
  tissue_cols <- c(
    "Normal colon"     = "#55A868",
    "Primary colon"    = "#4C72B0",
    "Liver metastasis" = "#C44E52",
    "Colon"            = "#8172B2",
    "Liver"            = "#F0AD4E"
  )
  
  other_tissues <- setdiff(
    levels(plot_df$Tissue),
    names(tissue_cols)
  )
  
  tissue_cols <- c(
    tissue_cols,
    setNames(
      scales::hue_pal()(length(other_tissues)),
      other_tissues
    )
  )
  
  p <- ggplot(
    plot_df,
    aes(
      Tissue,
      TPM_plot
    )
  ) +
    
    geom_jitter(
      aes(colour = Tissue),
      width = 0.12,
      alpha = 0.8,
      size = 2.5
    ) +
    
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.28,
      colour = "black",
      linewidth = 1
    ) +
    
    scale_y_continuous(
      trans = "log2",
      breaks = c(
        0.1, 0.5, 1, 2, 5,
        10, 20, 50, 100
      ),
      labels = c(
        "0", "0.5", "1", "2",
        "5", "10", "20", "50", "100"
      )
    ) +
    
    scale_x_discrete(
      labels = tissue_labels
    ) +
    
    scale_colour_manual(
      values = tissue_cols
    ) +
    
    annotate(
      "segment",
      x = 1,
      xend = 3,
      y = 0.06,
      yend = 0.06,
      colour = "steelblue4",
      linewidth = 1.5
    ) +
    
    annotate(
      "text",
      x = 2,
      y = 0.035,
      label = "Discovery cohort\n(GSE50760)",
      colour = "steelblue4",
      fontface = "bold",
      size = 4
    ) +
    
    annotate(
      "segment",
      x = 4,
      xend = length(tissue_order),
      y = 0.06,
      yend = 0.06,
      colour = "darkgreen",
      linewidth = 1.5
    ) +
    
    annotate(
      "text",
      x = (4 + length(tissue_order)) / 2,
      y = 0.035,
      label = "GTEx",
      colour = "darkgreen",
      fontface = "bold",
      size = 4.5
    ) +
    
    labs(
      title = te_of_interest,
      subtitle = "Tumor-associated candidate (<=1/18 normal samples, TPM > 1)",
      x = NULL,
      y = "TE expression (TPM, log2 scale)"
    ) +
    
    theme_bw(base_size = 18) +
    
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        size = 22,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        face = "italic",
        size = 14,
        hjust = 0.5
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 18
      ),
      axis.text.x = ggtext::element_markdown(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 14
      )
    )
  
  outfile <- file.path(
    outdir,
    paste0(
      gsub("[^A-Za-z0-9_]", "_", te_of_interest),
      "_GTEx_expression.tiff"
    )
  )
  
  ggsave(
    filename = outfile,
    plot = p,
    device = "tiff",
    width = 14,
    height = 8,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
  
  cat(
    "Saved:",
    basename(outfile),
    "\n"
  )
}

cat(
  "\nAll plots saved in folder:",
  normalizePath(outdir),
  "\n"
)



#23) Unsupervised heatmap of tumor-associated TE candidates with GTEx. Z-score, including class and age

# Discovery cohort (GSE50760) samples
discovery_samples <- c(
  normal_samples,
  primary_samples,
  met_samples
)


# GTEx info
gtex_info <- gtex %>%
  distinct(
    Sample,
    Tissue
  ) %>%
  filter(
    Sample %in% colnames(gtex_tpm)
  )


# Order GTEx:
# Liver -> Colon -> Others
gtex_liver_samples <- gtex_info %>%
  filter(
    Tissue == "Liver"
  ) %>%
  pull(Sample)

gtex_colon_samples <- gtex_info %>%
  filter(
    Tissue == "Colon"
  ) %>%
  pull(Sample)

gtex_other_samples <- gtex_info %>%
  filter(
    !Tissue %in% c(
      "Liver",
      "Colon"
    )
  ) %>%
  pull(Sample)

gtex_order <- c(
  gtex_liver_samples,
  gtex_colon_samples,
  gtex_other_samples
)


# Common TEs
common_tes <- intersect(
  TE_tumor_associated,
  intersect(
    rownames(tpm),
    rownames(gtex_tpm)
  )
)


# Expression matrix
tumor_mat <- cbind(
  
  tpm[
    common_tes,
    discovery_samples,
    drop = FALSE
  ],
  
  gtex_tpm[
    common_tes,
    gtex_order,
    drop = FALSE
  ]
  
)

# Replace NAs
tumor_mat[is.na(tumor_mat)] <- 0


# LOG2(TPM+1)
tumor_mat <- log2(
  tumor_mat + 1
)

#Row Z-score
tumor_mat <- t(
  scale(
    t(tumor_mat)
  )
)

tumor_mat[
  is.na(tumor_mat)
] <- 0


#TE class and age information
row_info <- te_annotation %>%
  
  filter(
    TE %in% rownames(tumor_mat)
  ) %>%
  
  distinct(
    TE,
    .keep_all = TRUE
  ) %>%
  
  select(
    TE,
    TE_class
  ) %>%
  
  left_join(
    
    te_age %>%
      
      select(
        TE,
        age_myr
      ),
    
    by = "TE"
    
  )

row_info <- row_info[
  match(
    rownames(tumor_mat),
    row_info$TE
  ),
]

#Class colors
class_cols <- c(
  
  "DNA"        = "#E41A1C",
  "LINE"       = "#377EB8",
  "LTR"        = "#4DAF4A",
  "SINE"       = "#FF7F00",
  "Satellite"  = "#A65628",
  "LTR?"       = "#F781BF",
  "Retroposon" = "#984EA3",
  "Unknown"    = "#999999"
  
)


# Row annotations 
ha_row <- rowAnnotation(
  
  Class = row_info$TE_class,
  
  `Age (Myr)` = row_info$age_myr,
  
  col = list(
    
    Class = class_cols,
    
    `Age (Myr)` = colorRamp2(
      
      c(
        min(row_info$age_myr, na.rm = TRUE),
        median(row_info$age_myr, na.rm = TRUE),
        max(row_info$age_myr, na.rm = TRUE)
      ),
      
      c(
        "#440154",
        "#21908C",
        "#F97306"
      )
      
    )
    
  ),
  
  annotation_width = unit(
    c(6, 6),
    "mm"
  )
  
)

# Dataset annotation
dataset_annotation <- c(
  
  rep(
    "Discovery cohort (GSE50760)",
    length(discovery_samples)
  ),
  
  rep(
    "GTEx",
    length(gtex_order)
  )
  
)

dataset_annotation <- factor(
  dataset_annotation,
  levels = c(
    "Discovery cohort (GSE50760)",
    "GTEx"
  )
)


# Tissue annotation
tissue_annotation <- c(
  
  rep(
    "Normal Colon (Discovery Cohort)",
    length(normal_samples)
  ),
  
  rep(
    "Primary Colon (Discovery Cohort)",
    length(primary_samples)
  ),
  
  rep(
    "CRLM (Discovery Cohort)",
    length(met_samples)
  ),
  
  gtex_info$Tissue[
    match(
      gtex_order,
      gtex_info$Sample
    )
  ]
  
)

tissue_annotation <- factor(
  tissue_annotation,
  levels = c(
    
    "Normal Colon (Discovery Cohort)",
    "Primary Colon (Discovery Cohort)",
    "CRLM (Discovery Cohort)",
    
    "Liver",
    "Colon",
    
    "Adipose_Tissue",
    "Adrenal_Gland",
    "Blood_Vessel",
    "Brain",
    "Esophagus",
    "Heart",
    "Lung",
    "Muscle",
    "Nerve",
    "Ovary",
    "Pancreas",
    "Pituitary",
    "Prostate",
    "Salivary_Gland",
    "Skin",
    "Small_Intestine",
    "Spleen",
    "Stomach",
    "Testis",
    "Thyroid",
    "Uterus",
    "Vagina"
  )
)

# Tissue colours
tissue_cols <- c(
  
  "Normal Colon (Discovery Cohort)" = "#55A868",
  "Primary Colon (Discovery Cohort)" = "#4C72B0",
  "CRLM (Discovery Cohort)" = "#C44E52",
  
  "Colon" = "#8172B2",
  "Liver" = "#F0AD4E",
  
  "Adipose_Tissue" = "#8DD3C7",
  "Adrenal_Gland" = "#FFFFB3",
  "Blood_Vessel" = "#BEBADA",
  "Brain" = "red",
  "Esophagus" = "darkorange",
  "Heart" = "#B3DE69",
  "Lung" = "#D9D9D9",
  "Muscle" = "#BC80BD",
  "Nerve" = "blue",
  "Ovary" = "#FFED6F",
  "Pancreas" = "pink",
  "Pituitary" = "#B3DE69",
  "Prostate" = "#FB9A99",
  "Salivary_Gland" = "#A50F15",
  "Skin" = "#A6CEE3",
  "Small_Intestine" = "#B2DF8A",
  "Spleen" = "yellow",
  "Stomach" = "green",
  "Testis" = "purple",
  "Thyroid" = "#6A3D9A",
  "Uterus" = "#B15928",
  "Vagina" = "#F781BF"
)


# Check missing colours
setdiff(
  unique(as.character(tissue_annotation)),
  names(tissue_cols)
)


# Top annotation
ha <- HeatmapAnnotation(
  
  Dataset = dataset_annotation,
  
  Tissue = tissue_annotation,
  
  col = list(
    
    Dataset = c(
      "Discovery cohort (GSE50760)" = "#6BAED6",
      "GTEx" = "#31A354"
    ),
    
    Tissue = tissue_cols
  )
)


# Z-score colour scale
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)


# Heatmap
tumor_allgtex <- Heatmap(
  
  tumor_mat,
  
  name = "Z-score",
  
  col = col_fun,
  
  top_annotation = ha,
  
  left_annotation = ha_row,
  
  cluster_rows = TRUE,
  
  cluster_columns = FALSE,
  
  show_row_names = FALSE,
  
  show_column_names = FALSE,
  
  use_raster = TRUE,
  
  raster_quality = 2,
  
  row_title = paste0(
    "Tumor-associated TEs (n=",
    nrow(tumor_mat),
    ")"
  ),
  
  column_title =
    "Tumor-associated TE Signature Across Discovery Cohort and GTEx",
  
  heatmap_legend_param = list(
    title = "Z-score on TPM",
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2")
  )
)


# Preview
draw(
  tumor_allgtex,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# ----------------------------------------------------------
# TIFF export
# ----------------------------------------------------------

tiff(
  "Tumor_Associated_TE_Signature_All_GTEx_Tissues_Zscore.tiff",
  width = 10000,
  height = 6000,
  res = 900,
  compression = "lzw"
)

draw(
  tumor_allgtex,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()



#And for the metastasis-associated. Discovery cohort and all gtex in zscore
# Common TEs
common_tes_meta <- intersect(
  TE_metastasis_associated,
  intersect(
    rownames(tpm),
    rownames(gtex_tpm)
  )
)

# Expression matrix
meta_mat <- cbind(
  
  tpm[
    common_tes_meta,
    discovery_samples,
    drop = FALSE
  ],
  
  gtex_tpm[
    common_tes_meta,
    gtex_order,
    drop = FALSE
  ]
  
)

# Replace NAs
meta_mat[is.na(meta_mat)] <- 0

# LOG2(TPM+1)
meta_mat <- log2(
  meta_mat + 1
)

# Row Z-score
meta_mat <- t(
  scale(
    t(meta_mat)
  )
)

meta_mat[
  is.na(meta_mat)
] <- 0

#TE class and age information
row_info_meta <- te_annotation %>%
  
  filter(
    TE %in% rownames(meta_mat)
  ) %>%
  
  distinct(
    TE,
    .keep_all = TRUE
  ) %>%
  
  select(
    TE,
    TE_class
  ) %>%
  
  left_join(
    
    te_age %>%
      
      select(
        TE,
        age_myr
      ),
    
    by = "TE"
    
  )

row_info_meta <- row_info_meta[
  match(
    rownames(meta_mat),
    row_info_meta$TE
  ),
]

#Row annotation
ha_row_meta <- rowAnnotation(
  
  Class = row_info_meta$TE_class,
  
  `Age (Myr)` = row_info_meta$age_myr,
  
  col = list(
    
    Class = class_cols,
    
    `Age (Myr)` = colorRamp2(
      
      c(
        min(row_info_meta$age_myr, na.rm = TRUE),
        median(row_info_meta$age_myr, na.rm = TRUE),
        max(row_info_meta$age_myr, na.rm = TRUE)
      ),
      
      c(
        "#440154",
        "#21908C",
        "#F97306"
      )
      
    )
    
  ),
  
  annotation_width = unit(
    c(6, 6),
    "mm"
  )
  
)
#Heatmap
meta_allgtex <- Heatmap(
  
  meta_mat,
  
  name = "Z-score on TPM",
  
  col = col_fun,
  
  top_annotation = ha,
  
  left_annotation = ha_row_meta,
  
  cluster_rows = TRUE,
  
  cluster_columns = FALSE,
  
  show_row_names = FALSE,
  
  show_column_names = FALSE,
  
  use_raster = TRUE,
  
  raster_quality = 2,
  
  row_title = paste0(
    "Metastasis-associated TEs (n=",
    nrow(meta_mat),
    ")"
  ),
  
  column_title =
    "Metastasis-associated TE Signature Across Discovery Cohort and GTEx",
  
  heatmap_legend_param = list(
    title = "Z-score on TPM",
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2")
  )
  
)

#Preview
draw(
  meta_allgtex,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

#TIFF export
tiff(
  "Metastasis_Associated_TE_Signature_All_GTEx_Tissues_Zscore.tiff",
  width = 10000,
  height = 6000,
  res = 900,
  compression = "lzw"
)

draw(
  meta_allgtex,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()




#24) The candidates that were more promising were: Tigger7_dup2343 andLTR66_dup140
#Make tables with their information
library(gt)

#Helper function
get_tpm <- function(te, samples){
  
  as.numeric(
    unlist(
      tpm[
        te,
        samples,
        drop = FALSE
      ]
    )
  )
  
}

#Candidate card
make_candidate_card <- function(te){
  
  anno <- candidate_annotation %>%
    filter(
      TE == te
    )
  
  
  # Discovery cohort prevalence
  normal_pct <- round(
    100 * mean(
      get_tpm(te, normal_samples) > 1
    ),
    1
  )
  
  primary_pct <- round(
    100 * mean(
      get_tpm(te, primary_samples) > 1
    ),
    1
  )
  
  met_pct <- round(
    100 * mean(
      get_tpm(te, met_samples) > 1
    ),
    1
  )
  
  # GTEx prevalence
  
  gtex_pct <- round(
    100 * mean(
      as.numeric(gtex_tpm[te, ]) > 1,
      na.rm = TRUE
    ),
    1
  )
  
  # Table
  
  card_df <- data.frame(
    
    Feature = c(
      
      "Class",
      "Family",
      "Chromosome",
      "Genomic Context",
      "Nearest Gene",
      "Gene Type",
      "Distance to TSS",
      "Evolutionary Age",
      
      "",
      
      "Normal Colon (>1 TPM)",
      "Primary CRC (>1 TPM)",
      "CRLM (>1 TPM)",
      "GTEx Samples (>1 TPM)"
      
    ),
    
    Value = c(
      
      anno$TE_class,
      anno$TE_family,
      anno$Chromosome,
      anno$Annotation,
      anno$Gene,
      anno$Gene_Type,
      
      paste0(
        format(
          anno$Distance_to_TSS,
          big.mark = ","
        ),
        " bp"
      ),
      
      paste0(
        round(
          anno$age_myr,
          1
        ),
        " Myr"
      ),
      
      "",
      
      paste0(normal_pct, "%"),
      paste0(primary_pct, "%"),
      paste0(met_pct, "%"),
      paste0(gtex_pct, "%")
      
    ),
    
    stringsAsFactors = FALSE
    
  )
  
  gt(card_df) %>%
    
    tab_header(
      title = md(
        paste0("**", te, "**")
      ),
      subtitle = md(
        "*Final Immunotherapy Candidate*"
      )
    ) %>%
    
    cols_label(
      Feature = "",
      Value = ""
    ) %>%
    
    tab_style(
      style = list(
        cell_fill(
          color = "#1F4E79"
        ),
        cell_text(
          color = "white",
          weight = "bold",
          size = px(24)
        )
      ),
      locations = cells_title()
    ) %>%
    
    tab_style(
      style = list(
        cell_fill("#EAF2F8"),
        cell_text(weight = "bold")
      ),
      locations = cells_body(
        rows = 1:8
      )
    ) %>%
    
    tab_style(
      style = list(
        cell_fill("#FDEDEC"),
        cell_text(weight = "bold")
      ),
      locations = cells_body(
        rows = 10:13
      )
    ) %>%
    
    opt_row_striping() %>%
    
    cols_align(
      align = "left",
      columns = everything()
    ) %>%
    
    tab_options(
      
      table.font.size = px(18),
      
      heading.title.font.size = px(30),
      
      heading.subtitle.font.size = px(16),
      
      data_row.padding = px(10),
      
      table.border.top.width = px(3),
      
      table.border.bottom.width = px(3)
      
    )
  
}

#Generate tables
tbl_tigger <- make_candidate_card(
  "Tigger7_dup2343"
)

tbl_ltr66 <- make_candidate_card(
  "LTR66_dup140"
)

# View
tbl_tigger
tbl_ltr66

#Save the tables
gtsave(
  tbl_tigger,
  "Tigger7_dup2343_Candidate_Card.png",
  zoom = 3
)

gtsave(
  tbl_ltr66,
  "LTR66_dup140_Candidate_Card.png",
  zoom = 3
)



##Jitter plots of the 2 final candidates

#Candidates
candidates <- c(
  "Tigger7_dup2343",
  "LTR66_dup140"
)

# GTEx annotation
gtex_info <- gtex %>%
  distinct(
    Sample,
    Tissue
  )

#Loop through the candidates
for(te_of_interest in candidates){
  
  # Discovery cohort
  
  
  cohort_df <- data.frame(
    
    TPM = c(
      as.numeric(unlist(
        tpm[
          te_of_interest,
          normal_samples,
          drop = FALSE
        ]
      )),
      as.numeric(unlist(
        tpm[
          te_of_interest,
          primary_samples,
          drop = FALSE
        ]
      )),
      as.numeric(unlist(
        tpm[
          te_of_interest,
          met_samples,
          drop = FALSE
        ]
      ))
    ),
    
    Tissue = c(
      rep(
        "Normal colon",
        length(normal_samples)
      ),
      rep(
        "Primary colon",
        length(primary_samples)
      ),
      rep(
        "Liver metastasis",
        length(met_samples)
      )
    )
    
  )
  
  
  # GTEx
  gtex_df <- data.frame(
    
    Sample = colnames(gtex_tpm),
    
    TPM = as.numeric(
      gtex_tpm[
        te_of_interest,
        colnames(gtex_tpm)
      ]
    )
    
  ) %>%
    
    left_join(
      gtex_info,
      by = "Sample"
    ) %>%
    
    select(
      TPM,
      Tissue
    )
  
  # Combine
  plot_df <- bind_rows(
    cohort_df,
    gtex_df
  )
  
  # Tissue order
  
  gtex_tissues <- c(
    "Colon",
    "Liver",
    sort(
      setdiff(
        unique(gtex_df$Tissue),
        c(
          "Colon",
          "Liver"
        )
      )
    )
  )
  
  tissue_order <- c(
    "Normal colon",
    "Primary colon",
    "Liver metastasis",
    gtex_tissues
  )
  
  plot_df$Tissue <- factor(
    plot_df$Tissue,
    levels = tissue_order
  )
  
  # Statistics shown below each tissue

  tissue_stats <- plot_df %>%
    group_by(Tissue) %>%
    summarise(
      N = n(),
      Pct_Expressed =
        round(
          100 * mean(TPM > 1),
          1
        ),
      .groups = "drop"
    )
  
  bold_tissues <- c(
    "Normal colon",
    "Primary colon",
    "Liver metastasis",
    "Colon",
    "Liver"
  )
  
  tissue_labels <- sapply(
    tissue_stats$Tissue,
    function(x){
      
      n_text <- paste0(
        "n=",
        tissue_stats$N[
          tissue_stats$Tissue == x
        ]
      )
      
      pct_text <- paste0(
        tissue_stats$Pct_Expressed[
          tissue_stats$Tissue == x
        ],
        "%"
      )
      
      if(x %in% bold_tissues){
        
        paste0(
          "<b>",
          x,
          "</b><br><b>(",
          n_text,
          " | ",
          pct_text,
          ")</b>"
        )
        
      } else {
        
        paste0(
          x,
          "<br>(",
          n_text,
          " | ",
          pct_text,
          ")"
        )
        
      }
      
    }
  )
  
  names(tissue_labels) <- tissue_stats$Tissue
  
  # Keep zeros visible on log scale
  
  plot_df$TPM_plot <- plot_df$TPM
  plot_df$TPM_plot[
    plot_df$TPM_plot == 0
  ] <- 0.1
  
  # Colours
  
  tissue_cols <- c(
    
    "Normal colon"     = "#55A868",
    "Primary colon"    = "#4C72B0",
    "Liver metastasis" = "#C44E52",
    
    "Colon"            = "#8172B2",
    "Liver"            = "#F0AD4E"
    
  )
  
  other_tissues <- setdiff(
    levels(plot_df$Tissue),
    names(tissue_cols)
  )
  
  tissue_cols <- c(
    tissue_cols,
    setNames(
      hue_pal()(
        length(other_tissues)
      ),
      other_tissues
    )
  )
  
  # Plot

  p <- ggplot(
    plot_df,
    aes(
      Tissue,
      TPM_plot
    )
  ) +
    
    geom_jitter(
      aes(
        colour = Tissue
      ),
      width = 0.12,
      alpha = 0.85,
      size = 3.8
    ) +
    
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.35,
      colour = "black",
      linewidth = 1.5
    ) +
    
    scale_y_continuous(
      trans = "log2",
      
      breaks = c(
        0.1,
        0.5,
        1,
        2,
        5,
        10,
        20,
        50,
        100
      ),
      
      labels = c(
        "0",
        "0.5",
        "1",
        "2",
        "5",
        "10",
        "20",
        "50",
        "100"
      )
    ) +
    
    scale_x_discrete(
      labels = tissue_labels
    ) +
    
    scale_colour_manual(
      values = tissue_cols
    ) +
    
    geom_hline(
      yintercept = 1,
      colour = "red",
      linewidth = 1,
      linetype = "dashed"
    ) +
    
    geom_hline(
      yintercept = 2,
      colour = "firebrick",
      linewidth = 1,
      linetype = "dotted"
    ) +
    
    annotate(
      "segment",
      x = 1,
      xend = 3,
      y = 0.06,
      yend = 0.06,
      colour = "steelblue4",
      linewidth = 2
    ) +
    
    annotate(
      "text",
      x = 2,
      y = 0.035,
      label = "Discovery cohort\n(GSE50760)",
      colour = "steelblue4",
      fontface = "bold",
      size = 7
    ) +
    
    annotate(
      "segment",
      x = 4,
      xend = length(tissue_order),
      y = 0.06,
      yend = 0.06,
      colour = "darkgreen",
      linewidth = 2
    ) +
    
    annotate(
      "text",
      x = (4 + length(tissue_order))/2,
      y = 0.035,
      label = "GTEx",
      colour = "darkgreen",
      fontface = "bold",
      size = 7
    ) +
    
    labs(
      title = te_of_interest,
      subtitle = "Top-ranked candidate",
      x = NULL,
      y = "TE expression (TPM, log2 scale)"
    ) +
    
    theme_bw(
      base_size = 26
    ) +
    
    theme(
      
      legend.position = "none",
      
      plot.title = element_text(
        face = "bold",
        size = 34,
        hjust = 0.5
      ),
      
      plot.subtitle = element_text(
        face = "italic",
        size = 20,
        hjust = 0.5
      ),
      
      axis.title.y = element_text(
        face = "bold",
        size = 28
      ),
      
      axis.text.y = element_text(
        size = 20
      ),
      
      axis.text.x = ggtext::element_markdown(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 20
      ),
      
      panel.grid.major = element_line(
        colour = "grey70"
      ),
      
      panel.grid.minor = element_line(
        colour = "grey85"
      )
    )
  
  print(p)
  
  ggsave(
    filename = paste0(
      te_of_interest,
      "_Final_Candidate.tiff"
    ),
    plot = p,
    device = "tiff",
    width = 20,
    height = 11,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
  
  cat(
    "Saved:",
    paste0(
      te_of_interest,
      "_Final_Candidate.tiff"
    ),
    "\n"
  )
  
}


##25) TE gene correlation: Tigger7_dup2343 <-> LINC02418 and LTR66_dup140  <-> CASC21

# Load gene TPM vectors
load("LINC02418_TPM_GSE50760.RData")
load("CASC21_TPM_GSE50760.RData")

# Correlation function
run_correlations <- function(
    te,
    gene_name,
    gene_tpm
){
  
  cat("\n")
  cat(rep("=", 80), sep = "")
  cat("\n")
  cat(te, "<->", gene_name, "\n")
  cat(rep("=", 80), sep = "")
  cat("\n\n")
  
  analyses <- list(
    
    All = rownames(metadata_qc),
    
    Normal = rownames(
      metadata_qc[
        metadata_qc$cond == "normal_colon",
        ,
        drop = FALSE
      ]
    ),
    
    Primary = rownames(
      metadata_qc[
        metadata_qc$cond == "primary_colon",
        ,
        drop = FALSE
      ]
    ),
    
    Metastasis = rownames(
      metadata_qc[
        metadata_qc$cond == "liver_metastasis",
        ,
        drop = FALSE
      ]
    )
    
  )
  
  results <- lapply(
    
    names(analyses),
    
    function(group){
      
      samples <- analyses[[group]]
      
      te_expr <- as.numeric(
        unlist(
          tpm[
            te,
            samples,
            drop = FALSE
          ]
        )
      )
      
      gene_expr <- gene_tpm[samples]
      
      keep <- complete.cases(
        te_expr,
        gene_expr
      )
      
      te_expr <- te_expr[keep]
      gene_expr <- gene_expr[keep]
      
      if(length(te_expr) < 3){
        
        return(
          data.frame(
            Group = group,
            N = length(te_expr),
            Spearman_rho = NA,
            Spearman_P = NA,
            Pearson_r = NA,
            Pearson_P = NA
          )
        )
        
      }
      
      spearman <- suppressWarnings(
        cor.test(
          te_expr,
          gene_expr,
          method = "spearman"
        )
      )
      
      pearson <- suppressWarnings(
        cor.test(
          te_expr,
          gene_expr,
          method = "pearson"
        )
      )
      
      data.frame(
        
        Group = group,
        
        N = length(te_expr),
        
        Spearman_rho =
          round(
            as.numeric(
              spearman$estimate
            ),
            3
          ),
        
        Spearman_P =
          signif(
            spearman$p.value,
            3
          ),
        
        Pearson_r =
          round(
            as.numeric(
              pearson$estimate
            ),
            3
          ),
        
        Pearson_P =
          signif(
            pearson$p.value,
            3
          )
        
      )
      
    }
    
  )
  
  results <- bind_rows(results)
  
  print(results)
  
  invisible(results)
  
}

# Tigger7_dup2343 <-> LINC02418
cor_tigger <- run_correlations(
  
  te = "Tigger7_dup2343",
  
  gene_name = "LINC02418",
  
  gene_tpm = linc02418_tpm
  
)

# LTR66_dup140 <-> CASC21
cor_ltr66 <- run_correlations(
  
  te = "LTR66_dup140",
  
  gene_name = "CASC21",
  
  gene_tpm = casc21_tpm
  
)



#Correlation plots: all 52 samples, overall cohort correlation and log2(TPM+1) transformed
#Function
make_corr_plot <- function(
    te,
    gene_name,
    gene_tpm,
    cor_results
){
  
  #############################################################
  # DATA
  #############################################################
  
  samples <- rownames(metadata_qc)
  
  plot_df <- data.frame(
    Sample = samples,
    TE = as.numeric(
      tpm[
        te,
        samples
      ]
    ),
    Gene = gene_tpm[samples]
  ) %>%
    
    left_join(
      metadata_qc %>%
        rownames_to_column("Sample"),
      by = "Sample"
    ) %>%
    
    mutate(
      
      Cohort = case_when(
        
        cond == "normal_colon" ~
          "Normal Colon",
        
        cond == "primary_colon" ~
          "Primary CRC",
        
        cond == "liver_metastasis" ~
          "Liver Metastasis"
        
      ),
      
      TE_log = log2(TE + 1),
      
      Gene_log = log2(Gene + 1)
      
    ) %>%
    
    filter(
      complete.cases(
        TE_log,
        Gene_log
      )
    )
  
  #############################################################
  # LEGEND ORDER
  #############################################################
  
  plot_df$Cohort <- factor(
    plot_df$Cohort,
    levels = c(
      "Normal Colon",
      "Primary CRC",
      "Liver Metastasis"
    )
  )
  
  #############################################################
  # OVERALL CORRELATION
  #############################################################
  
  all_stats <- cor_results %>%
    filter(Group == "All")
  
  rho <- all_stats$Spearman_rho
  
  pval <- format(
    all_stats$Spearman_P,
    scientific = TRUE,
    digits = 2
  )
  
  #############################################################
  # BRIGHTER COLOURS
  #############################################################
  
  cohort_cols <- c(
    "Normal Colon"     = "#4DAF4A",
    "Primary CRC"      = "#377EB8",
    "Liver Metastasis" = "#E41A1C")
  
  #############################################################
  # PLOT
  #############################################################
  
  p <- ggplot(
    plot_df,
    aes(
      x = TE_log,
      y = Gene_log
    )
  ) +
    
    geom_point(
      aes(
        colour = Cohort
      ),
      size = 4.5,
      alpha = 0.90
    ) +
    
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      colour = "black",
      linewidth = 1,
      se = FALSE
    ) +
    
    scale_colour_manual(
      values = cohort_cols,
      breaks = c(
        "Normal Colon",
        "Primary CRC",
        "Liver Metastasis"
      )
    ) +
    
    labs(
      
      title = paste0(
        "Correlation between ",
        te,
        " and ",
        gene_name
      ),
      
      subtitle = paste0(
        "All samples (n = ",
        nrow(plot_df),
        ") | log2(TPM + 1) expression | Spearman \u03c1 = ",
        rho,
        " | P = ",
        pval
      ),
      
      x = paste0(
        te,
        " expression"
      ),
      
      y = paste0(
        gene_name,
        " expression"
      ),
      
      colour = NULL
      
    ) +
    
    theme_classic(
      base_size = 16
    ) +
    
    theme(
      
      plot.title =
        element_text(
          face = "bold",
          size = 20,
          hjust = 0.5
        ),
      
      plot.subtitle =
        element_text(
          size = 12,
          hjust = 0.5
        ),
      
      axis.title =
        element_text(
          face = "bold",
          size = 16
        ),
      
      axis.text =
        element_text(
          size = 13
        ),
      
      legend.position =
        "bottom",
      
      legend.direction =
        "horizontal",
      
      legend.text =
        element_text(
          size = 12
        )
      
    )
  
  return(p)
  
}


# Tigger7_dup2343 ↔ LINC02418
library(tibble)
plot_tigger <- make_corr_plot(
  
  te = "Tigger7_dup2343",
  
  gene_name = "LINC02418",
  
  gene_tpm = linc02418_tpm,
  
  cor_results = cor_tigger
  
)

plot_tigger

ggsave(
  "Tigger7_dup2343_LINC02418_Correlation.tiff",
  plot = plot_tigger,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# LTR66_dup140 ↔ CASC21
plot_ltr66 <- make_corr_plot(
  
  te = "LTR66_dup140",
  
  gene_name = "CASC21",
  
  gene_tpm = casc21_tpm,
  
  cor_results = cor_ltr66
  
)

plot_ltr66

ggsave(
  "LTR66_dup140_CASC21_Correlation.tiff",
  plot = plot_ltr66,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)




#Create the table with the 2 candidates and the coordinates (for the supplementary files)
library(dplyr)
library(gt)
library(magick)

candidate_table <- repeatmasker %>%
  filter(
    `PeakID (cmd=annotatePeaks.pl TEs.bed hg38)` %in% c("LTR66_dup140", "Tigger7_dup2343")
  ) %>%
  transmute(
    TE = `PeakID (cmd=annotatePeaks.pl TEs.bed hg38)`,
    `Genomic Coordinates` = paste0(Chr, ":", Start, "-", End),
    Strand,
    `Length (bp)` = End - Start + 1
  )

gt_table <- candidate_table %>%
  gt() %>%
  tab_header(title = "Final Candidate TEs for Immunotherapy") %>%
  opt_stylize(style = 1, color = "blue")

gtsave(gt_table, "candidate_TE_table.png", vwidth = 1200, vheight = 400)

