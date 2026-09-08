#GSE144259

#Load the TE count matrix

#Read the TE count matrix from file (rows = TE subfamilies, columns = samples)
counts <- read.table(
  "GSE144259_TE_counts_matrix.final.txt",
  header      = TRUE,
  row.names   = 1,
  sep         = "\t",
  check.names = FALSE
)

#Check the dimensions of the matrix
dim(counts)  #1180 TEs x 9 samples


#2) Load and prepare the metadata

#Load the GEOquery library to download sample metadata from GEO
library(GEOquery)

#Download the GEO ExpressionSet for dataset GSE144259
gse  <- getGEO("GSE144259", GSEMatrix = TRUE, getGPL = FALSE)

#Extract the sample metadata table (rows = samples, columns = annotations)
meta <- pData(gse[[1]])

#Build a clean metadata dataframe with relevant columns
metadata <- data.frame(
  GSM        = meta$geo_accession,
  tissue     = gsub("tissue: ", "", meta$`tissue:ch1`),
  patient    = factor(sub(".*(CRC[0-9]+).*", "\\1", meta$title)),
  stringsAsFactors = FALSE
)

#Rename tissue labels to shorter, cleaner names
metadata$tissue <- c(
  "Primary colon cancer tissue"          = "primary_colon",
  "Matched normal colon tissue"          = "normal_colon",
  "Liver metastasis colon cancer tissue" = "liver_metastasis"
)[metadata$tissue]

#Set rownames of metadata to match column names of the count matrix
rownames(metadata) <- colnames(counts)

#Verify that metadata rows and count matrix columns are aligned
all(rownames(metadata) == colnames(counts))


#Create the score tables for the GSE144259, using the 2 signatures of the discovery cohort (GSE50760)
library(scuttle)

# TPM matrix
tpm_mat <- calculateTPM(
  counts,
  lengths = rep(1, nrow(counts))
)


#68 TEs- tumor-associated signature

#Use the signatures from the discovery cohort
tumor_associated <- read.table(
  "tumor_associated_68_TEs.txt",
  stringsAsFactors = FALSE
)[,1]

metastasis_associated <- read.table(
  "metastasis_associated_141_TEs.txt",
  stringsAsFactors = FALSE
)[,1]



#Compute the signature scores

#Tumor-associated
expr_val_68 <- tpm_mat[tumor_associated, ]

score_val_68 <- colMeans(expr_val_68, na.rm = TRUE)

df_val_68 <- data.frame(
  sample = names(score_val_68),
  score  = score_val_68,
  cond   = metadata$tissue
)

write.csv(
  df_val_68,
  "TE_signature68_scores_GSE144259.csv",
  row.names = FALSE
)


#Metastasis-associated

expr_val_141 <- tpm_mat[metastasis_associated, ]

score_val_141 <- colMeans(expr_val_141, na.rm = TRUE)

df_val_141 <- data.frame(
  sample = names(score_val_141),
  score  = score_val_141,
  cond   = metadata$tissue
)

write.csv(
  df_val_141,
  "TE_signature141_scores_GSE144259.csv",
  row.names = FALSE
)






