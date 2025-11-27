#load Libraries
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(mclust)
library(tximport)
library(biomaRt)

# Define file paths and sample names
samples <- c("SHSY5Y_EWSR1KO_1", "SHSY5Y_EWSR1KO_2", "SHSY5Y_EWSR1KO_3",
             "SHSY5Y_FUSKO_1", "SHSY5Y_FUSKO_2", "SHSY5Y_FUSKO_3",
             "SHSY5Y_TAF15KO_1", "SHSY5Y_TAF15KO_2", "SHSY5Y_TAF15KO_3",
             "SHSY5Y_WT_1", "SHSY5Y_WT_2", "SHSY5Y_WT_3")

files <- file.path("D:/FETseq/Kallisto/Kallisto_bootstrap", samples, "abundance.tsv")

# Assign sample names to the file paths
names(files) <- samples

# Load transcript-to-gene mapping file
tx2gene <- read.csv("D:/FETseq/Kallisto/transcript_to_gene_mapping_isoforms.csv", 
                    col.names = c("transcript", "gene"))

# Verify the mapping file
head(tx2gene)
# 
# # (Optional) Check if IDs match the Kallisto target_id
# abundance_sample <- read.csv("C:/Users/Roman/Desktop/Kallisto_out/SHSY5Y_EWSR1KO_1/abundance.tsv", 
#                              sep = "\t", header = TRUE)
# 
# head(abundance_sample$target_id)

# Import data
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)

# Inspect the resulting gene-level counts
head(txi$counts)

count_data <- as.data.frame(txi)

colnames(count_data) <- samples

count_data <- as.data.frame(txi$counts)
any_na_nan <- any(unlist(lapply(count_data, function(x) any(is.na(x) | is.nan(x)))))
print(paste("Are there any NA or NaN values in the DataFrame?", any_na_nan))

# Convert each column to integers
count_data[] <- lapply(count_data, function(x) as.integer(round(x)))

#count_data <- as.matrix(count_data)
colnames(count_data)
head(count_data)
# 
# #Make some nice colors to facilitate the visualisation of time-points
# mytime                 <- factor(as.character(info$DIV),levels=c(0,3,7,14,22,35))
# mycols_days            <- c("#CCFF00","#33CC33","#669999","#6699FF","#3300FF","#CC33CC")
# names(mycols_days)     <- c(0,3,7,14,22,35)
# mycols                 <- unlist(lapply(info$Cell_Line,function(Z)return(mycols_days[match(as.character(Z),names(mycols_days))])))

#Load the info
info <- read.csv('D:/FETseq/Kallisto/info.csv', header = TRUE, row.names=1)
info$DIV <- NULL
info$Treatment <- c('WT','WT','WT','KO','KO','KO','KO','KO','KO','KO','KO','KO')
info$Sequencing <- c('PE','PE','PE','PE','PE','PE','PE','PE','PE','PE','PE','PE')
colnames(info)
head(info)


#calculate mean and variance of the rows
row_avg<-apply(count_data,1,mean)
row_var <-apply(count_data,1,var)

#Log-transform your count data
count_datal <- log2(1+count_data)

#DESeq2 variance stabilisation
vsd    <- DESeq2::varianceStabilizingTransformation(as.matrix(count_data))


# 1 Fit bimodal distribution
bimdens <- lapply(c(1:ncol(count_datal)),function(IX)return(mclust::densityMclust(data=count_datal[,IX],G=2,plot=FALSE)))

# 2. Identify limit that discriminate foreground from background
# Lims <- unlist(lapply(bimdens,function(x)return(qnorm(0.99,mean=x$parameters$mean[1],sd=sqrt(x$parameters$variance$sigmasq[1])))))
Lims <- unlist(lapply(bimdens, function(x) {
  lower_index <- which.min(x$parameters$mean)  # Identify background component
  background_mean <- x$parameters$mean[lower_index]
  background_sd <- sqrt(x$parameters$variance$sigmasq[lower_index])
  
  # Avoid unrealistically high values
  if (background_sd > 5) background_sd <- 5  # Cap SD at 5
  return(qnorm(0.99, mean = background_mean, sd = background_sd))
}))


# 1. Select reliably expressed genes in each sample
is_expressed_samples <- do.call(lapply(c(1:ncol(count_datal)),function(IX)return(count_datal[,IX]>Lims[IX])),what=cbind)
no_reliably_expressed_genes_samples <- apply(is_expressed_samples,2,sum)

# 2. Select reliably expressed genes in each group
is_expressed_groups <- t(apply(is_expressed_samples,1,function(Z)return(tapply(Z,INDEX=factor(info$Cell_Line),FUN=function(W)return(sum(W)==length(W))))))
no_reliably_expressed_genes_group <- apply(is_expressed_groups,2,sum)

# 3. Select reliably expressed genes in at least one of the groups
is_expressed_global <- apply(is_expressed_groups,1,sum)>=1
num_reliably_expressed_genes <- sum(is_expressed_global)
print(num_reliably_expressed_genes)#15550 genes are reliably expressed in at least one group




# SPLICING
rownames(is_expressed_groups) <- rownames(count_datal)

# 3.2 Select reliably expressed genes in all of the groups
is_expressed_all_groups <- apply(is_expressed_groups,1,sum) == ncol(is_expressed_groups)
num_reliably_expressed_all_groups <- sum(is_expressed_all_groups)
print(num_reliably_expressed_all_groups) # Number of genes reliably expressed in all groups
# Add results as a new column in a dataframe
count_datal_reliablyExpressed <- count_datal
count_datal_reliablyExpressed$Reliably_Expressed_All_Groups <- is_expressed_all_groups


# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

options(timeout = 60)

# Remove version numbers if present
rownames(count_datal_reliablyExpressed) <- sub("\\..*$", "", rownames(count_datal_reliablyExpressed))
head(count_datal_reliablyExpressed)

# Extract gene symbols for the Ensembl IDs present in your dataframe
gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(count_datal_reliablyExpressed),
  mart = ensembl
)

# Merge gene symbols into the dataframe
count_datal_reliablyExpressed <- merge(gene_symbols, count_datal_reliablyExpressed, by.x = "ensembl_gene_id", by.y = "row.names", all.y = TRUE)

# Rename the merged column to make it clearer
colnames(count_datal_reliablyExpressed)[1] <- "Gene_ID"
colnames(count_datal_reliablyExpressed)[2] <- "Gene_Symbol"

# Load the data into a dataframe
# inclusion_levels_df <- read.delim("H:/Мой диск/Kallisto/INCLUSION_LEVELS_FULL-hg38-16.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
inclusion_levels_df <- read.delim("D:/FETseq/vast-tools_results/INCLUSION_LEVELS_FULL-hg38-12-v251.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Keep only genes that are reliably expressed in all samples
filtered_second_df <- inclusion_levels_df[inclusion_levels_df$GENE %in% count_datal_reliablyExpressed$Gene_Symbol[count_datal_reliablyExpressed$Reliably_Expressed_All_Groups], ]

filtered_second_df_naFiltered <- na.omit(filtered_second_df)
filtered_second_df_naFiltered_for_SVD <- filtered_second_df_naFiltered[, c("GENE", "EVENT",
                                                                           "SHSY5Y_EWSR1KO_1", "SHSY5Y_EWSR1KO_2","SHSY5Y_EWSR1KO_3",
                                                                           "SHSY5Y_FUSKO_1", "SHSY5Y_FUSKO_2","SHSY5Y_FUSKO_3",
                                                                           "SHSY5Y_TAF15KO_1", "SHSY5Y_TAF15KO_2","SHSY5Y_TAF15KO_3",
                                                                           "SHSY5Y_WT_1", "SHSY5Y_WT_2","SHSY5Y_WT_3")]
write.csv(filtered_second_df_naFiltered_for_SVD, "filtered_second_df_naFiltered_for_SVD.csv", row.names = TRUE)


split_by_splicing_events <- function(df) {
  # Define patterns for each splicing event type
  event_patterns <- list(
    "EX" = "^HsaEX",
    "INT" = "^HsaINT",
    "ALTD" = "^HsaALTD",
    "ALTA" = "^HsaALTA"
  )
  
  # Initialize a list to store dataframes for each event type
  split_data <- list()
  
  # Filter the dataframe for each event type based on the patterns
  for (event_type in names(event_patterns)) {
    pattern <- event_patterns[[event_type]]
    split_data[[event_type]] <- df %>%
      filter(grepl(pattern, EVENT))
  }
  
  # Capture any events that do not match the known patterns
  known_patterns <- paste(unlist(event_patterns), collapse = "|")
  split_data[["Other"]] <- df %>%
    filter(!grepl(known_patterns, EVENT))
  
  return(split_data)
}
  
  
split_splicing_data <- split_by_splicing_events(filtered_second_df_naFiltered_for_SVD)

split_splicing_data_EX <- t(split_splicing_data$EX[, c(1, 3: 14)])
colnames(split_splicing_data_EX) <- split_splicing_data_EX[1,]
split_splicing_data_EX <- cbind(rownames(split_splicing_data_EX),split_splicing_data_EX)
colnames(split_splicing_data_EX)[1] <- "#samples"
split_splicing_data_EX <- t(split_splicing_data_EX[-1,])
rownames(split_splicing_data_EX) <- make.unique(rownames(split_splicing_data_EX))

split_splicing_data_INT <- t(split_splicing_data$INT[, c(1, 3: 14)])
colnames(split_splicing_data_INT) <- split_splicing_data_INT[1,]
split_splicing_data_INT <- cbind(rownames(split_splicing_data_INT),split_splicing_data_INT)
colnames(split_splicing_data_INT)[1] <- "#samples"
split_splicing_data_INT <- t(split_splicing_data_INT[-1,])
rownames(split_splicing_data_INT) <- make.unique(rownames(split_splicing_data_INT))

split_splicing_data_ALTD <- t(split_splicing_data$ALTD[, c(1, 3: 14)])
colnames(split_splicing_data_ALTD) <- split_splicing_data_ALTD[1,]
split_splicing_data_ALTD <- cbind(rownames(split_splicing_data_ALTD),split_splicing_data_ALTD)
colnames(split_splicing_data_ALTD)[1] <- "#samples"
split_splicing_data_ALTD <- t(split_splicing_data_ALTD[-1,])
rownames(split_splicing_data_ALTD) <- make.unique(rownames(split_splicing_data_ALTD))

split_splicing_data_ALTA <- t(split_splicing_data$ALTA[, c(1, 3: 14)])
colnames(split_splicing_data_ALTA) <- split_splicing_data_ALTA[1,]
split_splicing_data_ALTA <- cbind(rownames(split_splicing_data_ALTA),split_splicing_data_ALTA)
colnames(split_splicing_data_ALTA)[1] <- "#samples"
split_splicing_data_ALTA <- t(split_splicing_data_ALTA[-1,])
rownames(split_splicing_data_ALTA) <- make.unique(rownames(split_splicing_data_ALTA))

write.table(split_splicing_data_EX, "split_splicing_data_EX.tsv", sep = "\t", row.names = T, quote = FALSE)
write.table(split_splicing_data_INT, "split_splicing_data_INT.tsv", sep = "\t", row.names = T, quote = FALSE)
write.table(split_splicing_data_ALTD, "split_splicing_data_ALTD.tsv", sep = "\t", row.names = T, quote = FALSE)
write.table(split_splicing_data_ALTA, "split_splicing_data_ALTA.tsv", sep = "\t", row.names = T, quote = FALSE)

