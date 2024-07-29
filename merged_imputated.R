#cat merged_imputated.R 
library(vegan)
library(edgeR)
library(tidyverse)
library(SYNCSA)
library(mice)

setwd("/hb/groups/kelley_lab/tina/mytilus/Rscripts/data/merged_cov")
targets <- read.delim("sample.txt", row.names = "sample", stringsAsFactors = FALSE)
meta<-targets %>%
  rownames_to_column(var = "sample_name")
Sample <- row.names(targets)

meta_data<-read.delim("sample_full.txt", row.names = "sample", stringsAsFactors = FALSE) %>%
  na.omit()

files <- paste0(Sample,".CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov")
yall <- readBismark2DGE(files, sample.names=Sample)

#pca on original samples

Methylation <- gl(2, 1, ncol(yall), labels = c("Me", "Un"))
Me <- yall$counts[ , Methylation == "Me" ]
Un <- yall$counts[ , Methylation == "Un" ]

#pca on original samples
TotalLibSize <- yall$samples$lib.size[Methylation=="Me"] +
  +                 yall$samples$lib.size[Methylation=="Un"]
yall$samples$lib.size <- rep(TotalLibSize, each=2)
yall$samples

Me <- yall$counts[ , Methylation == "Me" ]
Un <- yall$counts[ , Methylation == "Un" ]


prop_meth_matrix <- Me/(Me+Un)


# Ensure column names start with a letter
colnames(prop_meth_matrix) <- make.names(colnames(prop_meth_matrix), unique = TRUE)
colnames(prop_meth_matrix)
imputed_data <- mice(prop_meth_matrix, m = 1, method = 'pmm', seed = 123)

# Complete the data using the first imputed dataset
complete_data <- complete(imputed_data, 1)

# Perform PCA using prcomp
imp_pca_result <- prcomp(t(complete_data))

# Create df of pca loadings with imputation
imp_scores <- as.data.frame(imp_pca_result$x)

write.csv(imp_scores, "imp_scores_pca_merged.csv")


check_G_between <- function(x) {
  # Find '.G_'
  if (grepl("\\.G_", x)) {
    return("Gill")
  } else {
    return("Foot")
  }
}

# Add metadata
imp_scores$Tissue <- sapply(row.names(imp_scores), check_G_between)

# Plot the PCA impute results
pca_imputated_merged<- ggplot(imp_scores, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 2) +
  theme_minimal() +
  theme_classic() +
  labs(
    title = "prcomp() with imputated data merged",
    x = "PC1",
    y = "PC2"
  )

ggsave(filename = "pca_unfiltered_imputed_merged_rm_small_library.png",pca_imputated_merged)
