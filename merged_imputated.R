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



#after getting the csv of pca loading, load the csv into R and make graphs grouped by tissue & sites.


#with imputated data
imputated_data<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/imp_scores_pca_merged_rmoutlier.csv")
  imputated_data<-as.data.frame(imputated_data)
rownames(imputated_data) <-  imputated_data[,1]
# Assuming imputated_data is your dataframe
rownames(imputated_data) <- gsub("^X", "", rownames(imputated_data))


imputated_data<-imputated_data[,-1]

check_G_between <- function(x) {
  # Find '.G_'
  if (grepl("\\.G_", x)) {
    return("Gill")
  } else {
    return("Foot")
  }
}

# Add metadata
imputated_data$Tissue <- sapply(row.names(imputated_data), check_G_between)



# Plot the PCA impute results
pca_imputated_merged<- ggplot(imputated_data, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 2) +
  theme_minimal() +
  theme_classic() +
  labs(
    title = "prcomp() with imputated data merged",
    x = "PC1",
    y = "PC2"
  )

calculate_hull <- function(df) df[chull(df$PC1, df$PC2), ]


# Calculate the convex hulls for each tissue group
hulls <- imputated_data %>%
  group_by(Tissue) %>%
  do(calculate_hull(.))

# Create the ggplot
p <- ggplot(imputated_data, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 2) +
  geom_polygon(data = hulls, aes(x = PC1, y = PC2, fill = Tissue, group = Tissue), 
               alpha = 0.3, linetype = 2, size = 0.5) +
  theme_minimal() +
  theme_classic() +
  labs(
    title = "prcomp() no NA values_unfiltered data merged",
    x = "PC1",
    y = "PC2"
  )

print(p)

ggsave(filename = "pca_merged_rmoutlier.png", p)


#for site specific
imputated_data$sample_prefix <- sub("\\..*", "", rownames(imputated_data))


# Function to transform sample names
transform_sample_name <- function(name) {
  # Replace the first dot with a dash
  name <- sub("\\.", "-", name)
  # Remove everything after the first underscore
  name <- sub("\\.Me", "", name)
  return(name)
}

# Apply the function to the sample_name column
imputated_data$sample_name <- sapply(rownames(imputated_data), transform_sample_name)


merged_filtered<-merge(imputated_data, meta_data, by = "sample_name", all.x = TRUE)



merged_data$origin_site
  
# Plot the PCA impute results
pca_imputated_merged<- ggplot(imputated_data, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 2) +
  theme_minimal() +
  theme_classic() +
  labs(
    title = "prcomp() with imputated data merged",
    x = "PC1",
    y = "PC2"
  )



calculate_hull <- function(df) df[chull(df$PC1, df$PC2), ]


#filtered data
# Plot the PCA impute results
ggplot(imputated_data, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 2) +
  theme_minimal() +
  theme_classic() +
  labs(
    title = "prcomp() with imputated data",
    x = "PC1",
    y = "PC2"
  )

calculate_hull <- function(df) df[chull(df$PC1, df$PC2), ]

# Calculate the convex hulls for each tissue group
hulls <-merged_filtered %>%
  group_by(origin_site) %>%
  do(calculate_hull(.))
merged_filtered$origin_site
# Create the ggplot
p <- ggplot(imputated_data, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 2) +
  # geom_line(aes(group = sample_prefix), color = "black", linetype = 1) +
  #geom_text(aes(label = sample_prefix), hjust = 1.5, vjust = 1.5, size = 3, check_overlap = TRUE) +
  geom_polygon(data = hulls, aes(x = PC1, y = PC2, fill = origin_site, group = origin_site), 
               alpha = 0.3, linetype = 2, size = 0.5) +
  theme_minimal() +
  theme_classic() +
  labs(
    title = "pca_merged_rmoutlier_unfiltered_originsite",
    x = "PC1",
    y = "PC2"
  )


print(p)

ggsave(filename = "pca_merged_rmoutlier_origin_unfiltered.png", p)


#filtering 

meta_data<-read.delim("sample_full.txt", row.names = "sample", stringsAsFactors = FALSE) %>%
  na.omit()
meta_data<-meta_data %>% 
  rownames_to_column(var = "sample_name")




