# PCA analysis following the UC Davis Business Anayltics tutorial from this site:
# https://uc-r.github.io/pca

# Principle Component Analysis

library(tidyverse)  # data manipulation and visualization
library(gridExtra)  # plot arrangement

# Load data:----
# dataframe:
dat <- readRDS("data/flame3.supercleaned.rds")

# remove missing values:
which(is.na(dat))
dat_noNAs <- na.omit(dat) # so having to use a dataset that is roughly half the size bc of all the rows I have to omit bc of NAs- talk to andrew about this. This cuts the data in half 
which(is.na(dat_noNAs))
# keep only columns that you need:
dat_cols <- dat_noNAs[ , c(6:8, 11:15)]

#prep data:----
#. compute variance of each variable: i think this is just to visulaize the variance?
apply(dat_cols, 2, var)
# log transform:
dat_log <- apply(dat_cols, 2, log) 
# now rescale data:
dat_log_scaled<- apply(dat_log, 2, scale)

# calc PCA for Whole River DATA SET: not subclustered one:----
# this  function centers the variables to have mean zero.
pca_result <- prcomp(dat_log_scaled, scale = TRUE)
names(pca_result) # The center and scale components correspond to the means and standard deviations of the variables that were used for scaling prior to implementing PCA.
print(pca_result)
# means
pca_result$center
# standard deviations
pca_result$scale
# PCA loadings:----
# file:///Users/colbyhause/Downloads/AtchleyOct19%20(1).pdf this pdf give definition of loadings 
pca_result$rotation <- -pca_result$rotation # switch the sign bc by default eigenvectors point in the neg direction
pca_result$rotation

#principal components scores from our results:----
pca_result$x <- - pca_result$x # also make them point in the positive direction
head(pca_result$x)
# Plot PCA:----
# this by default only plots principle components 1 and 2
pdf(file = "figure_output/PCA_biplot_WholeRiver.pdf")
biplot(pca_result, scale = 0, main ="Biplot Whole River T3" )
dev.off()
# can make it plot other PCs with choices argument:
#biplot(pca_result, scale = 0, choices = 3:4) shows just 1 blob bc no/very little variace 

#The variance explained by each principal component is obtained by squaring these values:
  (VE <- pca_result$sdev^2)

#To compute the proportion of variance explained by each principal component, we simply divide the variance explained by each principal component by the total variance explained by all principal components:
PVE <- VE / sum(VE)
round(PVE, 2) 
#0.78 0.12 0.06 0.02 0.01 0.01 0.00 0.00
# so 78% of the variace is explained by PC1, 12% of variance explained by PC2

#calculations by hand:----
#these give same results as above, I am just doing them by hand to understand what is going on
# Calculate eigenvalues & eigenvectors ----
f3.cov <- cov(dat_log_scaled) # this is the covariance matrix
f3.eigen <- eigen(f3.cov) # 
str(f3.eigen)

#Extract the loadings
(phi <- f3.eigen$vectors[,1:2])
#swittch eigen vectors to a pos direction:
phi <- -phi
row.names(phi) <- c("temp", "spCond", "pH", "ODO", "turb", "fdom", "CHL", "NO3")
colnames(phi) <- c("PC1", "PC2")
phi

# To get the percent of variance in allthe variables
#This is the same as dividing the factor's eigenvalueby the number of variables.
f3.eigen$values
(f3.eigen$values)/8
#0.7848133483 0.1233119461 0.0552868682 0.0210671182 0.0087514945 0.0051203883 0.0012521998 0.0003966367
#this is the same output from the code above that comnputes the variance explained by each principle component

# PCA on sub clusters:----

# load data that is already seperated into clustered groups
f3.clust1 <-read_csv("data_output/f3.clust1group.csv")
f3.clust2 <-read_csv("data_output/f3.clust2group.csv")

# log and scale data:
# Cluster1 group:
# log transform:
f3.clust1.log <- apply(f3.clust1, 2, log) 
# now rescale data:
f3.clust1.log.scaled<- apply(f3.clust1.log, 2, scale)

#cluster 2 group:
# log transform:
f3.clust2.log <- apply(f3.clust2, 2, log) 
# now rescale data:
f3.clust2.log.scaled<- apply(f3.clust2.log, 2, scale)

# subclust PCA results:
pca_clust1_result <- prcomp(f3.clust1.log.scaled, scale = TRUE)
pca_clust2_result <- prcomp(f3.clust2.log.scaled, scale = TRUE)

# look at plot:
pdf(file = "figure_output/PCA_biplot_subclust1.pdf")
biplot(pca_clust1_result, scale = 0) # cluster group 1
dev.off()

pdf(file = "figure_output/PCA_biplot_subclust2.pdf")
biplot(pca_clust2_result, scale = 0) # cluster group 2 
dev.off()
# subcluster loadings:
# file:///Users/colbyhause/Downloads/AtchleyOct19%20(1).pdf this pdf gives definition of loadings, which are akso the correlation coefficients 
#pca_clust1_result$rotation <- -pca_clust1_result$rotation ....the website says to make multiply by -1 to switch values to the postive direction, but since I didnt do this with the PC scores I am going to stay consistent and opt out of doing it here 
pca_clust1_result$rotation 

#pca_clust2_result$rotation <- -pca_clust2_result$rotation
pca_clust2_result$rotation

#The variance explained by each principal component is obtained by squaring these values:
VE_clust1 <- pca_clust1_result$sdev^2
VE_clust2 <- pca_clust2_result$sdev^2

#compute the proportion of variance explained by each principal component:
PVE_clust1 <- VE_clust1 / sum(VE_clust1)
round(PVE_clust1, 2) 
#  0.69 0.19 0.07 0.02 0.01 0.01 0.00 0.00

PVE_clust2 <- VE_clust2 / sum(VE_clust2)
round(PVE_clust2, 2) 
# 0.53 0.24 0.12 0.05 0.03 0.01 0.01 0.00


# Summary of loadings table to show AR:----

# Whole River Cluster Analysis PCA derived loadings: 
wholeriver_loadings <- pca_result$rotation 
wholeriver_loadings <- as.data.frame(wholeriver_loadings)
biplot(pca_result, scale = 0, main= "Whole River") # plot
write_csv(wholeriver_loadings, "data_output/wholeriver_loadings.csv")

# Group1 (Upper River) Sub Cluster Analysis PCA-derived loadings: 
clust1_loadings <- pca_clust1_result$rotation 
clust1_loadings <- as.data.frame(clust1_loadings)
biplot(pca_clust1_result, scale = 0, main = "Upper River sub-cluster") # plot
write_csv(clust1_loadings, "data_output/clust1_loadings.csv")
# Group2 (Lower River/ Delta) Sub Cluster Analysis PCA-derived loadings: 
clust2_loadings <- pca_clust2_result$rotation
clust2_loadings <- as.data.frame(clust2_loadings)
biplot(pca_clust2_result, scale = 0, main = "Lower River/Delta sub-cluster") # plot
write_csv(clust2_loadings, "data_output/clust2_loadings.csv")


# Determining PC scores ----

#1.Calculate eigenvalues & eigenvectors
f3.cov <- cov(dat_log_scaled) # this is the covariance matrix
f3.eigen <- eigen(f3.cov) # 
str(f3.eigen)

#Extract the loadings
(phi <- f3.eigen$vectors[,1:8])
#swittch eigen vectors to a pos direction:
phi <- -phi
row.names(phi) <- c("temp", "spCond", "pH", "ODO", "turb", "fdom", "CHL", "NO3")
colnames(phi) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")
phi

# Calculate Principal Components scores
PC1 <- as.matrix(dat_log_scaled) %*% phi[,1]
PC2 <- as.matrix(dat_log_scaled) %*% phi[,2]
PC3 <- as.matrix(dat_log_scaled) %*% phi[,3]
PC4 <- as.matrix(dat_log_scaled) %*% phi[,4]
PC5 <- as.matrix(dat_log_scaled) %*% phi[,5]
PC6 <- as.matrix(dat_log_scaled) %*% phi[,6]
PC7 <- as.matrix(dat_log_scaled) %*% phi[,7]
PC8 <- as.matrix(dat_log_scaled) %*% phi[,8]

PCscores_df<- data.frame(PC1, PC2, PC3, PC4, PC5,PC6,PC7, PC8)

# Look at Just PC1 and PC2 on graph:
ggplot(PCscores_df[,1:2]) +
       geom_point(aes(x = PC1,y = PC2)) 
# write PC1 and PC2 to own dataframe to do a kmeans analysis on
PCscores_1_2 <- PCscores_df[,1:2] 

