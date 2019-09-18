# The goal was to simplify the number of variables within the cluster/ pca analysis to better interpret whats going on. Did this by removing certain parameters 

library(tidyverse)
if(!require(psych)){install.packages("car")}
if(!require(MASS)){install.packages("MASS")}
if(!require(rcompanion)){install.packages("rcompanion")}

library(rcompanion)
library(psych)
library(MASS)
library(RColorBrewer)

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization

# Read in data from transect 3 ----
dat <- readRDS("data/flame3.supercleaned.rds")

# transform and scale data ----
dat_transformed <- apply(dat[ , c(6:8, 10:15)], 2, log) 
dat_rescaled <- apply(dat_transformed, 2, scale)

# look at the correlation matrix, which is the correlation factor between each parameter.
cor(dat_rescaled) # looks like pH is the only parameter that is not correlated with any other. Most are highly correlated with eachother

# Try removing pH 
dat_no_pH <- dat_rescaled[ , -c(3)]

# Determine Clust groups:----
#Elbow Method:
set.seed(123)
fviz_nbclust(dat_no_pH, kmeans, method = "wss") # looks like 2 clusters are optimal 

#Silhousette method:
fviz_nbclust(dat_no_pH, kmeans, method = "silhouette") # says 2 clusters are optimal

# Gap statistic method:
set.seed(123)
gap_stat <- clusGap(dat_rescaled, FUN = kmeans, nstart = 25,
                    K.max = 30, B = 50)
print(gap_stat, method = "firstmax")
#We can visualize the results with fviz_gap_stat which suggests four clusters as the optimal number of clusters.

fviz_gap_stat(gap_stat) # shows that 10 iterations were not suffient (did not converge), so trying 50...shows that optimal is 12...WARNING: 50 iterations take a long time to run

library(NbClust)
es.nbclust <- NbClust(dat_no_pH, distance = "euclidean",
                      min.nc = 2, max.nc = 15, 
                      method = "complete", index ="all")
factoextra::fviz_nbclust(es.nbclust) + theme_minimal() + ggtitle("NbClust's optimal number of clusters") # says 4 clusters are best

# Plot kmean analysis----
set.seed(123)
dat_nopH_2 <- kmeans(dat_no_pH, 2, nstart = 25)
print(dat_nocond_2)

nopH_wholeriver_2 <- fviz_cluster(dat_nopH_2, data = dat_no_pH, labelsize = NA, main = " No cond Whole River Cluster Analysis", xlab = "PC1", ylab ="PC2",ggtheme =theme_minimal())
nopH_wholeriver_2


# PCA get factor loadings:
pca_result <- prcomp(dat_no_pH, scale = TRUE)
names(pca_result)
pca_result$rotation

# Calculate variance explained by each PC
#The variance explained by each principal component is obtained by squaring Sdev:
(VE <- pca_result$sdev^2)

#To compute the proportion of variance explained by each principal component, we simply divide the variance explained by each principal component by the total variance explained by all principal components:
PVE <- VE / sum(VE)
round(PVE, 2) 

# plot PCA: NEXT STEP: MAKE BETTER PCA PLOTS

pdf(file = "figure_output/PCA_biplot_WholeRiver_nopH.pdf")
biplot(pca_result, scale = 0, main ="Biplot Whole River T3" )
dev.off()
