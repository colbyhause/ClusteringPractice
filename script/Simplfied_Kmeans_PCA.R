# The goal was to simplify the number of variables within the cluster/ pca analysis to better interpret whats going on. Did this by removing certain parameters 

library(tidyverse)
if(!require(psych)){install.packages("car")}
if(!require(MASS)){install.packages("MASS")}
if(!require(rcompanion)){install.packages("rcompanion")}
library(devtools)
#install_github("vqv/ggbiplot")

library(ggbiplot)
library(rcompanion)
library(psych)
library(MASS)
library(RColorBrewer)
library(NbClust)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization

# Read in data from transect 3 ----
dat <- readRDS("data/flame3.supercleaned.rds")

# Remove missing values:
which(is.na(dat))
dat_noNAs <- na.omit(dat) # so having to use a dataset that is roughly half the size bc of all the rows I have to omit bc of NAs- talk to andrew about this. This cuts the data in half 
which(is.na(dat_noNAs))
# keep only columns that you need:
dat_cols <- dat_noNAs[ , c(6:8, 11:15)]


# transform and scale data ----
dat_transformed <- apply(dat_cols, 2, log) 
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
print(dat_nopH_2)

nopH_wholeriver_2 <- fviz_cluster(dat_nopH_2, data = dat_no_pH, labelsize = NA, main = " No pH Whole River Cluster Analysis", xlab = "PC1", ylab ="PC2",ggtheme =theme_minimal())
nopH_wholeriver_2


# PCA get factor loadings:
pca_result <- prcomp(dat_no_pH, scale = TRUE)
pca_result_all <- prcomp(dat_rescaled, scale = TRUE) 
names(pca_result)
pca_result$rotation

# Calculate variance explained by each PC
#The variance explained by each principal component is obtained by squaring Sdev:
(VE <- pca_result$sdev^2)

#To compute the proportion of variance explained by each principal component, we simply divide the variance explained by each principal component by the total variance explained by all principal components:
PVE <- VE / sum(VE)
round(PVE, 2) 

# plot PCA: NEXT STEP: MAKE BETTER PCA PLOTS

pdf(file = "figure_output/PCA_biplot_WholeRiver_nopH_092619.pdf")
biplot(pca_result, scale = 0, main ="Biplot Whole River T3" )
dev.off()

# bi plot using ggbiplot:
# No pH data included:
biplot_no_pH<- ggbiplot(pca_result, ellipse = T, groups = )+
  coord_equal(ratio = 0.4)

biplot_all<- ggbiplot(pca_result_all, ellipse = T)+
  coord_equal(ratio = 0.4)


########################################Look at sub cluster plots with the no_pH data:
# First do kmeans analysis on whole river data minus pH:
#1. Determine number of cluster groups:
#1a. Elbow Method: 
set.seed(123)
fviz_nbclust(dat_no_pH, kmeans, method = "wss") # looks like 2 clusters are optimal 

#1b.Silhousette method:
fviz_nbclust(dat_no_pH, kmeans, method = "silhouette") # says 2 clusters are optimal

#1c.Gap statistic method: this takes a long time to run
set.seed(123)
gap_stat <- clusGap(dat_rescaled, FUN = kmeans, nstart = 25,
                    K.max = 50, B = 50)
print(gap_stat, method = "firstmax")
#visualize the results with fviz_gap_stat which suggests four clusters as the optimal number of clusters.
fviz_gap_stat(gap_stat) 

#1e. Determine number of clusters with NbClust: warning : this also takes a long time to run 
es.nbclust <- NbClust(dat_no_pH, distance = "euclidean",
                      min.nc = 2, max.nc = 15, 
                      method = "complete", index ="all")
factoextra::fviz_nbclust(es.nbclust) + theme_minimal() + ggtitle("NbClust's optimal number of clusters")

#2. Kmeans analysis:
kmeans_no_pH_2clusts <- kmeans(dat_no_pH, 2, nstart = 25)
print(kmeans_no_pH_2clusts)

#plot Kmeans:
pdf("figure_output/Plots_no_pH/f3.wholeRiverClusterAnalysis_nopH.pdf")
kmeans_wholeriver_nopH <- fviz_cluster(kmeans_no_pH_2clusts, data = dat_no_pH, labelsize = NA, main = " Whole River Cluster Analysis (no pH)", xlab = "PC1", ylab ="PC2",ggtheme =theme_minimal())
kmeans_wholeriver_nopH
dev.off()


#3. Sub Cluster Analysis ###########################################################################################

#3a.Add cluster to data----
# make sure you switch it from dataframe to matrix
dat_no_pH <- data.frame(dat_no_pH)

#Add cluster groups:
wholeRiver_nopH_clustgrps<-  dat_no_pH%>%
  mutate(cluster = wholeriver_nopH_clusters) # this is found in the code above under First Cluster Analysis

#3b. Make dataframes of those clusters:
dat_nopH_clust1 <- wholeRiver_nopH_clustgrps %>% 
  filter(cluster == 1)
class(wholeRiver_nopH__clustgrps) #################### STOPPED HERE: NEED TO MAKE THIS WHOLE CODE MORE USER FRIENDLY
dat_clust1_params<- dat_clust1[ , c(6:8, 11:15)]
write_csv(x = dat_clust1_params, "data_output/No_pH/f3.clust1group.csv")


dat_clust2 <- wholeRiver_nopH_clustgrps %>% 
  filter(cluster == 2)
write_csv(x = dat_clust2_params, "data_output/f3.clust2group.csv")

class(wholeRiver_nopH__clustgrps)
