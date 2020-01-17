#Determining Clustering Algorithm

# The goal of this code is to quantitatively determine the clustering method to use and choosing the appropriate number of clusters.
# Resources include: 
# Numerical Ecology with R book (2nd edition) (NER)

#install.packages("labdsv")
#install.packages("gclus")
#install.packages("ade4")
library(labdsv)
library(ade4)
library(cluster)
library(tidyverse)
library(vegan)
library(gclus)
library(RColorBrewer)
library(dendextend)

# Means dataset----
# load data ---
means_dat <- read_csv("data/Mean_fromPredicted_AllParams_across_actualRKMS.csv")
#Log transform and scale:----
dat_transformed <- apply(means_dat[ , 1:8], 2, log) # transform the data
dat_rescaled <- apply(dat_transformed, 2, scale) # rescale the data
#dat_stand <- decostand(dat_transformed[ , 1:8], "standardize") # this is the same out put as the line above (NER p. 56)
#dat_norm <- decostand(dat_transformed[ , 1:8], "normalize") # this is what the book does, except without log transforming it before hand (NER p. 56)

# create distance matrix:----
dist_m <- dist(x = dat_rescaled, method = "euclidean")

# Running data through many different Hierarchical clustering methods:
# Complete linkage:----
means_complete <-  hclust(dist_m, method = "complete")
plot(means_complete, hang = -5)

# upgma (average) clustering----
means_upgma <- hclust(dist_m, method = "average")
plot(means_upgma, hang = -5)

# Single Clustering ----
means_single <- hclust(dist_m, method = "single")
plot(means_single, hang = -5)

# ward clustering----
means_ward <- hclust(dist_m, method = "ward.D")
plot(means_ward, hang = -5)

# Comparing Clustering Results ----
# Cophenetic Correlation:
# The cophenetic distance btwn 2 onjects in a dendrogram is the distance where 2 objects become memebers of the same group

#complete----
means_complete_coph <- cophenetic(means_complete)
cor(dist_m, means_complete_coph) #.8767 best 
#upgma----
means_upgma_coph <- cophenetic(means_upgma)
cor(dist_m, means_upgma_coph) #.8751 close 2nd
# single----
means_single_coph <- cophenetic(means_single)
cor(dist_m, means_single_coph) # .484 last
#ward ----
means_ward_coph <- cophenetic(means_ward)
cor(dist_m, means_ward_coph) #.8398 3rd

# can visualize this with shepard-like diagrams: plot original distances against cophenetic distances:
par(mfrow= c(2, 2))
plot(dist_m, means_complete_coph, xlab = "Eucl Distance", ylab = "Cophenetic Distance", asp = 1, main = "Complete Linkage", paste("Cophenetic correlation ", round(cor(dist_m, means_complete_coph), 3)))
abline(0, 1)
lines(lowess(dist_m, means_complete_coph), col = "red")
# after this one plots, do this for the rest of them^

#Gower Distance: Another way to compare cluster methods----
#Computed as the sum of squared differences between the original and cophenetic distances. Cluster method that produces the smallest gower number is the best method
gower.complete <- sum((dist_m - means_complete_coph)^2)
gower.complete/1000000 #2nd best               
gower.upgma <- sum((dist_m - means_upgma_coph)^2)                    
gower.upgma/1000000  # says this is the best method                     
gower.single<- sum((dist_m - means_single_coph)^2)
gower.single/1000000 # 3rd best
gower.ward <- sum((dist_m - means_ward_coph)^2)
gower.ward/1000000 # worst

# After running these 2 tests it looks like upgma would be the best method, with complete in 2nd

# Determining clusters:----
# Fusion Values----
# These help to determine where to cut the dendrogram. Fusion levels are the dissimilarity values where a fusion between 2 branches of a dendrogram occurs

par(mfrow = c(2, 2))
summary(means_complete)
# complete linkage fusion values
plot(means_complete$height, nrow(means_dat):2, type = "S", main = "Fusion Levels-Eucl Dist- Complete", ylab = "k (number of clusters", xlab = "h (node height)", col = "grey")
text(means_complete$height, nrow(means_dat):2, nrow(means_dat):2, col = "red", cex = .8)
# there is too much data for this method to be useful i think

# Silhouette widths:----
# The greater the value is, the better the object is clustered
# Create empty vector in which the asw (average sil widths) will be written
asw <- numeric(nrow(means_dat))
#Retrieve and write asw values into the vector
# complete linkage
for (k in 2:30) { # only doing uo to 30 clusters
  sil <- silhouette(cutree(means_complete, k = k), dist_m)
  asw[k] <- summary(sil)$avg.width
}

# best (largest) silhouette width
k.best<- which.max(asw)

plot(1:30, asw[1:30], type = "h", main = "Silhouette- optimal number of clusters, Complete", xlab = "k (number of groups)", ylab = "Average silhouette width")
axis(1, k.best, paste("optimum", k.best, sep = "\n"), col = "red", font = 2, col.axis = "red")
points(k.best, max(asw), pch = 16, col = "red", cex = 1.5)
cat("", "Silhouette-optimal number of clusters k = ", k.best, "\n", "wih an average silhouette width of", max(asw), "\n")

# Do this same thing but for the UPGMA method
# Create empty vector in which the asw (average sil widths) will be written
asw <- numeric(nrow(means_dat))
#Retrieve and write asw values into the vector
# complete linkage
for (k in 2:30) { # only doing uo to 30 clusters
  sil <- silhouette(cutree(means_upgma, k = k), dist_m)
  asw[k] <- summary(sil)$avg.width
}

# best (largest) silhouette width
k.best<- which.max(asw)

plot(1:30, asw[1:30], type = "h", main = "Silhouette- optimal number of clusters, UPGMA", xlab = "k (number of groups)", ylab = "Average silhouette width", cex = .5)
axis(1, k.best, paste("optimum", k.best, sep = "\n"), col = "red", font = 2, col.axis = "red")
points(k.best, max(asw), pch = 16, col = "red", cex = 1.5)
cat("", "Silhouette-optimal number of clusters k = ", k.best, "\n", "wih an average silhouette width of", max(asw), "\n")

#this returns 8 as the optimal number of clusters. Going to go with the UPGMA method because the gower method strongly chose this method, and it was basically tied for first using the cophenetic correlation method

# Mantel Coeff:
grpdist <- function(x)
{
  require(cluster)
  gr <- as.data.frame(as.factor(x))
  distgr <- daisy(gr, "gower")
  distgr
}

kt <- data.frame(k = 1:30, r = 0) # 

for (i in 2:30) { # only going to look at clusters up to 30
  gr <- cutree(means_upgma, i)
  distgr <- grpdist(gr)
  mt <- cor(dist_m, distgr, method = "pearson")
  kt[i, 2] <- mt
}
kt
k.best <- which.max(kt$r)
k.best # says 6 clusters are best 
#plot:
plot(kt$k, kt$r, type = "h", main = "Mantel-optimal clusters-upgma", xlab = "k", ylab = "Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep = "/n"), col = "red", font = 2, col.axis = "red")
points(k.best, max(kt$r), pch = 16, col = "red", cex = 1.5)

# Silhouette plot of the final partition:
k <- 8
cutg <- cutree(means_upgma, k = k)
sil <- silhouette(cutg, dist_m)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(means_dat)

# from book:
k <- 4
cutg <- cutree(spe.ch.ward, k = k)
sil <- silhouette(cutg, spe.ch)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(means_dat)

plot(silo)
plot(silo, main = "Silhouette plot", cex.names = .4, col = cutg+1, nmax.lab = 100)

# with flame data:
# for some reason data isnt showing up on plot, but still the data for the sil width for each cluster is useful
#k_i <- 6
k_i <- 8 # 8 has an avg sil width of .64, vs .56 for 6 clusters. so 8 is best
cutg_i <- cutree(means_upgma, k = k_i)
sil_i <- silhouette(cutg_i, dist_m)
silo_i <- sortSilhouette(sil_i)
rownames(silo) <- row.names(means_dat)

plot(silo_i)
plot(silo_i, main = "Silhouette plot", cex.names = .4, col = cutg+1, nmax.lab = 100)
min(sil_i)
max(sil_i)

# cut clusters :
means_ct<- cutree(means_upgma, k = 8) # chose 12 to break them into 2 clusters

# Dendextend:-----
means_dend_clust <- as.dendrogram(means_upgma)
means_color_dend <- color_branches(dend = means_dend_clust, k = 8)
plot(means_color_dend)

# CV DATASET:----
# load data ---
cv_dat <- read_csv("data/CV_fromPredicted_AllParams_across_actualRKMS.csv")
#Log transform and scale:----
dat_transformed <- apply(cv_dat[ , 1:8], 2, log) # transform the data
dat_rescaled <- apply(dat_transformed, 2, scale) # rescale the data
#dat_stand <- decostand(dat_transformed[ , 1:8], "standardize") # this is the same out put as the line above (NER p. 56)
#dat_norm <- decostand(dat_transformed[ , 1:8], "normalize") # this is what the book does, except without log transforming it before hand (NER p. 56)

# create distance matrix:----
dist_cv <- dist(x = dat_rescaled, method = "euclidean")

# After going through all this code, its cleaer the UPGMA is the best option. I have to remove the outliers that make the really small clusters (<10 dets) and then will compare with the rest of the other methods 

# Running data through many different Hierarchical clustering methods:
# Complete linkage:----
cv_complete <-  hclust(dist_cv, method = "complete") # still contains small group
plot(cv_complete, hang = -5)

# upgma (average) clustering----
pdf("figure_output/DeterminingOptimalClusters/upgma_dendro_original.pdf")
cv_upgma <- hclust(dist_cv, method = "average") # has small group
plot(cv_upgma, hang = -5)
cv_upgma
dev.off()

# Single Clustering ----
cv_single <- hclust(dist_cv, method = "single") 
plot(cv_single, hang = -5)

# ward clustering----
cv_ward <- hclust(dist_cv, method = "ward.D") # looks pretty good 
plot(cv_ward, hang = -5)


# Comparing Clustering Results ----
# Cophenetic Correlation:
# The cophenetic distance btwn 2 onjects in a dendrogram is the distance where 2 objects become memebers of the same group

#complete----
cv_complete_coph <- cophenetic(cv_complete)
cor(dist_cv, cv_complete_coph) #0.8983385 2nd
#upgma----
cv_upgma_coph <- cophenetic(cv_upgma)
cor(dist_cv, cv_upgma_coph) # 0.9448407 # best
# single----
cv_single_coph <- cophenetic(cv_single)
cor(dist_cv, cv_single_coph) # 0.7626776 Last
#ward ----
cv_ward_coph <- cophenetic(cv_ward)
cor(dist_cv, cv_ward_coph) # 0.8180965 #3rd

# can visualize this with shepard-like diagrams: plot original distances against cophenetic distances:
# RUNNING THIS PLOT FREEZES R UP- trying running over night
# Complete:
par(mfrow= c(1, 1))
plot(dist_cv, cv_complete_coph, xlab = "Eucl Distance", ylab = "Cophenetic Distance", asp = 1,xlim = c(0, sqrt(2)), ylim = c(0, sqrt(2)),
main = c("Complete Linkage", paste("Cophenetic correlation ", round(cor(dist_cv, cv_complete_coph),3))))
abline(0, 1)
lines(lowess(dist_cv, cv_complete_coph), col = "red")

# Upgma:
plot(dist_cv, cv_upgma_coph, xlab = "Eucl Distance", ylab = "Cophenetic Distance", asp = 1,xlim = c(0, sqrt(2)), ylim = c(0, sqrt(2)),
main = c("UPGMA Linkage", paste("Cophenetic correlation ", round(cor(dist_cv, cv_upgma_coph),3))))
abline(0, 1)
lines(lowess(dist_cv, cv_upgma_coph), col = "red")

#Gower Distance: Another way to compare cluster methods----
#Computed as the sum of squared differences between the original and cophenetic distances. Cluster method that produces the smallest gower number is the best method
gower.complete <- sum((dist_cv - cv_complete_coph)^2)
gower.complete/1000000 # 660.1258  2nd           
gower.upgma <- sum((dist_cv - cv_upgma_coph)^2)                    
gower.upgma/1000000  # 22.38617 best              
gower.single<- sum((dist_cv - cv_single_coph)^2)
gower.single/1000000 # 345.3784 3rd
gower.ward <- sum((dist_cv - cv_ward_coph)^2)
gower.ward/1000000 #4584179316 worst

# After running these 2 tests it looks like upgma would be the best method but it produces a weird small clusters, need to get rid of that

#------
cv_dat <- read_csv("data/CV_fromPredicted_AllParams_across_actualRKMS.csv")
#Log transform and scale:----
dat_transformed <- apply(cv_dat[ , 1:8], 2, log) # transform the data
dat_rescaled <- apply(dat_transformed, 2, scale) # rescale the data
#dat_stand <- decostand(dat_transformed[ , 1:8], "standardize") # this is the same out put as the line above (NER p. 56)
#dat_norm <- decostand(dat_transformed[ , 1:8], "normalize") # this is what the book does, except without log transforming it before hand (NER p. 56)

# create distance matrix:----
dist_cv <- dist(x = dat_rescaled, method = "euclidean")
dist_unscaled<- dist(x = cv_dat, method = "euclidean")


cv_scaled_upgma <- hclust(dist_cv, method = "average") # also has small group
plot(cv_scaled_upgma, hang = -5)

cv_unscaled_upgma <- hclust(dist_unscaled, method = "average") # also has small group
plot(cv_unscaled_upgma, hang = -5)

# ok: so what i was seeing before in the difference in dendrograms after removing the two records: i was actually not putting the scaled and transformed file into the dist function, so it was spitting out a completely different plot. so, i need to just keep trying this method but continue to scale and transform data 

#------
# cutree:-----
cv_ct_outliers <- cutree(cv_upgma, k = 2) #  break them into 2 clusters
unique(cv_ct_outliers)
length(which(cv_ct_outliers == 2)) # only 2 observations in cluster 2, remove these points
length(cv_ct_outliers)
# put the clusters into original df:
cv_hclust_clusters <- cv_dat %>% 
  mutate(cluster = cv_ct_outliers)
unique(cv_hclust_clusters$cluster)
length(which(cv_hclust_clusters$cluster == 2)) # 2 rows, correct

# Removing the 2 outlier points
outliers <- which(cv_hclust_clusters$cluster == 2)
cv_dat_edited <- cv_hclust_clusters[-outliers, ]
unique(cv_dat_edited$cluster) # only 1 unique cluster

#outliers <- which(cv_hclust_clusters$cluster == 2)
#cv_dat_edited <- cv_dat[-outliers, ]

# transform and scale new data:
#Log transform and scale:----
dat_transformed_edited <- apply(cv_dat_edited[ , 1:8], 2, log) # transform the data
dat_rescaled_edited <- apply(dat_transformed_edited, 2, scale) # rescale the data

# run dist function on edited data:
dist_cv_edited <- dist(x = dat_rescaled_edited, method = "euclidean")
# run the upgma method again:
#upgma edited-----
pdf("figure_output/DeterminingOptimalClusters/upgma_dendro_firstedits.pdf")
cv_upgma_edited_test <- hclust(dist_cv_edited, method = "average") # also has small group
plot(cv_upgma_edited_test, hang = -5) # still showing another weird cluster, so this process again:
cv_upgma_edited_test
dev.off()
# looks the exact same. try again:

# Remove the new small cluster that still exists:----
# cutree:
cv_ct_outliers2 <- cutree(cv_upgma_edited_test, k = 2) #  break them into 2 clusters
unique(cv_ct_outliers2)
length(which(cv_ct_outliers2 == 2)) # now 8 observations in cluster 2, remove these points
# put the clusters into original df:
cv_hclust_clusters2 <- cv_dat_edited %>% 
  mutate(cluster = cv_ct_outliers2)
unique(cv_hclust_clusters2$cluster)
length(which(cv_hclust_clusters2$cluster == 2)) # 8 rows, correct

# Removing the 2 outlier points
outliers2 <- which(cv_hclust_clusters2$cluster == 2)
cv_dat_edited2 <- cv_dat_edited[-outliers2, ]
nrow(cv_dat_edited2) # looks good

# transform and scale new data:
#Log transform and scale:----
dat_transformed_edited2 <- apply(cv_dat_edited2[ , 1:8], 2, log) # transform the data
dat_rescaled_edited2 <- apply(dat_transformed_edited2, 2, scale) # rescale the data

# run dist function on edited2 data:
dist_cv_edited2 <- dist(x = dat_rescaled_edited2, method = "euclidean")

# run the upgma method again:
#upgma edited-----
cv_upgma_edited2 <- hclust(dist_cv_edited2, method = "average") # also has small group
plot(cv_upgma_edited2, hang = -5) # still showing another weird cluster, so this process again:
# yay this finally fixed it!
# Now use cv_upgma_edited2 as your upgma cluster 

# Determining clusters:----
# Fusion Values----
# These help to determine where to cut the dendrogram. Fusion levels are the dissimilarity values where a fusion between 2 branches of a dendrogram occurs

par(mfrow = c(2, 2))
summary(cv_complete)
# complete linkage fusion values
plot(cv_complete$height, nrow(cv_dat):2, type = "S", main = "Fusion Levels-Eucl Dist- Complete", ylab = "k (number of clusters", xlab = "h (node height)", col = "grey")
text(cv_complete$height, nrow(cv_dat):2, nrow(cv_dat):2, col = "red", cex = .8)
# there is too much data for this method to be useful i think

# Silhouette widths:----
# The greater the value is, the better the object is clustered
# Create empty vector in which the asw (average sil widths) will be written
asw <- numeric(nrow(cv_dat))
#Retrieve and write asw values into the vector
# complete linkage- just to see what this one says
for (k in 2:30) { # only doing uo to 30 clusters
  sil <- silhouette(cutree(cv_complete, k = k), dist_cv)
  asw[k] <- summary(sil)$avg.width
}

# best (largest) silhouette width
k.best<- which.max(asw)

plot(1:30, asw[1:30], type = "h", main = "Silhouette- optimal number of clusters, Complete", xlab = "k (number of groups)", ylab = "Average silhouette width")
axis(1, k.best, paste("optimum", k.best, sep = "\n"), col = "red", font = 2, col.axis = "red")
points(k.best, max(asw), pch = 16, col = "red", cex = 1.5)
cat("", "Silhouette-optimal number of clusters k = ", k.best, "\n", "wih an average silhouette width of", max(asw), "\n")

#UPGMA method
# Create empty vector in which the asw (average sil widths) will be written
asw <- numeric(nrow(cv_dat))
#Retrieve and write asw values into the vector
# complete linkage
for (k in 2:30) { # only doing uo to 30 clusters
  sil <- silhouette(cutree(cv_upgma, k = k), dist_cv)
  asw[k] <- summary(sil)$avg.width
}

# best (largest) silhouette width
k.best<- which.max(asw)

plot(1:30, asw[1:30], type = "h", main = "Silhouette- optimal number of clusters, UPGMA", xlab = "k (number of groups)", ylab = "Average silhouette width", cex = .5)
axis(1, k.best, paste("optimum", k.best, sep = "\n"), col = "red", font = 2, col.axis = "red")
points(k.best, max(asw), pch = 16, col = "red", cex = 1.5)
cat("", "Silhouette-optimal number of clusters k = ", k.best, "\n", "wih an average silhouette width of", max(asw), "\n")

#Both the complete and the upgma returns 2 as the optimal number of clusters. Going to go with the UPGMA method because the both the cophenetic correlation and the gower method strongly chose this method

# Mantel Coeff:
grpdist <- function(x)
{
  require(cluster)
  gr <- as.data.frame(as.factor(x))
  distgr <- daisy(gr, "gower")
  distgr
}

kt <- data.frame(k = 1:30, r = 0) # 

for (i in 2:30) { # only going to look at clusters up to 30
  gr <- cutree(cv_upgma, i)
  distgr <- grpdist(gr)
  mt <- cor(dist_cv, distgr, method = "pearson")
  kt[i, 2] <- mt
}
kt
k.best <- which.max(kt$r)
k.best # says 9 clusters are best 
#plot:
plot(kt$k, kt$r, type = "h", main = "Mantel-optimal clusters-upgma", xlab = "k", ylab = "Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep = "/n"), col = "red", font = 2, col.axis = "red")
points(k.best, max(kt$r), pch = 16, col = "red", cex = 1.5)
# according to the plot, it looks like 7-11 clusters look pretty even

# Silhouette plot of the final partition:
# from book:
k <- 4
cutg <- cutree(spe.ch.ward, k = k)
sil <- silhouette(cutg, spe.ch)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(cv_dat)

plot(silo)
plot(silo, main = "Silhouette plot", cex.names = .4, col = cutg+1, nmax.lab = 100)

# with flame data:
# for some reason data isnt showing up on plot, but still the data for the sil width for each cluster is useful
k_i <- 2
k_i <- 9 # this is the result from the mantel coeff,  has an avg sil width of .64, vs .56 for 6 clusters. so 8 is best
cutg_i <- cutree(cv_upgma, k = k_i)
sil_i <- silhouette(cutg_i, dist_cv)
silo_i <- sortSilhouette(sil_i)
rownames(silo) <- row.names(cv_dat)

plot(silo_i)
plot(silo_i, main = "Silhouette plot", cex.names = .4, col = cutg_i+1, nmax.lab = 100)
min(sil_i)
max(sil_i)
# for 2 clusters, the avg sil width is .66, for 9 clusters, the average sil width is .62. So they are pretty close. 
# cut clusters :
cv_ct<- cutree(cv_upgma, k = 2) # chose 12 to break them into 2 clusters

# Dendextend:-----
cv_dend_clust <- as.dendrogram(cv_upgma)
cv_color_dend <- color_branches(dend = cv_dend_clust, k = 2)
plot(cv_color_dend) # this is no good, need to remove the small cluster. 
