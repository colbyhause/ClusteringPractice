# Fri Sep 27 09:29:29 2019 ------------------------------
# This code provide easy workflow to perform kmeans and PCA analysis on Flame data.

library(ggbiplot)
library(rcompanion)
library(psych)
library(MASS)
library(RColorBrewer)
library(NbClust)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(tidyverse)

#FOR SPATIAL DATA:
library(naniar)
library(tidyverse)
library(udunits2)
#install.packages("tmap")
library(tmap)
#install.packages("viridis")
library(viridis)
#gg_miss_var(as.data.frame(flame2.edited), show_pct= T)
#install.packages("sf",dependencies = T) 
library(sf)

# ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*> 
#1. Read in Data ----
#data_raw <- readRDS("data/flame3.supercleaned.rds")
data_raw <- readRDS("data/flame2.supercleaned.rds")

#2. Format data ----
which(is.na(data_raw)) # determine where there are NAs
dat_noNAs <- na.omit(data_raw) # remove NAs
which(is.na(dat_noNAs))

dat_cols <- dat_noNAs[ , c("temp", "spCond", "pH", "ODO", "turb",
                           "fDOM", "CHL", "NO3")] # keep only the columns of the parameters 
# NO pH:                        
#dat_cols <- dat_noNAs[ , c("temp", "spCond", "ODO", "turb",
                           #"fDOM", "CHL", "NO3")] # keep only the columns of the parameters ****TOOK pH OUT

#write_csv(dat_cols, "data_output/No_pH/WholeRiverDataRaw_NoNas.csv") #save this as a csv to pull in later

dat_transformed <- apply(dat_cols, 2, log) # transform the data
dat_rescaled <- apply(dat_transformed, 2, scale) # rescale the data


# ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*> 
#3.  Determine number of cluster groups (or centroids) ----
# Use multiple methods to confirm you are choosing the correct number

#3a. Elbow Method:
set.seed(123)
fviz_nbclust(dat_rescaled, kmeans, method = "wss") 

#3b. Silhousette method: 
fviz_nbclust(dat_rescaled, kmeans, method = "silhouette") # says 6 clusters for Transect 2

#3c. Gap Statistic Method: ***WARNING: this method takes a long time to run****
set.seed(123)
gap_stat <- clusGap(dat_rescaled, FUN = kmeans, nstart = 25, K.max = 15, B = 15) 
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat) # look at the results with fviz_gap_stat
# says optimal # of clusters for transect 2 is 9

#3d. NB Clust: ***WARNING: this method takes a long time to run****
# This package provides 30 indices for determining the relevant number of clusters and proposes to users the best clustering scheme from the different results obtained by varying all combinations of number of clusters, distance measures, and clustering methods.

results.nbclust <- NbClust(dat_rescaled, distance = "euclidean",
                      min.nc = 2, max.nc = 15, 
                      method = "complete", index ="all")
factoextra::fviz_nbclust(results.nbclust) + theme_minimal() + ggtitle("NbClust's optimal number of clusters")

# ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*> 
# 4. Whole River Cluster Analysis ----
kmeans_wholeriver <- kmeans(dat_rescaled, 2, nstart = 25) # decided on 2 cluster groups
#**MAKE SURE THAT CLUSTER 1 IS PINK (SMALLER CLUSTER, ON THE LEFT) AND 2 (BIGGER CLUSTER, ON THE RIGHT) SHOULD BE GREEN... OTHERWISE YOU WILL HAVE PROBLEMS WITH THE SUBCLUSTER GROUPS***

#4a. Plot Kmeans analysis
pdf("figure_output/Plots_no_pH/WholeRiverClusterAnalysis_nopH.pdf")
kmeans_wholeriver_plot <- fviz_cluster(kmeans_wholeriver, data = dat_rescaled, labelsize = NA, main = " Whole River Cluster Analysis", xlab = "PC1", ylab ="PC2",ggtheme =theme_minimal())
kmeans_wholeriver_plot
dev.off()

#4b. Put the cluster groupings into dataset----
dat_wholeRiver_clustgrps<- dat_noNAs %>%
  mutate(cluster = kmeans_wholeriver$cluster)

#Save this as its own csv
write_csv(dat_wholeRiver_clustgrps, "data_output/No_pH/WholeRiverData_w_2ClusterGroups_nopH.csv")

# ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>
# ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*> 
# ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*> 

#5. Subcluster Analysis:----
#5a. Read in data:
data_w_clusters<- read_csv("data_output/No_pH/WholeRiverData_w_2ClusterGroups_nopH.csv")

#5b. Make 2 dataframes from those clusters:----
dat_clust1 <- dat_wholeRiver_clustgrps %>% 
  filter(cluster == 1)
write_csv(dat_clust1, "data_output/No_pH/f3_clustgrp1_nopH.csv")

dat_clust2 <- dat_wholeRiver_clustgrps %>% 
  filter(cluster == 2)
write_csv(dat_clust2, "data_output/No_pH/f3_clustgrp2_nopH.csv")

#5c. Normalize and scale data for each cluster group: ----
#cluster1 df:
#dat_clust1.log <- apply(dat_clust1[ ,c("temp", "spCond", "pH", "ODO", 
                                       #"turb","fDOM", "CHL", "NO3")], 2, log) # log transform
dat_clust1.log <- apply(dat_clust1[ ,c("temp", "spCond", "ODO", 
                                       "turb","fDOM", "CHL", "NO3")], 2, log) # log transform **TOOK OUT pH****
dat_clust1.log.scale <- apply(dat_clust1.log, 2, scale) # scale data
#cluster2 df:
#dat_clust2.log <- apply(dat_clust2[ ,c("temp", "spCond", "pH", "ODO", 
                                      #"turb","fDOM", "CHL", "NO3")], 2, log) # log transform
dat_clust2.log <- apply(dat_clust2[ ,c("temp", "spCond", "ODO", 
                                       "turb","fDOM", "CHL", "NO3")], 2, log) # log transform **TOOK OUT pH****
dat_clust2.log.scale <- apply(dat_clust2.log, 2, scale)   # scale data

# ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>
#6. Get Optimal clusters----
#6a. Elbow Method:
set.seed(123)
elbow_clust1 <- fviz_nbclust(dat_clust1.log.scale, kmeans, method = "wss") 
elbow_clust1
elbow_clust2 <- fviz_nbclust(dat_clust2.log.scale, kmeans, method = "wss") 
elbow_clust2
#6b. Silhousette method: 
sil_clust1 <- fviz_nbclust(dat_clust1.log.scale, kmeans, method = "silhouette")
sil_clust1
sil_clust2 <-fviz_nbclust(dat_clust2.log.scale, kmeans, method = "silhouette")
sil_clust2

#6c. Gap Statistic Method: ***WARNING: this method takes a long time to run****
gap_stat_clust1 <- clusGap(dat_clust1.log.scale, FUN = kmeans, nstart = 25, K.max = 50, B = 50)
print(gap_stat_clust1, method = "firstmax")
fviz_gap_stat(gap_stat_clust1) # look at the results with fviz_gap_stat

gap_stat_clust2 <- clusGap(dat_clust1.log.scale, FUN = kmeans, nstart = 25, K.max = 50, B = 50)
print(gap_stat_clust2, method = "firstmax")
fviz_gap_stat(gap_stat_clust2) # look at the results with fviz_gap_stat

#6d. NB Clust: ***WARNING: this method takes a long time to run****
# This package provides 30 indices for determining the relevant number of clusters and proposes to users the best clustering scheme from the different results obtained by varying all combinations of number of clusters, distance measures, and clustering methods.

results.nbclust_clust1 <- NbClust(dat_rescaled, distance = "euclidean",
                           min.nc = 2, max.nc = 15, 
                           method = "complete", index ="all")
factoextra::fviz_nbclust(results.nbclust_clust1) + theme_minimal() + ggtitle("NbClust's optimal number of clusters")


results.nbclust_clust2 <- NbClust(dat_rescaled, distance = "euclidean",
                                min.nc = 2, max.nc = 15, 
                                method = "complete", index ="all")
factoextra::fviz_nbclust(results.nbclust_clust2) + theme_minimal() + ggtitle("NbClust's optimal number of clusters")

# ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*> 
# 7. Sub-Cluster Analysis ----

clust1_2groups <- kmeans(dat_clust1.log.scale, 2, nstart = 25) # 2 subclusters for cluster 1
clust2_4groups <- kmeans(dat_clust2.log.scale, 4, nstart = 25) # 4 subclusters for cluster 2

#7a. Plot Kmeans Subcluster analysis
# Cluster group 1:
pdf(file = "figure_output/Plots_no_pH/f3.grp1.subclusts.clusterPlot_nopH.pdf")
kmeans_subcluster1 <- fviz_cluster(clust1_2groups, data = dat_clust1.log.scale, labelsize = NA, main = "Upper River Cluster Analysis SubCluster (nopH)", ggtheme =theme_minimal(), xlab = "PC1", ylab ="PC2" )
kmeans_subcluster1
dev.off()

#Cluster group 2:
pdf(file = "figure_output/Plots_no_pH/f3.grp2.subclusts.clusterPlot_nopH.pdf")
kmeans_subcluster2<- fviz_cluster(clust2_4groups, data = dat_clust2.log.scale, labelsize = NA, main = "Lower River/ Delta Cluster Anlysis SubClusters (no pH)", ggtheme =theme_minimal(),  palette = "Dark2", xlab = "PC1", ylab ="PC2")
kmeans_subcluster2
dev.off()


# ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>
# ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*> 
# ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*>  ------------- ><)))))*> 

#8. Spatialize the Clusters:----

# 8a. Add subcluster groups to the orginal dataframes of raw data, write to csv:----
# cluster group 1:
dat_subClusts1 <- dat_clust1 %>% 
  mutate(Subcluster = clust1_2groups$cluster)
write_csv(x = dat_subClusts1, path = "data_output/No_pH/flame3_clust1_subclusts_noPH.csv")

# cluster group 2:
dat_subClusts2 <- dat_clust2 %>% 
  mutate(Subcluster = clust2_4groups$cluster) 
write_csv(x = dat_subClusts2, path = "data_output/No_pH/flame3_clust2_subclusts_nopH.csv")

# 8b.Make data Spatial----
# Read in dataframe csvs:----
f3.dat_subClusts1 <- read_csv(file = "data_output/No_pH/flame3_clust1_subclusts_noPH.csv")
f3.dat_subClusts2 <- read_csv(file = "data_output/No_pH/flame3_clust2_subclusts_nopH.csv")
data_w_clusters<- read_csv("data_output/No_pH/WholeRiverData_w_2ClusterGroups_nopH.csv")

f3_wholeRiver.geodata.sf<- st_as_sf(data_w_clusters, 
                                    coords = c("lon", "lat"), #specify lat lon cols
                                    remove = F, # don't remove these lat/lon cols from df
                                    crs = 4326) # add projection

f3_clust1.geodata.sf<- st_as_sf(f3.dat_subClusts1, 
                                coords = c("lon", "lat"), #specify lat lon cols
                                remove = F, # don't remove these lat/lon cols from df
                                crs = 4326) # add projection

f3_clust2.geodata.sf<- st_as_sf(f3.dat_subClusts2, 
                                coords = c("lon", "lat"), #specify lat lon cols
                                remove = F, # don't remove these lat/lon cols from df
                                crs = 4326) # add projection

# write spatial data to creat shp file----
st_write(f3_wholeRiver.geodata.sf, "data_output/No_pH/f3_wholeriver_clusts.shp", delete_dsn = T)
st_write(f3_clust1.geodata.sf, "data_output/No_pH/f3_clust1_subclusts.shp", delete_dsn = T)
st_write(f3_clust2.geodata.sf, "data_output/No_pH/f3_clust2_subclusts.shp", delete_dsn = T)

# 8c. Look at plots with the clusters as the variable:----
# i.plot spatial data: using base R
# Whole River :
pdf("figure_output/Plots_no_pH/f3_spatial_wholeriver_nopH_clust.pdf")
f3.wholeRiver.clustPlot<- plot(f3_wholeRiver.geodata.sf["cluster"],graticule = TRUE, axes = TRUE, main="Whole River Cluster Analysis", pal = c("rosybrown2", "turquoise3"))#, cex.lab=.6)
f3.wholeRiver.clustPlot
dev.off()

# Group1 cluster 
pdf("figure_output/Plots_no_pH/f3_spatial_clust1_nopH_subclusts.pdf")
f3.clust1.subclustsPlot<- plot(f3_clust1.geodata.sf["Subcluster"],graticule = TRUE, axes = TRUE, main="Cluster 1 sub clusters No pH")#, cex.lab=.6)
dev.off()
f3.clust1.subclustsPlot

# Group 2 Cluster
pdf("figure_output/Plots_no_pH/f3_spatial_clust2_nopH_subclusts.pdf")
f3.clust2.subclustsPlot<- plot(f3_clust2.geodata.sf["Subcluster"],graticule = TRUE, axes = TRUE, main="Cluster 2 sub clusters no pH")#, cex.lab=4)
dev.off()
f3.clust2.subclustsPlot

# ii. plot spatial data: using tmap----
# cluster group1:
f3.clust1.subclusts.tmap <- tm_shape(f3_clust1.geodata.sf["Subcluster"]) + 
  tm_layout(title = "Group 1 sub-clusters no pH", frame=F, inner.margin=0.1, title.size = 1.2, title.position = c("center", "top"), legend.position = c("left", "bottom"), legend.text.size = 1.2, legend.title.size = 1) +
  tm_compass(type = "8star", size = 2, position = c(0.8, 0.8)) +
  tm_scale_bar(breaks = c(0, 10, 20), position = c("left", "bottom"),
               size = 0.7) +
  tm_symbols(col="Subcluster",n = 2, title.col = "clusters", palette =  viridis(n=2, direction = -1), size = .1, border.lwd   = NA) 
f3.clust1.subclusts.tmap
tmap_save(f3.clust1.subclusts.tmap, "figure_output/Plots_no_pH/f3_ClusterGrp1_subclusters_tmap.png")

# cluster group2:
f3.clust2.subclusts.tmap <- tm_shape(f3_clust2.geodata.sf["Subcluster"]) + 
  tm_layout(title = "Group 2 sub-clusters no pH", frame=F, inner.margin=0.1, title.size = 1.2, title.position = c("center", "top"), legend.position = c("left", "bottom"), legend.text.size = 1.2, legend.title.size = 1) +
  tm_compass(type = "8star", size = 2, position = c(0.8, 0.8)) +
  tm_scale_bar(breaks = c(0, 10, 20), position = c("left", "bottom"),
               size = 0.7) +
  tm_symbols(col="Subcluster",n = 5, title.col = "clusters", palette =  viridis(n=4, direction = -1), size = .1, border.lwd   = NA) 
f3.clust2.subclusts.tmap
tmap_save(f3.clust1.subclusts.tmap, "figure_output/Plots_no_pH/f3_ClusterGrp2_subclusters_tmap.png")

# Merge the subclust dataframes:----
# Put group1 and subclust group2 dataframes into 1 so you can make a map of the whole river with difference subclust groups:
# first change the subclust group2 numbers so that 1 and 2 are different from the subclust group1. You need to change them to 5-8 instead 4-6 bc when you change 1 to 3 it will eventually change that 3 to 5 as you go down the line.
dat_subClusts2_edit1and2 <- dat_subClusts2

dat_subClusts2_edit1and2$Subcluster[dat_subClusts2_edit1and2$Subcluster == 1]  <- 5
dat_subClusts2_edit1and2$Subcluster[dat_subClusts2_edit1and2$Subcluster == 2]  <- 6
dat_subClusts2_edit1and2$Subcluster[dat_subClusts2_edit1and2$Subcluster == 3]  <- 7
dat_subClusts2_edit1and2$Subcluster[dat_subClusts2_edit1and2$Subcluster == 4]  <- 8

f3.wholeriver.subclusts <- rbind(dat_subClusts1, dat_subClusts2_edit1and2)
# just double checking to make sure subcluster numbers are correct

unique(f3.wholeriver.subclusts$Subcluster)
unique(dat_subClusts2_edit1and2$Subcluster)
unique(dat_subClusts2$Subcluster)
unique(dat_subClusts1$Subcluster)

# Make the whole river subcluster dataframew spatial:

f3_wholeRiver_subclusts.geodata.sf<- st_as_sf(f3.wholeriver.subclusts, 
                                              coords = c("lon", "lat"), #specify lat lon cols
                                              remove = F, # don't remove these lat/lon cols from df
                                              crs = 4326) # add projection

st_write(f3_wholeRiver_subclusts.geodata.sf, "data_output/No_pH/f3_wholeriver_subclusts.shp", delete_dsn = T)

# Plot Whole River with all subclusters ----
pdf("figure_output/Plots_no_pH/f3.spatial.wholeriver.subclusts.pdf")
f3.wholeRiver.subclustsPlot<- plot(f3_wholeRiver_subclusts.geodata.sf["Subcluster"],graticule = TRUE, axes = TRUE)
dev.off()


# Make a shapefile for each cluster group ----
#  This is to import into qgis so you can manipulate those colors
# first: split into dataframes by group  and cluster
grp1_subclust1 <- dat_subClusts1[dat_subClusts1$Subcluster == 1, ]
grp1_subclust2 <- dat_subClusts1[dat_subClusts1$Subcluster == 2, ]

grp2_subclust1 <- dat_subClusts2[dat_subClusts2$Subcluster == 1, ]
grp2_subclust2 <- dat_subClusts2[dat_subClusts2$Subcluster == 2, ]
grp2_subclust3 <- dat_subClusts2[dat_subClusts2$Subcluster == 3, ]
grp2_subclust4 <- dat_subClusts2[dat_subClusts2$Subcluster == 4, ]

# create shapefiles:
f3_grp1_subclust1.geodata.sf<- st_as_sf(grp1_subclust1, 
                                        coords = c("lon", "lat"), #specify lat lon cols
                                        remove = F, # don't remove these lat/lon cols from df
                                        crs = 4326) # add projection

f3_grp1_subclust2.geodata.sf<- st_as_sf(grp1_subclust2, 
                                        coords = c("lon", "lat"), #specify lat lon cols
                                        remove = F, # don't remove these lat/lon cols from df
                                        crs = 4326) # add projection

f3_grp2_subclust1.geodata.sf<- st_as_sf(grp2_subclust1, 
                                        coords = c("lon", "lat"), #specify lat lon cols
                                        remove = F, # don't remove these lat/lon cols from df
                                        crs = 4326) # add projection

f3_grp2_subclust2.geodata.sf<- st_as_sf(grp2_subclust2, 
                                        coords = c("lon", "lat"), #specify lat lon cols
                                        remove = F, # don't remove these lat/lon cols from df
                                        crs = 4326) # add projection

f3_grp2_subclust3.geodata.sf<- st_as_sf(grp2_subclust3, 
                                        coords = c("lon", "lat"), #specify lat lon cols
                                        remove = F, # don't remove these lat/lon cols from df
                                        crs = 4326) # add projection

f3_grp2_subclust4.geodata.sf<- st_as_sf(grp2_subclust4, 
                                        coords = c("lon", "lat"), #specify lat lon cols
                                        remove = F, # don't remove these lat/lon cols from df
                                        crs = 4326) # add projection

# wrte spatial files:
st_write(f3_grp1_subclust1.geodata.sf, "data_output/No_pH/f3_grp1_subclust1.shp", delete_dsn = T)
st_write(f3_grp1_subclust2.geodata.sf, "data_output/No_pH/f3_grp1_subclust2.shp", delete_dsn = T)

st_write(f3_grp2_subclust1.geodata.sf, "data_output/No_pH/f3_grp2_subclust1.shp", delete_dsn = T)
st_write(f3_grp2_subclust2.geodata.sf, "data_output/No_pH/f3_grp2_subclust2.shp", delete_dsn = T)
st_write(f3_grp2_subclust3.geodata.sf, "data_output/No_pH/f3_grp2_subclust3.shp", delete_dsn = T)
st_write(f3_grp2_subclust4.geodata.sf, "data_output/No_pH/f3_grp2_subclust4.shp", delete_dsn = T)


# Note to self:----
# Andrew was concerned about using the component loadings, or coefficient correlations, from the PCA to directly interpret the clusters. Have a day and a half of reseraching and googling this issue, I determined that it is totally fine to make a direct unterpretation of our kmeans cluster plots using the loadings for the PCA.
# Reason: The kmeans analysis only creates a massive distance matrix and calculates clusters based on the position of the centroids, and the number of centroids are determined by you (when you set the number of groups for the cluster). There is now way to plot this directly, as there are too many dimensions.You mnust then, as a step in the plotting process, prerform a PCA to reduce those dimensions down. If you looking at how the function fviz_cluster is working ( getAnywhere(fviz_cluster)) you will see that it calculates a PCA and what you are actually plotting are the PC1 and PC2 scores from those results. Below is also a link that basically says this same thing:
# https://www.datanovia.com/en/lessons/k-means-clustering-in-r-algorith-and-practical-examples/#k-means-algorithm
# So, to validate this is, as in, make sure that the PCA I calulated by hand is the same data as my cluster plot, I can look st the PC scores that I calculated and compare them with the data output from the fviz_cluster results for both the whole river analysis and the subcluster analysis:
# PC scores for whole river PCA analysis: 
# first need to make them point in the positive direction:
head(pca_result$x)
plot(x =  pca_result$x[ , 1], y = pca_result$x[, 2])
# PC scores being plotted from the whole river fviz_zluster function:
kmeans_wholeriver$data
plot(x = kmeans_wholeriver$data[ ,2], y = kmeans_wholeriver$data[ ,3]) # same plot as above

# Now verify for the subcluster results:
# PC scores from PCA on subcluster 1: ( first need to make them point in the positive direction:)
head(pca_clust1_result$x)
plot(x =  pca_clust1_result$x[ , 1], y = pca_clust1_result$x[, 2])

head(pca_clust2_result$x)
plot(x =  pca_clust2_result$x[ , 1], y = pca_clust2_result$x[, 2])
# PC scores being plotted from the subcluster 1 fviz_zluster function:
kmeans_subcluster1$data
plot(x = kmeans_subcluster1$data[ ,2], y = kmeans_subcluster1$data[ ,3]) # same plot as above

kmeans_subcluster2$data
plot(x = kmeans_subcluster2$data[ ,2], y = kmeans_subcluster2$data[ ,3]) # same plot as above
