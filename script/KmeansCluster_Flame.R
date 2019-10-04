# first attempt at cluster analysis with own data:

if(!require(psych)){install.packages("car")}
if(!require(MASS)){install.packages("MASS")}
if(!require(rcompanion)){install.packages("rcompanion")}

library(rcompanion)
library(psych)
library(MASS)
library(RColorBrewer)
dat <- readRDS("data/flame3.supercleaned.rds")

# Determining the right transformation to use:

#1. looking at untransformed data in histogram:
turbhist.3 <- ggplot(dat, aes(x = turb)) +
  geom_histogram(bins = 50) +
  xlim(c(0, 70)) +
  xlab(label = "Turbidity (FNU)") +
  ggtitle(label = "Flame3: Turbidity")+
  theme(plot.title = element_text(size = 9), 
        axis.title = element_text(size = 9), 
        axis.text.y = element_text(size = 7))
turbhist.3

# creasting vector of turbidity values to play with:
turb <- dat$turb

# Log-transforming turb values, plotting witht a normal curve superimposed:
T_log = log(turb)
plotNormalHistogram(T_log)

# Trying the square root transformation:
T_sqrt = sqrt(turb)
plotNormalHistogram(T_sqrt)

# going to go with log transformation- kind of an arbitrary decision...
# transform data:
dat_transformed <- apply(dat[ , c(6:8, 10:15)], 2, log) 
# now rescale data:
dat_rescaled <- apply(dat_transformed, 2, scale)

# first try clustering on just turbidity data:
turb_rescaled <-dat_rescaled[ ,6]
#remove missing values:
turb_rescaled <- turb_rescaled[!is.na(turb_rescaled)]
plot(turb_rescaled) # just checking 

# Now determine how many clusters: 
# Elbow method ----
kmean_withinss <- function(k) {
  cluster <- kmeans(turb_rescaled, k) # run the algorithm k times
  return (cluster$tot.withinss) # Store the total within clusters sum of squares
  
}

max_k <-20
# Run algorithm over a range of k 
wss <- sapply(2:max_k, kmean_withinss) # run function over range of 2- 20
wss

# create a data frame with the results of the algorithm
# # Create a data frame to plot the graph
elbow <-data.frame(2:max_k, wss)

# Plot the graph with ggplot
ggplot(elbow, aes(x = X2.max_k, y = wss)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1, max_k, by = 1))
# from looking at the elbow plot- maybe 5 groups would be best?

# Silhouette Method ---- NOT WORKING WITH DATA: i think it need to be a dataframe not a vector
## function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(turb_rescaled, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(df))
  mean(ss[, 3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")

fviz_nbclust(turb_rescaled, kmeans, method = "silhouette")

# Gap statstic method ---- ALSO NOT WORKING BC DATA NOT A DATAFRAME
set.seed(123)
gap_stat <- clusGap(turb_rescaled, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
print(gap_stat, method = "firstmax")

# Trying cluster analysis with 5:

turb_cluster <-kmeans(turb_rescaled, 5)
turb_cluster$cluster
turb_cluster$centers
turb_cluster$size 


# Try looking at clusters in a heat map:

#install.packages("RColorBrewer")
library(RColorBrewer)

center_t <-turb_cluster$centers 
center_t # rows refer to numeration of the cluster and the columns are the variables used by the algorithm.
#The values are the average score by each cluster for the interested column. Standardization makes the interpretation easier. Positive values indicate the z-score for a given cluster is above the overall mean. For instance, cluster 2 has the highest price average among all the clusters.

#1. Build a dataframe:
#create dataset with the cluster number
cluster <- c(1: 5)
center_df <- data.frame(cluster, center)

#2. Reshape the data
center_reshape <- gather(center_df, features, values, price_scal: trend_scal)# this just reshapes the data- cant do this with just turbidity 
head(center_reshape)

# Create the palette
hm.palette <-colorRampPalette(rev(brewer.pal(10, 'RdYlGn')),space='Lab')

#3. Plot the heat map
ggplot(data = center_df, aes(x = features, y = cluster, fill = values)) +
  scale_y_continuous(breaks = seq(1, 7, by = 1)) +
  geom_tile() +
  coord_equal() +
  scale_fill_gradientn(colours = hm.palette(90)) +
  theme_classic()

####### Trying again using methods from 3rd tutorial from ClusterExamples script page: ----

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization

# dataframe:
dat <- readRDS("data/flame3.supercleaned.rds")

# remove missing values:
which(is.na(dat))
dat_noNAs <- na.omit(dat) # so having to use a dataset that is roughly half the size bc of all the rows I have to omit bc of NAs- talk to andrew about this. This cuts the data in half 
which(is.na(dat_noNAs))
# keep only columns that you need:
dat_cols <- dat_noNAs[ , c(6:8, 11:15)]
class(dat_cols)
# log transform:
dat_transformed <- apply(dat_cols, 2, log) 
# now rescale data:
dat_rescaled <- apply(dat_transformed, 2, scale)

distance <- get_dist(dat_rescaled) # calculate distances
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) # visualize by creating a dustance matrix- This takes a long time to run bc of the large dataset!

# Determine Clust groups:

#Elbow Method:
wss <- function(k) {
  kmeans(df, k, nstart = 10 )$tot.withinss
}

k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

#OR better yet use this single function for elbow method:
set.seed(123)
fviz_nbclust(dat_rescaled, kmeans, method = "wss") # looks like 2 clusters are optimal 

#Silhousette method:
fviz_nbclust(dat_rescaled, kmeans, method = "silhouette") # says 2 clusters are optimal

# Gap statistic method:
set.seed(123)
gap_stat <- clusGap(dat_rescaled, FUN = kmeans, nstart = 25,
                    K.max = 50, B = 50)
print(gap_stat, method = "firstmax")
#We can visualize the results with fviz_gap_stat which suggests four clusters as the optimal number of clusters.

fviz_gap_stat(gap_stat) # shows that 10 iterations were not suffient (did not converge), so trying 50...shows that optimal is 12...WARNING: 50 iterations take a long time to run


# Determine number of clusters with NbClust package: this package provides 30 indices for determining the relevant number of clusters and proposes to users the best clustering scheme from the different results obtained by varying all combinations of number of clusters, distance measures, and clustering methods.
#install.packages("NbClust")
library(NbClust)
es.nbclust <- NbClust(dat_rescaled, distance = "euclidean",
                      min.nc = 2, max.nc = 15, 
                      method = "complete", index ="all")
factoextra::fviz_nbclust(es.nbclust) + theme_minimal() + ggtitle("NbClust's optimal number of clusters")

# comes out to optimal number of clusters is 2

# since 2 methods gave 2 optiminal groups and the other gave 12 optimal groups, try running k means with both:
# Whole River Cluster Analysis:----

set.seed(123)
final_2 <- kmeans(dat_rescaled, 2, nstart = 25)
print(final_2)

#final_12 <- kmeans(dat_rescaled, 12, nstart = 25)
#print(final_12)

#final_4 <- kmeans(dat_rescaled, 4, nstart = 25)
#print(final_12)
#visualize:

pdf("figure_output/f3.wholeRiverClusterAnalysis.pdf")
kmeans_wholeriver <- fviz_cluster(final_2, data = dat_rescaled, labelsize = NA, main = " Whole River Cluster Analysis", xlab = "PC1", ylab ="PC2",ggtheme =theme_minimal())
kmeans_wholeriver
#dat_from_fviz <-test$data
dev.off()
#fviz_cluster(final_12, data = dat_rescaled,  labelsize = NA)

# Put the cluster group back into original data 
dat_wholeRiver_clustgrps<- dat_noNAs %>%
  mutate(cluster = final_2$cluster)

# Sub Cluster Analysis:----

# Load and prep data:
# dataframe:
dat <- readRDS("data/flame3.supercleaned.rds")

# remove missing values:
which(is.na(dat))
dat_noNAs <- na.omit(dat) # so having to use a dataset that is roughly half the size bc of all the rows I have to omit bc of NAs- talk to andrew about this. This cuts the data in half 
which(is.na(dat_noNAs))
# keep only columns that you need:
dat_cols <- dat_noNAs[ , c(6:8, 10:15)]
#Add cluster to data----
#1. add the cluster group assignments to the dataset:

dat_wholeRiver_clustgrps<- dat_noNAs %>%
  mutate(cluster = final_2$cluster) # this is found in the code above under First Cluster Analysis , eventually rearrange all this to make better work flow
  class(dat_wholeRiver_clustgrps) 

#2. make 2 dataframe from those clusters:----
class(dat_wholeRiver_clustgrps)
dat_clust1 <- dat_wholeRiver_clustgrps %>% 
  filter(cluster == 1)
dat_clust1_params<- dat_clust1[ , c(6:8, 11:15)]
write_csv(x = dat_clust1_params, "data_output/f3.clust1group.csv")

dat_clust2 <- dat_wholeRiver_clustgrps %>% 
  filter(cluster == 2)
dat_clust2_params<- dat_clust2[ , c(6:8, 11:15)]
write_csv(x = dat_clust2_params, "data_output/f3.clust2group.csv")

#3. normalize and scale data: ----
#cluster1 df:
dat_clust1.log <- apply(dat_clust1_params, 2, log) # log transform
dat_clust1.log.scale <- apply(dat_clust1.log, 2, scale) # scale data
#cluster2 df:
dat_clust2.log <- apply(dat_clust2_params, 2, log) # log transform
dat_clust2.log.scale <- apply(dat_clust2.log, 2, scale)   # scale data

#4. Get distances to look at distance matrix:----
#cluster group1:
dist.clust1 <- get_dist(dat_clust1.log.scale) # calculate distances
fviz_dist(dist.clust1, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) # this takes a long time to run!
# cluster group2:
dist.clust2 <- get_dist(dat_clust2.log.scale) # calculate distances
fviz_dist(dist.clust2, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) # this takes a long time to run!


#5. Get Optimal clusters----
#Elbow Method:
#cluster group1:
wss1 <- function(k) {
  kmeans(dat_clust1.log.scale, k, nstart = 10 )$tot.withinss
}
k.values <- 1:15

# cluster group2:
wss2 <- function(k) {
  kmeans(dat_clust2.log.scale, k, nstart = 10 )$tot.withinss
}
k.values <- 1:15
# extract wss for 2-15 clusters
wss_clust1 <- map_dbl(k.values, wss1)
wss_clust2 <- map_dbl(k.values, wss2)

plot(k.values, wss_clust1,
     type="b", pch = 19, frame = FALSE,
     main = "cluster group1",
     xlab="Number of clusters K ",
     ylab="Total within-clusters sum of squares")

plot(k.values, wss_clust2,
     type="b", pch = 19, frame = FALSE,
     main = "cluster group2",
     xlab="Number of clusters K ",
     ylab="Total within-clusters sum of squares")

#OR better yet use this single function for elbow method: 
set.seed(123)
fviz_nbclust(dat_clust1.log.scale, kmeans, method = "wss") #maybe 3 clusters optimal?
fviz_nbclust(dat_clust2.log.scale, kmeans, method = "wss") # maybe 4 clusters?

#Silhousette method: shows the same plot for both, dont use this method
fviz_nbclust(dat_clust1.log.scale, kmeans, method = "silhouette") # says 2 clusters are optimal
fviz_nbclust(dat_clust2.log.scale, kmeans, method = "silhouette") # says 3 clusters are optimal

# Gap statistic method: not using this method to determine clusters
set.seed(123)
gap_stat1 <- clusGap(dat_clust1.log.scale, FUN = kmeans, nstart = 25,
                    K.max = 17, B = 50)
print(gap_stat1, method = "firstmax")

set.seed(123)
gap_stat2 <- clusGap(dat_clust2.log.scale, FUN = kmeans, nstart = 25,
                     K.max = 25, B = 50)
print(gap_stat2, method = "firstmax")

# Nb cluster method: provides 30 indices for determining the relevant number of clusters and proposes to users the best clustering scheme from the different results obtained by varying all combinations of number of clusters, distance measures, and clustering methods.

library(NbClust)
res.nbclust1 <- NbClust(dat_clust1.log.scale, distance = "euclidean",
                      min.nc = 2, max.nc = 15, 
                      method = "complete", index ="all")
factoextra::fviz_nbclust(res.nbclust1) + theme_minimal() + ggtitle("NbClust's optimal # of clusters Group1") # Says optimal cluster number  = 2 

res.nbclust2 <- NbClust(dat_clust2.log.scale, distance = "euclidean",
                        min.nc = 2, max.nc = 15, 
                        method = "complete", index ="all")
factoextra::fviz_nbclust(res.nbclust2) + theme_minimal() + ggtitle("NbClust's optimal # of clusters Group2")
# Says optimal cluster number = 4

#We can visualize the results with fviz_gap_stat which suggests four clusters as the optimal number of clusters.
fviz_gap_stat(gap_stat1)  # said didnt converge in 10 iterations...
fviz_gap_stat(gap_stat2) # said didnt converge in 10 iterations...

# SubCluster plots----
# Cluster group 1: trying a couple different cluster groups
set.seed(123)
cluster1_final2 <- kmeans(dat_clust1.log.scale, 2, nstart = 25) 
#cluster1_final3 <- kmeans(dat_clust1.log.scale, 3, nstart = 25)

# Cluster group 2: trying a couple different cluster groups
#cluster2_final2 <- kmeans(dat_clust2.log.scale, 2, nstart = 25)
#cluster2_final3 <- kmeans(dat_clust2.log.scale, 3, nstart = 25)
cluster2_final4 <- kmeans(dat_clust2.log.scale, 4, nstart = 25) 

#visualize Plots:----
#Cluster group 1:
pdf(file = "figure_output/f3.grp1.subclusts.clusterPlot.pdf")
kmeans_subcluster1 <- fviz_cluster(cluster1_final2, data = dat_clust1.log.scale, labelsize = NA, main = "Upper River Cluster Analysis", ggtheme =theme_minimal(), xlab = "PC1", ylab ="PC2" )
kmeans_subcluster1
dev.off()
#fviz_cluster(cluster1_final3, data = dat_clust1.log.scale, labelsize = NA, main = "Cluster Group 1 Plot")

#Cluster group 2:
#fviz_cluster(cluster2_final2, data = dat_clust2.log.scale, labelsize = NA, main = "Cluster Group 2 Plot")
#fviz_cluster(cluster2_final3, data = dat_clust2.log.scale, labelsize = NA, main = "Cluster Group 2 Plot")
pdf(file = "figure_output/f3.grp2.subclusts.clusterPlot.pdf")
kmeans_subcluster2<- fviz_cluster(cluster2_final4, data = dat_clust2.log.scale, labelsize = NA, main = "Lower River/ Delta Cluster Anlysis", ggtheme =theme_minimal(),  palette = "Dark2", xlab = "PC1", ylab ="PC2")
kmeans_subcluster2
dev.off()


#Spatialize Sub clusters: ----

# Add subcluster groups to dataframes, write to csv:----
# cluster group 1:
dat_subClusts1 <- dat_clust1 %>% 
  mutate(cluster = cluster1_final2$cluster)
write_csv(x = dat_subClusts1, path = "data_output/flame3_clust1_subclusts.csv")

# cluster group 2:
dat_subClusts2 <- dat_clust2 %>% 
  mutate(cluster = cluster2_final4$cluster) 
write_csv(x = dat_subClusts2, path = "data_output/flame3_clust2_subclusts.csv")

# read in dataframe csvs:----
f3.dat_subClusts1 <- read_csv(file = "data_output/flame3_clust1_subclusts.csv")
f3.dat_subClusts2 <- read_csv(file = "data_output/flame3_clust2_subclusts.csv")

# Make data Spatial----
#install.packages("naniar")
library(naniar)
library(tidyverse)
library(udunits2)
#install.packages("tmap")
library(tmap)
#install.packages("viridis")
library(viridis)
gg_miss_var(as.data.frame(flame2.edited), show_pct= T)
#install.packages("sf",dependencies = T) 
library(sf)
f3_wholeRiver.geodata.sf<- st_as_sf(dat_wholeRiver_clustgrps, 
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
st_write(f3_wholeRiver.geodata.sf, "data_output/f3_wholeriver_clusts.shp", delete_dsn = T)
st_write(f3_clust1.geodata.sf, "data_output/f3_clust1_subclusts.shp", delete_dsn = T)
st_write(f3_clust2.geodata.sf, "data_output/f3_clust2_subclusts.shp", delete_dsn = T)

# Look at plots with the clusters as the variable:----
# plot spatial data: using base R
# Whole River :
pdf("figure_output/f3.spatial.wholeriver.clust.pdf")
f3.wholeRiver.clustPlot<- plot(f3_wholeRiver.geodata.sf["cluster"],graticule = TRUE, axes = TRUE, main="Whole River Cluster Analysis", pal = c("rosybrown2", "turquoise3"))#, cex.lab=.6)
dev.off()

# Group1 cluster 
pdf("figure_output/f3.spatial.clust1.subclusts.pdf")
f3.clust1.subclustsPlot<- plot(f3_clust1.geodata.sf["cluster"],graticule = TRUE, axes = TRUE, main="Cluster 1 sub clusters")#, cex.lab=.6)
dev.off()
f3.clust1.subclustsPlot

# Group 2 Cluster
pdf("figure_output/f3.spatial.clust2.subclusts.pdf")
f3.clust2.subclustsPlot<- plot(f3_clust2.geodata.sf["cluster"],graticule = TRUE, axes = TRUE, main="Cluster 2 sub clusters")#, cex.lab=4)
dev.off()
f3.clust2.subclustsPlot

# plot spatial data: using tmap----
# cluster group1:
f3.clust1.subclusts.tmap <- tm_shape(f3_clust1.geodata.sf["cluster"]) + 
  tm_layout(title = "Cluster Group 1 sub-clusters", frame=F, inner.margin=0.1, title.size = 1.2, title.position = c("center", "top"), legend.position = c("left", "bottom"), legend.text.size = 1.2, legend.title.size = 1) +
  tm_compass(type = "8star", size = 2, position = c(0.8, 0.8)) +
  tm_scale_bar(breaks = c(0, 10, 20), position = c("left", "bottom"),
               size = 0.7) +
  tm_symbols(col="cluster",n = 2, title.col = "clusters", palette =  viridis(n=2, direction = -1), size = .1, border.lwd   = NA) 
f3.clust1.subclusts.tmap
tmap_save(f3.clust1.subclusts.tmap, "figure_output/f3.ClusterGrp1_subclusters.tmap.png")

# cluster group2:
f3.clust2.subclusts.tmap <- tm_shape(f3_clust2.geodata.sf["cluster"]) + 
  tm_layout(title = "Cluster Group 1 sub-clusters", frame=F, inner.margin=0.1, title.size = 1.2, title.position = c("center", "top"), legend.position = c("left", "bottom"), legend.text.size = 1.2, legend.title.size = 1) +
  tm_compass(type = "8star", size = 2, position = c(0.8, 0.8)) +
  tm_scale_bar(breaks = c(0, 10, 20), position = c("left", "bottom"),
               size = 0.7) +
  tm_symbols(col="cluster",n = 5, title.col = "clusters", palette =  viridis(n=4, direction = -1), size = .1, border.lwd   = NA) 
f3.clust2.subclusts.tmap
tmap_save(f3.clust1.subclusts.tmap, "figure_output/f3.ClusterGrp1_subclusters.tmap.png")

# Merge the subclust group1 and subclust group2 dataframes into 1 so you can make a map of the whole river with difference subclust groups:
# first change the subclust group2 numbers so that 1 and 2 are different from the subclust group1
dat_subClusts2_edit1and2 <- dat_subClusts2

dat_subClusts2_edit1and2$cluster[dat_subClusts2_edit1and2$cluster == 1]  <- 5
dat_subClusts2_edit1and2$cluster[dat_subClusts2_edit1and2$cluster == 2]  <- 6
dat_subClusts2_edit1and2$cluster[dat_subClusts2_edit1and2$cluster == 3]  <- 7
dat_subClusts2_edit1and2$cluster[dat_subClusts2_edit1and2$cluster == 4]  <- 8

f3.wholeriver.subclusts <- rbind(dat_subClusts1, dat_subClusts2_edit1and2)
unique(f3.wholeriver.subclusts$cluster)
unique(dat_subClusts2_edit1and2$cluster)
unique(dat_subClusts2$cluster)
unique(dat_subClusts1$cluster)

# need to fix issue that 

# make spatial:
f3_wholeRiver_subclusts.geodata.sf<- st_as_sf(f3.wholeriver.subclusts, 
                                    coords = c("lon", "lat"), #specify lat lon cols
                                    remove = F, # don't remove these lat/lon cols from df
                                    crs = 4326) # add projection

st_write(f3_wholeRiver_subclusts.geodata.sf, "data_output/f3_wholeriver_subclusts.shp", delete_dsn = T)

# plot:
pdf("figure_output/f3.spatial.wholeriver.subclusts.pdf")
f3.wholeRiver.subclustsPlot<- plot(f3_wholeRiver_subclusts.geodata.sf["Subcluster"],graticule = TRUE, axes = TRUE)
dev.off()

f3.clust2.subclustsPlot<- plot(f3_clust2.geodata.sf["cluster"],graticule = TRUE, axes = TRUE, main="Cluster 2 sub clusters")#, cex.lab=4)

# Make a shapefile for each cluster group to import into qgis so you can manipulate those colors
# first: split into dataframes by group  and cluster
grp1_subclust1 <- dat_subClusts1[dat_subClusts1$cluster == 1, ]
grp1_subclust2 <- dat_subClusts1[dat_subClusts1$cluster == 2, ]

grp2_subclust1 <- dat_subClusts2[dat_subClusts2$cluster == 1, ]
grp2_subclust2 <- dat_subClusts2[dat_subClusts2$cluster == 2, ]
grp2_subclust3 <- dat_subClusts2[dat_subClusts2$cluster == 3, ]
grp2_subclust4 <- dat_subClusts2[dat_subClusts2$cluster == 4, ]

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
st_write(f3_grp1_subclust1.geodata.sf, "data_output/f3_grp1_subclust1.shp", delete_dsn = T)
st_write(f3_grp1_subclust2.geodata.sf, "data_output/f3_grp1_subclust2.shp", delete_dsn = T)

st_write(f3_grp2_subclust1.geodata.sf, "data_output/f3_grp2_subclust1.shp", delete_dsn = T)
st_write(f3_grp2_subclust2.geodata.sf, "data_output/f3_grp2_subclust2.shp", delete_dsn = T)
st_write(f3_grp2_subclust3.geodata.sf, "data_output/f3_grp2_subclust3.shp", delete_dsn = T)
st_write(f3_grp2_subclust4.geodata.sf, "data_output/f3_grp2_subclust4.shp", delete_dsn = T)


# Note to self:
# Andrew was concerned about using the component loadings, or coefficient correlations, from the PCA to directly interpret the clusters. Have a day and a half of reseraching and googling this issue, I determined that it is totally fine to make a direct unterpretation of our kmeans cluster plots using the loadings for the PCA.
# Reason: The kmeans analysis only creates a massive distance matrix and calculates clusters based on the position of the centroids, and the number of centroids are determined by you (when you set the number of groups for the cluster). There is now way to plot this directly, as there are too many dimensions.You mnust then, as a step in the plotting process, prerform a PCA to reduce those dimensions down. If you looking at how the function fviz_cluster is working ( getAnywhere(fviz_cluster)) you will see that it calculates a PCA and what you are actually plotting are the PC1 and PC2 scores from those results. Below is also a link that basically says this same thing:
# https://www.datanovia.com/en/lessons/k-means-clustering-in-r-algorith-and-practical-examples/#k-means-algorithm
# Som to validate this is, as in, make sure that the PCA I calulated by hand is the same data as my cluster plot, I can look st the PC scores that I calculated and compare them with the data output from the fviz_cluster results for both the whole river analysis and the subcluster analysis:
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
