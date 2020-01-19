# The goal of this script is to plot the result of the hierarchical clustering method (UPGMA) spatially:
library(tidyverse)
library(sf)

# Read in files:
meansDat_upgma_8clusters <- read_csv("data_output/mean_fromPredicted_AllParams_across_actualRKMS_8clusts_upgma.csv")
CVdat_upgma_6clusters <- read_csv("data_output/CV_fromPredicted_AllParams_across_actualRKMS_NoOutliers_6clusts_upgma.csv")

# Run spatialize clusters function:
spatialize_clusters <- function(data, coords, filename) {
  spatial_data <- st_as_sf(data, coords = coords, #specify lat lon
                           remove = F, # don't remove these lat/lon cols from df
                           crs = 4326) # add projection
  # write spatial data to creat shp file----
  st_write(spatial_data, paste0("data_output/SpatialData/",filename,".shp"), delete_dsn = T)
}

# Spatilaize data:----
# MEANS----
meansDat_upgma_8clusters_spatial <- spatialize_clusters(meansDat_upgma_8clusters, c("geomean_lon", "geomean_lat"), "HierarchicalClustering/Means/meansDat_upgma_8clusters_spatial")
#Cluster 1 has only 12 points in it- need to determine where to keep or get rid of these points...

#CV----
CVdat_upgma_6clusters_spatial <- spatialize_clusters(CVdat_upgma_6clusters, c("geomean_lon", "geomean_lat"), "HierarchicalClustering/CV/CVdat_upgma_6clusters_spatial")
# This also has a clusters with only 6 points for clust 6, should probably get ride of these

# After removing the cluster 6 in the CV data---
cv_dat_edited3_upgma_5Clusters<- read_csv("data_output/CV_fromPredicted_AllParams_across_actualRKMS_edited3_5clusts_upgma.csv")

CVdat_upgma_5clusters_spatial <- spatialize_clusters(cv_dat_edited3_upgma_5Clusters, c("geomean_lon", "geomean_lat"), "HierarchicalClustering/CV/CVdat_upgma_5clusters_spatial")

# to discuss w/ Andrew:
# I looked at 4 different hierarchical clustering methods:
# complete, single, average (UPGMA), ward
# Following guideline from the Num Ecology with R book, I used 2 methods to choose the one of the methods above to cluster data
# For both means and CV data, UPGMA was best method
# I used silhouette method to determine the number of clusters, and check using the mantel coefficient. Means: 8, CV: 6
# I had to make edits to the CV data: first 2 round of edits were bc only a couple data points were creating a 3rd arm, removed those so there were 2 main brainches. Then one of the 6 clusters only had 6 records in it, so i removed those and did the analysis over again. Gave me 5 clusters this time, see qgis.
# see Qgis porject for the clusters shown spatially. 
# ASK: is log transforming and scaling data necessary for CV? why do we need to log transform it period?


