# first attempt at cluster analysis with own data:

if(!require(psych)){install.packages("car")}
if(!require(MASS)){install.packages("MASS")}
if(!require(rcompanion)){install.packages("rcompanion")}

library(rcompanion)
library(psych)
library(MASS)

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

# first try clustering on just turbidity data :
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
dat_cols <- dat_noNAs[ , c(6:8, 10:15)]

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


# since 2 methods gave 2 optiminal groups and the other gave 12 optimal groups, try running k means with both:

set.seed(123)
final_2 <- kmeans(dat_rescaled, 2, nstart = 25)
print(final_2)

final_12 <- kmeans(dat_rescaled, 12, nstart = 25)
print(final_12)

final_4 <- kmeans(dat_rescaled, 4, nstart = 25)
print(final_12)
#visualize:
fviz_cluster(final_2, data = dat_rescaled, labelsize = NA)
fviz_cluster(final_12, data = dat_rescaled,  labelsize = NA)


# trying to determine number of clusters with NbClust package: this package provides 30 indices for determining the relevant number of clusters and proposes to users the best clustering scheme from the different results obtained by varying all combinations of number of clusters, distance measures, and clustering methods.
#install.packages("NbClust")
library(NbClust)
es.nbclust <- NbClust(dat_rescaled, distance = "euclidean",
                      min.nc = 2, max.nc = 15, 
                      method = "complete", index ="all")
factoextra::fviz_nbclust(es.nbclust) + theme_minimal() + ggtitle("NbClust's optimal number of clusters")

# comes out to optimal number of clusters is 2


