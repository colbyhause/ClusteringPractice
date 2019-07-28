# clustering practice/ tutorial from this site: wasnt super helpful----
# https://www.r-bloggers.com/k-means-clustering-in-r/


library(ggplot2)

# we see from the data above that petal width and length is similar among species but varied btwn species :
ggplot(iris, aes(Petal.Length, Petal.Width, color = Species)) + geom_point()

# K-means clustering with 3 clusters of sizes 46, 54, 50


set.seed(20) 
irisCluster <- kmeans(iris[, 3:4], 3, nstart = 20) # since there are 3 species, tell R to make 3 groups, nstart means R will try 20 different random starting assignments and then select the one with the lowest within cluster variation 
irisCluster

table(irisCluster$cluster, iris$Species)

irisCluster$cluster <- as.factor(irisCluster$cluster)
ggplot(iris, aes(Petal.Length, Petal.Width, color = iris$cluster)) + geom_point()


mydata

#another practice set: ----

# practive from this site:
# https://www.guru99.com/r-k-means-clustering.html

library(dplyr)
library(tidyverse)
PATH <-"https://raw.githubusercontent.com/thomaspernet/data_csv_r/master/data/Computers.csv"
df <- read.csv(PATH) %>% 
  select(-c(X, cd, multi, premium))
glimpse(df)


# rescale:
rescale_df <- df %>%
  mutate(price_scal = scale(price),
         hd_scal = scale(hd),
         ram_scal = scale(ram),
         screen_scal = scale(screen),
         ads_scal = scale(ads),
         trend_scal = scale(trend)) %>%
  select(-c(price, speed, hd, ram, screen, ads, trend))


install.packages("animation")	


# example running the algorithm on only rescaled variables hd and ram with three clusters:
set.seed(2345)
library(animation)
kmeans.ani(rescale_df[2:3], 3)


pc_cluster <-kmeans(rescale_df, 5)


# determing optimal K- which means the number of groups you wantto cluster:
# create the function that runs the k-mean algorithm and store the total within clusters sum of squares:

kmean_withinss <- function(k) {
  cluster <- kmeans(rescale_df, k) # run the algorithm k times
  return (cluster$tot.withinss) # Store the total within clusters sum of squares

}

# trying just 2 times:
kmean_withinss(2)

# use the sapply() function to run the algorithm over a range of k. 
# This technique is faster than creating a loop and store the value:
# Set maximum cluster 
max_k <-20
# Run algorithm over a range of k 
wss <- sapply(2:max_k, kmean_withinss) # run function over range of 2- 20
wss


# create a data frame with the results of the algorithm
# # Create a data frame to plot the graph
elbow <-data.frame(2:max_k, wss)

# Plot the graph with gglop
ggplot(elbow, aes(x = X2.max_k, y = wss)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1, max_k, by = 1))

# From the graph, you can see the optimal k is seven, where the curve is starting to have a diminishing return.
# Once you have our optimal k, you re-run the algorithm with k equals to 7 and evaluate the clusters.

# Now examine the cluster
pc_cluster_2 <-kmeans(rescale_df, 7)
pc_cluster_2$cluster
pc_cluster_2$centers
pc_cluster_2$size 


#look at clusters in a heat map

install.packages("RColorBrewer")
library(RColorBrewer)

center <-pc_cluster_2$centers 
center # rows refer to numeration of the cluster and the columns are the variables used by the algorithm.
#The values are the average score by each cluster for the interested column. Standardization makes the interpretation easier. Positive values indicate the z-score for a given cluster is above the overall mean. For instance, cluster 2 has the highest price average among all the clusters.

#1. Build a dataframe:
#create dataset with the cluster number
cluster <- c(1: 7)
center_df <- data.frame(cluster, center)

#2. Reshape the data
center_reshape <- gather(center_df, features, values, price_scal: trend_scal)# this just reshapes the data 
head(center_reshape)

# Create the palette
hm.palette <-colorRampPalette(rev(brewer.pal(10, 'RdYlGn')),space='Lab')

#3. Plot the heat map
ggplot(data = center_reshape, aes(x = features, y = cluster, fill = values)) +
  scale_y_continuous(breaks = seq(1, 7, by = 1)) +
  geom_tile() +
  coord_equal() +
  scale_fill_gradientn(colours = hm.palette(90)) +
  theme_classic()


# Cluster tutorial sent from Christine:----
# from site: https://uc-r.github.io/kmeans_clustering#kmeans

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization

# dataframe to work with:
df <- USArrests

# remove missing data
df <- na.omit(df)
which(is.na(df))

# scale the data to standardize it:
df <- scale(df)
head(df)

# calculate dissimilarity by getting distances:
class(df)
distance <- get_dist(df) # calculate distances
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) # visualize by creating a dustance matrix

# K-means clustering:
# K-means algorithm can be summarized as follows:
#1. Specify the number of clusters (K) to be created (by the analyst)
#2. Select randomly k objects from the data set as the initial cluster centers or means
#3. Assigns each observation to their closest centroid, based on the Euclidean distance between the object and the centroid
#4. For each of the k clusters update the cluster centroid by calculating the new mean values of all the data points in the cluster. The centroid of a Kth cluster is a vector of length p containing the means of all variables for the observations in the kth cluster; p is the number of variables.
#5.Iteratively minimize the total within sum of square (Eq. 7). That is, iterate steps 3 and 4 until the cluster assignments stop changing or the maximum number of iterations is reached. By default, the R software uses 10 as the default value for the maximum number of iterations.

k2 <- kmeans(df, centers = 2, nstart = 25) # nstart means the function will attempt that many initial configurations and report on the best one 
str(k2)

# View clustering results with fviz_cluster--> provides visualization:
# fviz_cluster will perform principal component analysis (PCA) and plot the data points according to the first two principal components that explain the majority of the variance.

fviz_cluster(k2, data = df)

#you can also use pair-wise scatter plot to see the clustering pattern:
df %>%
  as_tibble() %>%
  mutate(cluster = k2$cluster,
         state = row.names(USArrests)) %>%
  ggplot(aes(UrbanPop, Murder, color = factor(cluster), label = state)) +
  geom_text()


# better to try many different grouping (k) options, not just set on 2:
k3 <- kmeans(df, centers = 3, nstart = 25)
k4 <- kmeans(df, centers = 4, nstart = 25)
k5 <- kmeans(df, centers = 5, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = df) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = df) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = df) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = df) + ggtitle("k = 5")

library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)

# Method above ^ does not provide information on the ideal number of k groups,use one of these methods:
# Elbow method
# Silhouette method
# Gap statistic

#1. Elbow method: ----
# want to determine how many clusters are needed so that the total within-cluster sum of square (wss) (which measures the compactness of the clustering) is as small as possible and does not get much smaller with more groups
# do this by calculating the wss over a range of group options (ie. 1-10) and then plot those values. Where the curve bends is the appropriate number of groups 

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
fviz_nbclust(df, kmeans, method = "wss")

#2. Silhouette method:----
#This method measure the quality of the clustering. So a high average silhouette width indicates good clustering, so when plotting this method the group that is maximized is the optimal group number for clustering:

# the long way:
## function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(df, centers = k, nstart = 25)
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

# wrapped up into one function:
fviz_nbclust(df, kmeans, method = "silhouette") # so says 2 is the optimal number of groups, but 4 groups is a close second option

#3. Gap statistic: ----
#this approach can be applied to ANY clustering method (kmeans, hierarchical, etc.). The gap statistic compares the total intracluster variation with different values of K with their expected values under null reference distribution of the data (i.e. a distribution with no obvious clustering).
#SO: That is, for each variable  in the data set we compute its range and generate values for the n points uniformly from the interval min to max.

# compute gap statistic
set.seed(123)
gap_stat <- clusGap(df, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
print(gap_stat, method = "firstmax")

#We can visualize the results with fviz_gap_stat which suggests four clusters as the optimal number of clusters.

fviz_gap_stat(gap_stat) # shows that 4 clusters are the best option

# overall, the results of the 3 methods mostly point towards 4 clusters, so use that many in the analysis:

set.seed(123)
final <- kmeans(df, 4, nstart = 25)
print(final)

#visualize:
fviz_cluster(final, data = df)

# do some descriptive statistics based onthe clustering:
USArrests %>%
  mutate(Cluster = final$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")
