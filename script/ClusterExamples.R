# clustering practice/ tutorial from this site: wasnt super helpful...
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

