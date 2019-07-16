# Tue Jul 16 10:28:44 2019 ------------------------------
# Hierarchical Clustering practice from this website:
# http://cc.oulu.fi/~jarioksa/opetus/metodi/sessio3.pdf

# load data
library(vegan)
library()
data(dune)
 
# Calculate dissimilarity by getting distances btwn the data 
d<- vegdist(dune)

# define 2 panels side by side for the plots:
par(mfrow=c(1,2))

# single linkage clustering :
csin <- hclust(d, method="single")
csin 
plot(csin) # plot the dendrogram
plot(csin, hang=-1) # another way to look at the plot with the branches hanging down

# complete linkage clustering:
par(mfrow=c(2,2))

ccom <- hclust(d, method="complete")
plot(ccom, hang=-1)

# average linkage clustering 
caver <- hclust(d, method="aver")
plot(caver, hang=-1)


vegemite(dune, caver)


# UC Business Analytics H Clustering Tutorial ----
# From this website: https://uc-r.github.io/hc_clustering
# Hierarchical clustering is an alternative approach to k-means clustering for identifying groups in the dataset. It does not require us to pre-specify the number of clusters to be generated as is required by the k-means approach. Furthermore, hierarchical clustering has an added advantage over K-means clustering in that it results in an attractive tree-based representation of the observations, called a dendrogram.

# load data set, remeber that : 
#1. Rows are observations (individuals) and columns are variables
#2. Any missing value in the data must be removed or estimated.
#3. The data must be standardized (i.e., scaled) to make variables comparable. Recall that, standardization consists of transforming the variables such that they have mean zero and standard deviation one.[^scale]

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

df <- USArrests

# remove missing values:
df <- na.omit(df)

#scale data:
df <- scale(df)
head(df)

# Agglomerative Hierarchical Clustering (using hclust):

# first compute the dissimilarity values (make a Dissimilarity matrix):
d <- dist(df, method = "euclidean")

# Feed the disim matric into hclust using Complete Linkage:
hc1 <- hclust(d, method = "complete" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

# Alternatively, use the agnes fucntion, which is the same but allows you toget the agglomerative coefficient, which measures the amount of clustering structure found (values closer to 1 suggest strong clustering structure).
hc2 <- agnes(df, method = "complete")

# Agglomerative coefficient
hc2$ac
## [1] 0.8531583

# Use the Aglomerative coeff. to assess stringer clustering structures:
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x) {
  agnes(df, method = x)$ac
}

map_dbl(m, ac)
##   average    single  complete      ward 
## 0.7379371 0.6276128 0.8531583 0.9346210

# visualize dendrogram using ward method:
hc3 <- agnes(df, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 

# Divise Hierachical Clustering ----
# compute divisive hierarchical clustering
hc4 <- diana(df)

# Divise coefficient; amount of clustering structure found
hc4$dc
## [1] 0.8514345

# plot dendrogram
pltree(hc4, cex = 0.6, hang = -1, main = "Dendrogram of diana")


# notes on Dendrograms ----
# Each leaf corresponds to one observation. As we move up the tree, observations that are similar to each other are combined into branches, which are themselves fused at a higher height.
# The height of the fusion, provided on the vertical axis, indicates the (dis)similarity between two observations. 
# So: The higher the height of the fusion, the less similar the observations are.


# In order to identify sub-groups (i.e. clusters), we can cut the dendrogram with cutree:
# Ward's method
hc5 <- hclust(d, method = "ward.D2" )

# Cut tree into 4 groups 
sub_grp <- cutree(hc5, k = 4)

    # to use cuttree with agnes and Diana do the folowing:
    # Cut agnes() tree into 4 groups
    hc_a <- agnes(df, method = "ward")
    cutree(as.hclust(hc_a), k = 4)

    # Cut diana() tree into 4 groups
    hc_d <- diana(df)
    cutree(as.hclust(hc_d), k = 4)

# Number of members in each cluster
table(sub_grp)
## sub_grp
##  1  2  3  4 
##  7 12 19 12

# use cuttree output to add the cliuster number to each observation in our data
USArrests %>%
mutate(cluster = sub_grp) %>%
  head

#You can also draw the dendrogram with a border around the 4 clusters. The argument border is used to specify the border colors for the rectangles:
  
plot(hc5, cex = 0.6)
rect.hclust(hc5, k = 4, border = 2:5)  

# can also use the fviz_cluster function from the factoextra package to visualize the result in a scatter plot.- like in the kmeans tutorial


fviz_cluster(list(data = df, cluster = sub_grp))  

# Comparing Dendrograms: allows you to compare dendrograms from different analytical methods side by side  with their labels connected by lines
# Compute distance matrix
res.dist <- dist(df, method = "euclidean")

# Compute 2 hierarchical clusterings
hc1 <- hclust(res.dist, method = "complete")
hc2 <- hclust(res.dist, method = "ward.D2")

# Create two dendrograms
dend1 <- as.dendrogram (hc1)
dend2 <- as.dendrogram (hc2)

tanglegram(dend1, dend2)
# the quality of the alignment can be measured with entanglement. Entanglement is a measure between 1 (full entanglement) and 0 (no entanglement). A lower entanglement coefficient corresponds to a good alignment. The output of tanglegram can be customized using many other options as follow:
dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2,
           highlight_distinct_edges = FALSE, # Turn-off dashed lines
           common_subtrees_color_lines = FALSE, # Turn-off line colors
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2))
)

# Determining Optimal Clusterings in Hierachical Clustering:
# use the same methods as in kmeans:
#1. Elbow method:
fviz_nbclust(df, FUN = hcut, method = "wss")

#2. Average Sihlouette Method:
fviz_nbclust(df, FUN = hcut, method = "silhouette")

#Gap statistic method:
gap_stat <- clusGap(df, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)
