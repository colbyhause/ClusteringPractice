# Tue Jul 16 11:34:35 2019 ------------------------------
# First attempt at Hierarchical Clustering Analysis on Flame 3 data
# following the UC Davis Business Analytics tutorial ( worked through with examnple data in the HclusterExamples script)

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

# load and prep data:
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

# first compute the dissimilarity values (make a Dissimilarity matrix):
d <- dist(dat_rescaled, method = "euclidean")

# Feed the disim matric into hclust using Complete Linkage:
hc1 <- hclust(d, method = "complete" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

# Alternatively, use the agnes fucntion, which is the same but allows you toget the agglomerative coefficient, which measures the amount of clustering structure found (values closer to 1 suggest strong clustering structure).
hc2 <- agnes(dat_rescaled, method = "complete")

# Agglomerative coefficient
hc2$ac  # warning: this takes a long time to run!!!!
## [1] 0.9972277

# Use the Aglomerative coeff. to assess stringer clustering structures:
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x) {
  agnes(dat_rescaled, method = x)$ac
}

map_dbl(m, ac)
##   average    single  complete      ward 
## 0.7379371 0.6276128 0.8531583 0.9346210

# visualize dendrogram using ward method:
hc3 <- agnes(df, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 

# Divise Hierachical Clustering ----
# compute divisive hierarchical clustering
hc4 <- diana(dat_rescaled)

# Divise coefficient; amount of clustering structure found
hc4$dc
## [1] 0.8514345

# plot dendrogram
pltree(hc4, cex = 0.6, hang = -1, main = "Dendrogram of diana")

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

# use cuttree output to add the cliuster number to each observation in our data
test <- dat_noNAs %>%
  mutate(cluster = sub_grp) 

#You can also draw the dendrogram with a border around the 4 clusters. The argument border is used to specify the border colors for the rectangles:

plot(hc5, cex = 0.6, labels = F)
rect.hclust(hc5, k = 20, border = 2:5)

# can also use the fviz_cluster function from the factoextra package to visualize the result in a scatter plot.- like in the kmeans tutorial


fviz_cluster(list(data = dat_rescaled, cluster = sub_grp), labelsize = NA)  

# Comparing Dendrograms: allows you to compare dendrograms from different analytical methods side by side  with their labels connected by lines
# Compute distance matrix
res.dist <- dist(dat_rescaled, method = "euclidean")

# Compute hierarchical clusterings
hc1 <- hclust(res.dist, method = "complete")
hc2 <- hclust(res.dist, method = "ward.D2")
hc3<- agnes(dat_rescaled, method = "average")
hc4<- agnes(dat_rescaled, method = "single")

# Create two dendrograms
dend1 <- as.dendrogram (hc1)
dend2 <- as.dendrogram (hc2)
dend3 <- as.dendrogram (hc3)
dend4 <- as.dendrogram (hc4)

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
fviz_nbclust(dat_rescaled, FUN = hcut, method = "wss")

#2. Average Sihlouette Method:
fviz_nbclust(dat_rescaled, FUN = hcut, method = "silhouette")

#Gap statistic method:
gap_stat <- clusGap(dat_rescaled, FUN = hcut, nstart = 25, K.max = 20, B = 50)
fviz_gap_stat(gap_stat) # says optimal is 20 clusters ???



