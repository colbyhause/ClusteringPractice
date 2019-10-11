###Trying our Density-based clustering:

#install.packages("dbscan")
library(fpc)
library(dbscan)

data_raw <- readRDS("data/flame2.supercleaned.rds")
class(data_raw)

which(is.na(data_raw)) # determine where there are NAs
dat_noNAs <- na.omit(data_raw) # remove NAs
which(is.na(dat_noNAs))

dat_cols <- dat_noNAs[ , c("temp", "spCond", "pH", "ODO", "turb",
                           "fDOM", "CHL", "NO3")] # keep only the columns of the parameters 

dat_transformed <- apply(dat_cols, 2, log) # transform the data
dat_rescaled <- apply(dat_transformed, 2, scale) # rescale the data

data_matrix<- as.matrix(dat_rescaled)
class(data_matrix)
kNNdistplot(data_matrix, k=15)
abline(h=0.4, col="red")

###### using dbscan from fpc package to determine number of clusters
df <- dat_rescaled

# Compute DBSCAN using fpc package: ISSUE: you have to set the parameters Minpts and Esp. MinPts should be you data dimension*2, but there is now good advice on how to set esp
library("fpc")
set.seed(123)
db <- fpc::dbscan(df, eps = .5, MinPts = 14)

# Plot DBSCAN results
library("factoextra")
fviz_cluster(db, data = df, stand = FALSE,
             ellipse = FALSE, show.clust.cent = FALSE,
             geom = "point",palette = "jco", ggtheme = theme_classic())





