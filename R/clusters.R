install.packages(c("compositions", "cluster", "NbClust", "factoextra"))
library(compositions)
library(cluster)
library(NbClust)
library(factoextra)

# Exclude the 'position..mm.' column
data_clr200 <- df.clr$raw.200[, !names(df.clr) %in% "position..mm."]

# Determine the optimal number of clusters

# Elbow method
fviz_nbclust(data_clr200, kmeans, method = "wss")

# Silhouette method
fviz_nbclust(data_clr200, kmeans, method = "silhouette")

# Gap statistic
set.seed(123)
gap_stat <- clusGap(data_clr200, FUN = kmeans, K.max = 10, B = 500)
fviz_gap_stat(gap_stat)

# Alternatively, use NbClust
nb <- NbClust(data_clr200, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans", diss=NULL)
fviz_nbclust(nb)

# Perform clustering with the chosen number of clusters
set.seed(123)  # Ensure reproducibility
k <- 3  # Replace with the chosen number of clusters
kmeans_result <- kmeans(data_clr200, centers = k, nstart = 25)

# Visualize clusters
fviz_cluster(kmeans_result, data = data_clr200)
