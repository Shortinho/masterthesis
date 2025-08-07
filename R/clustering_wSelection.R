# clustering 
# needs to have selected vars in global environment (var_selection scripts)

source('xrf_functions.R')
library(FactoMineR)
library(rioja)
library(scales)
library(ggplot2)
library(forcats)

set.seed(123)

selected_vars <- c(
  "Mn_Fe", 'Ba', 'Ti', 'K', 'Si', "S_Ti", "Fe_Ti", 'TCWt', 
  "Ca_Ti", "Ti_Al"
)
# combined contains all the variables in one data frame
combined <- combined_df_creator()

df_clust <- combined %>%
  select(any_of(c("position..mm.", 'facies_class', selected_vars))) %>%
  drop_na()

df_clust <- read.csv(file = 'data/generated/multivar/multivar_10', ) %>%
  select(-X)

combined <- read.csv('data/generated/combined/combined_AC_200_inc_coh.csv')
xrf_elem_sel <- c('Fe.y','Mn.y','Ba.y','S.y','Ti.y','Al.y','Si.y','K.y', 'Ca.y', 'Rb.y', 'Zr.y', 'Sr.y', 'P.y', 'Ba.y','position..mm.')

df_clust <- combined %>% 
  select(any_of(c(xrf_elem_sel, 'facies_class')))

# ---- fviz method ------
X <- scale(select(df_clust, -c(position..mm., facies_class)))

# Choose optimal number of clusters (here 2 is optimal)
fviz_nbclust(X, kmeans, method = 'silhouette') +
  theme_minimal() +
  labs(title = 'Optimal number of clusters \n 500µm - All Core')


# ---- Perform k-means -----
km <- kmeans(X, centers = 2, nstart = 25)

df_clust$kmeans_cluster <- as.factor(km$cluster)

# Visualize clusters (PCA plot)
fviz_cluster(km, data = X, geom = "point", ellipse.type = "norm",
             main = "K-means Clustering")
# Compare to facies plot on PCA
pca_res <- prcomp(X, center = TRUE, scale. = TRUE)

# Step 2: Extract PCA scores and add facies and kmeans info
pca_scores <- as.data.frame(pca_res$x)
pca_scores$facies_class <- as.factor(df_clust$facies_class)
pca_scores$kmeans <- as.factor(df_clust$kmeans_cluster)

# Step 3: Plot
ggplot(pca_scores, aes(x = PC1, y = PC2, color = facies_class)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA: Redox/Productivity Variables Colored by Facies Class",
       x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "%)"),
       y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "%)"),
       color = "Facies") +
  theme_minimal() +
  scale_color_manual(values = c("#1b9e77", "#d95f02")) + 
  stat_ellipse(type = "norm", level = 0.68)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = kmeans)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA: Redox/Productivity Variables Colored by K-means cluster",
       x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "%)"),
       y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "%)"),
       color = "Facies") +
  theme_minimal() +
  scale_color_manual(values = c("#1b9e77", "#d95f02")) + 
  stat_ellipse(type = "norm", level = 0.68)

# ward's tree method
dist_mat <- dist(X)
hc <- hclust(dist_mat, method = "ward.D2")

plot(hc, labels = FALSE, hang = -1)
rect.hclust(hc, k = 2, border = "red")

df_clust$hclust_cluster <- as.factor(cutree(hc, k = 2))

pca_scores$h_clusters <- as.factor(df_clust$hclust_cluster)

# plot ward's tree method on PCA plot
ggplot(pca_scores, aes(x = PC1, y = PC2, color = h_clusters)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA: Redox/Productivity Variables Colored by hclust",
       x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "%)"),
       y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "%)"),
       color = "Facies") +
  theme_minimal() +
  scale_color_manual(values = c("#1b9e77", "#d95f02")) + 
  stat_ellipse(type = "norm", level = 0.68)

# comparing clustering to facies 
table(Kmeans = df_clust$kmeans_cluster, Facies = combined$facies_class)
table(Hierarchical = df_clust$hclust_cluster, Facies = combined$facies_class)

# Adjusted Rand Index (needs mclust)
library(mclust)
adjustedRandIndex(df_clust$kmeans_cluster, combined$facies_class)


# supervised

# using random forest
# Prepare data
df_rf <- df_clust %>%
  select(facies_class, Mn_Fe, Ba, Ti, K, Si, S_Ti, Fe_Ti, TCWt, Ca_Ti, Ti_Al) %>%
  na.omit()

df_rf$facies_class <- as.factor(df_rf$facies_class)

# Fit Random Forest model
set.seed(123)
rf_model <- randomForest(facies_class ~ ., data = df_rf, importance = TRUE, ntree = 500)

# Print results
print(rf_model)

# Plot variable importance
varImpPlot(rf_model)

df_rf$predicted_facies <- predict(rf_model)
# Add PCA scores and predictions
pca_df <- as.data.frame(pca_res$x)
pca_df$Predicted <- df_rf$predicted_facies
pca_df$True <- df_rf$facies_class

ggplot(pca_df, aes(x = PC1, y = PC2, color = Predicted)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA Colored by Random Forest Predicted Facies")

ggplot(pca_df, aes(x = PC1, y = PC2, color = Predicted)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA: Redox/Productivity Variables Colored by RF predicted",
       x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "%)"),
       y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "%)"),
       color = "Facies") +
  theme_minimal() +
  scale_color_manual(values = c("#1b9e77", "#d95f02")) + 
  stat_ellipse(type = "norm", level = 0.68)
# visualize

ggplot(df_clust, aes(x = position..mm., y = kmeans_cluster)) +
  geom_point() +
  scale_y_discrete(name = "Cluster") +
  coord_flip() +
  labs(title = "K-means Clusters by Depth (mm)") +
  theme_minimal()

# create bar plot for k-means clusters
cluster_bar <- df_clust %>%
  arrange(position..mm.) %>%
  mutate(depth_top = position..mm.,
         depth_bottom = lead(position..mm.),
         group = kmeans_cluster) %>%
  filter(!is.na(depth_bottom))

# Optional: choose color palette
col_pal_cluster <- RColorBrewer::brewer.pal(5, "Set2")  # Adjust number to your cluster count

p_cluster <- ggplot(cluster_bar) +
  geom_rect(aes(ymin = depth_top, ymax = depth_bottom,
                xmin = 0, xmax = 1, fill = factor(group))) +
  scale_fill_manual(values = c('#af8dc3','#f7f7f7','#7fbf7b')) +
  scale_y_reverse(name = "Depth (mm)") +  # Add depth axis label
  scale_x_continuous(expand = c(0, 0)) +  # Fix spacing on x-axis
  theme_minimal() +  # Use minimal theme instead of void
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    aspect.ratio = 5
  ) +
  labs(title = "Cluster Classification Downcore \n 200µm All Core - KMeans", fill = "Cluster")

print(p_cluster)

output_dir <- 'plots/PCdowncore/'
filename <- file.path(output_dir, "cluster_200_AC.svg")
save_svg_plot(p_cluster, filename, width = 3, height = 8)



# ---- CONISS -----

diss <- df_clust %>%
  select(-c(position..mm., facies_class, kmeans_cluster)) %>%
  scale() %>%
  dist()  
clust <- chclust(diss, method = "coniss")

# Plot the constrained dendrogram
plot(clust)


# Optional: use broken stick to choose number of clusters
bstick(clust, 15)


# Step 1: Assign clusters (2 in this case)
clusters <- cutree(clust, k = 7)

# Step 2: Add cluster labels and calculate top/bottom depths
cluster_bar <- df_clust
cluster_bar$group <- clusters

# Rename for easier handling
cluster_bar$position_mm <- cluster_bar$`position..mm.`

# Sort by depth
cluster_bar <- cluster_bar[order(cluster_bar$position_mm), ]

# Compute sample interval (assumes roughly constant spacing)
interval <- median(diff(cluster_bar$position_mm))

# Calculate top and bottom positions around the center point
cluster_bar$depth_top <- cluster_bar$position_mm - interval / 2
cluster_bar$depth_bottom <- cluster_bar$position_mm + interval / 2

# Step 3: Plot
library(ggplot2)

p_cluster <- ggplot(cluster_bar) +
  geom_rect(aes(ymin = depth_top, ymax = depth_bottom,
                xmin = 0, xmax = 1, fill = factor(group))) +
  #scale_fill_manual(values = c('#8856a7','#b3cde3')) +
  scale_y_reverse(name = "Depth (mm)") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    aspect.ratio = 5
  ) +
  labs(title = "Cluster Classification Downcore \n 200µm CONISS 7 clusters", fill = "Cluster")

p_cluster


# ---- CONISS but 3 clusters -----

clusters_3 <- cutree(clust, k = 3)
cluster_bar_3 <- df_clust
cluster_bar_3$group <- clusters_3
cluster_bar_3$position_mm <- cluster_bar_3$`position..mm.`

# Ensure samples are ordered
cluster_bar_3 <- cluster_bar_3[order(cluster_bar_3$position_mm), ]

# Compute interval
interval <- median(diff(cluster_bar_3$position_mm))

# Calculate top and bottom
cluster_bar_3$depth_top <- cluster_bar_3$position_mm - interval / 2
cluster_bar_3$depth_bottom <- cluster_bar_3$position_mm + interval / 2

colors_3 <- c('#e41a1c', '#377eb8', '#4daf4a')

p_cluster_3 <- ggplot(cluster_bar_3) +
  geom_rect(aes(ymin = depth_top, ymax = depth_bottom,
                xmin = 0, xmax = 1, fill = factor(group))) +
  scale_fill_manual(values = colors_3) +
  scale_y_reverse(name = "Depth (mm)") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    aspect.ratio = 5
  ) +
  labs(title = "3-Cluster Classification Downcore", fill = "Cluster")

p_cluster_3



#----- Zander et al style cluster centers ----- 

# Combine cluster labels with X (scaled variables used in clustering)
X_df <- as.data.frame(X) %>%
  mutate(kmeans_cluster = factor(km$cluster))

# Compute cluster centers (means of each variable per cluster)
cluster_centers <- X_df %>%
  group_by(kmeans_cluster) %>%
  summarise(across(.cols = everything(), .fns = mean), .groups = "drop")

# Reshape to long format for plotting
centers_long <- cluster_centers %>%
  pivot_longer(cols = -kmeans_cluster, names_to = "variable", values_to = "center_value")

# Determine global x-axis range for comparability
x_limits <- range(centers_long$center_value)

# Plot for cluster 1
p1 <- ggplot(filter(centers_long, kmeans_cluster == 1),
             aes(x = center_value, y = fct_rev(variable))) +
  geom_col(fill = "#8856a7") +
  scale_x_continuous(limits = x_limits) +
  labs(title = "Cluster 1", x = "Cluster Center (Scaled)", y = "Variable") +
  theme_minimal()

# Plot for cluster 2
p2 <- ggplot(filter(centers_long, kmeans_cluster == 2),
             aes(x = center_value, y = fct_rev(variable))) +
  geom_col(fill = "#b3cde3") +
  scale_x_continuous(limits = x_limits) +
  labs(title = "Cluster 2", x = "Cluster Center (Scaled)", y = "Variable") +
  theme_minimal()

# Display side by side
p1 + p2



