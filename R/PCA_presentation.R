# clean PCA with selection for presentation

library(ggplot2)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(tibble)
library(ggrepel)  # for optional smart labels
library(svglite)
library(scales)


source('xrf_functions.R')

combined <- read.csv('/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/generated/combined/combined_200_fig1.csv')

xrf_elem_sel <- c("Fe","Mn",'Ba',"S","Ti","Al","Si","K", 'Ca', 'Rb', 'Zr', 'Sr', 'P', 'Ba','Ni','Zn','position..mm.')
non_xrf_sel <- c('position..mm.', 'TCWt')
all_sel <- append(xrf_elem_sel, non_xrf_sel)

res <- 0.5
half_res <- res/2
df.raw <- load_raw_data()
df.xrf <- clean_df(df.raw, sel = xrf_elem_sel)
df.xrf$raw.200 <- NULL
df.xrf <- tibble(df.xrf$raw.500) %>%
  select_top_varves()

df.cs <- closed_sum(df.xrf) 
df.clr <- clr_transform(df.xrf)



# ---- Compute ratios on closed-sum data ----
df.clr <- df.cs %>%
  compute_ratios('Fe', 'Ti') %>%
  compute_ratios('Ca', 'Ti') %>%
  compute_ratios('S', 'Ti') %>%
  compute_ratios('Fe', 'Mn') %>%
  compute_ratios('Ti', 'Al')



# Match TCWt from combined using position
df.clr <- df.clr %>%
  left_join(combined %>% select(position..mm., TCWt),
            by = c("position..mm." = "position..mm."))
# ---- Select final variables for PCA ----
df.pca <- df.clr %>%
  select(
    Fe_Ti, Ca_Ti, S_Ti, Fe_Mn, Ti_Al,  # ratios
    Si, Ba, K, Ti,          # clr-transformed elements
    TCWt
  ) %>%
  na.omit()
df.multivar <- df.pca %>%
  mutate(position..mm. = df.clr$position..mm.)

df.multivar.clust <- df.multivar %>%
  mutate(kmeans_cluster = df_clust$kmeans_cluster) %>%
  mutate(ward_hcluster = df_clust$hclust_cluster) %>%
  mutate(facies = combined$facies_class)
  

write.csv(df.multivar, file = 'data/generated/multivar/multivar_10')


# ---- Run PCA ----
pca_result <- PCA(df.pca, scale.unit = TRUE, graph = FALSE)

# ---- Plot PCA with fviz ---- ----
fviz_pca_biplot(pca_result,
                geom.ind = "point",       # show points only (no text)
                pointshape = 21,
                pointsize = 2,
                fill.ind = "gray70",      # light gray for sample points
                col.ind = "black",
                col.var = "black",        # variables in black
                arrowsize = 0.8,
                repel = TRUE) +           # repel labels for variables only
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ggtitle("PCA of Geochemical Proxies - 500µm dataset") +
  labs(x = paste0("PC1 (", round(pca_result$eig[1, 2], 1), "%)"),
       y = paste0("PC2 (", round(pca_result$eig[2, 2], 1), "%)"))






# ---- Run PCA with prcomp() ----
# Assumes df.pca is already preprocessed and standardized as needed

pca_result <- prcomp(df.pca, scale. = TRUE, center = TRUE)
pca_result <- fix_pca_signs(pca_result)
# ---- Extract scores and loadings ----
scores <- as_tibble(pca_result$x)
loadings <- as_tibble(pca_result$rotation, rownames = "variable")

# Rescale loadings to match sample space
loadings_scaled <- loadings %>%
  mutate(
    PC1 = PC1 * pca_result$sdev[1],
    PC2 = PC2 * pca_result$sdev[2]
  )

#  Add metadata for faceting
scores <- scores %>%
  mutate(depth = df.clr$position..mm.) %>%
  mutate(facies = combined$facies_class) %>%
  mutate(kcluster = df.multivar.clust$kmeans_cluster) %>%
  mutate(hcluster = df.multivar.clust$ward_hcluster)

var_expl <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# ---- plot pca biplot with ggplot----
# Assuming loadings_scaled has already been created by scaling with sdev

p.pca <- ggplot() +
  # Sample points colored by depth
  geom_point(data = scores, aes(x = PC1, y = PC2, colour = depth),
             shape = 16, size = 2.5) +
  
  # Loadings arrows
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.6) +
  
  # Loadings labels
  geom_text_repel(data = loadings_scaled,
                  aes(x = PC1 * 3.2, y = PC2 * 3.2, label = variable),
                  size = 6, color = "black") +
  
  # Center lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  
  # Axes and theme
  xlab(paste0("PC1 (", var_expl[1], "%)")) +
  ylab(paste0("PC2 (", var_expl[2], "%)")) +
  coord_equal() +
  scale_colour_gradient(
    low = "#c6dbef",
    high = "#08306b",
    name = "Depth (mm)",
    guide = guide_colorbar(reverse = TRUE)  # must go here
  ) +  # blue scale
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank()
  ) +
  ggtitle("PCA Biplot Colored by Depth")

print(p.pca)
output_dir <- 'plots/fig_pres/'
filename <- file.path(output_dir, "PCA_depth.svg")
save_svg_plot(p.pca, filename, width = 15, height = 10)


# ---- pca biplot with color by cluster -----
p.pca.kclust <- ggplot() +
  # Sample points colored by depth
  geom_point(data = scores, aes(x = PC1, y = PC2, color = kcluster),
             shape = 16, size = 2.5) +
  
  # ellipse over points
  stat_ellipse(data = scores, aes(x = PC1, y = PC2, color = kcluster), 
               type = "norm", level = 0.68) +
  # Loadings arrows
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.6) +
  
  # Loadings labels
  geom_text_repel(data = loadings_scaled,
                  aes(x = PC1 * 3.2, y = PC2 * 3.2, label = variable),
                  size = 5, color = "black") +
  
  # Center lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  
  # Axes and theme
  xlab(paste0("PC1 (", var_expl[1], "%)")) +
  ylab(paste0("PC2 (", var_expl[2], "%)")) +
  coord_equal() +
  scale_colour_manual(values = c("#b3cde3", "#8856a7")) + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank()
  ) +
  ggtitle("PCA Biplot Colored by K-means Cluster") 

print(p.pca.kclust)



# ---- downcore plotting ---- 

#prepare data for plotting function
pc.dc <- perform_pca2(df.multivar)
pc.dc <- fix_pca_signs(pc.dc)
scores_df <- as.data.frame(pc.dc$x)
scores_df$position <- pc.dc$position..mm.
pca_downcore_plot(pc.dc, pc_num = 'PC1', df_name = 'multivar')


pc_num <- 'PC2'
scores_df[[pc_num]] <- pc.dc$x[, pc_num]

# plot
p.pc2 <- ggplot(scores_df, aes(x = position, y = .data[[pc_num]])) +
  # PCA line
  geom_line(color = "black", linewidth = 0.3) +
  
  # Positive ribbon (blue)
  geom_ribbon(aes(
    ymin = 0,
    ymax = ifelse(.data[[pc_num]] > 0, .data[[pc_num]], 0)
  ),
  fill = "#0072B2", alpha = 0.6) +  # Colorblind-safe blue
  
  # Negative ribbon (red)
  geom_ribbon(aes(
    ymin = ifelse(.data[[pc_num]] < 0, .data[[pc_num]], 0),
    ymax = 0
  ),
  fill = "#D55E00", alpha = 0.6) +  # Colorblind-safe red
  
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.2) +
  
  # Labels and scales
  labs(
    title = paste0(pc_num, " Downcore Profile"),
    x = "Depth (mm)",
    y = paste0(pc_num, " Score")
  ) +
  scale_x_reverse(expand = expansion(mult = c(0.01, 0.01))) +
  
  # Theme and layout
  coord_flip() +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 0),
    axis.text.y = element_text(angle = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(5, 5, 5, 5),
    plot.background = element_blank(),
    aspect.ratio = 3
  )


print(p.pc2)
output_dir <- 'plots/PCdowncore/'
filename <- file.path(output_dir, "PC2_downcore.svg")
save_svg_plot(p.pc2, filename, width = 3, height = 8)

ggsave("plots/PCdowncore/PC2_downcore_profile.png", plot = p, width = 80, height = 120, units = "mm", dpi = 600)



# ---- PCA using only elements, no ratios, 500µm ----

combined.elem <- read.csv('/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/generated/combined/combined_500_inc_coh.csv')

xrf_elem_sel <- c("Fe","Mn",'Ba',"S","Ti","Al","Si","K", 'Ca', 'Rb', 'Zr', 'Sr', 'P', 'Ba','position..mm.')
non_xrf_sel <- c('position..mm.', 'TCWt')
all_sel <- append(xrf_elem_sel, non_xrf_sel)

res <- 0.5
half_res <- res/2
df.raw <- load_raw_data()
df.xrf <- clean_df(df.raw, sel = xrf_elem_sel)
df.xrf$raw.200 <- NULL
df.xrf <- tibble(df.xrf$raw.500) %>%
  select_top_varves()

df.cs <- closed_sum(df.xrf) 
df.clr <- clr_transform(df.xrf)
df.pca <- df.clr %>%
  select(-position..mm.)

pca_result <- prcomp(df.pca, scale. = TRUE, center = TRUE)
pca_result <- fix_pca_signs(pca_result)

# prepare cluster data
X <- scale(df.pca)

# Choose optimal number of clusters (here 2 is optimal)
fviz_nbclust(X, kmeans, method = "silhouette") +
  theme_minimal()


# ---- Perform k-means -----
km <- kmeans(X, centers = 2, nstart = 25)

combined.elem$kmeans_cluster <- as.factor(km$cluster)

# ---- Extract scores and loadings ----
scores <- as_tibble(pca_result$x)
loadings <- as_tibble(pca_result$rotation, rownames = "variable")

# Rescale loadings to match sample space
loadings_scaled <- loadings %>%
  mutate(
    PC1 = PC1 * pca_result$sdev[1],
    PC2 = PC2 * pca_result$sdev[2],
    PC3 = PC3 * pca_result$sdev[3]
  )

#  Add metadata for faceting
scores <- scores %>%
  mutate(depth = combined.elem$position..mm.) %>%
  mutate(facies = combined.elem$facies_class) %>%
  mutate(kcluster = combined.elem$kmeans_cluster)

var_expl <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# ---- PLOTTING ----
# ---- plotting by depth ----

# by depth
p.pca <- ggplot() +
  # Sample points colored by depth
  geom_point(data = scores, aes(x = PC1, y = PC2, colour = depth),
             shape = 16, size = 2.5) +
  
  # Loadings arrows
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.6) +
  
  # Loadings labels
  geom_text_repel(data = loadings_scaled,
                  aes(x = PC1 * 3.2, y = PC2 * 3.2, label = variable),
                  size = 6, color = "black") +
  
  # Center lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  
  # Axes and theme
  xlab(paste0("PC1 (", var_expl[1], "%)")) +
  ylab(paste0("PC2 (", var_expl[2], "%)")) +
  coord_equal() +
  scale_colour_gradient(
    low = "#c6dbef",
    high = "#08306b",
    name = "Depth (mm)",
    guide = guide_colorbar(reverse = TRUE)  # must go here
  ) +  # blue scale
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank()
  ) +
  ggtitle("PCA Biplot Colored by Depth - 500µm")

print(p.pca)
output_dir <- 'plots/fig_pres/'
filename <- file.path(output_dir, "PCA_all_elem_depth.svg")
save_svg_plot(p.pca, filename, width = 15, height = 10)

# ---- by cluster ----
p.pca.kclust <- ggplot() +
  # Sample points colored by depth
  geom_point(data = scores, aes(x = PC1, y = PC2, color = kcluster),
             shape = 16, size = 2.5) +
  
  # ellipse over points
  stat_ellipse(data = scores, aes(x = PC1, y = PC2, color = kcluster), 
               type = "norm", level = 0.68) +
  # Loadings arrows
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.6) +
  
  # Loadings labels
  geom_text_repel(data = loadings_scaled,
                  aes(x = PC1 * 3.2, y = PC2 * 3.2, label = variable),
                  size = 5, color = "black") +
  
  # Center lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  
  # Axes and theme
  xlab(paste0("PC1 (", var_expl[1], "%)")) +
  ylab(paste0("PC2 (", var_expl[2], "%)")) +
  coord_equal() +
  scale_colour_manual(values = c("#8856a7","#b3cde3")) + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank()
  ) +
  ggtitle("PCA Biplot Colored by K-means Cluster") 

print(p.pca.kclust)
output_dir <- 'plots/fig_pres/'
filename <- file.path(output_dir, "PCA_all_elem_cluster.svg")
save_svg_plot(p.pca.kclust, filename, width = 15, height = 10)




# ---- PCA using only elements, no ratios, 200µm ----

combined.elem.200 <- read.csv('/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/generated/combined/combined_200_inc_coh.csv')


xrf_elem_sel <- c("Fe","Mn",'Ba',"S","Ti","Al","Si","K", 'Ca', 'Rb', 'Zr', 'Sr', 'P', 'Ba','position..mm.')
non_xrf_sel <- c('position..mm.', 'TCWt')
all_sel <- append(xrf_elem_sel, non_xrf_sel)

res <- 0.2
half_res <- res/2
df.raw <- load_raw_data()
df.xrf <- clean_df(df.raw, sel = xrf_elem_sel)
df.xrf$raw.500 <- NULL
df.xrf <- tibble(df.xrf$raw.200) %>%
  select_top_varves()

df.cs <- closed_sum(df.xrf) 
df.clr <- clr_transform(df.xrf)
df.pca <- df.clr %>%
  select(-position..mm.)

pca_result <- prcomp(df.pca, scale. = TRUE, center = TRUE)
pca_result <- fix_pca_signs(pca_result)

# prepare cluster data
X <- scale(df.pca)

# Choose optimal number of clusters (here 2 is optimal)
fviz_nbclust(X, kmeans, method = "silhouette") +
  theme_minimal()


# ---- Perform k-means -----
km <- kmeans(X, centers = 2, nstart = 25)

combined.elem.200$kmeans_cluster <- as.factor(km$cluster)

# ---- Extract scores and loadings ----
scores <- as_tibble(pca_result$x)
loadings <- as_tibble(pca_result$rotation, rownames = "variable")

# Rescale loadings to match sample space
loadings_scaled <- loadings %>%
  mutate(
    PC1 = PC1 * pca_result$sdev[1],
    PC2 = PC2 * pca_result$sdev[2],
    PC3 = PC3 * pca_result$sdev[3]
  )

#  Add metadata for faceting
scores <- scores %>%
  mutate(depth = combined.elem.200$position..mm.) %>%
  mutate(facies = combined.elem.200$facies_class) %>%
  mutate(kcluster = combined.elem.200$kmeans_cluster)

var_expl <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# ---- PLOTTING ----
# ---- plotting by depth ----

# by depth
p.pca <- ggplot() +
  # Sample points colored by depth
  geom_point(data = scores, aes(x = PC1, y = PC3, colour = depth),
             shape = 16, size = 2.5) +
  
  # Loadings arrows
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC3 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.6) +
  
  # Loadings labels
  geom_text_repel(data = loadings_scaled,
                  aes(x = PC1 * 3.2, y = PC3 * 3.2, label = variable),
                  size = 6, color = "black") +
  
  # Center lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  
  # Axes and theme
  xlab(paste0("PC1 (", var_expl[1], "%)")) +
  ylab(paste0("PC3 (", var_expl[3], "%)")) +
  coord_equal() +
  scale_colour_gradient(
    low = "#c6dbef",
    high = "#08306b",
    name = "Depth (mm)",
    guide = guide_colorbar(reverse = TRUE)  # must go here
  ) +  # blue scale
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank()
  ) +
  ggtitle("PCA Biplot Colored by Depth - 200µm")

print(p.pca)
output_dir <- 'plots/fig_pres/'
filename <- file.path(output_dir, "PCA_1_3_all_elem_depth.svg")
save_svg_plot(p.pca, filename, width = 15, height = 10)

# ---- by cluster ----
p.pca.kclust <- ggplot() +
  # Sample points colored by depth
  geom_point(data = scores, aes(x = PC1, y = PC3, color = kcluster),
             shape = 16, size = 2.5) +
  
  # ellipse over points
  stat_ellipse(data = scores, aes(x = PC1, y = PC3, color = kcluster), 
               type = "norm", level = 0.68) +
  # Loadings arrows
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC3 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.6) +
  
  # Loadings labels
  geom_text_repel(data = loadings_scaled,
                  aes(x = PC1 * 3.2, y = PC3 * 3.2, label = variable),
                  size = 5, color = "black") +
  
  # Center lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  
  # Axes and theme
  xlab(paste0("PC1 (", var_expl[1], "%)")) +
  ylab(paste0("PC3 (", var_expl[3], "%)")) +
  coord_equal() +
  scale_colour_manual(values = c("#8856a7","#b3cde3")) + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank()
  ) +
  ggtitle("PCA Biplot Colored by K-means Cluster") 

print(p.pca.kclust)
output_dir <- 'plots/fig_pres/'
filename <- file.path(output_dir, "PCA_1_3_all_elem_cluster.svg")
save_svg_plot(p.pca.kclust, filename, width = 15, height = 10)


# ---- PCA using only elements, no ratios, 200µm ALL CORE ----

combined.elem.200.AC <- read.csv('/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/generated/combined/combined_AC_200_inc_coh.csv')


xrf_elem_sel <- c("Fe","Mn",'Ba',"S","Ti","Al","Si","K", 'Ca', 'Rb', 'Zr', 'Sr', 'P', 'Ba','position..mm.')
non_xrf_sel <- c('position..mm.', 'TCWt')
all_sel <- append(xrf_elem_sel, non_xrf_sel)

res <- 0.2
half_res <- res/2
df.raw <- load_raw_data()
df.xrf <- clean_df(df.raw, sel = xrf_elem_sel)
df.xrf$raw.500 <- NULL
df.xrf <- tibble(df.xrf$raw.200)

df.cs <- closed_sum(df.xrf) 
df.clr <- clr_transform(df.xrf)
df.pca <- df.clr %>%
  select(-position..mm.)

pca_result <- prcomp(df.pca, scale. = TRUE, center = TRUE)
pca_result <- fix_pca_signs(pca_result)

# prepare cluster data
X <- scale(df.pca)

# Choose optimal number of clusters (here 3 is optimal)
fviz_nbclust(X, kmeans, method = "silhouette") +
  theme_minimal()


# ---- Perform k-means -----
km <- kmeans(X, centers = 3, nstart = 25)

combined.elem.200.AC$kmeans_cluster <- as.factor(km$cluster)

# ---- Extract scores and loadings ----
scores <- as_tibble(pca_result$x)
loadings <- as_tibble(pca_result$rotation, rownames = "variable")

# Rescale loadings to match sample space
loadings_scaled <- loadings %>%
  mutate(
    PC1 = PC1 * pca_result$sdev[1],
    PC2 = PC2 * pca_result$sdev[2],
    PC3 = PC3 * pca_result$sdev[3]
  )

#  Add metadata for faceting
scores <- scores %>%
  mutate(depth = combined.elem.200.AC$position..mm.) %>%
  mutate(facies = combined.elem.200.AC$facies_class) %>%
  mutate(kcluster = combined.elem.200.AC$kmeans_cluster)

var_expl <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# ---- PLOTTING ----
# ---- plotting by depth ----

# by depth
p.pca <- ggplot() +
  # Sample points colored by depth
  geom_point(data = scores, aes(x = PC1, y = PC2, colour = depth),
             shape = 16, size = 2.5) +
  
  # Loadings arrows
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.6) +
  
  # Loadings labels
  geom_text_repel(data = loadings_scaled,
                  aes(x = PC1 * 3.2, y = PC2 * 3.2, label = variable),
                  size = 6, color = "black") +
  
  # Center lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  
  # Axes and theme
  xlab(paste0("PC1 (", var_expl[1], "%)")) +
  ylab(paste0("PC2 (", var_expl[2], "%)")) +
  coord_equal() +
  scale_colour_gradient(
    low = "#c6dbef",
    high = "#08306b",
    name = "Depth (mm)",
    guide = guide_colorbar(reverse = TRUE)  # must go here
  ) +  # blue scale
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank()
  ) +
  ggtitle("PCA Biplot Colored by Depth - 200µm ALL CORE")

print(p.pca)
output_dir <- 'plots/fig_pres/'
filename <- file.path(output_dir, "PCA_AC_1_2_200_all_elem_depth_.svg")
save_svg_plot(p.pca, filename, width = 15, height = 10)

# ---- by cluster ----
p.pca.kclust <- ggplot() +
  # Sample points colored by depth
  geom_point(data = scores, aes(x = PC1, y = PC2, color = kcluster),
             shape = 16, size = 2.5) +
  
  # ellipse over points
  stat_ellipse(data = scores, aes(x = PC1, y = PC2, color = kcluster), 
               type = "norm", level = 0.68) +
  # Loadings arrows
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.6) +
  
  # Loadings labels
  geom_text_repel(data = loadings_scaled,
                  aes(x = PC1 * 3.2, y = PC2 * 3.2, label = variable),
                  size = 5, color = "black") +
  
  # Center lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  
  # Axes and theme
  xlab(paste0("PC1 (", var_expl[1], "%)")) +
  ylab(paste0("PC2 (", var_expl[2], "%)")) +
  coord_equal() +
  scale_colour_manual(values = c('#af8dc3','#f7f7f7','#7fbf7b')) + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank()
  ) +
  ggtitle("PCA Biplot Colored by K-means Cluster") 

print(p.pca.kclust)
output_dir <- 'plots/fig_pres/'
filename <- file.path(output_dir, "PCA_AC_1_2_200_all_elem_cluster.svg")
save_svg_plot(p.pca.kclust, filename, width = 15, height = 10)


# ---- PCA using only elements, no ratios, 500µm ALL CORE ----

combined.elem.500.AC <- read.csv('/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/generated/combined/combined_AC_500_inc_coh.csv')


xrf_elem_sel <- c("Fe","Mn",'Ba',"S","Ti","Al","Si","K", 'Ca', 'Rb', 'Zr', 'Sr', 'P', 'Ba','position..mm.')
non_xrf_sel <- c('position..mm.', 'TCWt')
all_sel <- append(xrf_elem_sel, non_xrf_sel)

res <- 0.5
half_res <- res/2
df.raw <- load_raw_data()
df.xrf <- clean_df(df.raw, sel = xrf_elem_sel)
df.xrf$raw.200 <- NULL
df.xrf <- tibble(df.xrf$raw.500)

df.cs <- closed_sum(df.xrf) 
df.clr <- clr_transform(df.xrf)
df.pca <- df.clr %>%
  select(-position..mm.)

pca_result <- prcomp(df.pca, scale. = TRUE, center = TRUE)
pca_result <- fix_pca_signs(pca_result)

# prepare cluster data
X <- scale(df.pca)

# Choose optimal number of clusters (here 2 is optimal)
fviz_nbclust(X, kmeans, method = "silhouette") +
  theme_minimal()


# ---- Perform k-means -----
km <- kmeans(X, centers = 2, nstart = 25)

combined.elem.500.AC$kmeans_cluster <- as.factor(km$cluster)

# ---- Extract scores and loadings ----
scores <- as_tibble(pca_result$x)
loadings <- as_tibble(pca_result$rotation, rownames = "variable")

# Rescale loadings to match sample space
loadings_scaled <- loadings %>%
  mutate(
    PC1 = PC1 * pca_result$sdev[1],
    PC2 = PC2 * pca_result$sdev[2],
    PC3 = PC3 * pca_result$sdev[3]
  )

#  Add metadata for faceting
scores <- scores %>%
  mutate(depth = combined.elem.500.AC$position..mm.) %>%
  mutate(facies = combined.elem.500.AC$facies_class) %>%
  mutate(kcluster = combined.elem.500.AC$kmeans_cluster)

var_expl <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# ---- PLOTTING ----
# ---- plotting by depth ----

# by depth
p.pca <- ggplot() +
  # Sample points colored by depth
  geom_point(data = scores, aes(x = PC1, y = PC2, colour = depth),
             shape = 16, size = 2.5) +
  
  # Loadings arrows
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.6) +
  
  # Loadings labels
  geom_text_repel(data = loadings_scaled,
                  aes(x = PC1 * 3.2, y = PC2 * 3.2, label = variable),
                  size = 6, color = "black") +
  
  # Center lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  
  # Axes and theme
  xlab(paste0("PC1 (", var_expl[1], "%)")) +
  ylab(paste0("PC2 (", var_expl[2], "%)")) +
  coord_equal() +
  scale_colour_gradient(
    low = "#c6dbef",
    high = "#08306b",
    name = "Depth (mm)",
    guide = guide_colorbar(reverse = TRUE)  # must go here
  ) +  # blue scale
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank()
  ) +
  ggtitle("PCA Biplot Colored by Depth - 500µm ALL CORE")

print(p.pca)
output_dir <- 'plots/fig_pres/'
filename <- file.path(output_dir, "PCA_AC_1_2_500_all_elem_depth_.svg")
save_svg_plot(p.pca, filename, width = 15, height = 10)

# ---- by cluster ----
p.pca.kclust <- ggplot() +
  # Sample points colored by depth
  geom_point(data = scores, aes(x = PC1, y = PC2, color = kcluster),
             shape = 16, size = 2.5) +
  
  # ellipse over points
  stat_ellipse(data = scores, aes(x = PC1, y = PC2, color = kcluster), 
               type = "norm", level = 0.68) +
  # Loadings arrows
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.6) +
  
  # Loadings labels
  geom_text_repel(data = loadings_scaled,
                  aes(x = PC1 * 3.2, y = PC2 * 3.2, label = variable),
                  size = 5, color = "black") +
  
  # Center lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  
  # Axes and theme
  xlab(paste0("PC1 (", var_expl[1], "%)")) +
  ylab(paste0("PC2 (", var_expl[2], "%)")) +
  coord_equal() +
  scale_colour_manual(values = c('#8856a7','#e0ecf4')) + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank()
  ) +
  ggtitle("PCA Biplot Colored by K-means Cluster - 500µm All Core") 

print(p.pca.kclust)
output_dir <- 'plots/fig_pres/'
filename <- file.path(output_dir, "PCA_AC_1_2_500_all_elem_cluster.svg")
save_svg_plot(p.pca.kclust, filename, width = 15, height = 10)
