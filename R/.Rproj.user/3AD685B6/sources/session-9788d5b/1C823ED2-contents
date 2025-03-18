# Top varves XRF data analysis

# Source the functions
source("xrf_functions.R")

#### Import data ####

df_list <- load_data()

#### Execution ####

# Clean dataframes
df_l <- clean_df(df_list)
rm(df_list)

# Filter out the top portion
df_l <- map(df_l, select_top_varves)

# Create and save correlation maps
map(names(df_l), function(df_name) {
  pdf(paste0("plots/corr_maps_top/correlation_map_toppart", df_name, ".pdf"))
  create_corr_map(df_l[[df_name]])
  dev.off()
})

# Create and save boxplots
map(names(df_l), function(df_name) {
  pdf(paste0("boxplot_", df_name, ".pdf"))
  print(create_boxplots(df_l[[df_name]]))
  dev.off()
})

# Compute ratios
ratios <- map(df_l, compute_ratios)

# Perform PCA
# Perform PCA and generate plots (no need to scale data beforehand)
pca_resutls <- map2(df_l, names(df_l), ~perform_pca(.x, df_name = .y, 
                                                    plot_scores = F, 
                                                    plot_loadings = F, 
                                                    print_results = F,
                                                    output_dir = "plots/PCA_plots_top"))

element_plots <- map2(df_l, names(df_l), ~plot_elements(.x, df_name = .y,
                                                        individual_plots = T,
                                                        output_dir = "plots/indiv_element_plots_top"))

# individual element plots 
indiv_plotting_top <- map2(df_l, names(df_l), ~plot_individual_top(.x, .y, "plots/individual_top"))

indiv_plotting_clustered <- map2(df_cluster, names(df_cluster), ~plot_individual_clust(.x, .y, "plots/individual_cluster_top"))
#clustering
df_cluster <- map2(df_l, names(df_l), ~perform_clustering(.x, .y, k = 2, method = "kmeans"))

map2(df_cluster, names(df_cluster), ~visualize_clusters(.x, .y, output_dir = 'plots/clusterViz_top'))

