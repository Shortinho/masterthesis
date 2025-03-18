#XRF preliminary data analysis


# Source the functions
source("xrf_functions.R")

#### Import data ####

df_list <- load_data()

#### Execution ####

# Clean dataframes
df_l <- clean_df(df_list)
rm(df_list)

# Create and save correlation maps
map(names(df_l), function(df_name) {
  pdf(paste0("correlation_map_", df_name, ".pdf"))
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
# Perform PCA and generate plots
pca_resutls <- map2(df_l, names(df_l), ~perform_pca(.x, df_name = .y, 
                                                    plot_scores = T, 
                                                    plot_loadings = T, 
                                                    print_results = F,
                                                    output_dir = "plots/PCA_plots"))

element_plots <- map2(df_l, names(df_l), ~plot_elements(.x, df_name = .y,
                                                        individual_plots = T,
                                                        output_dir = "plots/indiv_element_plots"))

# individual element plots 
indiv_plotting <- map2(df_l, names(df_l), ~plot_individual(.x, .y, "plots/individual"))

indiv_plotting_clustered <- map2(df_cluster, names(df_cluster), ~plot_individual_clust(.x, .y, "plots/individual_cluster"))
#clustering
df_cluster <- map2(df_l, names(df_l), ~perform_clustering(.x, .y, k = 3, method = "kmeans"))

map2(df_cluster, names(df_cluster), ~visualize_clusters(.x, .y, output_dir = 'plots/clusterViz'))

