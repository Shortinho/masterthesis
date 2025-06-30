# FILE WORKING WITH IMPORTED RAW DATA
source('xrf_functions.R')

#### data prep and transformation ####

# import
df.raw <- load_raw_data()
sel_pc <- c('position..mm.','Al','Si','P','S','K','Ca','Ti','Fe','Ni','Sr','Zr')
df.clean <- clean_df(df.raw, sel = sel_pc)
df.clean.top <- map(df.clean, select_top_varves)
df.clr <- map(df.clean, clr_transform)
df.clr.top <- map(df.clean.top, clr_transform)
pca.500.2 <- perform_pca2(df.clr$raw.500, df_name = 'clr500_all_core2_sel', plot_scores = F, plot_loadings = T, output_dir = 'plots/pca_20_6')
pca.500.2.top <- perform_pca2(df.clr.top$raw.500, df_name = 'clr500_top_core2_sel', plot_scores = F, plot_loadings = T, output_dir = 'plots/pca_20_6')
pca_downcore_plot(pca.500.2, df_name = 'clr500_sel_2', 'PC2', 'plots/downcore_PCA')
pca_downcore_plot(pca.500.2.top, df_name = 'clr500_top_sel_2', 'PC2', 'plots/downcore_PCA')

# Selection of elements
all_elements <- c("Al", "Si", "P", "S", "K", "Ca", "Ti", "Mn", "Fe", "Ni", "Zn", "Rb", "Sr", "Zr", "Ba")
selection <- c("Si", "P", "S", "K", "Ca", "Ti", "Mn", "Fe", "Zn", "Rb", "Sr", "Ba")
selection2 <- c("P", "S", "K", "Ca", "Ti", "Fe", "Zn", "Rb", "Ba")
selection <- c('position..mm.','Al','Ar','Si','P','S','K','Ca','Ti','Mn','Fe','Ni','Zn','Rb','Sr','Zr','Ba')
sel_pc <- c('position..mm.','Al','Si','P','S','K','Ca','Ti','Fe','Ni','Sr','Zr')

# remove unwanted variables (keep only position and element counts)
df.clean <- clean_df(df.raw, sel = sel_pc)
df.clean.top <- map(df.clean, select_top_varves)



df.clean.partial <- map(df.clean, ~select_elements(.x, elements = sel_pc))
df.clean.partial.top <- map(df.clean.top, ~select_elements(.x, elements = selection))

df.clean.partial.top2 <- map(df.clean.top, ~select_elements(.x, elements = selection2))

#### transforms ####
# closed sum
df.cs <- map(df.clean, closed_sum)

# centered log ratio
df.clr <- map(df.clean, clr_transform)
df.clr.top <- map(df.clean.top, clr_transform)
df.clr.partial <- map(df.clean.partial, clr_transform)
df.clr.partial.top <- map(df.clean.partial.top, clr_transform)

df.clr.partial.top2 <- map(df.clean.partial.top2, clr_transform)

# data is now ready for multivariate analysis!

#### PCA for all core ####

# all elements
pca.200 <- perform_pca(df.clr$raw.200, df_name = 'clr200_all_core', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots')

pca.500 <- perform_pca(df.clr$raw.500, df_name = 'clr500_all_core', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots')

pca.200.2 <- perform_pca2(df.clr$raw.200, df_name = 'clr200_all_core2', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots')

pca.500.2 <- perform_pca2(df.clr$raw.500, df_name = 'clr500_all_core2_sel', plot_scores = F, plot_loadings = T, output_dir = 'plots/pca_20_6')

pca.sel2 <- perform_pca2(df.clr.partial.top2$raw.200, df_name = 'clr200_top_core_selection2', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots_top')
# selection of elements
pca.200.partial <- perform_pca(df.clr.partial$raw.200, df_name = 'clr200_all_core_partial', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots')

pca.500.partial <- perform_pca(df.clr.partial$raw.500, df_name = 'clr500_all_core_partial', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots')

pca.200.partial2 <- perform_pca2(df.clr.partial$raw.200, df_name = 'clr200_all_core_partial2', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots')

pca.500.partial <- perform_pca2(df.clr.partial$raw.500, df_name = 'clr500_all_core_partial2', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots')

#### PCA for top portion of core ####

pca.200.top <- perform_pca(df.clr.top$raw.200, df_name = 'clr200_top_core', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots_top')

pca.500.top <- perform_pca(df.clr.top$raw.500, df_name = 'clr500_top_core', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots_top')

pca.200.top2 <- perform_pca2(df.clr.top$raw.200, df_name = 'clr200_top_core2', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots_top')

pca.500.top2 <- perform_pca2(df.clr.top$raw.500, df_name = 'clr500_top_core2', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots_top')

# selection of elements
pca.200.partial <- perform_pca(df.clr.partial.top$raw.200, df_name = 'clr200_top_core_partial', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots_top')

pca.500.partial <- perform_pca(df.clr.partial.top$raw.500, df_name = 'clr500_top_core_partial', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots_top')

pca.200.partial2 <- perform_pca2(df.clr.partial.top$raw.200, df_name = 'clr200_top_core_partial2', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots_top')

pca.500.partial <- perform_pca2(df.clr.partial.top$raw.500, df_name = 'clr500_top_core_partial2', plot_scores = T, plot_loadings = T, output_dir = 'plots/PCA_plots_top')

#### Fe over Mn ####
df.ratio <- df.raw %>%
  clean_df()
df.ratio <- df.ratio$raw.200 %>% 
  compute_ratios()

selectionFe_mn <- c("P", "S", "K", "Ca", "Ti", "Zn", "Rb", "Ba", "Fe_Mn")
df.ratio.slct <- df.ratio %>% select_elements(elements = selectionFe_mn)
df.ratio.slct.top <- df.ratio.slct %>% select_top_varves()

df.ratio.clr <- df.ratio.slct %>% clr_transform()
df.ratio.clr.top <- df.ratio.slct.top %>% clr_transform()

pca.FeMn <- perform_pca2(df.ratio.clr, df_name = 'Fe_Mn_ratio', plot_loadings = T, output_dir = 'plots/PCA_plots')

pca.FeMn.top <- perform_pca2(df.ratio.clr.top, df_name = 'Fe_Mn_ratio_top', plot_loadings = T, output_dir = 'plots/PCA_plots_top')


#### remove Mn ####

df.ratio <- df.raw %>%
  clean_df()
df.ratio <- df.ratio$raw.200 %>% 
  compute_ratios()

selectionFe_mn <- c("P", "S", "K", "Ca", "Ti", "Zn", "Rb", "Ba", "Fe_Mn")
df.ratio.slct <- df.ratio %>% select_elements(elements = selectionFe_mn)
df.ratio.slct.top <- df.ratio.slct %>% select_top_varves()

df.ratio.clr <- df.ratio.slct %>% clr_transform()
df.ratio.clr.top <- df.ratio.slct.top %>% clr_transform()

pca.FeMn <- perform_pca2(df.ratio.clr, df_name = 'Fe_Mn_ratio', plot_loadings = T, output_dir = 'plots/PCA_plots')

pca.FeMn.top <- perform_pca2(df.ratio.clr.top, df_name = 'Fe_Mn_ratio_top', plot_loadings = T, output_dir = 'plots/PCA_plots_top')


dc_plot_pc1 <- pca_downcore_plot(pca.200, df_name = 'clr200', 'PC1')

pca_downcore_plot(pca.500.2, df_name = 'clr500_sel_2', 'PC1', 'plots/downcore_PCA')
