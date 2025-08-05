# synchronize all datasets for top 35cm

# load all datasets

source('xrf_functions.R')

library(readxl)

sel <- c('position..mm.','Fe', 'Mn', 'S', 'Ti', 'Al', 'Si', 'K', 'Rb', 'Zr', 'Ca', 'Sr', 'P', 'Ba')


#### XRF 0.5mm dataset ####
res <- 0.5
half_res <- res/2
df.raw <- load_raw_data()
df.xrf <- clean_df(df.raw, sel = sel)
df.xrf$raw.200 <- NULL
df.xrf <- tibble(df.xrf$raw.500) %>%
  select_top_varves()

df.cs <- closed_sum(df.xrf) 
df.clr <- clr_transform(df.xrf)

# extract Cr inc and Cr coh 500µm

inc.coh <- df.raw$raw.500 %>%
  select('position (mm)', 'Cr inc', 'Cr coh')
names <- names(inc.coh)
tidy.names <- names %>% make.names(unique = T)
names(inc.coh) <- tidy.names
inc.coh <- filter(inc.coh, `position..mm.`> 25& `position..mm.` < 1231.6) %>%
  select_top_varves()


#### HSI dataset ####
path_HSI <- '/Users/maxshore/Documents/Unibe/MasterThesis/POS-22-20_DATA/POS_22_20_HSI/RABD673_Chla/POS22-20_230802-145515_refl_sub_RABD673_spl_645_674.csv'

raw.HSI <- read_excel('/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/POS_22_20_HSI/HSI_calibrated.xlsx') %>%
  tibble()
names(raw.HSI) <- make.names(names(raw.HSI), unique = T)
hsi <- raw.HSI %>% select(c(Core.Depth..mm., RABD673, Moving.Average, TChl.ug.g))
hsi <- rename(hsi, position..mm.= Core.Depth..mm.) %>%
  select_top_varves()

#### CNS dataset ####


path_CNS <- '/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/POS-22-20_CNS.xlsx'

raw.cns <- read_excel(path_CNS) %>%
  tibble()
# convert depth from cm to mm
raw.cns$Depth <- raw.cns$Depth*10
cns <- rename(raw.cns, position..mm. = Depth) %>%
  select_top_varves()
cns <- cns %>%
  rename_with(~ gsub('[%/\" ]', '', .x), .cols = matches('TN|TC|TS'))

#### grayscale dataset on 20 bandwidth ####
path_gs <- '/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/generated/grayscale/gray_profile_bw50groups_2filter_0.2mm.csv'

gs <- read.csv(path_gs)


#### synchronize ###############################################################

# Separate facies from grayscale
gs_fac <- gs %>% select(position..mm., facies_class)
gs_gray <- gs %>% select(position..mm., mean_gray)

gray_sync <- synchronize_to_reference(df.xrf, gs_gray,
                                      ref_depth_col = "position..mm.",
                                      target_depth_col = "position..mm.",
                                      target_vars = c("mean_gray"))

library('fuzzyjoin')
facies_sync_all <- fuzzyjoin::difference_left_join(
  df.xrf, gs_fac, 
  by = "position..mm.",
  max_dist = half_res,
  distance_col = "dist"
)

facies_sync <- facies_sync_all %>%
  group_by(position..mm..x) %>%
  slice_min(order_by = dist, n = 1) %>%
  ungroup() %>%
  select(position..mm. = position..mm..x, facies_class)


hsi_sync <- synchronize_to_reference(df.xrf, hsi,
                                     target_vars = c("RABD673", 'Moving.Average', 'TChl.ug.g'))

cns_sync <- synchronize_to_reference(df.xrf, cns,
                                     target_vars = c("TNWt", "TCWt", "TSWt"))


combined_df <- df.xrf %>%
  left_join(gray_sync, by = "position..mm.") %>%
  left_join(hsi_sync, by = "position..mm.") %>%
  left_join(cns_sync, by = "position..mm.") %>%
  left_join(df.clr, by = 'position..mm.') %>% #.y
  left_join(df.cs, by = 'position..mm.') %>% # normal no suffix
  left_join(facies_sync, by = "position..mm.") %>% 
  left_join(inc.coh, by = 'position..mm.')

# compute ratios (make sure to perform on the correct xrf normalisation (cs))
combined_df <- compute_ratios(combined_df, 'Mn', 'Fe')
combined_df <- compute_ratios(combined_df, 'S', 'Ti')
combined_df <- compute_ratios(combined_df, 'Fe', 'Ti')
combined_df <- compute_ratios(combined_df, 'Si', 'Ti')
combined_df <- compute_ratios(combined_df, 'Ba', 'Ti')
combined_df <- compute_ratios(combined_df, 'Ca', 'Ti')
combined_df <- compute_ratios(combined_df, 'Ti', 'Al')
combined_df <- compute_ratios(combined_df, 'Si', 'Al')
combined_df <- compute_ratios(combined_df, 'Zr', 'Rb')
combined_df <- compute_ratios(combined_df, 'TCWt', 'TNWt')
combined_df <- compute_ratios(combined_df, 'Cr.inc', 'Cr.coh')

# save combined df
write.csv(combined_df, file = 'data/generated/combined/combined_500_inc_coh.csv')

#### descriptive stats by facies ###############################################
vars_to_summarize <- c(
  # closed-sum elements
  'Fe', 'Mn', 'S', 'Ti', 'Al', 'Si', 'K', 'Rb', 'Zr', 'Ca', 'Sr', 'P', 'Ba',
  # centered log ratio elements
  'Fe.y', 'Mn.y', 'S.y', 'Ti.y', 'Al.y', 'Si.y', 'K.y', 'Rb.y', 'Zr.y', 'Ca.y', 'Sr.y', 'P.y', 'Ba.y',
  # raw counts
  'Fe.x', 'Mn.x', 'S.x', 'Ti.x', 'Al.x', 'Si.x', 'K.x', 'Rb.x', 'Zr.x', 'Ca.x', 'Sr.x', 'P.x', 'Ba.x',
  # ratios
  'Mn_Fe', 'S_Ti', 'Fe_Ti', 'Si_Ti', 'Ti_Al', 'Ba_Ti', 'Ca_Ti',
  # CNS
  'TNWt', 'TCWt', 'TSWt',
  # HSI
  'RABD673', 'Moving.Average',
  # grayscale
  'mean_gray'
)

summary_all <- combined_df %>%
  select(all_of(vars_to_summarize)) %>%
  summarise(across(everything(),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE),
                        min  = ~min(.x, na.rm = TRUE),
                        max  = ~max(.x, na.rm = TRUE),
                        median = ~median(.x, na.rm = TRUE)
                   ),
                   .names = "{.col}_{.fn}"
  ))

summary_by_facies <- combined_df %>%
  group_by(facies_class) %>%
  summarise(across(all_of(vars_to_summarize),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE),
                        min  = ~min(.x, na.rm = TRUE),
                        max  = ~max(.x, na.rm = TRUE),
                        median = ~median(.x, na.rm = TRUE)
                   ),
                   .names = "{.col}_{.fn}"
  ),
  .groups = "drop"
  )

vars_to_test <- vars_to_summarize


test_results <- lapply(vars_to_test, function(var) {
  formula <- as.formula(paste(var, "~ facies_class"))
  result <- wilcox.test(formula, data = combined_df)
  tibble(variable = var, p_value = result$p.value, method = result$method)
})

test.df <- bind_rows(test_results) # small p-values indicate difference in distributions

write.csv(test.df, file = "data/generated/tables/facies_test.csv", row.names = FALSE)
write.csv(summary_by_facies, file = "data/generated/tables/facies_basic_stats.csv", row.names = FALSE)
write.csv(summary_all, file = "data/generated/tables/basic_stats.csv", row.names = FALSE)

### plotting of basic stats ###################################################
# Variables to plot
plot_vars <- vars_to_summarize

# Output folder
output_dir <- "plots/facies_comparison"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Plotting function
plot_facies_comparison <- function(var) {
  p_box <- ggplot(combined_df, aes(x = factor(facies_class), y = .data[[var]])) +
    geom_boxplot() +
    labs(x = "Facies", y = var, title = paste("Boxplot of", var, "by Facies")) +
    theme_bw()
  
  p_violin <- ggplot(combined_df, aes(x = factor(facies_class), y = .data[[var]])) +
    geom_violin(fill = "gray", alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.size = 0.5) +
    labs(x = "Facies", y = var, title = paste("Violin Plot of", var, "by Facies")) +
    theme_bw()
  
  p_density <- ggplot(combined_df, aes(x = .data[[var]], fill = factor(facies_class))) +
    geom_density(alpha = 0.5) +
    labs(x = var, fill = "Facies", title = paste("Density Plot of", var, "by Facies")) +
    theme_bw()
  
  ggsave(filename = file.path(output_dir, paste0(var, "_boxplot.png")), plot = p_box, width = 5, height = 5)
  ggsave(filename = file.path(output_dir, paste0(var, "_violin.png")), plot = p_violin, width = 5, height = 5)
  ggsave(filename = file.path(output_dir, paste0(var, "_density.png")), plot = p_density, width = 5, height = 5)
}

# Loop over all variables
#walk(plot_vars, plot_facies_comparison)




#### correlation btw datasets ##################################################

variance_sel <- c("Fe","Mn","S","Ti","Al","Si","K","Rb","Zr","Ca","Sr","P","Ba","Mn_Fe","S_Ti","Fe_Ti","Si_Ti","Ba_Ti", 'Ca_Ti',"Ti_Al",'Si_Al','Zr_Rb', "mean_gray", 'TChl.ug.g',"TNWt","TCWt","TSWt",'TCWt_TNWt', 'position..mm.')

combined_df %>% 
  select(any_of(variance_sel)) %>%
  create_corr_map()

#### PCA #########################################################

combined_df %>%
  select(any_of(variance_sel)) %>%
  perform_pca2(plot_loadings = T,
               output_dir = 'plots/PCA_variable_selection')

pca_sel1 <- c("Ti","K","Rb","Zr","Ca","Sr","P","Ba","Mn_Fe","S_Ti","Fe_Ti","Si_Ti","Ba_Ti", 'Ca_Ti',"Ti_Al",'Si_Al','Zr_Rb', "mean_gray",'Moving.Average',"TCWt", 'position..mm.')

combined_df %>%
  select(any_of(pca_sel1)) %>%
  perform_pca2(df_name = paste(pca_sel1, collapse = ', '),
               plot_loadings = T,
               output_dir = 'plots/PCA_variable_selection')

pca_sel2 <- c("Ti","K","Rb","Zr","Ca","Mn_Fe","S_Ti","Fe_Ti","S_Ti","Ba_Ti", 'Ca_Ti',"Ti_Al",'Zr_Rb', "mean_gray",'Moving.Average',"TCWt", 'position..mm.')

combined_df %>%
  select(any_of(pca_sel2)) %>%
  perform_pca2(df_name = paste(pca_sel2, collapse = ', '),
               plot_loadings = T,
               output_dir = 'plots/PCA_variable_selection')



