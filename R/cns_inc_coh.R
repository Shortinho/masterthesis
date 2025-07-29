# CNS plotting to compare with Cr inc and coh

#### interpolated CNS dataset xrf 500µm resolution ############################

# have combined_df loaded in global env (either from sync_datasets script or file)

om.df <- combined_df %>% select(c('position..mm.','TCWt', 'TSWt','Cr.inc_Cr.coh'))

plot_individual_ratio(om.df, 'oc_comparison', output_dir = 'plots/om.compare')

cor(om.df$TCWt, om.df$Cr.inc_Cr.coh) # low correlation (r=0.255)


#### use non interpolated dataset #############################################

carbon <- raw.cns %>% 
  select(any_of(c('Depth', 'TC%/Wt'))) %>%
  rename(position..mm. = Depth)

inc.coh <- combined_df %>% 
  select(c(position..mm., Cr.inc_Cr.coh))

library(fuzzyjoin)

inc.coh_matched <- difference_join(
  carbon, inc.coh,
  by = "position..mm.",
  max_dist = 0.25,      # or more if needed (units are mm)
  distance_col = NULL
) %>%
  rename(position..mm. = position..mm..x) %>%
  rename('TCWt' = 'TC%/Wt')

cor(inc.coh_matched$TCWt, inc.coh_matched$Cr.inc_Cr.coh) # low pearson correlation (r = 0.26)
cor(inc.coh_matched$TCWt, inc.coh_matched$Cr.inc_Cr.coh, method = 'spearman')

plot_individual_ratio(inc.coh_matched, 'oc_comparison_lowres', output_dir = 'plots/om.compare')

#### non interpolated - compare sulfur and CNS ################################

sulfur <- combined_df %>% 
  select(c(position..mm., S_Ti))

cns.sulf.joined <- difference_join(
  carbon, sulfur,
  by = 'position..mm.',
  max_dist = 0.25,
  distance_col = NULL) %>%
  rename(position..mm. = position..mm..x) %>%
  rename('TCWt' = 'TC%/Wt')

cor(cns.sulf.joined$TCWt, cns.sulf.joined$S_Ti) # 0.42

plot_individual_ratio(cns.sulf.joined, 'oc_comparison_lowres_sulf', output_dir = 'plots/om.compare')

cor(cns.sulf.joined$S_Ti, inc.coh_matched$Cr.inc_Cr.coh) # r = 0.7530153
