# CNS plotting to compare with Cr inc and coh

#### interpolated CNS dataset xrf 500µm resolution ############################

# have combined_df loaded in global env (either from sync_datasets script or file)
combined_df <- read.csv('data/generated/combined/combined_AC_500_inc_coh.csv')

path_CNS <- 'data/POS-22-20_CNS.xlsx'
raw.cns <- read_excel(path_CNS) %>%
  tibble()
# convert depth from cm to mm
raw.cns$Depth <- raw.cns$Depth*10
cns <- rename(raw.cns, position..mm. = Depth)
cns <- cns %>%
  rename_with(~ gsub('[%/\" ]', '', .x), .cols = matches('TN|TC|TS'))

om.df <- combined_df %>% select(c('position..mm.','TCWt', 'TSWt','Cr.inc_Cr.coh'))

plot_individual_ratio(om.df, 'oc_comparison', output_dir = 'plots/om.compare.AC')

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

# prep data for two lines side by side plotting
range_TCWt <- range(combined_df$TCWt, na.rm = TRUE)
range_TNWt <- range(combined_df$Cr.inc_Cr.coh, na.rm = TRUE)

rescale_TNWt <- function(x) {
  scales::rescale(x, to = range_TCWt, from = range_TNWt)
}

# Inverse function for secondary axis
inv_rescale_TNWt <- function(x) {
  scales::rescale(x, to = range_TNWt, from = range_TCWt)
}

inc.coh_matched$inc.coh.rescaled <- rescale_TNWt(inc.coh_matched$Cr.inc_Cr.coh)

library(ggplot2)

# Reshape the data to long format for easier legend control
library(tidyr)
long_data <- tidyr::pivot_longer(
  inc.coh_matched,
  cols = c(TCWt, inc.coh.rescaled),
  names_to = "Variable",
  values_to = "Value"
)

# Assign desired colors
line_colors <- c(
  TCWt = "black",
  inc.coh.rescaled = "steelblue"
)

# Plot
p.compare <- ggplot(data = long_data, aes(x = position..mm., y = Value, color = Variable)) +
  geom_line() +
  scale_color_manual(
    values = line_colors,
    labels = c(
      TCWt = "Total carbon (%)",
      inc.coh.rescaled = "CrInc/CrCoh"
    )
  ) +
  scale_y_continuous(
    name = "TC (wt%)",
    sec.axis = sec_axis(~inv_rescale_TNWt(.), name = "CrInc/CrCoh")
  ) +
  scale_x_reverse() +
  coord_flip() +
  theme_minimal() +
  theme(
    aspect.ratio = 3,
    axis.text.y.left = element_text(size = 7, color = "black"),
    axis.text.y.right = element_text(size = 7, color = "steelblue"),
    axis.title.y.left = element_text(size = 8, color = "black"),
    axis.title.y.right = element_text(size = 8, color = "steelblue"),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    axis.ticks.x = element_line(),
    axis.ticks.length = unit(.15, "cm"),
    axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )


print(p.compare)

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
