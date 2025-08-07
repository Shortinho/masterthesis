# plotting of downcore for figure showing characteristics of dark and light layers (presentation 31st of july)

# need to take XRF data only relevant for the plot, and then normalize
# do this by using the sync_datasets.R script and modifying 'sel' variable

sel <- c('position..mm.','Fe', 'Mn', 'S', 'Ti', 'Rb', 'Zr', 'Ca') # sel variable used on this run (file saved in data/generated/combined_500_fig1.csv)

source('xrf_functions.R')

library(dplyr)
library(ggplot2)
library(tidyr)
library(rlang)
library(svglite)

df <- df.multivar
df <- df.combined
var_plot <- 'Fe.x'
ratio_plot_text <- 'Fe (cps)'
element_plot <- ggplot(df, aes(x = position..mm., y = !!sym(var_plot)/1000)) +
  geom_path(lineend = "round", linejoin = "round", color = "black", linewidth = 0.3) +
  labs(
    title = NULL,
    # y = ratio_plot_text,
    y = expression(paste("Fe (", 10^{3},'cps)')),
    x = NULL
  ) +
  theme_minimal(base_size = 8) +
  theme(
    aspect.ratio = 5,
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    axis.ticks.x = element_line(),
    axis.ticks.length=unit(.15, "cm"),
    axis.line.x = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "none"
  ) +
  coord_flip() +
  scale_x_reverse()

print(element_plot)
output_dir <- 'plots/fig_pres/'
filename <- file.path(output_dir, paste0(var_plot, "_line.svg"))
save_svg_plot(element_plot, filename, width = 3, height = 8)




## for CNS on same plot ####

# Reshape to long format for plotting both variables
df_long <- df %>%
  pivot_longer(cols = c(TCWt, TNWt, TSWt), names_to = "variable", values_to = "value")

# Optional: Rename variable labels for clarity in legend
variable_labels <- c(TCWt = "TC (wt%)", TNWt = "TN (wt%)", TSWt = 'TS (wt%)')

# Create plot
library(ggplot2)
library(scales)

df <- combined_df

# Rescale TNWt to align with TCWt (for plotting on same y-scale)
range_TCWt <- range(df$TCWt, na.rm = TRUE)
range_TNWt <- range(df$TNWt, na.rm = TRUE)

# Function to scale TNWt to TCWt scale
rescale_TNWt <- function(x) {
  scales::rescale(x, to = range_TCWt, from = range_TNWt)
}

# Inverse function for secondary axis
inv_rescale_TNWt <- function(x) {
  scales::rescale(x, to = range_TNWt, from = range_TCWt)
}

df$TNWt_rescaled <- rescale_TNWt(df$TNWt)

element_plot <- ggplot(df, aes(x = position..mm.)) +
  geom_path(aes(y = TCWt), color = "black", linewidth = 0.3, lineend = "round", linejoin = "round") +
  geom_path(aes(y = TNWt_rescaled), color = "steelblue", linewidth = 0.3, lineend = "round", linejoin = "round") +
  scale_y_continuous(
    name = "TC (wt%)", 
    sec.axis = sec_axis(~inv_rescale_TNWt(.), name = "TN (wt%)")
  ) +
  labs(
    x = NULL,
    title = NULL
  ) +
  theme_minimal(base_size = 8) +
  theme(
    aspect.ratio = 5,
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
    plot.background = element_blank(),
    legend.position = "none"
  ) +
  coord_flip() +
  scale_x_reverse()

# Save
print(element_plot)
output_dir <- 'plots/fig_pres/'
filename <- file.path(output_dir, "TC_TN_dual_axis.svg")
save_svg_plot(element_plot, filename, width = 3, height = 8)




