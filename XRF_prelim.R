#XRF preliminary data analysis

# Load necessary libraries
library(readxl)
library(ggplot2)
library(purrr)

# Source the functions
source("xrf_functions.R")

# Import data
f_raw200mu <- read_excel("path/to/POS-22-20_0.2mm_Cr_slct.xlsx", skip = 2)
f_norm200mu <- read_excel("path/to/POS-22-20_0.2mm_Cr_norm.xlsx", skip = 2)
f_raw500mu <- read_excel("path/to/POS-22-20_0.5mm_Cr_slct.xlsx", skip = 2)
f_norm500mu <- read_excel("path/to/POS-22-20_0.5mm_Cr_norm.xlsx", skip = 2)

# Create a list of dataframes
df_list <- list(
  fraw200 = f_raw200mu,
  fnorm200 = f_norm200mu,
  fraw500 = f_raw500mu,
  fnorm500 = f_norm500mu
)

# Clean dataframes
df_l <- clean_df(df_list)

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
