### HSI

source('xrf_functions.R')

library(readxl)

path_HSI <- '/Users/maxshore/Documents/Unibe/MasterThesis/POS-22-20_DATA/POS_22_20_HSI/RABD673_Chla/POS22-20_230802-145515_refl_sub_RABD673_spl_645_674.csv'

raw.HSI <- read.csv(file = '/Users/maxshore/Documents/Unibe/MasterThesis/POS-22-20_DATA/POS_22_20_HSI/RABD673_Chla/POS22-20_230802-145515_refl_sub_RABD673_spl_645_674.csv') %>%
  tibble()

HSI <- raw.HSI %>% select(c(Core.Depth..mm., RABD673))
HSI <- rename(HSI, position..mm.= Core.Depth..mm.)

plot_individual_ratio(HSI, df_name = 'HSI', output_dir = 'plots/HSI')
