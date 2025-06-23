# ELEMENTAL RATIOS

source('xrf_functions.R')

#### data prep and transformation ####

# import
df.raw <- load_raw_data()
# remove unwanted variables (keep only position and element counts)
df.clean <- clean_df(df.raw)
# closed sum
df.cs <- map(df.clean, closed_sum)
# centered log ratio
df.clr <- map(df.clean, clr_transform)

# compute ratios 
# done on intensities (counts) following Bertrand et al. 2024 of 200µ data
ratio.df <- compute_ratios(df.clean$raw.200, 'S', 'Ti')
ratio.df <- compute_ratios(ratio.df, 'Fe', 'Ti')
ratio.df <- compute_ratios(ratio.df, 'Mn', 'Ti')
ratio.df <- compute_ratios(ratio.df, 'Mn', 'Fe')
ratio.df <- compute_ratios(ratio.df, 'Si', 'Al')
ratio.df <- compute_ratios(ratio.df, 'Ba', 'Ti')
ratio.df <- compute_ratios(ratio.df, 'Ca', 'K')

# on 500µ data
ratio.df.500 <- compute_ratios(df.clean$raw.500, 'S', 'Ti')
ratio.df.500 <- compute_ratios(ratio.df.500, 'Fe', 'Ti')
ratio.df.500 <- compute_ratios(ratio.df.500, 'Mn', 'Ti')
ratio.df.500 <- compute_ratios(ratio.df.500, 'Mn', 'Fe')
ratio.df.500 <- compute_ratios(ratio.df.500, 'Si', 'Al')
ratio.df.500 <- compute_ratios(ratio.df.500, 'Ba', 'Ti')
ratio.df.500 <- compute_ratios(ratio.df.500, 'Ca', 'K')

# plot the ratios

#remove elements
ratio.df <- select(ratio.df, -c(2:16))
ratio.df.500 <- select(ratio.df.500, -c(2:16))
plot_individual_ratio(ratio.df, df_name = "ratios", output_dir = "plots/ratios/")
plot_individual_ratio(ratio.df.500, df_name = "ratios500", output_dir = "plots/ratios500/")


### ratios without log

r.df <- mutate(df.cs$raw.200, "S_Ti" = S/Ti)
r.raw.df <- mutate(df.clean$raw.200, 'S_Ti' = S/Ti)

plot_individual(r.df, df_name = 'no_log', output_dir = 'plots/ratio_nolog')
plot_individual(r.raw.df, df_name = 'no_log_raw', output_dir = 'plots/ratio_nolog_raw')



