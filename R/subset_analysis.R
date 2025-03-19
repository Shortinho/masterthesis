# analysis with subset of elements

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

# Select elements
all_elements <- c("Al", "Si", "P", "S", "K", "Ca", "Ti", "Mn", "Fe", "Ni", "Zn", "Rb", "Sr", "Zr", "Ba")
selection <- c("Al", "Si", "P", "S", "K", "Ca", "Ti", "Mn", "Fe", "Ni", "Zn", "Rb", "Zr")
df_l_slct_elmnt <- map(df_l, ~select_elements(.x, all_elements))