# plot all elements downcore 

source('xrf_functions.R')
# prepare the data

# 1. Raw counts data on all elements (Zn and Ni included)
xrf_elem_sel_ZnNi <- c("Fe","Mn",'Ba',"S","Ti","Al","Si","K", 'Ca', 'Rb', 'Zr', 'Sr', 'P', 'Ba','Ni','Zn','position..mm.')
xrf_elem_sel <- c("Fe","Mn",'Ba',"S","Ti","Al","Si","K", 'Ca', 'Rb', 'Zr', 'Sr', 'P', 'Ba','position..mm.')

df.raw <- load_raw_data()
df.list.clean <- clean_df(df.raw, sel = xrf_elem_sel_ZnNi)
df.counts.200 <- df.list.clean$raw.200
df.counts.500 <- df.list.clean$raw.500

plot_elements(df.counts.200, df_name = 'AC_200_raw', output_dir = 'plots/downcore_raw')
