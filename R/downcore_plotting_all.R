# plot all elements downcore 

source('xrf_functions.R')
# prepare the data

# 1. Raw counts data on all elements (Zn and Ni included)
xrf_elem_sel_ZnNi <- c("Fe","Mn",'Ba',"S","Ti","Al","Si","K", 'Ca', 'Rb', 'Zr', 'Sr', 'P', 'Ba','Ni','Zn','position..mm.')
xrf_elem_sel <- c("Fe","Mn",'Ba',"S","Ti","Al","Si","K", 'Ca', 'Rb', 'Zr', 'Sr', 'P', 'Ba','position..mm.')
argon_sel <- c("Fe","Mn",'Ba',"S","Ti","Al","Si","K", 'Ca', 'Rb', 'Zr', 'Sr', 'P', 'Ba','Ni','Zn', 'Ar', 'position..mm.')

df.raw <- load_raw_data()
df.list.clean <- clean_df(df.raw, sel = xrf_elem_sel_ZnNi)
df.counts.200 <- df.list.clean$raw.200
df.counts.500 <- df.list.clean$raw.500

df.cs.200 <- closed_sum(df.counts.200)
df.cs.500 <- closed_sum(df.counts.500)

df.clr.200 <- clr_transform(df.counts.200)
df.clr.500 <- clr_transform(df.counts.500)

plot_elements(df.counts.200, df_name = '200_AC_raw', output_dir = 'plots/downcore_raw', xlab = 'cps')
plot_elements(df.counts.500, df_name = '500_AC_raw', output_dir = 'plots/downcore_raw', xlab = 'cps')
plot_elements(df.cs.200, df_name = '200_AC_CS', output_dir = 'plots/downcore_CS', xlab = 'closed-sum normalized intensity (%)')
plot_elements(df.cs.500, df_name = '500_AC_CS', output_dir = 'plots/downcore_CS', xlab = 'closed-sum normalized intensity (%)')
plot_elements(df.clr.200, df_name = '200_AC_CLR', output_dir = 'plots/downcore_CLR', xlab = 'Centered log-ratio')
plot_elements(df.clr.500, df_name = '500_AC_CLR', output_dir = 'plots/downcore_CLR', xlab = 'Centered log-ratio')


