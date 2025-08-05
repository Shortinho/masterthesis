### Element selection for multivariable analysis script

source('xrf_functions.R')

# combined contains all the variables in one data frame
combined <- combined_df_creator()

# next step is to make a df that contains all variables except the xrf elements and ratios
vars_not_xrf <- c('position..mm.', 'TCWt')

no_xrf_combined <- combined %>%
  select(any_of(vars_not_xrf))

# select variables for first PCA: here, all potential variables are included, in the following steps, we filter out variables that are not relevant i.e. keeping the smallest amount of variables that explain the environmental process

variance_sel <- c("Fe","Mn","S","Ti","Al","Si","K","Rb","Zr","Ca","Sr","P","Ba","Mn_Fe","S_Ti","Fe_Ti","Si_Ti","Ba_Ti", 'Ca_Ti',"Ti_Al",'Si_Al','Zr_Rb', "mean_gray","RABD673",'Moving.Average',"TNWt","TCWt","TSWt",'TCWt_TNWt', 'position..mm.', 'facies_class')




combined %>%
  select(any_of(variance_sel)) %>%
  perform_pca2(plot_loadings = T,
               df_name = paste('n_',length(variance_sel),
                               collapse = ', '),
               output_dir = 'plots/PCA_variable_selection')

# from this first loadings plot, we keep only TCWt, S, Fe, S_Ti, Fe_Ti, Ca_Ti, Ba_Ti, Mn_Fe, Ti_Al, Moving.Average (HSI), Ba, Zr_Rb, Ti, K, Si, Al, Si_Ti

xrf_vars_kept <- c("Fe","Mn","S","Ti","Al","Si","K","Rb","Ca","Ba",'Zr',"Mn_Fe","S_Ti","Fe_Ti","Si_Ti","Ba_Ti", 'Ca_Ti',"Ti_Al",'Si_Al','Zr_Rb')
others_kept <- c('Moving.Average', 'TCWt', 'position..mm.')

xrf_elem1 <- c("Fe","Mn","S","Ti","Al","Si","K","Rb","Ca","Ba",'Zr', 'position..mm.')

# load xrf df with updated element selection
df.raw <- load_raw_data()
df.xrf <- clean_df(df.raw, sel = xrf_elem1)
df.xrf$raw.200 <- NULL
df.xrf <- tibble(df.xrf$raw.500) %>%
  select_top_varves()
df.cs <- closed_sum(df.xrf) 

# join to other data
combined1 <- df.cs %>%
  left_join(no_xrf_combined, by = "position..mm.")

combined1 <- compute_ratios(combined1, 'Mn', 'Fe')
combined1 <- compute_ratios(combined1, 'S', 'Ti')
combined1 <- compute_ratios(combined1, 'Fe', 'Ti')
combined1 <- compute_ratios(combined1, 'Si', 'Ti')
combined1 <- compute_ratios(combined1, 'Ba', 'Ti')
combined1 <- compute_ratios(combined1, 'Ca', 'Ti')
combined1 <- compute_ratios(combined1, 'Ti', 'Al')
combined1 <- compute_ratios(combined1, 'Si', 'Al')
combined1 <- compute_ratios(combined1, 'Zr', 'Rb')


new.sel <- c(xrf_vars_kept, others_kept)
combined1 %>%
  select(any_of(new.sel)) %>%
  perform_pca2(plot_loadings = T,
               df_name = paste('n_',length(new.sel),
                               collapse = ', '),
               output_dir = 'plots/PCA_variable_selection')

# second iteration
xrf_vars_kept <- c("Ti","Ba",'Zr',"Mn_Fe","S_Ti","Fe_Ti","Si_Ti","Ba_Ti", 'Ca_Ti',"Ti_Al",'Zr_Rb')
others_kept <- c('TCWt', 'position..mm.')
xrf_elem2 <- c("Fe","Mn","S","Ti","Al","Si","K","Rb","Ca","Ba",'Zr', 'position..mm.')
df.xrf <- clean_df(df.raw, sel = xrf_elem2)
df.xrf$raw.200 <- NULL
df.xrf <- tibble(df.xrf$raw.500) %>%
  select_top_varves()
df.cs <- closed_sum(df.xrf) 

# join to other data
combined2 <- df.cs %>%
  left_join(no_xrf_combined, by = "position..mm.")

combined2 <- compute_ratios(combined2, 'Mn', 'Fe')
combined2 <- compute_ratios(combined2, 'S', 'Ti')
combined2 <- compute_ratios(combined2, 'Fe', 'Ti')
combined2 <- compute_ratios(combined2, 'Si', 'Ti')
combined2 <- compute_ratios(combined2, 'Ba', 'Ti')
combined2 <- compute_ratios(combined2, 'Ca', 'Ti')
combined2 <- compute_ratios(combined2, 'Ti', 'Al')
combined2 <- compute_ratios(combined2, 'Zr', 'Rb')

new.sel2 <- c(xrf_vars_kept, others_kept)
combined2 %>%
  select(any_of(new.sel2)) %>%
  perform_pca2(plot_loadings = T,
               df_name = paste('n_',length(new.sel2),
                               collapse = ', '),
               output_dir = 'plots/PCA_variable_selection')


# third iteration
xrf_vars_kept <- c("Mn_Fe","S_Ti","Fe_Ti","Si_Ti", 'Ca_Ti',"Ti_Al")
others_kept <- c('position..mm.')
xrf_elem2 <- c("Fe","Mn","S","Ti","Al","Si","Ca", 'position..mm.')
df.xrf <- clean_df(df.raw, sel = xrf_elem2)
df.xrf$raw.200 <- NULL
df.xrf <- tibble(df.xrf$raw.500) %>%
  select_top_varves()
df.cs <- closed_sum(df.xrf) 

# join to other data
combined2 <- df.cs %>%
  left_join(no_xrf_combined, by = "position..mm.")

combined2 <- compute_ratios(combined2, 'Mn', 'Fe')
combined2 <- compute_ratios(combined2, 'S', 'Ti')
combined2 <- compute_ratios(combined2, 'Fe', 'Ti')
combined2 <- compute_ratios(combined2, 'Si', 'Ti')
combined2 <- compute_ratios(combined2, 'Ba', 'Ti')
combined2 <- compute_ratios(combined2, 'Ca', 'Ti')
combined2 <- compute_ratios(combined2, 'Ti', 'Al')
combined2 <- compute_ratios(combined2, 'Zr', 'Rb')

new.sel2 <- c(xrf_vars_kept, others_kept)
combined2 %>%
  select(any_of(new.sel2)) %>%
  perform_pca2(plot_loadings = T,
               df_name = paste('n_',length(new.sel2),
                               collapse = ', '),
               output_dir = 'plots/PCA_variable_selection')
