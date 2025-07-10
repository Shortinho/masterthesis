source('xrf_functions.R')

# from correlation map I already remove highly correlated variables (r>0.9)
# S and S/Ti --> keep S/Ti
# Fe and Fe/Ti --> keep Fe/Ti
# RABD and Moving average --> keep Moving.Average
# TNWt and TCWt --> keep TCWt


vars_not_xrf <- c('position..mm.', 'mean_gray', 'RABD673', 'Moving.Average', 
                  'TNWt', 'TCWt', 'TSWt', 'facies_class', 'TCWt_TNWt')




# Iteration 1
no_xrf_combined <- create_combined_df()

xrf_vars_kept_1 <- c("Fe","Mn","S","Ti","Al","Si","K","Rb","Ca","Ba",'Zr',
                     "Mn_Fe","S_Ti","Fe_Ti","Si_Ti","Ba_Ti", 'Ca_Ti',"Ti_Al",'Si_Al','Zr_Rb')
other_vars_1 <- c('Moving.Average', 'TCWt', 'position..mm.')

ratios_1 <- list(
  c('Mn', 'Fe'), c('S', 'Ti'), c('Fe', 'Ti'),
  c('Si', 'Ti'), c('Ba', 'Ti'), c('Ca', 'Ti'),
  c('Ti', 'Al'), c('Si', 'Al'), c('Zr', 'Rb')
)

run_pca_selection_iteration(xrf_vars_kept_1, ratios_1, other_vars_1, no_xrf_combined, "1")


# iteration 2
xrf_vars_kept_2 <- c(
  "Ti", "Ba", "Zr", 
  "Mn_Fe", "S_Ti", "Fe_Ti", "Si_Ti", "Ba_Ti", "Ca_Ti", 
  "Ti_Al", "Zr_Rb"
)

other_vars_2 <- c("TCWt", "position..mm.")

ratios_2 <- list(
  c("Mn", "Fe"), c("S", "Ti"), c("Fe", "Ti"),
  c("Si", "Ti"), c("Ba", "Ti"), c("Ca", "Ti"),
  c("Ti", "Al"), c("Zr", "Rb")
)

run_pca_selection_iteration(
  xrf_vars_kept = xrf_vars_kept_2,
  ratio_list = ratios_2,
  other_vars = other_vars_2,
  no_xrf_combined = no_xrf_combined,
  iteration_name = "2"
)

# iteration 3
xrf_vars_kept_3 <- c(
  "Mn_Fe", "S_Ti", "Fe_Ti", "Si_Ti", 
  "Ca_Ti", "Ti_Al"
)

other_vars_3 <- c("position..mm.")

ratios_3 <- list(
  c("Mn", "Fe"), c("S", "Ti"), c("Fe", "Ti"),
  c("Si", "Ti"), c("Ca", "Ti"), c("Ti", "Al")
)

run_pca_selection_iteration(
  xrf_vars_kept = xrf_vars_kept_3,
  ratio_list = ratios_3,
  other_vars = other_vars_3,
  no_xrf_combined = no_xrf_combined,
  iteration_name = "3"
)


# iteration 4 (based on Boruta importance)
xrf_vars_kept_4 <- c(
  "Mn_Fe", 'Ba', 'Ti', 'Al', 'K', 'Si', 'Fe', "S_Ti", "Fe_Ti", 'S', 'TCWt', 
  "Ca_Ti", "Ti_Al"
)

other_vars_4 <- c("position..mm.")

ratios_4 <- list(
  c("Mn", "Fe"), c("S", "Ti"), c("Fe", "Ti"), c("Ca", "Ti"), c("Ti", "Al")
)

run_pca_selection_iteration(
  xrf_vars_kept = xrf_vars_kept_4,
  ratio_list = ratios_4,
  other_vars = other_vars_4,
  no_xrf_combined = no_xrf_combined,
  iteration_name = "4"
)

# iteration 5 (based on 4)
xrf_vars_kept_5 <- c(
  "Mn_Fe", 'Ba', 'Ti', 'K', 'Si', "S_Ti", "Fe_Ti", 'TCWt', 
  "Ca_Ti", "Ti_Al"
)

other_vars_5 <- c("position..mm.")

ratios_5 <- list(
  c("Mn", "Fe"), c("S", "Ti"), c("Fe", "Ti"), c("Ca", "Ti"), c("Ti", "Al")
)

run_pca_selection_iteration(
  xrf_vars_kept = xrf_vars_kept_5,
  ratio_list = ratios_5,
  other_vars = other_vars_5,
  no_xrf_combined = no_xrf_combined,
  iteration_name = "5"
)



