# statistics on combined dataset
source('xrf_functions.R')

library(ggplot2)
library(RColorBrewer)
library(dplyr)

combined <- read.csv('/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/generated/combined/combined_500_inc_coh.csv')

variance_sel <- c("Fe","Mn","S","Ti","Al","Si","K","Rb","Zr","Ca","Sr","P","Ba","Mn_Fe","S_Ti","Fe_Ti","Si_Ti","Ba_Ti", 'Ca_Ti',"Ti_Al",'Si_Al','Zr_Rb', "mean_gray", 'TChl.ug.g',"TNWt","TCWt","TSWt",'TCWt_TNWt', 'position..mm.')

# lagged correlation with grayscale
sel_vars <- c("Fe", "Mn", "S", "Ti", "Al", "K", "Si", "Ca", "Ba", "Zr", 
              "Mn_Fe", "Fe_Ti", "S_Ti", 'Si_Ti', "Ba_Ti", "Ca_Ti", "Ti_Al", "Zr_Rb", 
              "TCWt", 'TNWt', 'TSWt', "RABD673", "Moving.Average")

#sulfur/titanium
l.gs.STi <- ccf(combined$mean_gray, combined$S_Ti)
plot(l.gs.STi, main = 'Mean gray & S/Ti autocorrelation')
idx.max.autocorr.gs.STi <- l.gs.STi$lag[which.max(abs(l.gs.STi$acf))]

#sulfur
l.gs.S <- ccf(combined$mean_gray, combined$S)
plot(l.gs.S, main = 'Mean gray & S autocorrelation')
idx.max.autocorr.gs.S <- l.gs.S$lag[which.max(abs(l.gs.S$acf))]
#calcium
l.gs.Ca <- ccf(combined$mean_gray, combined$Ca)
plot(l.gs.Ca, main = 'Mean gray & Ca Autocorrelation')
idx.max.autocorr.gs.Ca <- l.gs.Ca$lag[which.max(abs(l.gs.Ca$acf))]
#titanium
l.gs.Ti <- ccf(combined$mean_gray, combined$Ti)
plot(l.gs.Ti, main = 'Mean gray & Ti Autocorrelation')
idx.max.autocorr.gs.Ti <- l.gs.Ti$lag[which.max(abs(l.gs.Ti$acf))]
# Iron
l.gs.Fe <- ccf(combined$mean_gray, combined$Fe)
plot(l.gs.Fe, main = 'Mean gray & Fe Autocorrelation')

idx.max.autocorr.gs.Fe <- l.gs.Fe$lag[which.max(abs(l.gs.Fe$acf))]

# Grayscale correlation w other vars ####
# Compute 10th and 90th percentiles
lower_decile <- quantile(combined$mean_gray, 0.10, na.rm = TRUE)
upper_decile <- quantile(combined$mean_gray, 0.90, na.rm = TRUE)

decile_df <- combined[combined$mean_gray <= lower_decile | 
                           combined$mean_gray >= upper_decile, ]
decile_df %>% 
       select(any_of(variance_sel)) %>%
       create_corr_map(title = 'Deciles of Mean_gray values')

# Filter for rows in the bottom or top decile
upper_df <- combined[combined$mean_gray >= upper_decile, ]
lower_df <- combined[combined$mean_gray <= lower_decile, ]





combined <- combined %>%
  mutate(decile_group = case_when(
    mean_gray <= lower_decile ~ "Lower decile",
    mean_gray >= upper_decile ~ "Upper decile",
    TRUE ~ "Middle"
  ))


# Step 2: Create segments between consecutive points
combined_segments <- combined %>%
  mutate(next_x = lead(position..mm.),
         next_y = lead(mean_gray),
         next_group = lead(decile_group)) %>%
  filter(!is.na(next_x))  # Remove last row (no next point)

# Step 3: Plot
col_pal <- brewer.pal(11, 'BrBG')

p <- ggplot(combined, aes(x = position..mm., y = mean_gray, color = decile_group)) +
  geom_point(size = 1.2, shape = 16) + 
  geom_line(aes(group = 1), color = "lightgrey", linewidth = 0.3)+
  scale_color_manual(values = c("Lower decile" = col_pal[1], 
                                "Middle" = "grey70", 
                                "Upper decile" = col_pal[4])) +
  labs(x = "Position (mm)", y = "Mean gray value", color = "Decile group") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8),
    strip.text = element_text(size = 8),
    legend.position = "none",
    aspect.ratio = 5,
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  coord_flip() +
  scale_x_reverse() +
  scale_y_reverse()

save_svg_plot(plot = p, filename = 'plots/grayscale_top/decile_grayscale', width = 4, height = 10)


# same but for 0.05 and 0.95

lower_quintile <- quantile(combined$mean_gray, 0.05, na.rm = TRUE)
upper_quintile <- quantile(combined$mean_gray, 0.95, na.rm = TRUE)

quintile_df <- combined[combined$mean_gray <= lower_quintile | 
                        combined$mean_gray >= upper_quintile, ]

quintile_df %>% 
       select(any_of(variance_sel)) %>%
       create_corr_map(title = '5th centile of Mean_gray values')

# Filter for rows in the bottom or top quintile
upper_df_5 <- combined[combined$mean_gray >= upper_quintile, ]
lower_df_5 <- combined[combined$mean_gray <= lower_quintile, ]

variance_sel <- c("Fe","Mn","S","Ti","Al","Si","K","Rb","Zr","Ca","Sr","P","Ba","Mn_Fe","S_Ti","Fe_Ti","Si_Ti","Ba_Ti", 'Ca_Ti',"Ti_Al",'Si_Al','Zr_Rb', "mean_gray", 'TChl.ug.g',"TNWt","TCWt","TSWt",'TCWt_TNWt', 'position..mm.')

upper_df_5 %>% 
  select(any_of(variance_sel)) %>%
  create_corr_map(title = 'Upper quintile of Mean_gray values')

lower_df_5 %>% 
  select(any_of(variance_sel)) %>%
  create_corr_map(title = 'Lower quintile of Mean_gray values')


combined <- combined %>%
  mutate(quintile_group = case_when(
    mean_gray <= lower_quintile ~ "Lower quintile",
    mean_gray >= upper_quintile ~ "Upper quintile",
    TRUE ~ "Middle"
  ))


# Step 2: Create segments between consecutive points
combined_segments <- combined %>%
  mutate(next_x = lead(position..mm.),
         next_y = lead(mean_gray),
         next_group = lead(quintile_group)) %>%
  filter(!is.na(next_x))  # Remove last row (no next point)

# Step 3: Plot
col_pal <- brewer.pal(11, 'BrBG')

p5 <- ggplot(combined, aes(x = position..mm., y = mean_gray, color = quintile_group)) +
  geom_point(size = 1.2, shape = 16) + 
  geom_line(aes(group = 1), color = "lightgrey", linewidth = 0.3)+
  scale_color_manual(values = c("Lower quintile" = col_pal[1], 
                                "Middle" = "grey70", 
                                "Upper quintile" = col_pal[4])) +
  labs(x = "Position (mm)", y = "Mean gray value", color = "quintile group") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8),
    strip.text = element_text(size = 8),
    legend.position = "none",
    aspect.ratio = 5,
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  coord_flip() +
  scale_x_reverse() +
  scale_y_reverse()

save_svg_plot(plot = p5, filename = 'plots/grayscale_top/quintile_grayscale', width = 4, height = 10)


# Compute density of the variable
dens <- density(combined$mean_gray)
densdf <- data.frame(x = dens$x, y = dens$y)

# Compute 5th and 95th percentiles of the raw data
quantiles <- quantile(combined$mean_gray, probs = c(0.1, 0.9))

# Add group info based on x (i.e., where in the distribution it is)
densdf$quantile <- cut(
  densdf$x,
  breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
  labels = c("bottom 10%", "middle", "top 10%"),
  include.lowest = TRUE
)

# Plot
ggplot(densdf, aes(x = x, y = y)) +
  geom_line() +
  geom_ribbon(aes(ymin = 0, ymax = y, fill = quantile)) +
  scale_fill_manual(values = c("bottom 10%" = col_pal[1], "middle" = "grey90", "top 10%" = col_pal[4])) +
  xlab('mean gray') + 
  ylab('density') +
  theme_minimal()


#### Check for differences in variance between top and bottom section
nopos_bottom <- df.clr %>% filter(position..mm.>326) %>% select(-position..mm.) 
nopos_top <- df.clr %>% filter(position..mm.<326) %>% select(-position..mm.) 

library(car)

# Run Levene's test for each numeric column
results_levene <- lapply(names(nopos), function(col) {
  if (is.numeric(nopos_bottom[[col]]) && is.numeric(nopos_top[[col]])) {
    df <- data.frame(
      value = c(nopos_bottom[[col]], nopos_top[[col]]),
      group = rep(c("nopos", "nopos_top"), c(nrow(nopos_bottom), nrow(nopos_top)))
    )
    test <- leveneTest(value ~ group, data = df, center = median) # Brown-Forsythe (more robust)
    data.frame(
      column = col,
      p.value = test$`Pr(>F)`[1]
    )
  } else {
    NULL
  }
})

# Combine results
results_levene_df <- do.call(rbind, results_levene)
print(results_levene_df)



#### Centers of facies -  ####

upper_decile <- quantile(combined$mean_gray, upper_q, na.rm = TRUE)
lower_decile <- quantile(combined$mean_gray, lower_q, na.rm = TRUE)

gs <- combined$mean_gray

mask_upper <- gs >= upper_decile
mask_lower <- gs <= lower_decile

clr_sel <- c('Fe.y', 'Mn.y', 'S.y', 'Ti.y', 'Al.y', 'Si.y', 'K.y', 'Rb.y', 'Zr.y', 'Ca.y', 'Sr.y', 'P.y', 'Ba.y')

combined.scaled <- scale(combined) %>% data.frame()

# --- Compute means ---
centers.upper <- combined.scaled %>%
  filter(mask_upper) %>%
  summarise(across(all_of(clr_sel), mean, na.rm = TRUE)) %>%
  mutate(group = "Upper")

centers.lower <- combined.scaled %>%
  filter(mask_lower) %>%
  summarise(across(all_of(clr_sel), mean, na.rm = TRUE)) %>%
  mutate(group = "Lower")

# --- Pivot to long format ---
centers.upper.long <- centers.upper %>%
  pivot_longer(-group, names_to = "variable", values_to = "center_value")

centers.lower.long <- centers.lower %>%
  pivot_longer(-group, names_to = "variable", values_to = "center_value")

# --- Define x limits (optional: same for both plots) ---
x_limits <- range(c(centers.upper.long$center_value, centers.lower.long$center_value))

# --- Plot Upper ---
p_upper <- ggplot(centers.upper.long, aes(x = center_value, y = fct_rev(variable))) +
  geom_col(fill = "#DFC27D") +
  scale_x_continuous(limits = x_limits) +
  labs(title = "Upper Decile Means", x = "Mean (Scaled)", y = "Variable") +
  theme_minimal()

# --- Plot Lower ---
p_lower <- ggplot(centers.lower.long, aes(x = center_value, y = fct_rev(variable))) +
  geom_col(fill = "#543005") +
  scale_x_continuous(limits = x_limits) +
  labs(title = "Lower Decile Means", x = "Mean (Scaled)", y = "Variable") +
  theme_minimal()

p_upper
p_lower


