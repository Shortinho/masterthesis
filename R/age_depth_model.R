#Plotting the age--depth model data provided by Adrianus Damanik

# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)

# Read the data (automatically handles tabs or multiple spaces)
data <- read_table("data/age_depth/POS-22-20_with selection_124_ages.txt", col_names = TRUE)

# Check the structure
str(data)

# Plot whole core
ggplot(data, aes(x = mean, y = depth)) +
  geom_ribbon(aes(xmin = min, xmax = max), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue", size = 1) +
  scale_y_reverse(breaks = seq(0, 130, by = 10)) +  # Depth increases downward
  scale_x_reverse() +
  
  labs(
    x = "Age (years)",
    y = "Depth (cm)",
    title = "Age-Depth Model",
    subtitle = "Mean age with min-max confidence interval"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 2)

data <- tibble(data)

data_top <- data %>% 
  tibble() %>%
  filter(depth <= 35)

# Plot top 35cm

ggplot(data_top, aes(x = mean, y = depth)) +
  geom_ribbon(aes(xmin = min, xmax = max), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue", size = 1) +
  scale_y_reverse(breaks = seq(0, 35, by = 5)) +  # Depth increases downward
  scale_x_reverse() +
  labs(
    x = "Age (years)",
    y = "Depth (cm)",
    title = "Age-Depth Model",
    subtitle = "Mean age with min-max confidence interval"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 2)


