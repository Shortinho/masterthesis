# GRAYSCALE ANALYSIS

library(magick)
library(ggplot2)
library(dplyr)

source('xrf_functions.R')

# load image
current_dir <- getwd()
img_path <- paste(current_dir, '/data/POS-22-20_highresruler copy-min copy.jpg', 
                  sep = "")
img <- image_read(img_path)
# Convert to grayscale
img_gray <- image_convert(img, colorspace = "gray")
img_array <- image_data(img_gray, channels = "gray")

# Get image dimensions
width <- dim(img_array)[2]
height <- dim(img_array)[3]

# Convert to grayscale matrix and data frame
gray_vals <- as.integer(img_array[1, , ]) / 255  # Normalize from 0–255 to 0–1
gray_df <- expand.grid(x = 1:width, y = 1:height)
gray_df$value <- as.vector(t(gray_vals))  # transpose to match orientation


# gray_profile <- gray_df %>%
#   group_by(y) %>%
#   summarize(mean_gray = mean(value))
# Define center band (±10 pixels around center)
center_band <- (width / 2 - 10):(width / 2 + 10)

# Filter to center band and compute mean grayscale at each depth (y)
gray_profile <- gray_df %>%
  filter(x %in% center_band) %>%
  group_by(y) %>%
  summarize(mean_gray = mean(value)) %>%
  ungroup()



core_length_mm <- 1231.6-25

gray_profile <- gray_profile %>% 
  mutate(depth_mm = ((y-1) / max(y)) * core_length_mm)

# find where to cut
gray_select <- gray_profile[560:26200,]

# gray_select <- gray_select %>% 
#   mutate(depth_mm = (y / max(y)) * core_length_mm +2.5)
# 
gray_select <- rename(gray_select, position..mm. = depth_mm)
plot_individual_ratio(gray_select, df_name = 'GrayScale_clean', output_dir = 'plots/grayscale') #remember to scale_y_reverse

gray_profile <- rename(gray_profile, position..mm. = depth_mm)
plot_individual_ratio(gray_profile, df_name = 'GrayScale', output_dir = 'plots/grayscale') #remember to scale_y_reverse

ggplot(gray_profile, aes(x = depth_mm, y = mean_gray)) +
  geom_line(linewidth = 0.25) +
  coord_flip() +
  scale_x_reverse() +  # Top of the core = top of the plot
  labs(title = "Grayscale Profile of Sediment Core",
       x = "Vertical Pixel Position",
       y = "Mean Grayscale Intensity") +
  ylim(c(0.05, 0.4))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, hjust = 0),
        legend.position = "none",
        aspect.ratio = 5,
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave("grayscale_profile_7_pixels.png", width = 3, height = 10, dpi = 300)
