# GRAYSCALE ANALYSIS

library(magick)
library(ggplot2)
library(dplyr)

source('xrf_functions.R')

# load image
current_dir <- getwd()
img_path <- paste(current_dir, '/data/bitmap.png', 
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


gray_profile <- gray_df %>%
  group_by(y) %>%
  summarize(mean_gray = mean(value))

core_length_cm <- 125.57-2.5  # Replace with your actual core length
gray_profile <- gray_profile %>%
  mutate(depth_cm = (y / max(y)) * core_length_cm +2.5)


ggplot(gray_profile, aes(x = depth_cm + 2.5, y = mean_gray)) +
  geom_line(linewidth = 0.05) +
  coord_flip() +
  scale_x_reverse() +  # Top of the core = top of the plot
  labs(title = "Grayscale Profile of Sediment Core",
       x = "Vertical Pixel Position",
       y = "Mean Grayscale Intensity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.position = "none",
        aspect.ratio = 5,
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave("grayscale_profile.png", width = 8, height = 10, dpi = 300)
