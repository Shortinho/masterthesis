# grayscale but better

library(magick)
library(ggplot2)
library(dplyr)

source('xrf_functions.R')

# Load image
current_dir <- getwd()
img_path <- paste(current_dir, '/data/POS-22-20_highresruler copy-min copy.jpg', sep = "")
img <- image_read(img_path)

# Convert to grayscale and extract pixel data
img_gray <- image_convert(img, colorspace = "gray")
img_array <- image_data(img_gray, channels = "gray")

# Get image dimensions
width <- dim(img_array)[2]
height <- dim(img_array)[3]

# Convert to grayscale matrix and data frame
gray_vals <- as.integer(img_array[1, , ]) / 255  # Normalize
gray_df <- expand.grid(x = 1:width, y = 1:height)
gray_df$value <- as.vector(t(gray_vals))  # transpose to match

# Extract grayscale profile from center vertical band (±10 px)
center_band <- (width / 2 - 10):(width / 2 + 10)

gray_profile <- gray_df %>%
  filter(x %in% center_band) %>%
  group_by(y) %>%
  summarize(mean_gray = mean(value)) %>%
  ungroup()

# --- STEP 1: Identify white band y-positions (value == 1)
white_band_rows <- gray_profile %>%
  filter(mean_gray > 0.99) %>%
  pull(y)

# Group consecutive white rows together
white_marker_positions <- split(white_band_rows, cumsum(c(1, diff(white_band_rows) != 1)))

# Extract start (bottom of first white bar) and end (top of last white bar)
y_start <- max(white_marker_positions[[1]])  # last row of first white bar
y_end <- min(white_marker_positions[[length(white_marker_positions)]])  # first row of second white bar

# These image y positions correspond to 25 mm and 1231.6 mm in real depth
image_y <- c(y_start, y_end)
known_depth_mm <- c(25, 1231.6)

# --- STEP 2: Fit linear model to calibrate image y → depth_mm
depth_model <- lm(known_depth_mm ~ image_y)

# Predict depth for all y in gray_profile
gray_profile <- gray_profile %>%
  mutate(position..mm. = predict(depth_model, newdata = data.frame(image_y = y)))

# Keep only calibrated range
gray_profile <- gray_profile %>%
  filter(y > y_start+2, y < y_end-2)

# Plot calibrated grayscale profile
plot_individual_ratio(gray_profile, df_name = 'GrayScale_calibrated', output_dir = 'plots/grayscale')


# facies classification
library(mclust)

# check for gaussianity
ggplot(gray_profile, aes(x = mean_gray)) +
  geom_histogram(aes(y = ..density..), bins = 100, fill = "lightgray", color = "black") +
  geom_density(color = "blue", size = 1) +
  labs(title = "Histogram and Density of Grayscale Values", x = "Grayscale", y = "Density")

qqnorm(gray_profile$mean_gray)
qqline(gray_profile$mean_gray, col = "red")

shapiro.test(sample(gray_profile$mean_gray, 5000))  # limit to 5000 points to avoid errors
# result rejets single normal distribution


# Assume gray_profile is already trimmed and has 'mean_gray' and 'position..mm.'
gmm <- Mclust(gray_profile$mean_gray, G = 2)

# Assign class based on maximum probability
gray_profile$facies_class <- gmm$classification  # 1 or 2
tapply(gray_profile$mean_gray, gray_profile$facies_class, mean) # check means of both groups

