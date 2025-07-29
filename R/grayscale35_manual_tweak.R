# grayscale analysis for top portion based on grayscale_TOP35 script
# in this script in take one of the classifications from the GMM and tweak it in order to manually refine the selection

library(magick)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

source('xrf_functions.R')

# Load image
current_dir <- getwd()
img_path <- paste(current_dir, '/data/POS-22-20_highresruler copy-min_35mm.jpg', sep = "")
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