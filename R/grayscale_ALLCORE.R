# grayscale analysis for all core portion based on grayscale_TOP35


library(magick)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(mclust)

source('xrf_functions.R')

# Load image
current_dir <- getwd()
img_path <- paste(current_dir, '/data/POS-22-20_highresruler copy-min GS.jpg', sep = "") # image is white before 25mm and after 1231.6mm
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

b_wdth_vector <- c(5, 10, 15, 20, 30, 40, 50, 60, 80, 100)
# b_wdth_vector <- c(50)
for (band_width in b_wdth_vector) {
  
  print(paste('bandwitdh = ', band_width))
  # Extract grayscale profile from center vertical band (±20 px)
  center_band <- floor(width / 2 - band_width / 2):ceiling(width / 2 + band_width / 2)  
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
  y_end <- min(white_marker_positions[[2]])  # first row of second white bar
  
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
  # plot_individual_ratio(gray_profile, df_name = 'GrayScale_calibrated_top35', output_dir = 'plots/grayscale_top')
  
  
  # facies classification
  ### GMM method ###
  
  # # check for gaussianity #####
  # ggplot(gray_profile, aes(x = mean_gray)) +
  #   geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "lightgray", color = "black") +
  #   geom_density(color = "blue", linewidth = 1) +
  #   labs(title = "Histogram and Density of Grayscale Values", x = "Grayscale", y = "Density")
  # 
  # qqnorm(gray_profile$mean_gray)
  # qqline(gray_profile$mean_gray, col = "red")
  # 
  # shapiro.test(sample(gray_profile$mean_gray, 5000))  # limit to 5000 points to avoid errors
  # # result rejets single normal distribution
  
  # --- GMM classification --- ####
  groups <- 2
  gmm <- Mclust(gray_profile$mean_gray, G = groups)
  gray_profile$facies_class <- gmm$classification
  gray_profile$facies_prob <- gmm$z[, 1]
  
  # Sort by depth just in case
  gray_profile <- gray_profile %>%
    arrange(position..mm.)
  
  # Store original facies
  facies_class <- gray_profile$facies_class
  depth_vec <- gray_profile$position..mm.
  
  # RLE to identify continuous segments
  rle_facies <- rle(facies_class)
  lengths <- rle_facies$lengths
  values <- rle_facies$values
  
  # Find start and end indices of each run
  end_idx <- cumsum(lengths)
  start_idx <- c(1, head(end_idx, -1) + 1)
  
  # Reassign small zones (< 0.2 mm i.e target thickness)
  for (i in seq_along(values)) {
    start_i <- start_idx[i]
    end_i <- end_idx[i]
    thickness <- abs(depth_vec[end_i] - depth_vec[start_i])
    
    target_thickness <- 0.2
    
    if (thickness < target_thickness) {
      # Reassign to neighbor class
      if (i == 1) {
        # At beginning: assign to next block
        facies_class[start_i:end_i] <- values[i + 1]
      } else if (i == length(values)) {
        # At end: assign to previous block
        facies_class[start_i:end_i] <- values[i - 1]
      } else {
        # Between two blocks: choose the neighbor with thicker zone
        prev_thick <- abs(depth_vec[end_idx[i - 1]] - depth_vec[start_idx[i - 1]])
        next_thick <- abs(depth_vec[end_idx[i + 1]] - depth_vec[start_idx[i + 1]])
        if (next_thick >= prev_thick) {
          facies_class[start_i:end_i] <- values[i + 1]
        } else {
          facies_class[start_i:end_i] <- values[i - 1]
        }
      }
    }
  }
  
  # Update gray_profile
  gray_profile$facies_class_cleaned <- facies_class
  
  tapply(gray_profile$mean_gray, gray_profile$facies_class, mean) # check means of both groups
  
  ### Couplet count ####
  facies_clean <- gray_profile$facies_class_cleaned
  
  # Get transitions between successive values
  facies_transitions <- diff(facies_clean)
  
  # Identify where a change occurred (±1)
  transition_points <- which(facies_transitions != 0)
  
  # Number of transitions
  num_transitions <- length(transition_points)
  
  # Number of varve couplets = transitions / 2
  num_couplets <- floor(num_transitions / 2)
  
  cat("Estimated number of varve couplets:", num_couplets, "\n")
  
  ### Otsu method ####
  # library(autothresholdr)
  # 
  # threshold <- auto_thresh(gray_profile$mean_gray, method = "Otsu")
  # 
  # gray_profile$facies_class_otsu <- ifelse(gray_profile$mean_gray < threshold, "dark", "light")
  # 
  # 
  # ### K-means method ####
  # set.seed(42)
  # kmeans_result <- kmeans(gray_profile$mean_gray, centers = 2)
  # 
  # gray_profile$facies_class_kmeans <- kmeans_result$cluster
  
  # write df ####
  write.csv(gray_profile, file = paste0("data/generated/grayscale/gray_profile_ALL_bw", band_width, 'groups_', groups,'filter_',target_thickness,'mm', ".csv"), row.names = FALSE)
  # plot ####  
  p <- ggplot(gray_profile, aes(x = position..mm., y = mean_gray)) +
    geom_line() +
    geom_vline(xintercept = c(25, 1231.6), linetype = "dashed") +
    coord_flip() +
    scale_x_reverse() +
    scale_y_reverse() +
    labs(title = "Grayscale Profile with GMM-based Facies", x = "Depth (mm)", y = "Grayscale")
  
  f_name <- paste0('plots/grayscale/grayscale_w_facies_class', band_width, '.svg')
  
  #save_svg_plot(p, f_name, width = 4, height = 10)
  col_pal <- brewer.pal(11, 'BrBG')
  p2 <- ggplot(gray_profile, aes(x = position..mm., y = mean_gray, color = factor(facies_class_cleaned))) +
    geom_point(size = 0.1) +
    geom_vline(xintercept = c(25, 1231.6), linetype = "dashed") +
    coord_flip() +
    scale_x_reverse() +
    scale_y_reverse() +
    scale_color_manual(values = c('1' = col_pal[1], '2' = col_pal[5])) +
    labs(title = paste("Grayscale Profile with GMM-based Facies, bandwidth = ", band_width), x = "Depth (mm)", y = "Grayscale") +
    theme_bw()
  
  f_name2 <- paste0('plots/grayscale/grayscale_w_facies_class_', band_width, 'BndWdth.svg')
  save_svg_plot(p2, f_name2, width = 4, height = 10)
  
  
  gray_segments <- gray_profile %>%
    arrange(position..mm.) %>%
    mutate(xend = lead(position..mm.),
           yend = lead(mean_gray),
           group = facies_class_cleaned) %>%
    filter(!is.na(xend))  # Remove last row (incomplete segment)
  
  p_line <- ggplot(gray_segments) +
    geom_segment(aes(x = position..mm., y = mean_gray,
                     xend = xend, yend = yend,
                     color = factor(group)),
                 linewidth = 0.3) +
    scale_color_manual(values = c('1' = col_pal[1], '2' = col_pal[5])) +
    coord_flip() +
    scale_x_reverse() +
    scale_y_reverse() +
    labs(title = paste("Grayscale Line with GMM-based Facies, bandwidth = ", band_width),
         x = "Depth (mm)", y = "Grayscale") +
    theme_bw()
  
  # Save
  f_name <- paste0('plots/grayscale/grayscale_line_w_facies_class_', band_width, 'BndWdth_groups', groups, '.svg')
  save_svg_plot(p_line, f_name, width = 4, height = 10)
  
  
  facies_bar <- gray_profile %>%
    arrange(position..mm.) %>%
    mutate(depth_top = position..mm.,
           depth_bottom = lead(position..mm.),
           group = facies_class_cleaned) %>%
    filter(!is.na(depth_bottom))  # Remove last row
  
  p_bar <- ggplot(facies_bar) +
    geom_rect(aes(ymin = depth_top, ymax = depth_bottom,
                  xmin = 0, xmax = 1, fill = factor(group))) +
    scale_fill_manual(values = c('1' = col_pal[1], '2' = col_pal[5])) +
    scale_y_reverse() +
    labs(
      title = paste("Facies Classification, bandwidth =", band_width),
      fill = "Facies",
      y = "Depth (mm)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      aspect.ratio = 3
    )
  
  f_name <- paste0('plots/grayscale/grayscale_bar', band_width, 'BndWdth_2mmfilter_wAxis.svg')
  save_svg_plot(p_bar, f_name, width = 4, height = 10)
  
}

