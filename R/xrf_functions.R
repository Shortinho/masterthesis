# preliminary data exploration functions

library(dplyr)
library(ggplot2)
library(tidyr)
library(corrplot)
library(purrr)
library(svglite)

load_data <- function(){
  library(readxl)
  #### Import data ####
  script_path <- getwd()
  xrf_data_path <- paste(script_path, "/data/XRF Data", sep="")
  # xrf_data_path <- "/Users/maxshore/Documents/Unibe/MasterThesis/R/data/XRF Data"
  
  f_norm200mu <- read_excel(paste(xrf_data_path,
                                  "/POS-22-20_0.2mm_Cr_Results/POS-22-20_0.2mm_Cr_normZnNi.xlsx", 
                                  sep = ""),
                            skip=2) %>%
    tibble()
  f_norm500mu <- read_excel(paste(xrf_data_path,
                                  "/POS-22-20_0.5mm_Cr_Results/POS-22-20_0.5mm_Cr_normZnNi.xlsx",
                                  sep = ""),
                            skip=2) %>%
    tibble()
  f_clr200mu <- read_excel(paste(xrf_data_path,
                                 "/POS-22-20_0.2mm_Cr_Results/POS-22-20_0.2mm_Cr_clrZnNi.xlsx",
                                 sep = ""),
                           skip=2) %>%
    tibble()
  f_clr500mu <- read_excel(paste(xrf_data_path,
                                 "/POS-22-20_0.5mm_Cr_Results/POS-22-20_0.5mm_Cr_clrZnNi.xlsx",
                                 sep = ""),
                           skip=2) %>%
    tibble()
  
  # Create a list of dataframes
  df_list <- list(norm200 = f_norm200mu,
                  norm500 = f_norm500mu,
                  clr200 = f_clr200mu,
                  clr500 = f_clr500mu)
}
load_raw_data <-function(){
  library(readxl)
  #### Import data ####
  script_path <- getwd()
  xrf_data_path <- paste(script_path, "/data/XRF Data", sep="")
  
  raw_200 <- read_excel(paste(xrf_data_path,
                              "/POS-22-20_0.2mm_Cr_Results/POS-22-20_0.2mm_Cr_slct_ZnNi.xlsx", 
                              sep = ""),
                        skip=2) %>%
    tibble()
  
  raw_500 <- read_excel(paste(xrf_data_path,
                              "/POS-22-20_0.5mm_Cr_Results/POS-22-20_0.5mm_Cr_slct_ZnNi.xlsx", 
                              sep = ""),
                        skip=2) %>%
    tibble()
  df.r <- list(raw.200 = raw_200,
               raw.500 = raw_500)
  return(df.r)
}

# Function to clean dataframes
clean_df <- function(df_list, sel = NULL) {
  library(purrr)
  library(dplyr)
  if (is.null(sel)){
    sel = c('position..mm.','Al','Si','P','S','K','Ca','Ti','Mn','Fe','Ni','Zn','Rb','Sr','Zr','Ba')}
  # make all column names universal (remove spaces)
  df_list <- purrr::map(df_list, ~{
    names(.x) <- make.names(names(.x), unique = TRUE)
    .x
  })
  # keep only desired data points (remove non-sediment measurements)
  df_pos <- purrr::map(df_list, ~dplyr::filter(.x, `position..mm.` > 25 & `position..mm.` < 1231.6))
  # remove unwanted variables
  cleaned_df <- purrr::map(df_pos, ~select(.x, any_of(sel)))
  return(cleaned_df)
}

clean_argon_artifact <- function(df){
  
  df_clean <- df %>%
    filter(
      !(
        (position..mm. >= 66 & position..mm. <= 91) |
          (position..mm. >= 290 & position..mm. <= 316) |
          (position..mm. >= 950 & position..mm. <= 976)
      )
    )
  return(df_clean)
}
# convert to long with only elements
make_long <- function(df){
  elements <- df %>%
    select(-c(position..mm.)) %>%
    names()
  df_long <- df %>%
    pivot_longer(cols = all_of(elements), names_to = "Variable", values_to = "Value")
  return(df_long)
}

# combine cluster to a long dataset
make_long_clust <- function(df){
  elements <- df %>%
    select(-c(position..mm., Cluster)) %>%
    names()
  df_long <- df %>%
    pivot_longer(cols = all_of(elements), names_to = "Variable", values_to = "Value")
  return(df_long)
}

#
closed_sum <- function(df_count) {
  # Separate the metadata (e.g., position) and numeric data
  meta_cols <- df_count %>% select(position..mm.)
  numeric_cols <- df_count %>% select(-position..mm.)
  
  # Convert to matrix for faster row-wise operations
  mat <- as.matrix(numeric_cols)
  row_sums <- rowSums(mat)
  
  # Avoid division by zero
  row_sums[row_sums == 0] <- NA
  
  mat_norm <- sweep(mat, 1, row_sums, FUN = "/")
  df_norm <- as.data.frame(mat_norm)
  
  # Recombine with metadata
  df_result <- bind_cols(meta_cols, df_norm)
  
  return(df_result)
}
# centered log ratio normalization. Input must be from df with count data
clr_transform <- function(df_count) {
  # Select only element values
  df_no_pos <- df_count %>% 
    select(-c('position..mm.'))
  
  # Compute geometric mean for each row using vectorized operations
  geometric_means <- exp(rowMeans(log(df_no_pos)))
  
  # Apply CLR transformation in a vectorized manner
  df_clr <- as.data.frame(log(df_no_pos / geometric_means))
  
  df_clr$position..mm. <- df_count$position..mm.
  
  # Reorder so 'position..mm.' is the first column
  df_clr <- df_clr %>% relocate(position..mm.)
  return(df_clr)
}

# function to only select top varves
select_top_varves <- function(df){
  library(dplyr)
  short_df <- df %>%
    dplyr::filter(`position..mm.` < 350)
}

select_bottom <- function(df){
  library(dplyr)
  short_df <- df %>%
    dplyr::filter(`position..mm.` > 350 & `position..mm.` < 1231.6)
}

# give a vector with elements desired and the df
select_elements <- function(df, elements) {
  print('The following elements are being kept:')
  cat(elements, sep=',')
  cat('\n')
  cat('\n')
  df %>%
    select("position..mm.", all_of(elements))
  
}


# Function to create downcore lineplots of all columns of df
plot_elements <- function(df, df_name = NULL, output_dir = NULL, xlab = 'value') {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(jpeg)
  
    elements <- df %>%
      select(-c(position..mm.)) %>%
      names()
    df_long <- df %>%
      pivot_longer(cols = all_of(elements), names_to = "Variable", values_to = "Value")
    
    # All elements in one figure
    p <- ggplot(df_long, aes(y = Value, x = position..mm.)) +
          geom_path(aes(group = Variable), linewidth=0.15) +
          facet_wrap(~ Variable, scales = "free_x", ncol = (length(elements))) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
                axis.text.y = element_text(size = 6),
                axis.title = element_text(size = 8),
                strip.text = element_text(size = 8),
                legend.position = "none",
                panel.background = element_rect(fill = 'white'),
                aspect.ratio = 8,
                panel.spacing = unit(0.1, 'pt'),
                panel.border = element_blank(),
                strip.background = element_rect(
                  color="white", fill="#c6dbef", linewidth=0.5, linetype="solid"
                )) +
          labs(title = paste("XRFnorm Analysis -", df_name),
               y = xlab,
               x = "Sample ID (depth)") +
          coord_flip() +
          scale_x_reverse()
    print(p)
    if(!is.null(output_dir)){
      ggsave(
        filename = file.path(output_dir, paste0("all_", df_name, ".png")),
        plot = p,
        width = 10, height = 8)}
}

# Plot individual elements
plot_individual <- function(df, df_name = NULL, output_dir = NULL){
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(jpeg)
  
  elements <- df %>%
    select(-c(position..mm.)) %>%
    names()
  df_long <- make_long(df)
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE)  # Create output directory if not exists
  }
  
  for (element in elements) {
    element_data <- df_long %>% dplyr::filter(Variable == element)
    
    element_plot <- ggplot(df_long %>% dplyr::filter(Variable == element), aes(x = position..mm., y = Value)) +
      geom_line(color = "black", linewidth = 0.2, alpha = 0.9) +
      labs(
        title = paste("Line Plot for", element, "in", df_name),
        x = "Position",
        y = element
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6),
            axis.title = element_text(size = 8),
            strip.text = element_text(size = 8),
            legend.position = "none",
            aspect.ratio = 5,
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)) +
      coord_flip() +
      scale_x_reverse()
    path_prfx <- getwd()
    path <- "/Users/maxshore/Documents/Unibe/MasterThesis/R/data/POS-22-20_highresruler_copy_min_crpd.jpg"
    
    # read the jpef file from device
    # img <- readJPEG(path, native = TRUE)
    
    # Extract axis limits from the built plot
    # plot_build <- ggplot_build(element_plot)
    # x_range <- plot_build$layout$panel_scales_x[[1]]$range$range
    # y_range <- plot_build$layout$panel_scales_y[[1]]$range$range
    # print(x_range)
    # print(paste(min(element_data$SampleID),max(element_data$SampleID)))
    # print(y_range)
    # print(paste(min(element_data$Value),max(element_data$Value)))
    # left <- (y_range[2]*1.05)/diff(y_range)
    # bottom <- -x_range[1]/diff(x_range)
    # right <- (y_range[2]*1.3)/diff(y_range)
    # top <- -x_range[2]/diff(x_range)
    # left <- 0.9
    # bottom <- 0.044
    # right <- 1.5
    # top <- 0.956
    # 
    # # Add inset image
    # img_element_plot <- element_plot +
    #   inset_element(
    #     p = img,
    #     left = left,
    #     bottom = bottom,
    #     right = right,
    #     top = top
    #   )
    
    # print(img_element_plot)
    
    #print(img_element_plot)
    #print(element_plot)
    # Save individual plots if output_dir is provided
    if (!is.null(output_dir)) {
      ggsave(
        filename = file.path(output_dir, paste0(df_name, "_", element, "_line_only.svg")),
        plot = element_plot,
        width = 4, height = 10,
        dpi = 800
        
      )
      save_svg_plot(element_plot,df_name,width = 4, height = 10)
    }
    
  }
}
plot_individual_top <- function(df, df_name = NULL, output_dir = NULL){
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(jpeg)
  library(png)
  library(cowplot)
  library(magick)
  library(svglite)
  
  elements <- df %>%
    select(-c(position..mm.)) %>%
    names()
  df_long <- make_long(df)
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE)  # Create output directory if not exists
  }
  
  for (element in elements) {
    element_data <- df_long %>% dplyr::filter(Variable == element)
    print(df_long%>%head())
    element_plot <- ggplot(df_long %>% dplyr::filter(Variable == element), aes(x = position..mm., y = Value)) +
      geom_line(color = "black", linewidth = 0.1, alpha = 0.9) +
      labs(
        title = paste("Line Plot for", element, "in", df_name),
        x = "Position",
        y = element
      ) +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6),
            axis.title = element_text(size = 8),
            strip.text = element_text(size = 8),
            legend.position = "none",
            aspect.ratio = 5) +
      coord_flip() +
      scale_x_reverse()
    path <- "/Users/maxshore/Documents/Unibe/MasterThesis/R/data/POS-22-20_highresruler_copy_min_crpd_top.jpg"
    

    # read the jpef file from device
    # img <- readJPEG(path, native = TRUE)
    # 
    # #Extract axis limits from the built plot
    # plot_build <- ggplot_build(element_plot)
    # x_range <- plot_build$layout$panel_scales_x[[1]]$range$range
    # y_range <- plot_build$layout$panel_scales_y[[1]]$range$range
    # print(x_range)
    # print(y_range)
    # print(paste(min(element_data$Value),max(element_data$Value)))
    # 
    # # top <- (y_range[2]*1.05)/diff(y_range)
    # # left <- -x_range[1]/diff(x_range)
    # # bottom <- (y_range[2]*1.3)/diff(y_range)
    # # right <- -x_range[2]/diff(x_range)
    # left <- 0.9
    # bottom <- 0
    # right <- 1.5
    # top <- 1
    # 
    # # Add inset image
    # img_element_plot <- element_plot +
    #   inset_element(
    #     p = img,
    #     left = left,
    #     bottom = bottom,
    #     right = right,
    #     top = top,
    #     align_to = "panel"
    #   )
    # 
    # print(img_element_plot)
    
    #print(img_element_plot)
    #print(element_plot)
    # Save individual plots if output_dir is provided
    if (!is.null(output_dir)) {
      ggsave(
        filename = file.path(output_dir, paste0(df_name, "_", "top", element,'.png')),
        plot = element_plot,
        width = 4, height = 10,
        dpi = 800
        
      )
    }
    
  }
}
plot_individual_clust <- function(df, df_name = NULL, output_dir = NULL){
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(jpeg)
    
    elements <- df %>%
      dplyr::select(-c(position..mm., Cluster)) %>%
      names()
    df_long <- make_long_clust(df)
    df_long$Cluster <- as.factor(df_long$Cluster)
    if (!is.null(output_dir)) {
      dir.create(output_dir, showWarnings = FALSE)  # Create output directory if not exists
    }
    
    for (element in elements) {
      element_data <- df_long %>% dplyr::filter(Variable == element)
      print(df_long%>%head())
      element_plot <- ggplot(df_long %>% dplyr::filter(Variable == element), aes(x = position..mm., y = Value, group = Variable, color = Cluster)) +
        geom_line(linewidth = 0.5, alpha = 0.9) +
        labs(
          title = paste("Line Plot for", element, "in", df_name),
          x = "Position",
          y = element
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
              axis.text.y = element_text(size = 6),
              axis.title = element_text(size = 8),
              strip.text = element_text(size = 8),
              legend.position = "none",
              aspect.ratio = 5) +
        coord_flip() +
        scale_x_reverse()
      path <- "/Users/maxshore/Documents/Unibe/MasterThesis/R/data/POS-22-20_highresruler_copy_min_crpd.jpg"
      
      # read the jpef file from device
      # img <- readJPEG(path, native = TRUE)
      
      # Extract axis limits from the built plot
      # plot_build <- ggplot_build(element_plot)
      # x_range <- plot_build$layout$panel_scales_x[[1]]$range$range
      # y_range <- plot_build$layout$panel_scales_y[[1]]$range$range
      # print(x_range)
      # print(paste(min(element_data$SampleID),max(element_data$SampleID)))
      # print(y_range)
      # print(paste(min(element_data$Value),max(element_data$Value)))
      # left <- (y_range[2]*1.05)/diff(y_range)
      # bottom <- -x_range[1]/diff(x_range)
      # right <- (y_range[2]*1.3)/diff(y_range)
      # top <- -x_range[2]/diff(x_range)
      # left <- 0.9
      # bottom <- 0.044
      # right <- 1.5
      # top <- 0.956
      # 
      # # Add inset image
      # img_element_plot <- element_plot +
      #   inset_element(
      #     p = img,
      #     left = left,
      #     bottom = bottom,
      #     right = right,
      #     top = top
      #   )
      
      # print(img_element_plot)
      
      #print(img_element_plot)
      #print(element_plot)
      # Save individual plots if output_dir is provided
      if (!is.null(output_dir)) {
        ggsave(
          filename = file.path(output_dir, paste0(df_name, "_", element, "clustered", ".png")),
          plot = element_plot,
          width = 4, height = 10,
          dpi = 800
          
        )
      }
      
    }
  }
save_svg_plot <- function(plot, filename, width = 6, height = 4) {
  svglite::svglite(filename, width = width, height = height)
  print(plot)
  dev.off()
}

plot_individual_ratio <- function(df, df_name = NULL, output_dir = NULL){
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  
  elements <- df %>%
    select(-c(position..mm.)) %>%
    names()
  df_long <- make_long(df)
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE)  # Create output directory if not exists
  }
  
  for (element in elements) {
    element_data <- df_long %>% filter(Variable == element)
    
    element_plot <- ggplot(element_data, aes(x = position..mm., y = Value)) +
      geom_line(color = "black", linewidth = 0.2, alpha = 0.9) +
      labs(
        title = paste("Line Plot for", element, "in", df_name),
        x = "Position",
        y = element
      ) +
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
      scale_x_reverse() 
      #scale_y_reverse()
    
    if (!is.null(output_dir)) {
      filename <- file.path(output_dir, paste0(df_name, "_", element, "_line.svg"))
      save_svg_plot(element_plot, filename, width = 4, height = 10)
    }
  }
}


# Function to create correlation map
create_corr_map_basic <- function(df) {
  if('position..mm.' %in% names(df)){
    df_numeric <- df %>% select(-position..mm.)}
  else{
    df_numeric <- df
  }
  corr_matrix <- cor(df_numeric, use = "pairwise.complete.obs")
  corrplot(corr_matrix, 
           method = "color", 
           type = "upper", 
           tl.col = "black", 
           tl.srt = 45,
           addCoef.col = "black",
           number.cex = 0.7,
           col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200))
}

create_corr_map <- function(df, title = NULL) {
  if (!requireNamespace("corrplot", quietly = TRUE)) stop("Please install 'corrplot'.")
  
  if ("position..mm." %in% names(df)) {
    df_numeric <- df %>% dplyr::select(-position..mm.)
  } else {
    df_numeric <- df
  }
  
  corr_matrix <- cor(df_numeric, use = "pairwise.complete.obs")
  
  # Define a diverging, color-blind-friendly palette (Blue-White-Orange)
  col_palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(200)
  
  # Set plotting parameters for high-quality output
  old_par <- par(no.readonly = TRUE)  # store current graphics settings
  on.exit(par(old_par))               # restore on function exit
  
  corrplot::corrplot(
    corr_matrix,
    method = "color",
    type = "upper",
    tl.col = "black",        # variable name color
    tl.srt = 45,             # rotate text labels
    tl.cex = 0.9,            # axis label size
    cl.cex = 0.8,            # color legend text size
    addCoef.col = "black",   # show correlation values
    number.cex = 0.7,        # size of numbers
    col = col_palette,
    mar = c(0, 0, 2, 0),     # top margin for optional title
    diag = FALSE             # hide diagonal (optional aesthetic)
  )
  
  if (!is.null(title)) {
    title(main = title, font.main = 2, cex.main = 1.2, line = 0.5)
  }
}




# Function to create boxplots
create_boxplots <- function(df) {
  df_long <- df %>%
    select(-SampleID) %>%
    pivot_longer(cols = everything(), 
                 names_to = "element", 
                 values_to = "value")
  ggplot(df_long, aes(x = element, y = value)) +
    geom_boxplot() +
    facet_wrap(~element, scales = "free") +
    theme_minimal()
}

# Function to compute ratios
compute_ratios <- function(df, el1, el2) {
  # Create a name for the new ratio column
  ratio_name <- paste(el1, el2, sep = "_")
  
  # Add the new column to the dataframe
  df <- df %>%
    mutate(!!ratio_name := log(.data[[el1]] / .data[[el2]]))
 
  return(df)
}


# Function to perform PCA
perform_pca <- function(df, df_name = NULL, plot_scores = FALSE, plot_loadings = FALSE, print_results = FALSE, output_dir = NULL) {
  library(ggfortify)
  library(ggplot2)
  
  # if (!is.null(output_dir)) {
  #   dir.create(output_dir, showWarnings = FALSE)  # Create output directory if not exists
  # }
  # Ensure the input dataframe is numeric
  df_numeric <- df %>% select_if(is.numeric)
  df_numeric <- df_numeric %>% select(-c(position..mm.))
  
  # Scale the data
  df_scaled <- scale(df_numeric)
  
  # Perform PCA
  pca_result.no.sign <- prcomp(df_scaled, center = TRUE, scale. = TRUE)
  pca_result <- fix_pca_signs(pca_result.no.sign)
  
  # Plot PCA scores (biplot) if plot_scores = TRUE
  if (plot_scores) {
    scores_plot <- autoplot(
      pca_result,
      data = df,
      colour = 'position..mm.'
    ) +
      labs(
        title = paste("PCA Scores Plot for", df_name),
        x = "Principal Component 1",
        y = "Principal Component 2"
      ) +
      guides(colour = guide_colorbar(reverse = T))+
      theme_light()
    print(scores_plot)
    
  
    if (!is.null(output_dir)) {
      ggsave(
        filename = file.path(output_dir, paste0(df_name, "_scores",".png")),
        plot = scores_plot,
        width = 10, height = 8
      )
    }
  }
  # Plot PCA loadings if plot_loadings = TRUE
  if (plot_loadings) {
    loadings_plot <- autoplot(
      pca_result,
      loadings = TRUE,
      loadings.label = TRUE,
      loadings.label.size = 3
    ) +
      labs(
        title = paste("PCA Loadings Plot for", df_name),
        x = "Principal Component 1",
        y = "Principal Component 2"
      ) +
      theme_light()
    
    
    print(loadings_plot)
  
    if (!is.null(output_dir)) {
      ggsave(
        filename = file.path(output_dir, paste0(df_name, "_loadings",".png")),
        plot = loadings_plot,
        width = 10, height = 8
      )
    }
  }
  # Print PCA results with dataframe name
  if (print_results) {
    message <- if (!is.null(df_name)) {
      paste("PCA results for dataframe:", df_name)
    } else {
      "PCA results:"
    }
    cat(message, "\n")
    print(summary(pca_result))
    cat("\n \n")
  }
  
  # Return the PCA result for further use
  
  return(pca_result)
}

perform_pca2 <- function(df, df_name = NULL, plot_scores = FALSE, plot_loadings = FALSE, print_results = FALSE, output_dir = NULL) {
  library(ggfortify)
  library(ggplot2)
  
  # Ensure the input dataframe is numeric
  df_numeric <- df %>% select_if(is.numeric) %>% select(-c(position..mm.))
  
  
  # Perform PCA
  pca_result.no.sign <- prcomp(df_numeric, center = TRUE, scale. = TRUE)
  pca_result <- fix_pca_signs(pca_result.no.sign)
  pca_result$position..mm. <- df$position..mm.
  
  # Plot PCA scores separately, colored by position
  if (plot_scores) {
    scores_plot <- ggplot(data.frame(pca_result$x), aes(x = PC1, y = PC2, color = df$position..mm.)) +
      geom_point(size = 3) +
      labs(
        title = paste("PCA Scores Plot for", df_name),
        x = "Principal Component 1",
        y = "Principal Component 2",
        color = "Position (mm)"
      ) +
      scale_color_viridis_c() + # Uses a perceptually uniform color scale
      guides(colour = guide_colorbar(reverse = TRUE)) +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA))
    
    print(scores_plot)
    
    if (!is.null(output_dir)) {
      ggsave(
        filename = file.path(output_dir, paste0(df_name, "_scores.png")),
        plot = scores_plot,
        width = 10, height = 8
      )
    }
  }
  
  # Plot PCA loadings separately
  if (plot_loadings) {
    loadings_data <- data.frame(pca_result$rotation[, 1:2], Variable = rownames(pca_result$rotation))
    
    loadings_plot <- ggplot(loadings_data, aes(x = 0, y = 0, xend = PC1, yend = PC2, label = Variable)) +
      geom_segment(arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "black") +
      geom_text(aes(x = PC1, y = PC2), size = 4, hjust = 1.2, vjust = 1.2) +
      labs(
        title = paste("PCA Loadings Plot for", df_name),
        x = "Principal Component 1",
        y = "Principal Component 2"
      ) +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA))
    
    print(loadings_plot)
    
    if (!is.null(output_dir)) {
      ggsave(
        filename = file.path(output_dir, paste0(df_name, "_loadings.png")),
        plot = loadings_plot,
        width = 10, height = 8
      )
    }
  }
  
  # Print PCA results
  if (print_results) {
    message <- if (!is.null(df_name)) {
      paste("PCA results for dataframe:", df_name)
    } else {
      "PCA results:"
    }
    cat(message, "\n")
    print(summary(pca_result))
    cat("\n \n")
  }

  return(pca_result)
}

fix_pca_signs <- function(pca_result) {
  if (pca_result$rotation[1, 1] < 0) {
    pca_result$x[, 1] <- -pca_result$x[, 1]
    pca_result$rotation[, 1] <- -pca_result$rotation[, 1]
  }
  if (pca_result$rotation[1, 2] < 0) {
    pca_result$x[, 2] <- -pca_result$x[, 2]
    pca_result$rotation[, 2] <- -pca_result$rotation[, 2]
  }
  if (pca_result$rotation[1, 3] < 0) {
    pca_result$x[, 3] <- -pca_result$x[, 3]
    pca_result$rotation[, 3] <- -pca_result$rotation[, 3]
  }
  return(pca_result)
}



pca_downcore_plot <- function(pc_df, df_name, pc_num, output_dir = NULL) {
  library(ggplot2)
  library(dplyr)
  library(svglite)
  
  # Create a clean data frame from prcomp output
  scores_df <- as.data.frame(pc_df$x)
  
  if (!"position..mm." %in% names(pc_df)) {
    stop("Missing 'position..mm.' column in pc_df")
  }

  scores_df$position <- pc_df$position..mm.

  
  # Ensure PC exists
  if (!pc_num %in% names(scores_df)) {
    stop(paste("Principal component", pc_num, "not found."))
  }
  
  # Plot
  p <- ggplot(scores_df, aes(x = position, y = .data[[pc_num]])) +
    geom_line(linewidth = 0.1) +
    geom_ribbon(aes(ymin = 0, ymax = ifelse(.data[[pc_num]] > 0, .data[[pc_num]], 0)),
                fill = "blue", alpha = 0.5) +
    geom_ribbon(aes(ymin = ifelse(.data[[pc_num]] < 0, .data[[pc_num]], 0), ymax = 0),
                fill = "red", alpha = 0.5) +
    labs(title = paste(pc_num, " Downcore Variation"),
         x = "Position (mm)", y = pc_num) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.position = "none",
          aspect.ratio = 5,
          plot.background = element_rect(fill = 'white'),
          ) +
    coord_flip() +
    scale_x_reverse()
  print(p)
  # Save using ggsave and svglite device (stable workaround)
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    filename <- file.path(output_dir, paste0(df_name, "_", pc_num, "_DownCore.pdf"))
    
    ggsave(
      filename = filename,
      plot = p,
      device = pdf,
      width = 4,
      height = 10,
      dpi = 300
    )
  }
  
}



# Function to perform clustering 

perform_clustering <- function(df, df_name, k = 3, method = "kmeans") {
  library(factoextra)
  library(cluster)
  library(NbClust)
  # Ensure numeric columns are used for clustering
  
  df_numeric <- df %>% select_if(is.numeric) %>% select(-c('position..mm.')) %>% scale()
  
  
  # Perform clustering based on the chosen method
  if (method == "kmeans") {
    clustering <- kmeans(df_numeric, centers = k)
  } else if (method == "hclust") {
    distance_matrix <- dist(df_numeric)
    hierarchical <- hclust(distance_matrix)
    clustering <- cutree(hierarchical, k = k)
  } else {
    stop("Unsupported clustering method. Use 'kmeans' or 'hclust'.")
  }
  
  # Add cluster labels to the dataframe
  df$Cluster <- if (method == "kmeans") clustering$cluster else clustering
  
  # Print a summary of the clustering
  cat("\nClustering results for dataframe:", df_name, "\n")
  print(table(df$Cluster))

  # Visualize clusters and optimal cluster number
  p <- fviz_cluster(clustering, data = df_numeric, geom = "point",
               stand = FALSE, frame.type = "norm")
  print(p)
  k.max <- 15
  wss <- sapply(1:k.max,
                function(k){kmeans(df_numeric, k, nstart=10 )$tot.withinss})
  plot(1:k.max, wss,
       type="b", pch = 19, frame = FALSE,
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
  fviz_nbclust(df_numeric, kmeans, method = "silhouette")

  # Return the dataframe with cluster labels
  return(df)
  
}

visualize_clusters <- function(df, df_name, output_dir = NULL, method = "PCA") {
  library(ggplot2)
  library(ggfortify)
  # Ensure numeric columns are used
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE)  # Create output directory if not exists
  }
  
  df_num <- df %>% select_if(is.numeric) %>% select(-c('position..mm.', 'Cluster'))
  df_numeric <- scale(df_num)
  
  # Perform dimensionality reduction (PCA)
  if (method == "PCA") {
    pca_result <- prcomp(df_numeric, center = TRUE, scale. = F)
    reduced_data <- as.data.frame(pca_result$x[, 1:2])
    colnames(reduced_data) <- c("Dim1", "Dim2")
  } else {
    stop("Currently, only PCA is supported for dimensionality reduction.")
  }
  
  # Add cluster labels to the reduced data
  reduced_data$Cluster <- factor(df$Cluster)  # Assuming 'Cluster' column exists
  
  # Get the number of clusters for labeling
  n_clust <- reduced_data$Cluster %>% unique() %>% as.numeric() %>% max()
  
  # Plotting the scores and loadings
  scores_plot <- ggplot(reduced_data, aes(x = Dim1, y = Dim2, color = Cluster)) +
    geom_point(size = 2, alpha = 0.8, na.rm = T) +
    labs(
      title = paste("PCA Cluster Visualization for", df_name),
      x = "Dimension 1",
      y = "Dimension 2",
      color = "Cluster"
    ) +
    theme_bw()
  
  # Loadings (arrows) for PCA
  loadings_data <- data.frame(pca_result$rotation[, 1:2], Variable = rownames(pca_result$rotation))
  
  # Scores plot with loadings overlaid
  scores_with_loadings_plot <- scores_plot + 
    geom_segment(data = loadings_data, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "blue") +
    geom_text(data = loadings_data, aes(x = PC1, y = PC2, label = Variable), size = 4, hjust = 1.2, vjust = 1.2, color = "blue") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
  
  # Save the plot with scores and loadings to the output directory
  plot_file <- file.path(output_dir, paste0("Cluster_PCA_Scores_and_Loadings_", n_clust,'clust_', df_name, ".png"))
  ggsave(plot_file, plot = scores_with_loadings_plot, width = 8, height = 6, dpi = 300)
  
  # Return the plot (optional)
  return(scores_with_loadings_plot)
}

perform_periodicity_analysis <- function(df, element, df_name = NULL, output_dir = NULL) {
  
  library(ggplot2)
  library(dplyr)
  library(signal)  # For detrending
  library(pracma)  # For Lomb-Scargle periodogram

  # Ensure the element exists
  if (!(element %in% colnames(df))) {
    stop(paste("Element", element, "not found in dataframe."))
  }
  
  # Extract signal and positions
  signal <- df[[element]]
  position <- df$position..mm.  # Assuming position column exists
  
  # Remove NA values
  valid_indices <- complete.cases(position, signal)
  position <- position[valid_indices]
  signal <- signal[valid_indices]
  
  # Ensure uniform spacing (if needed)
  if (any(diff(position) != mean(diff(position)))) {
    cat("Warning: Uneven spacing detected, interpolating data...\n")
    interp_func <- approxfun(position, signal, method = "linear")
    new_position <- seq(min(position), max(position), length.out = length(position))
    signal <- interp_func(new_position)
    position <- new_position
  }
  
  # Remove linear trend (helps eliminate artificial zero-frequency dominance)
  signal <- detrend(signal, tt = "linear")  # Removes trend without distorting periodicity
  
  # Compute sampling rate (spatial frequency in mm)
  dx <- mean(diff(position))
  sampling_rate <- 1 / dx
  
  # Apply FFT
  n <- length(signal)
  fft_result <- fft(signal)
  
  # Compute frequencies (only keep the positive half, ensuring correct length)
  freq <- seq(0, sampling_rate / 2, length.out = floor(n / 2) + 1)
  power_spectrum <- Mod(fft_result[1:length(freq)])^2  # Match lengths correctly
  
  # Find dominant frequency (ignore zero frequency)
  dominant_index <- which.max(power_spectrum[-1]) + 1  # Avoid index 1 (zero freq)
  dominant_freq <- freq[dominant_index]
  
  # Convert dominant frequency to period (avoid Inf)
  dominant_period <- ifelse(dominant_freq == 0, NA, 1 / dominant_freq)
  
  # Create power spectrum plot
  spectrum_df <- data.frame(Frequency = freq, Power = power_spectrum)  # Ensured equal lengths
  spectrum_plot <- ggplot(spectrum_df, aes(x = Frequency, y = Power)) +
    geom_line(color = "blue") +
    geom_vline(xintercept = dominant_freq, linetype = "dashed", color = "red") +
    labs(
      title = paste("Periodicity Analysis for", element, "in", df_name),
      subtitle = paste("Dominant Period:", ifelse(is.na(dominant_period), "N/A", round(dominant_period, 2))),
      x = "Frequency (1/mm)",
      y = "Power"
    ) +
    theme_minimal()
  print(spectrum_plot)
  # Save plot if needed
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE)
    plot_filename <- file.path(output_dir, paste0("Periodicity_", element, "_", df_name, ".png"))
    ggsave(plot_filename, plot = spectrum_plot, width = 8, height = 6, dpi = 300)
    cat("Periodicity analysis plot saved at:", plot_filename, "\n")
  }
  
  # Print dominant period
  cat("Dominant period for", element, ":", ifelse(is.na(dominant_period), "N/A", round(dominant_period, 2)), "mm\n")
  
  return(spectrum_plot)
}

plot_individual_2 <- function(df, df_name = NULL, output_dir = NULL){
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(jpeg)
  
  elements <- df %>%
    select(-c(position..mm.)) %>%
    names()
  df_long <- make_long(df)
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE)  # Create output directory if not exists
  }
  
  for (element in elements) {
    element_data <- df_long %>% dplyr::filter(Variable == element)
    print(df_long%>%head())
    element_plot <- ggplot(df_long %>% dplyr::filter(Variable == element), aes(x = position..mm., y = Value)) +
      geom_point(color = "white", size=0.00001, shape='.') +
      labs(
        title = paste("Line Plot for", element, "in", df_name),
        x = "Position",
        y = element
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6),
            axis.title = element_text(size = 8),
            strip.text = element_text(size = 8),
            legend.position = "none",
            aspect.ratio = 5) +
      coord_flip() +
      scale_x_reverse()
    path <- "/Users/maxshore/Documents/Unibe/MasterThesis/R/data/POS-22-20_highresruler_copy_min_crpd copy.jpg"
    
    # read the jpef file from device
    img <- readJPEG(path, native = TRUE)
    
    # Extract axis limits from the built plot
    # plot_build <- ggplot_build(element_plot)
    # x_range <- plot_build$layout$panel_scales_x[[1]]$range$range
    # y_range <- plot_build$layout$panel_scales_y[[1]]$range$range
    # print(x_range)
    # print(paste(min(element_data$SampleID),max(element_data$SampleID)))
    # print(y_range)
    # print(paste(min(element_data$Value),max(element_data$Value)))
    # left <- (y_range[2]*1.05)/diff(y_range)
    # bottom <- -x_range[1]/diff(x_range)
    # right <- (y_range[2]*1.3)/diff(y_range)
    # top <- -x_range[2]/diff(x_range)
    left <- 0.3
    bottom <- 0.044
    right <- 1.3
    top <- 0.956
    
    # Add inset image
    img_element_plot <- element_plot +
      inset_element(
        p = img,
        left = left,
        bottom = bottom,
        right = right,
        top = top,
        on_top = F
      )
    
    # print(img_element_plot)
    
    #print(img_element_plot)
    #print(element_plot)
    # Save individual plots if output_dir is provided
    if (!is.null(output_dir)) {
      ggsave(
        filename = file.path(output_dir, paste0(df_name, "_", element, ".png")),
        plot = img_element_plot,
        width = 4, height = 10,
        dpi = 800
        
      )
    }
    
  }
}

save_svg_plot <- function(plot, filename, width = 6, height = 4) {
  svglite::svglite(filename, width = width, height = height)
  print(plot)
  dev.off()
}

synchronize_to_reference <- function(ref_df, target_df,
                                     ref_depth_col = "position..mm.",
                                     target_depth_col = "position..mm.",
                                     target_vars = NULL,
                                     method = "linear") {
  # Ensure depth columns are numeric
  ref_depths <- ref_df[[ref_depth_col]]
  target_depths <- target_df[[target_depth_col]]
  
  # Choose which variables to interpolate
  if (is.null(target_vars)) {
    target_vars <- setdiff(names(target_df), target_depth_col)
  }
  
  # Initialize output data frame with reference positions
  interp_df <- data.frame(ref_df[[ref_depth_col]])
  colnames(interp_df) <- ref_depth_col
  
  # Interpolate each variable
  for (var in target_vars) {
    if (!is.numeric(target_df[[var]])) next  # skip non-numeric vars
    
    interp_vals <- approx(
      x = target_depths,
      y = target_df[[var]],
      xout = ref_depths,
      method = method,
      rule = 2  # allows extrapolation with flat values at ends
    )$y
    
    interp_df[[var]] <- interp_vals
  }
  
  return(interp_df)
}
combined_df_creator <- function(sel = NULL){
  
  if(is.null(sel)){
    sel <- c('position..mm.','Fe', 'Mn', 'S', 'Ti', 'Al', 'Si', 'K', 'Rb', 'Zr', 'Ca', 'Sr', 'P', 'Ba')
  }
  
  #### XRF 0.5mm dataset ####
  df.raw <- load_raw_data()
  df.xrf <- clean_df(df.raw, sel = sel)
  df.xrf$raw.200 <- NULL
  df.xrf <- as_tibble(df.xrf$raw.500) %>%
    select_top_varves()
  
  df.cs <- closed_sum(df.xrf) 
  df.clr <- clr_transform(df.xrf)
  
  
  #### HSI dataset ####
  path_HSI <- '/Users/maxshore/Documents/Unibe/MasterThesis/POS-22-20_DATA/POS_22_20_HSI/RABD673_Chla/POS22-20_230802-145515_refl_sub_RABD673_spl_645_674.csv'
  
  raw.HSI <- read.csv(file = '/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/POS_22_20_HSI/HSI_calibrated.xlsx') %>%
    tibble()
  
  hsi <- raw.HSI %>% select(c(Core.Depth..mm., RABD673, Moving.Average, TChl.ug.g))
  hsi <- rename(hsi, position..mm.= Core.Depth..mm.) %>%
    select_top_varves()
  
  #### CNS dataset ####
  library(readxl)
  
  path_CNS <- '/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/POS-22-20_CNS.xlsx'
  
  raw.cns <- read_excel(path_CNS) %>%
    tibble()
  # convert depth from cm to mm
  raw.cns$Depth <- raw.cns$Depth*10
  cns <- rename(raw.cns, position..mm. = Depth) %>%
    select_top_varves()
  cns <- cns %>%
    rename_with(~ gsub('[%/\" ]', '', .x), .cols = matches('TN|TC|TS'))
  
  #### grayscale dataset on 20 bandwidth ####
  path_gs <- '/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/generated/grayscale/gray_profile_bw20.csv'
  
  gs <- read.csv(path_gs)
  
  
  #### synchronize ###############################################################
  
  # Separate facies from grayscale
  gs_fac <- gs %>% select(position..mm., facies_class)
  gs_gray <- gs %>% select(position..mm., mean_gray)
  
  gray_sync <- synchronize_to_reference(df.xrf, gs_gray,
                                        ref_depth_col = "position..mm.",
                                        target_depth_col = "position..mm.",
                                        target_vars = c("mean_gray"))
  
  library('fuzzyjoin')
  facies_sync_all <- fuzzyjoin::difference_left_join(
    df.xrf, gs_fac, 
    by = "position..mm.",
    max_dist = 0.25,
    distance_col = "dist"
  )
  
  facies_sync <- facies_sync_all %>%
    group_by(position..mm..x) %>%
    slice_min(order_by = dist, n = 1) %>%
    ungroup() %>%
    select(position..mm. = position..mm..x, facies_class)
  
  
  hsi_sync <- synchronize_to_reference(df.xrf, hsi,
                                       target_vars = c("RABD673", 'Moving.Average'))
  
  cns_sync <- synchronize_to_reference(df.xrf, cns,
                                       target_vars = c("TNWt", "TCWt", "TSWt"))
  
  
  combined_df <- df.xrf %>%
    left_join(gray_sync, by = "position..mm.") %>%
    left_join(hsi_sync, by = "position..mm.") %>%
    left_join(cns_sync, by = "position..mm.") %>%
    left_join(df.clr, by = 'position..mm.') %>% #.y
    left_join(df.cs, by = 'position..mm.') %>% # normal no suffix
    left_join(facies_sync, by = "position..mm.")
  
  #### compute ratios ########################################################
  # (make sure to perform on the correct xrf normalisation (cs))
  
  combined_df <- compute_ratios(combined_df, 'Mn', 'Fe')
  combined_df <- compute_ratios(combined_df, 'S', 'Ti')
  combined_df <- compute_ratios(combined_df, 'Fe', 'Ti')
  combined_df <- compute_ratios(combined_df, 'Si', 'Ti')
  combined_df <- compute_ratios(combined_df, 'Ba', 'Ti')
  combined_df <- compute_ratios(combined_df, 'Ca', 'Ti')
  combined_df <- compute_ratios(combined_df, 'Ti', 'Al')
  combined_df <- compute_ratios(combined_df, 'Si', 'Al')
  combined_df <- compute_ratios(combined_df, 'Zr', 'Rb')
  combined_df <- compute_ratios(combined_df, 'TCWt', 'TNWt')
  
  # save combined df
  
  # write.csv(combined_df, file = 'data/generated/combined/combined_500.csv')
  return(combined_df)
}

create_combined_df <- function() {
  combined <- combined_df_creator()
  no_xrf_combined <- combined %>% select(any_of(vars_not_xrf))
  return(no_xrf_combined)
}

load_and_clean_xrf <- function(raw_df, elements) {
  df <- clean_df(raw_df, sel = elements)
  df$raw.200 <- NULL
  df <- tibble(df$raw.500) %>% select_top_varves()
  df <- closed_sum(df)
  return(df)
}

add_ratios <- function(df, ratio_list) {
  for (pair in ratio_list) {
    df <- compute_ratios(df, pair[1], pair[2])
  }
  return(df)
}

run_pca_selection_iteration <- function(xrf_vars_kept, ratio_list, other_vars, no_xrf_combined, iteration_name) {
  # Get all raw elements needed for computing ratios
  raw_elements_for_ratios <- get_unique_elements_from_ratios(ratio_list)
  
  # Ensure we load all variables required for computing ratios
  xrf_elements <- unique(c(xrf_vars_kept, raw_elements_for_ratios, "position..mm."))
  
  df.raw <- load_raw_data()
  df.xrf <- load_and_clean_xrf(df.raw, elements = xrf_elements)
  
  combined <- df.xrf %>%
    left_join(no_xrf_combined, by = "position..mm.") %>%
    add_ratios(ratio_list)
  
  selected_vars <- c(xrf_vars_kept, other_vars)
  
  combined %>%
    select(any_of(selected_vars)) %>%
    perform_pca2(
      plot_loadings = TRUE,
      df_name = paste0('iter_', iteration_name, '_n_', length(selected_vars)),
      output_dir = 'plots/PCA_variable_selection'
    )
}

get_unique_elements_from_ratios <- function(ratios) {
  unique(unlist(ratios))
}