# PCA element selector script

# ==========================================================
# Automated PCA over multiple element selections, resolutions and core sections
# ==========================================================

library(ggplot2)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(tibble)
library(ggrepel)
library(svglite)
library(scales)
library(readr)

source("xrf_functions.R")
elem_sets <- list(
  
  # ------------------------------------------------------------------------
  # 1. Minimal & robust sets (low noise, ideal for structure)
  # ------------------------------------------------------------------------
  mineral_core = c("Fe", "Ti", "Al", "Si", "Ca", "K", "Zr", "position..mm."),
  quartz_clay = c("Si", "Al", "K", "Ti", "Zr", "position..mm."),
  detrital_basic = c("Ti", "Al", "K", "Zr", "Ca", "position..mm."),
  
  # ------------------------------------------------------------------------
  # 2. Lithogenic + carbonate interaction sets
  # ------------------------------------------------------------------------
  litho_carb = c("Fe", "Ti", "Al", "Si", "Ca", "K", "position..mm."),
  broadened_litho = c("Fe", "Ti", "Al", "Si", "Ca", "K", "Rb", "Ba", "position..mm."),
  
  # ------------------------------------------------------------------------
  # 3. Redox-related selections (Fe–S–Ba axis)
  # ------------------------------------------------------------------------
  redox_basic = c("Fe", "S", "Ba", "Ti", "Al", "position..mm."),
  redox_extended = c("Fe", "S", "Ba", "Ti", "Al", "Si", "K", "Zr", "position..mm."),
  
  # ------------------------------------------------------------------------
  # 4. Productivity / carbonate / detrital mixes
  # ------------------------------------------------------------------------
  prod_carb_mix = c("Ba", "S", "Ca", "Fe", "Ti", "Al", "position..mm."),
  prod_detrital = c("Ba", "S", "Ca", "Ti", "Al", "K", "Zr", "position..mm."),
  weathering_signal = c("K", "Rb", "Zr", "Ti", "Al", "position..mm."),
  
  # ------------------------------------------------------------------------
  # 5. Full stable set (all relevant, excluding noisy Mn/P/Sr)
  # ------------------------------------------------------------------------
  full_clean = c("Fe","S","Ba","Ti","Al","Si","K","Ca","Rb","Zr","position..mm."),
  
  # ------------------------------------------------------------------------
  # 6. Compare impact of trace metals (Ni, Zn)
  # ------------------------------------------------------------------------
  trace_metals_light = c("Fe","Ti","Al","Si","K","Zr","Ni","Zn","position..mm."),
  trace_metals_extended = c("Fe","S","Ba","Ti","Al","Si","K","Ca","Zr","Ni","Zn","position..mm.")
)

# ----------------------------------------------------------
# 2. Resolution & core section parameters
# ----------------------------------------------------------

resolutions <- c(0.2, 0.5)
sections <- c("AC_", "top_", "bottom_")

# ----------------------------------------------------------
# 3. Output base folder
# ----------------------------------------------------------
output_dir <- "/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/plots/PCA_selector/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ----------------------------------------------------------
# 4. MAIN LOOP
# ----------------------------------------------------------

for (sel_name in names(elem_sets)) {
  
  sel_vec <- elem_sets[[sel_name]]
  
  for (res in resolutions) {
    
    res_micro <- res * 1000
    
    # make resolution subfolder
    res_folder <- file.path(output_dir, paste0(res_micro, "um"))
    dir.create(res_folder, showWarnings = FALSE)
    
    for (section in sections) {
      
      # ---------------------------
      # Label for text
      # ---------------------------
      l_text <- if (section == "AC_") {
        "all core"
      } else if (section == "top_") {
        "top 35cm"
      } else {
        "bottom 35-123.1cm"
      }
      
      # ---------------------------
      # Prepare file naming
      # ---------------------------
      section_folder <- file.path(res_folder, section)
      dir.create(section_folder, showWarnings = FALSE)
      
      selection_size <- length(sel_vec[sel_vec != "position..mm."])
      
      plot_file <- file.path(
        section_folder,
        paste0("PCA_", sel_name, "_", selection_size, "elem.svg")
      )
      
      message("----------")
      message("Running ", sel_name, 
              " | ", res_micro, " µm | ", section, 
              " | ", selection_size, " elements")
      message("Saving → ", plot_file)
      
      # -------------------------------------------------
      # 5. Load data
      # -------------------------------------------------
      path_combined_folder <- "/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/generated/combined/"
      combined_file <- paste0(path_combined_folder, 
                              "combined_", section, res_micro, "_inc_coh.csv")
      
      combined <- read.csv(combined_file)
      
      df.raw <- load_raw_data()
      df.xrf <- clean_df(df.raw, sel = sel_vec)
      
      # Select resolution field
      if (res_micro == 500) {
        df.xrf$raw.200 <- NULL
        df.xrf <- tibble(df.xrf$raw.500)
      } else {
        df.xrf$raw.500 <- NULL
        df.xrf <- tibble(df.xrf$raw.200)
      }
      
      # Select subsets
      if(section == "top_"){
        df.xrf <- df.xrf %>% select_top_varves()
      }
      if(section == "bottom_"){
        df.xrf <- df.xrf %>% select_bottom()
      }
      
      # -------------------------------------------------
      # 6. Transformations
      # -------------------------------------------------
      df.clr <- clr_transform(df.xrf)
      df.pca <- df.clr %>% select(-position..mm.)
      
      pca_result <- prcomp(df.pca, scale. = TRUE, center = TRUE)
      pca_result <- fix_pca_signs(pca_result)
      
      # Scores & loadings
      scores <- as_tibble(pca_result$x)
      loadings <- as_tibble(pca_result$rotation, rownames = "variable")
      
      loadings_scaled <- loadings %>%
        mutate(
          PC1 = PC1 * pca_result$sdev[1],
          PC2 = PC2 * pca_result$sdev[2],
          PC3 = PC3 * pca_result$sdev[3]
        )
      
      scores <- scores %>%
        mutate(depth = combined$position..mm.) %>%
        mutate(facies = combined$facies_class)
      
      var_expl <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
      
      # -------------------------------------------------
      # 7. Biplot
      # -------------------------------------------------
      p_plot <- ggplot() +
        geom_point(data = scores, aes(x = PC1, y = PC2), shape = 16, size = 0.5) +
        geom_segment(data = loadings_scaled,
                     aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
                     arrow = arrow(length = unit(0.2, "cm")),
                     color = "red", linewidth = 0.6) +
        geom_text_repel(data = loadings_scaled,
                        aes(x = PC1 * 3.2, y = PC2 * 3.2, label = variable),
                        size = 5, color = "red") +
        geom_hline(yintercept = 0, linetype = "dotted") +
        geom_vline(xintercept = 0, linetype = "dotted") +
        xlab(paste0("PC1 (", var_expl[1], "%)")) +
        ylab(paste0("PC2 (", var_expl[2], "%)")) +
        coord_equal() +
        theme_minimal(base_size = 14) +
        ggtitle(paste0(
          "PCA Biplot – ", res_micro, "µm – ", l_text,
          " – ", selection_size, " elements"
        ))
      
      # -------------------------------------------------
      # 8. Save
      # -------------------------------------------------
      svg_file <- file.path(
        section_folder,
        paste0("PCA_", sel_name, "_", selection_size, "elem.svg")
      )
      
      png_file <- file.path(
        section_folder,
        paste0("PCA_", sel_name, "_", selection_size, "elem.png")
      )
      
      # --- Save SVG ---
      svglite(svg_file, width = 8, height = 7)
      print(p_plot)
      dev.off()
      
      # --- Save PNG ---
      png(filename = png_file, width = 8, height = 7, units = "in", res = 300)
      print(p_plot)
      dev.off()
      
    } # sections
  } # resolutions
} # selections
