# autocorrelation analysis

source('xrf_functions.R')

library(ggplot2)
library(RColorBrewer)
library(dplyr)

combined <- read.csv('/Users/maxshore/Documents/Unibe/MasterThesis/masterthesis/R/data/generated/combined/combined_500_inc_coh.csv')

library(dplyr)
library(tidyr)
# Optional: install.packages(c("zoo"))
library(zoo)

# --- Params to tune ---
upper_q <- 0.90
lower_q <- 0.10
win_L   <- 20   # half-window length in samples (use what makes sense for your depth/time resolution)
min_sep <- 5    # minimum separation (in samples) between detected events to avoid duplicates

# --- Thresholds ---
upper_decile <- quantile(combined$mean_gray, upper_q, na.rm = TRUE)
lower_decile <- quantile(combined$mean_gray, lower_q, na.rm = TRUE)

gs <- combined$mean_gray

# --- Helper to find local extrema (peaks above 90th, troughs below 10th) ---
is_local_max <- function(x) {
  # internal points: strictly greater than neighbors
  nx <- length(x)
  res <- rep(FALSE, nx)
  idx <- 2:(nx-1)
  res[idx] <- (x[idx] > x[idx-1]) & (x[idx] > x[idx+1])
  res
}
is_local_min <- function(x) {
  nx <- length(x)
  res <- rep(FALSE, nx)
  idx <- 2:(nx-1)
  res[idx] <- (x[idx] < x[idx-1]) & (x[idx] < x[idx+1])
  res
}

peaks_idx  <- which(is_local_max(gs) & gs >= upper_decile)
troughs_idx <- which(is_local_min(gs) & gs <= lower_decile)
events_idx <- sort(c(peaks_idx, troughs_idx))

# --- Enforce minimum separation between events ---
if (length(events_idx) > 1) {
  keep <- logical(length(events_idx))
  last_kept <- -Inf
  for (i in seq_along(events_idx)) {
    if (events_idx[i] - last_kept >= min_sep) {
      keep[i] <- TRUE
      last_kept <- events_idx[i]
    }
  }
  events_idx <- events_idx[keep]
}

# --- Build windows around each event (drop those that would exceed series bounds) ---
n <- nrow(combined)
valid <- events_idx[(events_idx - win_L >= 1) & (events_idx + win_L <= n)]
if (length(valid) == 0) {
  stop("No valid event windows: adjust win_L/min_sep or thresholds.")
}

# --- Function: lag per event for one variable ---
lag_from_window_ccf <- function(x, y, lag.max = NULL) {
  # x and y are numeric vectors of equal length (one window)
  # optional: standardize to focus on shape
  xz <- scale(x, center = TRUE, scale = TRUE)[,1]
  yz <- scale(y, center = TRUE, scale = TRUE)[,1]
  res <- ccf(xz, yz, plot = FALSE, lag.max = lag.max)
  lag <- res$lag[which.max(abs(res$acf))]
  val <- max(abs(res$acf))
  c(lag = lag, max_abs_acf = val)
}

# --- Variables to analyze ---
vars_to_test <- c("S_Ti", "S", "Ca", "Ti", "Fe")

# --- Collect event-wise lags and also compute full-series lag for comparison ---
results_list <- lapply(vars_to_test, function(vn) {
  y <- combined[[vn]]
  # full-series CCF (for reference)
  ccf_full <- ccf(gs, y, plot = FALSE)
  lag_full <- ccf_full$lag[which.max(abs(ccf_full$acf))]
  max_full <- max(abs(ccf_full$acf))
  
  # per-event lags
  ev_lags <- lapply(valid, function(i0) {
    idx <- (i0 - win_L):(i0 + win_L)
    out <- lag_from_window_ccf(gs[idx], y[idx])
    tibble(event_center = i0,
           lag_event = as.integer(out["lag"]),
           max_abs_acf_event = as.numeric(out["max_abs_acf"]))
  }) %>% bind_rows()
  
  summary <- ev_lags %>%
    summarise(
      n_events = n(),
      median_lag_event = median(lag_event, na.rm = TRUE),
      iqr_lag_event_low = quantile(lag_event, 0.25, na.rm = TRUE),
      iqr_lag_event_high = quantile(lag_event, 0.75, na.rm = TRUE),
      mean_abs_acf_event = mean(max_abs_acf_event, na.rm = TRUE)
    ) %>%
    mutate(
      variable = vn,
      lag_full = lag_full,
      max_abs_acf_full = max_full
    ) %>%
    relocate(variable)
  
  list(summary = summary, events = ev_lags)
})

# --- Summary table across variables ---
lag_summary <- bind_rows(lapply(results_list, `[[`, "summary"))
lag_summary





#### second version for CCF ####

library(dplyr)
library(tidyr)

# --- Params ---
upper_q <- 0.90
lower_q <- 0.10
win_L   <- 20   # half-window length in samples
min_sep <- 5    # minimum separation between detected events

# --- Thresholds ---
upper_decile <- quantile(combined$mean_gray, upper_q, na.rm = TRUE)
lower_decile <- quantile(combined$mean_gray, lower_q, na.rm = TRUE)
gs <- combined$mean_gray

# --- Local max/min helpers ---
is_local_max <- function(x) {
  nx <- length(x); res <- rep(FALSE, nx)
  idx <- 2:(nx-1); res[idx] <- (x[idx] > x[idx-1]) & (x[idx] > x[idx+1])
  res
}
is_local_min <- function(x) {
  nx <- length(x); res <- rep(FALSE, nx)
  idx <- 2:(nx-1); res[idx] <- (x[idx] < x[idx-1]) & (x[idx] < x[idx+1])
  res
}

# --- Detect events ---
peaks_idx   <- which(is_local_max(gs) & gs >= upper_decile)
troughs_idx <- which(is_local_min(gs) & gs <= lower_decile)
events_idx  <- sort(c(peaks_idx, troughs_idx))

# --- Enforce minimum separation ---
if (length(events_idx) > 1) {
  keep <- logical(length(events_idx)); last_kept <- -Inf
  for (i in seq_along(events_idx)) {
    if (events_idx[i] - last_kept >= min_sep) {
      keep[i] <- TRUE; last_kept <- events_idx[i]
    }
  }
  events_idx <- events_idx[keep]
}

# --- Valid events (respecting window) ---
n <- nrow(combined)
valid <- events_idx[(events_idx - win_L >= 1) & (events_idx + win_L <= n)]

# --- Function: lag per event for one variable ---
lag_from_window_ccf <- function(x, y, lag.max = NULL) {
  xz <- scale(x, center = TRUE, scale = TRUE)[,1]
  yz <- scale(y, center = TRUE, scale = TRUE)[,1]
  res <- ccf(xz, yz, plot = FALSE, lag.max = lag.max)
  lag <- res$lag[which.max(abs(res$acf))]
  val <- max(abs(res$acf))
  c(lag = lag, max_abs_acf = val)
}

vars_to_test <- c("S_Ti", "S", "Ca",'Ca_Ti', "Ti", "Fe_Ti",'Fe','Mn_Ti', 'Mn', 'Ba_Ti', 'Ba')

# --- Main loop ---
results_list <- lapply(vars_to_test, function(vn) {
  y <- combined_df[[vn]]
  
  # full-series CCF
  ccf_full <- ccf(gs, y, plot = FALSE)
  lag_full <- ccf_full$lag[which.max(abs(ccf_full$acf))]
  max_full <- max(abs(ccf_full$acf))
  
  # event-based CCFs
  ev_lags <- lapply(valid, function(i0) {
    idx <- (i0 - win_L):(i0 + win_L)
    out <- lag_from_window_ccf(gs[idx], y[idx])
    tibble(event_center = i0,
           lag_event = as.integer(out["lag"]),
           max_abs_acf_event = as.numeric(out["max_abs_acf"]))
  }) %>% bind_rows()
  
  summary <- ev_lags %>%
    summarise(
      n_events = n(),
      median_lag_event = median(lag_event, na.rm = TRUE),
      iqr_lag_event_low = quantile(lag_event, 0.25, na.rm = TRUE),
      iqr_lag_event_high = quantile(lag_event, 0.75, na.rm = TRUE),
      mean_abs_acf_event = mean(max_abs_acf_event, na.rm = TRUE)
    ) %>%
    mutate(
      variable = vn,
      lag_full = lag_full,
      max_abs_acf_full = max_full
    ) %>%
    relocate(variable)
  
  list(summary = summary, events = ev_lags, ccf_full = ccf_full)
})

# --- Summary table ---
lag_summary <- bind_rows(lapply(results_list, `[[`, "summary"))
lag_summary

# --- Plots ---
par(mfrow = c(2, 1)) # two plots per variable
for (vn in vars_to_test) {
  res <- results_list[[which(vars_to_test == vn)]]
  
  # Full-series CCF plot
  plot(res$ccf_full, main = paste("Full-series CCF: mean_gray &", vn))
  
  # Event lag histogram
  hist(res$events$lag_event,
       breaks = 20,
       main = paste("Event-wise lag distribution:", vn),
       xlab = "Lag (samples)", col = "gray")
}
par(mfrow = c(1,1))