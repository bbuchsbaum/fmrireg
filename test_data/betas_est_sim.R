# -------------------------------------------------------------------
#     0) Load libraries (or source them) and define needed objects
# -------------------------------------------------------------------
library(tidyverse)    # for data wrangling & plotting
#library(fmrireg)      # your main package that has matrix_dataset, estimate_betas, etc.
library(foreach)      # for loops or parallel
library(doParallel)   # optional, if you want parallel

# (Ensure your "simulate_fmri_timecourse" and all "estimate_betas" code is loaded.)

# Example: If you have the above definitions in e.g. "my_fmrireg_extensions.R", then do:
# source("my_fmrireg_extensions.R")

# -------------------------------------------------------------------
#     1) Define a parameter grid for simulation
# -------------------------------------------------------------------
# We'll define some "ISI modes" for demonstration:
#   - "even":  evenly spaced events (default)
#   - "dense":  1 event every 1s or 2s (simulate with uniform w/ isi_min=1, isi_max=1)
#   - "sparse": uniform w/ isi_min=4, isi_max=6
# Adjust as you like.

sim_grid <- expand.grid(
  amplitude_sd = c(.1, .2, .4, .6),
  noise_type   = c("white", "ar1", "ar2"),
  noise_sd     = c(.2, 0.5, 1.0, 2),
  isi_mode     = c("even", "dense", "sparse"),
  replicate    = 1:3,  # do multiple replicates
  stringsAsFactors = FALSE
)

# We'll define a helper to map "isi_mode" to the actual arguments
# for "simulate_fmri_timecourse" (isi_dist + isi_min/max).
# For "dense" => uniform with (isi_min=1, isi_max=1).
# For "sparse" => uniform with (isi_min=4, isi_max=6).
# For "even" => just use isi_dist="even".

get_isi_args <- function(mode) {
  if (mode == "dense") {
    list(isi_dist="uniform", isi_min=1, isi_max=1)
  } else if (mode == "sparse") {
    list(isi_dist="uniform", isi_min=2, isi_max=6)
  } else {
    # "even"
    list(isi_dist="even")
  }
}

# -------------------------------------------------------------------
#     2) Define the set of methods we want to test
# -------------------------------------------------------------------
# Single-trial methods we test for each dataset:
#   - "fracridge" (with fracs=0.5 as example)
#   - "lss"
#   - "r1" (Rank-1)
#   - "pls" (with ncomp=1..6)
#   - "ols"
#   - "mixed"
#   - "mixed_cpp"
#
# If any method references a function not in your code, remove it or
# define it similarly to the others.

methods_to_test <- c("lss", "pls", "ols", "mixed", "mixed_cpp")
pls_components  <- 1:6  # for method="pls" only

# Define default HRF basis and reference for rank-1 methods
# Simple canonical HRF + derivatives
create_default_hrf <- function() {
  t <- seq(0, 24, by=0.5)
  canonical <- dgamma(t, shape=6, rate=1) - 0.1 * dgamma(t, shape=16, rate=1)
  canonical <- canonical / max(canonical)
  
  # First derivative
  derivative1 <- c(diff(canonical), 0) 
  derivative1 <- derivative1 / max(abs(derivative1))
  
  # Second derivative  
  derivative2 <- c(diff(diff(canonical)), 0, 0)
  derivative2 <- derivative2 / max(abs(derivative2))
  
  # Return both basis and reference
  list(
    basis = cbind(canonical, derivative1, derivative2),
    reference = canonical
  )
}

# Create default HRF objects
default_hrf <- create_default_hrf()

# -------------------------------------------------------------------
#     3) Loop over the simulation grid
# -------------------------------------------------------------------
# We'll accumulate all results in a data.frame "sim_results".
# Each row = one combination of simulation + replicate + method
# and the resulting correlation with true amplitudes.

res_list <- list()
counter  <- 1

for (rowi in seq_len(nrow(sim_grid))) {
  # Extract parameters
  parrow <- sim_grid[rowi, ]
  amp_sd   <- parrow$amplitude_sd
  ntype    <- parrow$noise_type
  nsd      <- parrow$noise_sd
  imode    <- parrow$isi_mode
  
  isiargs  <- get_isi_args(imode)
  
  # Calculate total scan time with buffer at the end to ensure HRF returns to baseline
  total_time <- 240
  hrf_buffer <- 20  # Allow for HRF to return to baseline (~20s is usually sufficient)
  effective_time <- total_time - hrf_buffer
  
  # Adjust number of events based on ISI mode and available time
  if (imode == "dense") {
    # For dense ISI, use available time / average ISI
    n_events <- floor(effective_time / 1.5)  # ~100 events with 1-1.5s ISI
  } else if (imode == "sparse") {
    # For sparse ISI, use available time / average ISI
    n_events <- floor(effective_time / 3)  # ~40 events with 4-6s ISI
  } else {
    # For even spacing, use available time / n_events to get even spacing
    n_events <- 50
  }
  
  # Ensure there's at least minimal noise and amplitude variation for meaningful correlation
  if (amp_sd == 0) {
    # Add minimal amplitude variation to ensure correlations are meaningful
    message("Adding minimal amplitude variation (SD=0.2) to ensure meaningful correlations")
    amp_sd <- max(amp_sd, 0.2)
  }
  
  # -----------------------------------------------------------------
  # 3A) Simulate data
  # -----------------------------------------------------------------
  sim_out <- simulate_fmri(
    n = 3,                   # say we do 2 columns for demonstration
    total_time = total_time, # Now using variable with buffer
    TR = 1,                  # let's do 1s TR so "dense" is truly event each second
    n_events = n_events,     # Adjusted based on ISI mode and buffer
    amplitude_sd = amp_sd,
    amplitude_dist = "gaussian",
    noise_type = ntype,
    noise_sd = nsd,
    # pass the isiargs
    isi_dist = isiargs$isi_dist,
    isi_min = isiargs$isi_min %||% 2,
    isi_max = isiargs$isi_max %||% 2,
    # keep durations=0 for single events
    single_trial = FALSE,
    buffer=50,
    random_seed = 1 + 1000*rowi  # for reproducibility
  )
  
  # After simulation, check if any events are too close to the end and remove them
  # This ensures all HRFs have time to return to baseline
  dset <- sim_out$time_series
  event_table <- dset$event_table
 
  
  # sim_out contains:
  #   - time_series  = a matrix_dataset with dims 240 x 2
  #   - ampmat       = n_events x 2
  #   - durmat       = n_events x 2
  # The first column in "ampmat" matches the event_table for the dataset. 
  # But we want "true" single-trial amplitudes for each column if we do a single-trial approach.
  
  # However, we set single_trial=FALSE above. If we truly want single-trial "true" amplitudes,
  # we can do single_trial=TRUE. Then ampmat is still the same. We'll proceed anyway:
  # We'll treat ampmat as the ground truth for each event & column.
  
  # The dataset:
  dset <- sim_out$time_series
  # We want to do single-trial estimation => "onset ~ trialwise()"
  # so that "betas_ran" yields one beta per event.
  
  # We'll define:
  run_form <- onset ~ trialwise()
  
  bmod <- baseline_model(basis="poly", degree=3, sframe=dset$sampling_frame)
  
  # Build a "block" formula. If we had multiple runs, we might do "~ run", but we have 1 run only.
  # so block can be ~1
  block_form <- ~1
  
  # Add a constant regressor for fixed effects
  dset$event_table$constant <- factor(rep(1, nrow(dset$event_table)))
  
  # Define fixed effect formula - this is the standard approach taken in test_betas.R
  fixed_form <- onset ~ hrf(constant)
  
  # Keep track of original event count (for validation)
  original_event_count <- nrow(dset$event_table)
  
  # -----------------------------------------------------------------
  # 3B) For each method, estimate single-trial betas & compute correlation
  # -----------------------------------------------------------------
  for (mm in methods_to_test) {
    # Determine whether to use fixed effects based on method
    use_fixed <- !(mm %in% c("ols", "lss"))
    
    # If "pls", we vary ncomp from 1..6
    if (mm == "pls") {
      for (nc in pls_components) {
        # Estimate betas
        bet_est <- tryCatch({
          estimate_betas(
            dset,
            fixed = if(use_fixed) fixed_form else NULL, 
            ran   = run_form,
            block = block_form,
            method = mm,
            basemod=bmod,
            ncomp = nc
          )
        }, error = function(e) {
          message("Error in beta estimation for method ", mm, ", ncomp=", nc, ": ", e$message)
          return(NULL)
        })
        
        # Skip if estimation failed
        if (is.null(bet_est)) {
          res_list[[counter]] <- data.frame(
            amplitude_sd = amp_sd,
            noise_type   = ntype,
            noise_sd     = nsd,
            isi_mode     = imode,
            replicate    = parrow$replicate,
            method       = mm,
            ncomp        = nc,
            correlation  = NA
          )
          counter <- counter + 1
          next
        }
        
        # Verify dimensions match before calculating correlation
        est_events_count <- nrow(bet_est$betas_ran)
        true_events_count <- nrow(sim_out$ampmat)
        
        if (est_events_count != true_events_count) {
          message("Dimension mismatch: estimated betas has ", est_events_count, 
                  " events, but true amplitudes has ", true_events_count, " events")
          # Use only the number of events that match
          min_events <- min(est_events_count, true_events_count)
          
          # Calculate correlation using available events
          cvals <- numeric(ncol(bet_est$betas_ran))
          for (col_j in seq_len(ncol(bet_est$betas_ran))) {
            xtrue <- sim_out$ampmat[1:min_events, col_j, drop=FALSE]
            xest  <- bet_est$betas_ran[1:min_events, col_j, drop=FALSE]
            cvals[col_j] <- cor(xtrue, xest)
          }
        } else {
          # Calculate correlation normally
          cvals <- numeric(ncol(bet_est$betas_ran))
          for (col_j in seq_len(ncol(bet_est$betas_ran))) {
            xtrue <- sim_out$ampmat[, col_j]
            xest  <- bet_est$betas_ran[, col_j]
            cvals[col_j] <- cor(xtrue, xest)
          }
        }
        
        mean_corr <- mean(cvals, na.rm=TRUE)
        
        # store in result
        res_list[[counter]] <- data.frame(
          amplitude_sd = amp_sd,
          noise_type   = ntype,
          noise_sd     = nsd,
          isi_mode     = imode,
          replicate    = parrow$replicate,
          method       = mm,
          ncomp        = nc,
          correlation  = mean_corr
        )
        counter <- counter + 1
      }
    } else {
      # For other methods
      bet_est <- estimate_betas(
        dset,
        fixed = if(use_fixed) fixed_form else NULL,
        ran   = run_form,
        block = block_form,
        basemod=bmod,
        method = mm,
        hrf_basis = default_hrf$basis,
        hrf_ref = default_hrf$reference
      )
      
      cvals <- numeric(ncol(bet_est$betas_ran))
      for (col_j in seq_len(ncol(bet_est$betas_ran))) {
        xtrue <- sim_out$ampmat[, col_j]
        xest  <- bet_est$betas_ran[, col_j]
        cvals[col_j] <- cor(xtrue, xest)
      }
      mean_corr <- mean(cvals)
      
      res_list[[counter]] <- data.frame(
        amplitude_sd = amp_sd,
        noise_type   = ntype,
        noise_sd     = nsd,
        isi_mode     = imode,
        replicate    = parrow$replicate,
        method       = mm,
        ncomp        = NA,
        correlation  = mean_corr
      )
      counter <- counter + 1
    }
  } # end for mm
} # end for rowi

# -------------------------------------------------------------------
#     4) Combine results and do some plotting
# -------------------------------------------------------------------
sim_results <- dplyr::bind_rows(res_list)

# For convenience, let's define a method label that includes "pls_ncomp" for PLS
sim_results <- sim_results %>%
  mutate(method_label = if_else(
    method == "pls" & !is.na(ncomp),
    paste0("pls_", ncomp),
    method
  ))

# Let's do a quick summarise of the average correlation across replicates
summary_results <- sim_results %>%
  #group_by(amplitude_sd, noise_type, noise_sd, isi_mode, method_label) %>%
  group_by(amplitude_sd, noise_type, isi_mode, method_label) %>%
  summarize(
    mean_corr = mean(correlation, na.rm=TRUE),
    sd_corr   = sd(correlation, na.rm=TRUE),
    .groups   = "drop"
  )

totsum <- summary_results %>%
  group_by(method_label) %>%
  summarize(
    mean_corr = mean(mean_corr, na.rm=TRUE),
    sd_corr   = mean(sd_corr, na.rm=TRUE),
    .groups   = "drop"
  ) %>% arrange(desc(mean_corr))

# Identify the "winner" per group (the method with highest mean_corr)
winners <- summary_results %>%
  #group_by(amplitude_sd, noise_type, noise_sd, isi_mode) %>%
  group_by(amplitude_sd, noise_type,  isi_mode) %>%
  slice_max(mean_corr, n=1, with_ties = FALSE) %>%
  ungroup()

cat("\n=== Winners per parameter setting ===\n")
print(winners)

# Optionally, we can see the overall best method across *all* settings
overall_winner <- summary_results %>%
  group_by(method_label) %>%
  summarize(all_mean_corr = mean(mean_corr), .groups="drop") %>%
  slice_max(all_mean_corr, n=1)
cat("\n=== Overall best method across all param combos ===\n")
print(overall_winner)

# -------------------------------------------------------------------
#     5) Plot the results
# -------------------------------------------------------------------
# Example: facet by amplitude_sd, color by method_label, x-axis = noise_type
library(ggplot2)

ggplot(summary_results, aes(x=noise_type, y=mean_corr, color=method_label, group=method_label)) +
  geom_point() + geom_line() +
  facet_grid(amplitude_sd ~ isi_mode, labeller=label_both) +
  labs(
    title="Single-trial Beta Estimation Correlations",
    y="Mean Correlation (true vs estimated betas)",
    x="Noise Type"
  ) +
  theme_bw()