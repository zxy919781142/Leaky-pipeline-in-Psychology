# ============================================================
# Title: Fig.1 : Statistical description of factors associated with 
# academic performance, collaboration, and institutional affiliations 
# Project: Leaky Pipeline in Psychology
# Author: Xinyi Zhao  |  Max Planck Institute for Human Development
# Date: 2025-09-22
# R version: >= 4.3
# Description: 
# Inputs: survival_data.csv (time_varying variables for each researcher)
# Outputs: fig_3.pdf and fig_4.pdf
# ============================================================

# ---------- Packages ----------
install_and_load <- function(packages, quietly = TRUE) {
  # install if missing, then load; suppress startup noise by default
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    if (quietly) {
      suppressPackageStartupMessages(
        library(pkg, character.only = TRUE)
      )
    } else {
      library(pkg, character.only = TRUE)
    }
  }
  invisible(NULL)
}

core_data <- c(
  "dplyr","tidyr","tibble","readr","stringr","janitor","purrr","here","httr2","fastDummies"
)

viz <- c(
  "ggplot2","ggrepel","scales","patchwork","viridisLite","gridExtra","GGally"
)


# Inference, tables, EDA
inference_reporting <- c(
  "margins","marginaleffects","lmtest","sandwich","gtsummary","sjPlot","skimr","DataExplorer","forestmodel"
)

# ---- Install & load (dedupe) ----
pkgs <- unique(c(core_data, viz, inference_reporting ))
install_and_load(pkgs)

# Helpful defaults for Stan (safe no-ops if rstan missing)
if ("rstan" %in% rownames(installed.packages())) {
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
}

# ----   TVEM-related functions ----
# =============================================================================
# tvem_tiv_result(): Tidy time-varying effects from a TVEM model
# =============================================================================

tvem_tv_result <- function(model, exponentiate = TRUE) {
  # Prefer names from the model frame (if present); otherwise, list names:
  tv_var_names <- try({
    unique(model.frame(model)$varying_effects_names)
  }, silent = TRUE)
  
  if (inherits(tv_var_names, "try-error") || is.null(tv_var_names)) {
    tv_var_names <- names(model$grid_fitted_coefficients)
  }
  
  if (is.null(tv_var_names) || length(tv_var_names) == 0L) {
    stop("No time-varying effects found in 'model'.")
  }
  
  # Build rows for each time-varying term
  out_list <- lapply(tv_var_names, function(v1) {
    comp <- model$grid_fitted_coefficients[[v1]]
    if (is.null(comp)) return(NULL)
    
    # Raw (log-scale)
    est_raw   <- comp$estimate
    se_raw    <- comp$standard_error
    up_raw    <- comp$upper
    lo_raw    <- comp$lower
    
    # Transformed (if requested)
    est <- if (isTRUE(exponentiate)) exp(est_raw) else est_raw
    up  <- if (isTRUE(exponentiate)) exp(up_raw)  else up_raw
    lo  <- if (isTRUE(exponentiate)) exp(lo_raw)  else lo_raw
    
    data.frame(
      estimate     = est,
      upper        = up,
      lower        = lo,
      estimate_raw = est_raw,
      std_raw      = se_raw,
      upper_raw    = up_raw,
      lower_raw    = lo_raw,
      tv_var       = v1,
      time_grid    = model$time_grid,
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, Filter(Negate(is.null), out_list))
  rownames(out) <- NULL
  out
}


# =============================================================================
# tvem_tiv_result(): Tidy time-invariant effects from a TVEM model.
# =============================================================================
tvem_tiv_result <- function(model, exponentiate = TRUE, conf_level = 0.95) {
  mat <- model$invar_effects_estimates
  if (is.null(mat)) stop("'model$invar_effects_estimates' not found.")
  
  # Determine term order
  tiv_var_names <- try({
    unique(model.frame(model)$invar_effects_names)
  }, silent = TRUE)
  
  if (inherits(tiv_var_names, "try-error") || is.null(tiv_var_names)) {
    tiv_var_names <- rownames(mat)
  }
  
  if (is.null(tiv_var_names) || length(tiv_var_names) == 0L) {
    stop("No time-invariant effects found in 'model'.")
  }
  
  # z-quantile for CI
  alpha <- 1 - conf_level
  z_q   <- stats::qnorm(1 - alpha / 2)
  
  out_list <- lapply(tiv_var_names, function(v1) {
    est <- as.numeric(mat[v1, "estimate"])
    se  <- as.numeric(mat[v1, "standard_error"])
    if (is.na(est) || is.na(se)) return(NULL)
    
    # Raw CI (log-scale)
    up_raw <- est + z_q * se
    lo_raw <- est - z_q * se
    
    # Transform if requested
    est_tr <- if (isTRUE(exponentiate)) exp(est)   else est
    up_tr  <- if (isTRUE(exponentiate)) exp(up_raw) else up_raw
    lo_tr  <- if (isTRUE(exponentiate)) exp(lo_raw) else lo_raw
    
    # Wald z and (two-sided) p-value on raw scale
    z_val <- est / se
    p_val <- 2 * stats::pnorm(abs(z_val), lower.tail = FALSE)
    
    data.frame(
      estimate     = est_tr,
      upper        = up_tr,
      lower        = lo_tr,
      estimate_raw = est,
      std_raw      = se,
      z            = z_val,
      p            = p_val,
      tiv_var      = v1,
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, Filter(Negate(is.null), out_list))
  rownames(out) <- NULL
  out
}


# =============================================================================
# tvem_cof_importance(): build a wide matrix of TV (time-varying) and TIV
# (time-invariant) coefficients across the model's time grid, with optional
# per-time (row-wise) standardization across variables.
#
# Notes:
# - Row-wise standardization (standard = TRUE) rescales each row (time point)
#   by subtracting the row mean and dividing by the row SD across variables.
#   This keeps relative *importance within a time point* comparable.
# - Time-invariant effects are repeated across all time points.
# - Returns a data.frame with 'time_grid' + one column per coefficient.
# =============================================================================
tvem_cof_importance <- function(model, standard = TRUE) {
  # Prefer names from model.frame (if available), else fall back to model slots
  mf <- try(model.frame(model), silent = TRUE)
  tv_var_names <- try(unique(mf$varying_effects_names), silent = TRUE)
  tiv_var_names <- try(unique(mf$invar_effects_names),  silent = TRUE)
  
  if (inherits(tv_var_names, "try-error") || all(is.na(tv_var_names))) {
    tv_var_names <- names(model$grid_fitted_coefficients)
    # Remove intercept if present in grid_fitted_coefficients
    tv_var_names <- setdiff(tv_var_names, "(Intercept)")
  }
  if (inherits(tiv_var_names, "try-error")) {
    tiv_var_names <- rownames(model$invar_effects_estimates)
  }
  
  time_grid <- model$time_grid
  if (is.null(time_grid)) stop("'model$time_grid' is missing.")
  
  # ---- Build TV part (one column per time-varying covariate) ----
  tv_df <- NULL
  if (length(tv_var_names) > 0) {
    tv_cols <- lapply(tv_var_names, function(v) {
      comp <- model$grid_fitted_coefficients[[v]]
      if (is.null(comp) || is.null(comp$estimate)) {
        stop("Missing grid_fitted_coefficients for TV var: ", v)
      }
      as.numeric(comp$estimate)
    })
    tv_df <- as.data.frame(tv_cols, stringsAsFactors = FALSE)
    colnames(tv_df) <- tv_var_names
  }
  
  # ---- Attach TIV part (repeat each constant across time) ----
  tiv_df <- NULL
  if (length(tiv_var_names) > 0 && any(!is.na(tiv_var_names))) {
    tiv_cols <- lapply(tiv_var_names, function(v) {
      est <- model$invar_effects_estimates[v, "estimate"]
      if (is.na(est)) est <- NA_real_
      rep(est, length(time_grid))
    })
    tiv_df <- as.data.frame(tiv_cols, stringsAsFactors = FALSE)
    colnames(tiv_df) <- tiv_var_names
  }
  
  # ---- Combine pieces ----
  if (is.null(tv_df) && is.null(tiv_df)) {
    stop("No coefficients found: neither TV nor TIV effects appear to be present.")
  }
  coef_df <- if (!is.null(tv_df) && !is.null(tiv_df)) {
    cbind(tv_df, tiv_df)
  } else if (!is.null(tv_df)) {
    tv_df
  } else {
    tiv_df
  }
  coef_df$time_grid <- time_grid
  # Put time_grid first
  coef_df <- coef_df[, c("time_grid", setdiff(names(coef_df), "time_grid")), drop = FALSE]
  
  # ---- Optional row-wise standardization across variables ----
  if (isTRUE(standard)) {
    # Requires dplyr ≥ 1.0
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      warning("Package 'dplyr' not available; returning unstandardized coefficients.")
      return(coef_df)
    }
    coef_df <- dplyr::rowwise(coef_df) |>
      dplyr::mutate(
        dplyr::across(
          -time_grid,
          ~ {
            row_vals <- dplyr::c_across(-time_grid)
            mu <- mean(row_vals, na.rm = TRUE)
            s  <- stats::sd(row_vals, na.rm = TRUE)
            if (!is.finite(s) || s == 0) 0 else ((.x - mu) / s)
          }
        )
      ) |>
      dplyr::ungroup()
  }
  
  coef_df
}


# =============================================================================
# tvem_predict(): cohort-level predictions over academic ages  (aggregated
# across individuals at the same age).
# Notes:
# - We align coefficients to time via model$time_grid. We **assume**
#   time_grid indexes the ages 1,2,...,max(time_grid). If your grid is
#   non-integer or not equal to academic_age, add a mapping step.
# - Missing predictors in `data` are treated as 0 (common with one-hot vars).
# - Intercept is taken from model$grid_fitted_coefficients$`(Intercept)`$estimate
# =============================================================================
tvem_predict <- function(model, data, maxage, var = NULL, format = "exp") {
  # --- Collect names and time grid ---
  mf <- try(model.frame(model), silent = TRUE)
  tv_names  <- try(unique(mf$varying_effects_names), silent = TRUE)
  tiv_names <- try(unique(mf$invar_effects_names),  silent = TRUE)
  
  if (inherits(tv_names, "try-error") || all(is.na(tv_names))) {
    tv_names <- names(model$grid_fitted_coefficients)
    tv_names <- setdiff(tv_names, "(Intercept)")
  }
  if (inherits(tiv_names, "try-error")) {
    tiv_names <- rownames(model$invar_effects_estimates)
  }
  
  time_grid <- model$time_grid
  if (is.null(time_grid)) stop("'model$time_grid' is missing.")
  
  # --- Build coefficient matrix at each time row ---
  # TV part:
  tv_mat <- NULL
  if (length(tv_names) > 0) {
    tv_mat <- sapply(tv_names, function(v) {
      comp <- model$grid_fitted_coefficients[[v]]
      if (is.null(comp) || is.null(comp$estimate)) {
        stop("Missing grid_fitted_coefficients for TV var: ", v)
      }
      as.numeric(comp$estimate)
    })
    # Ensure matrix shape even for single column
    tv_mat <- as.matrix(tv_mat)
    colnames(tv_mat) <- tv_names
  }
  
  # TIV part:
  tiv_mat <- NULL
  if (length(tiv_names) > 0 && any(!is.na(tiv_names))) {
    tiv_vec <- sapply(tiv_names, function(v) {
      model$invar_effects_estimates[v, "estimate"]
    })
    tiv_mat <- matrix(
      rep(tiv_vec, each = length(time_grid)),
      nrow = length(time_grid), byrow = FALSE,
      dimnames = list(NULL, tiv_names)
    )
  }
  
  # Combine coefficient matrix (no intercept)
  coef_mat <- NULL
  if (!is.null(tv_mat) && !is.null(tiv_mat)) {
    coef_mat <- cbind(tv_mat, tiv_mat)
  } else if (!is.null(tv_mat)) {
    coef_mat <- tv_mat
  } else if (!is.null(tiv_mat)) {
    coef_mat <- tiv_mat
  } else {
    stop("No coefficients found for prediction (neither TV nor TIV).")
  }
  
  # Intercept over time
  intercept <- model$grid_fitted_coefficients[["(Intercept)"]][["estimate"]]
  if (is.null(intercept)) stop("Intercept not found in model$grid_fitted_coefficients.")
  
  # --- Determine ages to predict (keep your original starting at 2) ---
  # Use ages present in data, <= maxage, and present in the time_grid.
  ages_in_data <- sort(unique(data$academic_age))
  ages <- ages_in_data[ages_in_data <= maxage & ages_in_data >= 2]
  if (length(ages) == 0) {
    stop("No academic ages >= 2 and <= maxage found in `data`.")
  }
  
  # Ensure predictors in `data` cover all coefficient columns
  needed <- colnames(coef_mat)
  missing_cols <- setdiff(needed, colnames(data))
  if (length(missing_cols) > 0) {
    # fill missing predictors with 0
    for (m in missing_cols) data[[m]] <- 0
  }
  
  # Pre-allocate results
  out <- vector("list", length(ages))
  
  # --- Predict loop over ages ---
  for (k in seq_along(ages)) {
    a <- ages[k]
    
    # Find the row in time_grid corresponding to age 'a'
    ri <- which(time_grid == a)
    if (length(ri) != 1) {
      # If no exact match, we could nearest-neighbor; for now, skip safely.
      warning("Age ", a, " not found exactly in model$time_grid; skipping.")
      next
    }
    
    beta_vec <- coef_mat[ri, , drop = TRUE]              # coefficients at age a
    X <- as.matrix(data[data$academic_age == a, needed, drop = FALSE])
    # Order columns to match beta order
    X <- X[, needed, drop = FALSE]
    
    eta <- as.numeric(X %*% beta_vec) + intercept[ri]
    
    if (identical(format, "origin")) {
      mean_pred <- mean(eta, na.rm = TRUE)
      sd_pred   <- stats::sd(eta, na.rm = TRUE)
      if (!is.finite(sd_pred) || is.na(sd_pred)) sd_pred <- 0
      lower_ci <- mean_pred - 1.96 * sd_pred
      upper_ci <- mean_pred + 1.96 * sd_pred
    } else {
      p <- 1 / (1 + exp(-eta))  # logistic
      mean_pred <- mean(p, na.rm = TRUE)
      # Match your original code: use IQR-style bounds as a robust interval
      qs <- stats::quantile(p, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
      lower_ci <- qs[1]
      upper_ci <- qs[2]
    }
    
    out[[k]] <- data.frame(
      academic_age = a,
      Mean        = mean_pred,
      CI_Lower    = lower_ci,
      CI_Upper    = upper_ci,
      row.names   = NULL
    )
  }
  
  # Bind and drop NULLs from skipped ages (if any)
  out <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(out) || nrow(out) == 0) {
    stop("No predictions could be computed; check time_grid vs data$academic_age.")
  }
  out
}

# ---------- Global options (safe defaults for all your figures) ----------
set.seed(42)
options(stringsAsFactors = FALSE, scipen = 999)
theme_set(theme_minimal(base_size = 14))





# ----   Read & basic cleanup ----
csv_path <- here::here("data", "survival_data_modified_filtered.csv")
survival_data_modified_filtered <- readr::read_csv(csv_path, show_col_types = FALSE) %>%
  janitor::clean_names()  # snake_case, ASCII-safe names


# ============================================================
# 1. TVEM by cohorts 
# ============================================================


# --- Parameters ---------------------------------------------------------------
end_year  <- 2018   # truncation year for feasible academic age
num_knots <- 4L     # spline knots for TVEM
min_rows  <- 10     # skip tiny cohorts (edit as needed)
min_ids   <- 5

# --- Quick sanity check -------------------------------------------------------
required_cols <- c(
  "cohort_1", "cohort", "academic_age", "scopus_author_id", "leave",
  "gender_female",
  "sum_fncr_5_years","mean_mncs_mean_yearly","affili_total_mncs_WT",
  "nunique_coauthors_yearly","firstauthor_num","lastauthor_num",
  "coauthor_female_prop","unique_affi",
  "continent_new_North_America","continent_new_Europe",
  "early_field_1_Applied_Psychology","early_field_1_Clinical_Psychology",
  "early_field_1_Developmental_and_Educational_Psychology",
  "early_field_1_Experimental_and_Cognitive_Psychology",
  "early_field_1_Neuropsychology_and_Physiological_Psychology",
  "early_field_1_Social_Psychology","early_field_1_Psychology_all"
)
missing_req <- setdiff(required_cols, colnames(survival_data_modified_filtered))
if (length(missing_req) > 0) {
  stop("Missing required columns: ", paste(missing_req, collapse = ", "))
}

# --- Prepare containers  ------------------------
model_tvem_list            <- list()
.model_tv_list             <- list()
.model_tiv_list            <- list()
.imp_std_list              <- list()
.imp_raw_list              <- list()
.pred_gender_list          <- list()
.pred_gender_ctr_list      <- list()

# --- loop over cohorts ---------------------------------------------------
for (c in sort(unique(survival_data_modified_filtered$cohort_1))) {
  
  data <- survival_data_modified_filtered %>%
    filter(cohort_1 == c) %>%
    # Only keep feasible ages given end_year truncation
    filter(academic_age <= (end_year - cohort + 1))
  
  # Skip undersized cohorts to avoid unstable fits
  if (nrow(data) < min_rows || length(unique(data$scopus_author_id)) < min_ids) {
    message("Skipping cohort ", c, " (n=", nrow(data), ", ids=", length(unique(data$scopus_author_id)), ").")
    next
  }
  
  max_age <- max(data$academic_age, na.rm = TRUE)
  
  # --- Fit TVEM ---------------------------------------------------------------
  model <- tvem(
    data = data,
    formula = leave ~ gender_female +
      sum_fncr_5_years + mean_mncs_mean_yearly + affili_total_mncs_WT +
      nunique_coauthors_yearly + firstauthor_num + lastauthor_num +
      coauthor_female_prop + unique_affi,
    invar_effect = ~ continent_new_North_America + continent_new_Europe +
      early_field_1_Applied_Psychology + early_field_1_Clinical_Psychology +
      early_field_1_Developmental_and_Educational_Psychology +
      early_field_1_Experimental_and_Cognitive_Psychology +
      early_field_1_Neuropsychology_and_Physiological_Psychology +
      early_field_1_Social_Psychology + early_field_1_Psychology_all,
    family    = binomial(),
    id        = scopus_author_id,
    time      = academic_age,
    num_knots = num_knots,
    grid      = max_age
  )
  
  # Keep the fitted model
  model_tvem_list[[as.character(c)]] <- model
  
  # --- Summaries --------------------------------------------------------------
  model_tv_result_cb  <- tvem_tv_result(model)  %>% mutate(cohort = c)
  model_tiv_result_cb <- tvem_tiv_result(model) %>% mutate(cohort = c)
  
  tvem_importance_cb          <- tvem_cof_importance(model, TRUE)  %>% mutate(cohort = c)
  tvem_importance_cb_nonstand <- tvem_cof_importance(model, FALSE) %>% mutate(cohort = c)
  
  # --- Predictions: overall by gender ----------------------------------------
  max_pred <- end_year - as.numeric(c) + 1
  
  tvem_predict_female <- tvem_predict(
    model,
    data %>% mutate(gender_female = 1),
    max_pred, "exp"                       # keep your original call signature
  ) %>% drop_na() %>% mutate(cohort = c, gender = "female")
  
  tvem_predict_male <- tvem_predict(
    model,
    data %>% mutate(gender_female = 0),
    max_pred, "exp"
  ) %>% drop_na() %>% mutate(cohort = c, gender = "male")
  
  tvem_predict_bothgender <- bind_rows(tvem_predict_female, tvem_predict_male)
  
  # --- Predictions: gender × continent ----------------------------------------
  tvem_predict_female_USA <- tvem_predict(
    model, data %>% mutate(gender_female = 1, continent_new_North_America = 1, continent_new_Europe = 0),
    max_pred, "exp"
  ) %>% drop_na() %>% mutate(cohort = c, gender = "female", ctr = "North_America")
  
  tvem_predict_male_USA <- tvem_predict(
    model, data %>% mutate(gender_female = 0, continent_new_North_America = 1, continent_new_Europe = 0),
    max_pred, "exp"
  ) %>% drop_na() %>% mutate(cohort = c, gender = "male", ctr = "North_America")
  
  tvem_predict_female_Euro <- tvem_predict(
    model, data %>% mutate(gender_female = 1, continent_new_North_America = 0, continent_new_Europe = 1),
    max_pred, "exp"
  ) %>% drop_na() %>% mutate(cohort = c, gender = "female", ctr = "Europe")
  
  # fixed the stray ",," in your original line
  tvem_predict_male_Euro <- tvem_predict(
    model, data %>% mutate(gender_female = 0, continent_new_North_America = 0, continent_new_Europe = 1),
    max_pred, "exp"
  ) %>% drop_na() %>% mutate(cohort = c, gender = "male", ctr = "Europe")
  
  tvem_predict_female_other <- tvem_predict(
    model, data %>% mutate(gender_female = 1, continent_new_North_America = 0, continent_new_Europe = 0),
    max_pred, "exp"
  ) %>% drop_na() %>% mutate(cohort = c, gender = "female", ctr = "Others")
  
  tvem_predict_male_other <- tvem_predict(
    model, data %>% mutate(gender_female = 0, continent_new_North_America = 0, continent_new_Europe = 0),
    max_pred, "exp"
  ) %>% drop_na() %>% mutate(cohort = c, gender = "male", ctr = "Others")
  
  tvem_predict_bothgender_ctr <- bind_rows(
    tvem_predict_female_USA, tvem_predict_male_USA,
    tvem_predict_female_Euro, tvem_predict_male_Euro,
    tvem_predict_female_other, tvem_predict_male_other
  )
  
  # --- Collect pieces (append to lists) ---------------------------------------
  .model_tv_list[[length(.model_tv_list) + 1]]        <- model_tv_result_cb
  .model_tiv_list[[length(.model_tiv_list) + 1]]      <- model_tiv_result_cb
  .imp_std_list[[length(.imp_std_list) + 1]]          <- tvem_importance_cb
  .imp_raw_list[[length(.imp_raw_list) + 1]]          <- tvem_importance_cb_nonstand
  .pred_gender_list[[length(.pred_gender_list) + 1]]  <- tvem_predict_bothgender
  .pred_gender_ctr_list[[length(.pred_gender_ctr_list) + 1]] <- tvem_predict_bothgender_ctr
}

# --- Bind everything into your original object names --------------------------
model_tv_result            <- bind_rows(.model_tv_list)
model_tiv_result           <- bind_rows(.model_tiv_list)
tvem_importance            <- bind_rows(.imp_std_list)
tvem_importance_nonstand   <- bind_rows(.imp_raw_list)
tvem_predict_result        <- bind_rows(.pred_gender_list)
tvem_predict_result_ctr    <- bind_rows(.pred_gender_ctr_list)




# ============================================================
# #2.  Importance summaries + bars/heatmaps by cohort
# Inputs expected in env:
#   - model_tv_result          (from tvem_tv_result + bind)
#   - model_tiv_result         (from tvem_tiv_result + bind)
#   - tvem_importance_nonstand (from tvem_cof_importance(..., standard = FALSE))
#   - variable_map: data.frame with columns:
#         variable (model term name), label (pretty label for plotting)
# Outputs:
#   - importance_summary  (cohort-level mean importance per variable)
#   - importance_standardize (long panel for heatmaps)
#   - plot (list of ggplot objects: bars then heatmaps per cohort)
# ============================================================
# ------------------------------------------------------------
# 1) Time-varying (TV) effects: average across ages 2..19
# ------------------------------------------------------------
tvem_importance_1 <- model_tv_result %>%
  filter(time_grid >= 2, time_grid <= 19) %>%
  select(cohort, time_grid, tv_var, estimate_raw) %>%
  rename(var = tv_var)

# Summarize per cohort × variable on the RAW scale (log-odds/hazard scale)
importance_summary <- tvem_importance_1 %>%
  group_by(cohort, var) %>%
  summarise(
    mean  = mean(estimate_raw, na.rm = TRUE),
    sd    = sd(estimate_raw,   na.rm = TRUE),
    n     = sum(!is.na(estimate_raw)),
    se    = sd / sqrt(pmax(n, 1)),                # guard n=0
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  ) %>%
  arrange(cohort, desc(mean))

# ------------------------------------------------------------
# 2) Time-invariant (TIV) effects: append with CIs from SE
# ------------------------------------------------------------
# tvem_tiv_result had: estimate_raw, std_raw, lower, upper already,
# but we recompute to ensure consistent columns.
tiv_append <- model_tiv_result %>%
  transmute(
    cohort,
    var   = tiv_var,
    mean  = estimate_raw,
    sd    = std_raw,
    n     = 1L,
    se    = sd,
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se
  )

importance_summary <- bind_rows(importance_summary, tiv_append)

# ------------------------------------------------------------
# 3) Add pretty labels and normalize absolute importance
# ------------------------------------------------------------
# If variable_map is missing or lacks columns, fall back to var as label
if (!exists("variable_map") ||
    !all(c("variable", "label") %in% names(variable_map))) {
  warning("variable_map not found or missing required columns; using raw variable names as labels.")
  variable_map <- data.frame(variable = unique(importance_summary$var),
                             label = unique(importance_summary$var),
                             stringsAsFactors = FALSE)
}

importance_summary <- importance_summary %>%
  left_join(variable_map, by = c("var" = "variable")) %>%
  mutate(
    label = if_else(is.na(label), var, label),
    mean_abs = abs(mean),
    dir = if_else(mean > 0, "p", "n")
  ) %>%
  group_by(cohort) %>%
  mutate(
    # Min-max on |mean| within cohort; safe if flat (denom=0)
    .den = max(mean_abs, na.rm = TRUE) - min(mean_abs, na.rm = TRUE),
    importance_normalized = if_else(.den > 0,
                                    (mean_abs - min(mean_abs, na.rm = TRUE)) / .den,
                                    0),
    .den = NULL
  ) %>%
  arrange(cohort, desc(importance_normalized)) %>%
  mutate(label = factor(label, levels = rev(unique(label)))) %>%
  ungroup()

# ------------------------------------------------------------
# 4) Color mapping (identity scale) for bars
#    Blue for negative, Orange/Red for positive; intensity = normalized
# ------------------------------------------------------------
pal_neg <- gradient_n_pal(c("#d1e5f0", "#2166ac"))  # light → dark blue
pal_pos <- gradient_n_pal(c("#fee8c8", "#e34a33"))  # light → dark orange/red

importance_summary <- importance_summary %>%
  mutate(
    A_rescaled = rescale(importance_normalized, to = c(0.3, 1)),   # optional opacity/intensity
    fill_color = if_else(dir == "n",
                         pal_neg(importance_normalized),
                         pal_pos(importance_normalized))
  )

# ------------------------------------------------------------
# 5) BAR PLOTS for selected cohorts
# ------------------------------------------------------------
target_cohorts <- c(2000, 2005, 2010)
plot <- vector("list", length = length(target_cohorts) * 2)  # bars + heatmaps

for (k in seq_along(target_cohorts)) {
  i <- target_cohorts[k]
  
  df_bar <- importance_summary %>% filter(cohort == i)
  
  plot[[2*k - 1]] <- ggplot(df_bar,
                            aes(x = importance_normalized,
                                y = reorder(label, importance_normalized))) +
    geom_col(aes(fill = fill_color), width = 0.6) +
    scale_fill_identity(guide = "none") +
    scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1.01)) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x  = element_text(size = 18, color = "black"),
      axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
      axis.title.x = element_text(size = 20, face = "bold", color = "black"),
      axis.title.y = element_text(size = 20, face = "bold", color = "black"),
      strip.text   = element_text(size = 22, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray80"),
      panel.grid.minor   = element_blank(),
      panel.background   = element_blank(),
      plot.background    = element_blank(),
      axis.line.x        = element_line(color = "black"),
      legend.position    = "none"
    )
}

# ------------------------------------------------------------
# 6) HEATMAP data: scale within (cohort × label) across time
# ------------------------------------------------------------
rescale_to_minus1_1 <- function(x) x / max(abs(x), na.rm = TRUE)

importance_wide <- tvem_importance_nonstand %>%
  # keep ages that are observable given end-year truncation per cohort
  filter(time_grid >= 2, time_grid <= (2018 - as.numeric(cohort) + 1))

# LONG format (modern replacement for melt)
importance_standardize <- importance_wide %>%
  pivot_longer(
    cols = -c(cohort, time_grid),
    names_to = "var",
    values_to = "importance"
  ) %>%
  left_join(variable_map, by = c("var" = "variable")) %>%
  mutate(label = if_else(is.na(label), var, label)) %>%
  group_by(cohort, label) %>%
  mutate(
    importance_scaled   = rescale_to_minus1_1(importance),
    importance_standard = as.numeric(scale(importance)),  # z over time
    importance_exp      = exp(importance)
  ) %>%
  ungroup()

# numeric cohorts for consistent joins/orderings
importance_standardize <- importance_standardize %>%
  mutate(cohort = as.numeric(cohort))
importance_summary <- importance_summary %>%
  mutate(cohort = as.numeric(cohort))

# ------------------------------------------------------------
# 7) HEATMAPS (one per target cohort, ordered by bar-plot ranking)
# ------------------------------------------------------------
for (k in seq_along(target_cohorts)) {
  i <- target_cohorts[k]
  
  data_hm <- importance_standardize %>% filter(cohort == i)
  
  # order labels using cohort’s bar-plot importance ranks
  label_order <- importance_summary %>%
    filter(cohort == i) %>%
    arrange(desc(importance_normalized)) %>%
    pull(label) %>%
    unique()
  
  data_hm <- data_hm %>%
    mutate(label = factor(label, levels = rev(label_order)))
  
  plot[[2*k]] <- ggplot(data_hm,
                        aes(x = time_grid, y = label, fill = importance_scaled)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("#2166ac", "#d1e5f0", "#fee8c8", "#e34a33"),
      values = rescale(c(-1, -0.01, 0.01, 1)),
      limits = c(-1, 1),
      name = "Importance\nacross time"
    ) +
    labs(title = NULL, x = "Academic age", y = "Variable") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x  = element_text(size = 18, color = "black"),
      axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
      axis.title.x = element_text(size = 20, face = "bold", color = "black"),
      axis.title.y = element_text(size = 20, face = "bold", color = "black"),
      strip.text   = element_text(size = 22, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray80"),
      panel.grid.minor   = element_blank(),
      panel.background   = element_blank(),
      plot.background    = element_blank(),
      axis.line          = element_line(color = "black"),
      legend.position    = "none"  # change to "right" if you want to show it
    )
}


# ------------------------------------------------------------
# 8) COmbine plots of different cohorts
# ------------------------------------------------------------

plot[[1]] <- plot[[1]] + ggtitle("Cohort: 2000–2004")
plot[[3]] <- plot[[3]] + ggtitle("Cohort: 2005–2009")
plot[[5]] <- plot[[5]] + ggtitle("Cohort: 2010–2014")
theme_update(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# build rows (bar | heatmap), then stack rows
row1 <- plot[[1]] | plot[[2]]
row2 <- plot[[3]] | plot[[4]]
row3 <- plot[[5]] | plot[[6]]

# combine; tweak relative column widths if needed
fig_all <- (row1 / row2 / row3) +
  plot_layout(widths = c(1, 1.15), heights = c(1, 1, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(10, 10, 10, 10)))


# save once as PNG and PDF
ggsave(here::here("result","fig_3.png"), fig_all, width = 12, height = 16, dpi = 300)
ggsave(here::here("result","fig_3.pdf"), fig_all, width = 12, height = 16, device = grDevices::cairo_pdf)



# =============================================================================
# #3. Peak ages + prediction plots (overall & by continent)
# Inputs:
#   tvem_predict_result      : cols [cohort, academic_age, gender, Mean, CI_Lower, CI_Upper]
#   tvem_predict_result_ctr  : cols above + [ctr]
# Outputs:
#   peak_ages, peak_ages_ctr, gender_gap_cohort, plot_predict, plot_predict_ctr
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# --- palettes -----------------------------------------------------------------
pal_gender <- c(male = "#3594cc", female = "#ea801c")
pal_gender_ctr <- c(
  "male.North_America"   = "#003366",
  "female.North_America" = "#a33d00",
  "male.Europe"          = "#3594cc",
  "female.Europe"        = "#ea801c",
  "male.Others"          = "#7bb8d9",
  "female.Others"        = "#ffc773"
)

# --- peak ages (overall by gender within cohort) ------------------------------
peak_ages <- tvem_predict_result %>%
  group_by(cohort, gender) %>%
  slice_max(order_by = Mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cohort, gender, academic_age)


# ------------------------------------------------------------
# 1) peak ages (by gender x continent within cohort)
# ------------------------------------------------------------
peak_ages_ctr <- tvem_predict_result_ctr %>%
  group_by(cohort, gender, ctr) %>%
  slice_max(order_by = Mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cohort, gender, ctr, academic_age)


# ------------------------------------------------------------
# 2) overall plot (female vs male)
# ------------------------------------------------------------
plot_predict <- ggplot(
  tvem_predict_result,
  aes(x = academic_age, y = Mean, group = gender, color = gender, fill = gender)
) +
  facet_wrap(~ cohort, ncol = 3) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.30, color = NA) + # band first
  geom_line(size = 1.2) +
  # vertical line at cohort- & gender-specific peak age
  geom_vline(
    data = peak_ages,
    aes(xintercept = academic_age),
    linetype = "dotted", size = 1.2, color = "black", show.legend = FALSE
  ) +
  scale_color_manual(values = pal_gender) +
  scale_fill_manual(values  = pal_gender) +
  scale_x_continuous(breaks = c(2, 5, 10, 15)) +
  labs(x = "Academic age", y = "Predicted risk of leaving academia", title = NULL) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    strip.text   = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80"),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_blank(),
    axis.line          = element_line(color = "black"),
    legend.position    = "none"
  )


# ------------------------------------------------------------
# 3) quantify gender gap (female vs male)
# ------------------------------------------------------------
# Safer join keyed by cohort + age; use pivot_wider for clarity
gender_gap_cohort <- tvem_predict_result %>%
  select(cohort, academic_age, gender, Mean) %>%
  tidyr::pivot_wider(names_from = gender, values_from = Mean, names_prefix = "mean_") %>%
  mutate(gender_gap = mean_male - mean_female) %>%
  arrange(cohort, academic_age)

# --- continent-specific plot --------------------------------------------------
# NOTE: The vertical lines below use overall peaks (from `peak_ages`).
# If you want continent-specific peaks, swap to `data = peak_ages_ctr` and
# map aes(color/fill) to interaction(gender, ctr) accordingly.
plot_predict_ctr <- ggplot(
  tvem_predict_result_ctr,
  aes(
    x = academic_age, y = Mean,
    group = interaction(gender, ctr),
    color = interaction(gender, ctr),
    fill  = interaction(gender, ctr)
  )
) +
  facet_wrap(~ cohort, ncol = 3) +
  geom_line(size = 1.2) +
  geom_vline(
    data = peak_ages_ctr,
    aes(xintercept = academic_age),
    linetype = "dotted", size = 1.2, color = "black", show.legend = FALSE
  ) +
  scale_color_manual(values = pal_gender_ctr) +
  scale_fill_manual(values  = pal_gender_ctr) +
  scale_x_continuous(breaks = c(2, 5, 10, 15)) +
  labs(x = "Academic age", y = "Predicted risk of leaving academia", title = NULL) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    strip.text   = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80"),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_blank(),
    axis.line          = element_line(color = "black"),
    legend.position    = c(0.97, 0.22)
  )

# Preview
plot_predict
plot_predict_ctr

# save once as PDF
ggsave(here::here("result","fig_4.png"), plot_predict, width = 12, height = 16, device = grDevices::cairo_pdf)
ggsave(here::here("result","fig_4_ctr.pdf"), plot_predict_ctr, width = 12, height = 45, device = grDevices::cairo_pdf)

