# ============================================================
# Title: Fig.1 : Statistical description of factors associated with 
# academic performance, collaboration, and institutional affiliations 
# Project: Leaky Pipeline in Psychology
# Author: Xinyi Zhao  |  Max Planck Institute for Human Development
# Date: 2025-09-22
# R version: >= 4.3
# Description: 
# Inputs: survival_data.csv (time_varying variables for each researcher)
# Outputs: fig_2.pdf
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


# ---------- Global options (safe defaults for all your figures) ----------
set.seed(42)
options(stringsAsFactors = FALSE, scipen = 999)
theme_set(theme_minimal(base_size = 14))





# ----   Read & basic cleanup ----
csv_path <- here::here("data", "2_earlycareer_character_0503.csv")
survival_data <- readr::read_csv(csv_path, show_col_types = FALSE) %>%
  janitor::clean_names()  # snake_case, ASCII-safe names



# ============================================================
# Data preprocessing 
# ============================================================


# ---- 1) Sanity checks ----
req_cols <- c("scopus_author_id","career_length","survive","cohort","academic_age")
stopifnot(all(req_cols %in% names(survival_data)))

# ---- 2) Event indicator & hard NA drop (only on key cols) ----
survival_data <- survival_data %>%
  mutate(
    leave = 1 - survive    # event indicator (1 = left academia), consistent with your code
  )

# keep rows with all key fields present (avoid dropping for ancillary NAs)
survival_data <- survival_data %>%
  tidyr::drop_na(any_of(req_cols), leave)

# ---- 3) Right-censor “still active” authors at end of window
# If last record of an author has leave==0 AND academic_age + cohort - 1 > 2018,
# set leave/survive to NA for that last row (treated as censored).
# -------------------------------------------------------------
survival_data_modified <- survival_data %>%
  group_by(scopus_author_id) %>%
  mutate(
    is_last = row_number() == n(),
    leave   = ifelse(is_last & leave == 0 & (academic_age + cohort - 1) > 2018, NA_real_, leave),
    survive = ifelse(!is.na(leave), 1 - leave, NA_real_)
  ) %>%
  ungroup() %>%
  select(-is_last)

# ---- 4) Reference levels before dummy-coding
# Set factor bases so dummy_cols(remove_first_dummy=TRUE) drops the desired references.
# Adjust levels only if the columns exist.
relevel_if_has <- function(df, col, lvls) {
  if (col %in% names(df)) df %>% mutate(!!col := factor(.data[[col]], levels = lvls)) else df
}

survival_data_modified <- survival_data_modified %>%
  relevel_if_has("gender",            c("male","female")) %>%
  relevel_if_has("country_us",        c("non_USA","USA")) %>%
  relevel_if_has("continent_new",     c("other","Europe","North America")) %>%
  relevel_if_has("early_field_1",     c("Psychology_miscellaneous","Applied_Psychology",
                                        "Clinical_Psychology","Developmental_and_Educational_Psychology",
                                        "Experimental_and_Cognitive_Psychology",
                                        "Neuropsychology_and_Physiological_Psychology",
                                        "Social_Psychology","Psychology_all")) %>%
  relevel_if_has("cohort_1",          c("2000","2005","2010")) %>%
  relevel_if_has("migration",         levels(survival_data_modified$migration)) %>%         # keep as-is if unknown
  relevel_if_has("instit_change",     levels(survival_data_modified$instit_change))



# ---- 5) Dummy coding (drop reference; drop original columns)
dummy_vars <- c("gender","country_us","continent_new","early_field_1","cohort_1",
                "migration","instit_change")
have_vars  <- intersect(dummy_vars, names(survival_data_modified))

survival_data_modified <- fastDummies::dummy_cols(
  survival_data_modified,
  select_columns            = have_vars,
  remove_first_dummy        = TRUE,   # avoids perfect collinearity
  remove_selected_columns   = FALSE,   # drop original factors
  ignore_na                 = TRUE
)

# Column names are already snake_case from clean_names(); no need for gsub() fixes.




# ---- 6) Outlier filtering (99.5th percentile) + minimal tenure filter
outlier_vars <- c(
  "sum_fncr_5_years","mean_mncs_max_yearly","affili_total_mncs_wt",
  "firstauthor_num","corresauthor_num","singleauthor_num","lastauthor_num"
)
outlier_vars <- intersect(outlier_vars, names(survival_data_modified))

# Compute upper bounds safely
upper_bounds <- survival_data_modified %>%
  summarise(across(all_of(outlier_vars), ~ quantile(.x, 0.995, na.rm = TRUE)))

# Apply filters sequentially for readability
survival_data_modified_filtered <- survival_data_modified
for (v in outlier_vars) {
  thr <- upper_bounds[[v]][1]
  survival_data_modified_filtered <- survival_data_modified_filtered %>%
    filter(.data[[v]] <= thr | is.na(.data[[v]]))
}

survival_data_ready <- survival_data_modified_filtered %>%
  filter(career_length > 1)

# ---- 7) Variable -> label mapping (for tables/figures)
# NOTE: these names must match your *post-cleaning* column names.
variable_map <- tibble::tibble(
  variable = c(
    "corresauthor_num","affili_total_mncs_wt","sum_fncr_5_years","mean_fncr_5_years",
    "firstauthor_num","mean_mncs_max_yearly","mean_mncs_median_yearly","mean_mncs_mean_yearly",
    "nunique_coauthors_yearly",
    "country_us_usa","continent_new_usa","continent_new_europe","continent_new_north_america",
    "cohort_1_2010","cohort_1_2005",
    "early_field_1_applied_psychology","early_field_1_clinical_psychology",
    "early_field_1_developmental_and_educational_psychology",
    "early_field_1_experimental_and_cognitive_psychology",
    "early_field_1_neuropsychology_and_physiological_psychology",
    "early_field_1_social_psychology","early_field_1_psychology_all",
    "gender_female",
    "singleauthor_num","lastauthor_num","coauthor_female_prop","unique_affi"
  ),
  label = c(
    "No. corresponding-author papers","Affiliation scores","Citation scores","Average citation scores",
    "No. first-authored papers","Max coauthor impact","Median coauthor impact","Mean coauthors' impact",
    "No. coauthors",
    "Origin: USA","Origin Ctr: USA","Origin Ctr: Europe","Origin Ctr: North America",
    "Cohort: 2010","Cohort: 2005",
    "Applied Psy","Clinical Psy","Develop & Edu Psy","Exp & Cog Psy","Neuro & Physio Psy","Social Psy","General Psy",
    "Gender: female",
    "No. single-author papers","No. last-authored papers","Female-coauthor proportion","Unique affiliations"
  )
)

variable_map_groups <- tibble::tibble(
  variable = c("country_us","continent_new","cohort_1","early_field_1"),
  label    = c("Origin country","Origin region","Cohort","Subfield")
)

# ---- 8) Quick report
message("Rows (raw / filtered): ",
        nrow(survival_data), " / ", nrow(survival_data_ready))
message("Columns (final): ", ncol(survival_data_ready))







# ============================================================
# Pick categorical variables to summarize
# ============================================================
# You can put multiple here, e.g., c("continent_new", "cohort_1", "early_field_1")
if (exists("survival_data_modified_filtered")) {
  df_cat <- survival_data_modified_filtered
} else if (exists("survival_data_modified")) {
  df_cat <- survival_data_modified
} else {
  stop("No suitable data frame found. Expected survival_data_modified_filtered or survival_data_modified.")
}

cate_vars_to_summarize <- c("continent_new")

# ---- 1)sanity check
stopifnot(all(c("gender", cate_vars_to_summarize) %in% names(df_cat)))

# ---- 2) Long format → counts → within-gender proportions (by variable) ----
category_summary_gender <- df_cat %>%
  select(gender, all_of(cate_vars_to_summarize)) %>%
  pivot_longer(
    cols = -gender, names_to = "variable", values_to = "category"
  ) %>%
  group_by(variable, gender, category) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(proportion = count / sum(count)) %>%  # proportion within each (variable, gender)
  ungroup()

# ---- 3) Nice labels for variables + tidy category names ----
# If you defined variable_map_1 earlier, join it to get human-friendly labels.
# variable_map_1 should have columns: variable, label
if (exists("variable_map_1")) {
  category_summary_gender <- category_summary_gender %>%
    left_join(variable_map_1, by = "variable")
} else {
  category_summary_gender <- category_summary_gender %>%
    mutate(label = variable)
}

# Harmonize some category strings
category_summary_gender <- category_summary_gender %>%
  mutate(category = dplyr::case_when(
    tolower(category) == "non_usa"           ~ "Non-USA",
    category == "Psychology (all)"           ~ "General Psychology",
    tolower(category) == "other"             ~ "Other countries",
    TRUE                                     ~ as.character(category)
  ))

# Optional: set a meaningful order for continent categories if present
if ("continent_new" %in% cate_vars_to_summarize) {
  category_summary_gender <- category_summary_gender %>%
    mutate(category = factor(
      category,
      levels = c("North America", "Europe", "Other countries")
    ))
}

# Gender order
category_summary_gender <- category_summary_gender %>%
  mutate(gender = factor(gender, levels = c("female", "male")))

# ---- 4) Plot ----
p_cate <- ggplot(category_summary_gender,
                 aes(x = category, y = proportion, fill = gender)) +
  geom_col(position = position_dodge(width = 0.4), width = 0.4) +
  # If you later summarize multiple variables, uncomment to facet:
  # facet_wrap(~ label, scales = "free_x") +
  scale_fill_manual(values = c("female" = "#ea801c", "male" = "#3594cc")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Origin region of institute",
    x = "",
    y = "",
    fill = "Gender"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black"),
    strip.text   = element_text(size = 22),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80"),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_blank(),
    axis.line          = element_line(color = "black"),
    legend.position    = "none"  # was "None"; lower-case is required
  )

p_cate



# ============================================================
# Pick continous variables to summarize
# ============================================================
df_cont <- if (exists("survival_data_modified_filtered")) {
  survival_data_modified_filtered
} else if (exists("survival_data_ready")) {
  survival_data_ready
} else if (exists("survival_data_modified")) {
  survival_data_modified
} else stop("No suitable data frame found for continuous summaries.")

# ---- 1) Variables to summarize (use the names you expect) ----
conti_vars_to_summarize <- c(
  "sum_fncr_5_years","firstauthor_num","lastauthor_num","mean_mncs_mean_yearly",
  "nunique_coauthors_yearly","coauthor_female_prop","mean_fncr_5_years",
  "affili_total_mncs_wt","unique_affi"
)

# Handle case-mismatch (e.g., WT vs wt)
to_available <- function(vars, df) {
  v_df <- names(df)
  out  <- sapply(vars, function(v) {
    if (v %in% v_df) return(v)
    # try case-insensitive match
    hit <- v_df[tolower(v_df) == tolower(v)]
    if (length(hit)) return(hit[1])
    NA_character_
  }, USE.NAMES = FALSE)
  stats::na.omit(out)
}
conti_vars_present <- to_available(conti_vars_to_summarize, df_cont)
if (length(conti_vars_present) == 0) stop("None of the requested continuous variables were found.")

# ---- 2) Long form + labels ----
survival_data_age_conti_long <- df_cont %>%
  select(gender, scopus_author_id, academic_age, all_of(conti_vars_present)) %>%
  pivot_longer(cols = all_of(conti_vars_present),
               names_to = "variable", values_to = "value")

# Attach pretty labels if variable_map exists; otherwise keep variable names
if (exists("variable_map")) {
  survival_data_age_conti_long <- survival_data_age_conti_long %>%
    left_join(variable_map, by = "variable")
} else {
  survival_data_age_conti_long <- survival_data_age_conti_long %>%
    mutate(label = variable)
}

# ---- 3) Mean + 95% CI by (variable, gender, academic_age) ----
mean_ci_values <- survival_data_age_conti_long %>%
  group_by(variable, label, gender, academic_age) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    se         = sd(value,   na.rm = TRUE) / sqrt(sum(!is.na(value))),
    lower_ci   = mean_value - 1.96 * se,
    upper_ci   = mean_value + 1.96 * se,
    .groups    = "drop"
  )

# ---- 4) Facet order (only keep those that exist) ----
custom_order <- c(
  "Citation scores",
  "No. first-authored papers",
  "No. last-authored papers",
  "Mean coauthors' impact",
  "No. coauthors",
  "Female-coauthor proportion",
  "Affiliation scores",
  "Average citation scores",
  "Unique affiliations"
)
present_labels <- intersect(custom_order, unique(mean_ci_values$label))
mean_ci_values <- mean_ci_values %>%
  mutate(
    label  = factor(label, levels = present_labels),
    gender = factor(gender, levels = c("female","male"))
  )

# ---- 5) Plot: all selected continuous variables ----
p_conti <- ggplot(mean_ci_values,
                  aes(x = academic_age, y = mean_value,
                      group = gender, color = gender, fill = gender)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, color = NA) +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  scale_color_manual(values = c("male" = "#3594cc", "female" = "#ea801c")) +
  scale_fill_manual(values  = c("male" = "#3594cc", "female" = "#ea801c")) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, 5, 10, 15, 20)) +
  labs(x = "Academic age", y = "Mean value", title = "", color = "Gender", fill = "Gender") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black"),
    strip.text   = element_text(size = 24, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80"),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_blank(),
    axis.line          = element_line(color = "black"),
    legend.position    = "none"   # (lowercase)
  )

p_conti

# ============================================================
# patch the plots
# ============================================================
# --- 1) Helper: function that draws ONE panel for a given facet label

lbl<-custom_order
panel_plot <- function(lbl) {
  ggplot(filter(mean_ci_values, label == lbl),
         aes(x = academic_age, y = mean_value,
             group = gender, color = gender, fill = gender)) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, color = NA) +
    scale_color_manual(values = c("male" = "#3594cc", "female" = "#ea801c")) +
    scale_fill_manual(values  = c("male" = "#3594cc", "female" = "#ea801c")) +
    scale_x_continuous(expand = c(0, 0), breaks = c(1, 5, 10, 15, 20)) +
    labs(x = "", y = "") +
    ggtitle(lbl) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray80"),
      panel.grid.minor   = element_blank(),
      axis.text.x  = element_text(size = 18, color = "black"),
      axis.text.y  = element_text(size = 20, color = "black"),
      axis.title.x = element_text(size = 20, color = "black"),
      axis.title.y = element_text(size = 20, color = "black"),
      axis.line    = element_line(color = "black")
    )
}



# --- 2) Build the 3×3 with plot1 in the 8th position
# Positions (ncol = 3): 1 2 3 / 4 5 6 / 7 **8** 9

label_df <- tibble(
  label = c( "Female", "Male"),
  x     = c( 10+3, 10-3),        # small horizontal separation
  y     = c( 2.4, 3),  # small vertical separation
)

p1<-panel_plot(lbl[1])+
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 7
  )

plots_list <- c(
  list(p1),
  lapply(lbl[2:7], panel_plot), # slots 1..7 from your data
  list(p_cate),                            # slot 8 = your custom plot
  list(panel_plot(lbl[9]))       # slot 9 from your data
)

final_3x3 <- wrap_plots(plots_list, ncol = 3)
final_3x3


ggsave(here::here("plots","fig_2.pdf"), final_3x3,
       width = 18, height = 12, dpi = 300,limitsize = FALSE)  # adjust size as needed
