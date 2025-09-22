# ============================================================
# Title: Fig.1 : Gender differences in academic career trajectories.
# Project: Leaky Pipeline in Psychology
# Author: Xinyi Zhao  |  Max Planck Institute for Human Development
# Date: 2025-09-22
# R version: >= 4.3
# Description: 
# Inputs: <list CSV/RDS/objects used>
# Outputs: <file(s) saved, e.g., figure_a.pdf>
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

pkgs <- c(
  # data handling
  "dplyr", "tidyr", "tibble", "readr", "stringr", "janitor","here",
  # plotting
  "ggplot2", "ggrepel", "scales", "patchwork", "viridisLite",
  # survival (for KM)
  "survival", "survminer"
)
install_and_load(pkgs)

# ---------- Global options (safe defaults for all your figures) ----------
set.seed(42)
options(stringsAsFactors = FALSE, scipen = 999)
theme_set(theme_minimal(base_size = 14))



# ---------- Data Import ----------

#df_plot1<-read.csv('//mpib-berlin.mpg.de/FB-ARC/Users/ZhaoXinyi/1_Project/2_Bibliometric_Data/2_Scopus_Data/3_result_forVIZ/1_female_field_cohort.csv')
#df_KM<-read.csv('//mpib-berlin.mpg.de/FB-ARC/Users/ZhaoXinyi/1_Project/2_Bibliometric_Data/2_Scopus_Data/2_analysis//2_earlycareer_character_0503.csv')
#df_plot1_3_yearly_cg<-read.csv('//mpib-berlin.mpg.de/FB-ARC/Users/ZhaoXinyi/1_Project/2_Bibliometric_Data/2_Scopus_Data/3_result_forVIZ/4_stage_trans_gap_yearly_censored_0529.csv')

df_trend<-read.csv(here::here("data", "1_female_field_cohort.csv")) #data for plot(a)
df_KM<-read.csv(here::here("data", "2_earlycareer_character_0503.csv"))
df_career<-read.csv(here::here("data", "4_stage_trans_gap_yearly_censored_0529.csv"))


# ---------- Analysis and result visualization ----------
# ---------- subplot (a)  Proportion of women among psychology entrants, by subfield and cohort----------


# ---- Quick data checks (fail early if something is missing) ----
req <- c("cohort", "subfield", "female_proportion")
stopifnot(all(req %in% names(df_trend)))

# ---- Tidy types / ordering ----
df_trend <- df_trend %>%
  mutate(
    # keep cohorts in numeric order on a discrete axis
    cohort = as.character(cohort),
    cohort = factor(cohort, levels = sort(unique(cohort), na.last = NA)),
    subfield = as.character(subfield)
  )

# ---- Subfield display names (avoid relying on 'unique()' order) ----
# If your raw subfield names are already what you want to show, you can skip this.
subfield_key <- tibble::tribble(
  ~subfield,          ~subfield_abbr,
  "Clinical Psy",     "Clinical Psy",
  "Social Psy",       "Social Psy",
  "Neu & Phy Psy",    "Neu & Phy Psy",
  "Applied Psy",      "Applied Psy",
  "General Psy",      "General Psy",
  "Exp & Cog Psy",    "Exp & Cog Psy",
  "DEV & Edu Psy",    "DEV & Edu Psy",
  "Miscellaneous",    "Miscellaneous",
  "Overall",          "Overall"
)

# join (safe even if some keys are missing—labels fall back to subfield)
df_trend <- df_trend %>%
  left_join(subfield_key, by = "subfield") %>%
  mutate(subfield_abbr = ifelse(is.na(subfield_abbr), subfield, subfield_abbr))

# ---- Colors: one per subfield, "Overall" forced to grey40 ----
all_subfields <- sort(unique(df_trend$subfield))
subfield_colors_1 <- setNames(
  viridisLite::viridis(length(all_subfields), begin = 0.2, end = 0.98),
  all_subfields
)
subfield_colors_1["Overall"] <- "grey40"

# ---- Label positions: last cohort per subfield ----
df_labels <- df_trend %>%
  group_by(subfield) %>%
  filter(cohort == max(cohort)) %>%
  ungroup()

# ---- Plot ----
p_0left <- ggplot(
  df_trend,
  aes(x = cohort, y = female_proportion, group = subfield, color = subfield)
) +
  # all subfields except Overall
  geom_line(
    data = dplyr::filter(df_trend, subfield != "Overall"),
    linewidth = 1.2, alpha = 0.9, show.legend = FALSE
  ) +
  # "Overall" highlighted
  geom_line(
    data = dplyr::filter(df_trend, subfield == "Overall"),
    color = "grey40", linewidth = 2.5, show.legend = FALSE
  ) +
  scale_color_manual(values = subfield_colors_1) +
  scale_y_continuous(limits = c(0.4, 0.8)) +
  scale_x_discrete(
    expand = c(0.05, 0.40),                            # extra right space for labels
    breaks = as.character(seq(2000, 2014, by = 2))
  ) +
  labs(
    x = "Cohort", y = "Proportion of female researchers", color = "Subfield",
    title = ""
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40", linewidth = 1.2) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(5.5, 80, 5.5, 5.5),  # room for right-hand labels
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80"),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_blank(),
    axis.line          = element_line(color = "black"),
    legend.position    = "none",
    plot.title         = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  # right-side end labels (spaced and not touching lines)
  +geom_text_repel(data = df_labels,
                  aes(label = subfield_abbr), color = "black",
                  hjust =0, nudge_x =1.8, size = 5, 
                  show.legend = FALSE, direction = "y",
                  force = 15,                  # push labels apart
                  segment.color = NA,          # no connecting lines
                  box.padding = grid::unit(0.2, "lines"),
                  point.padding = grid::unit(0.04, "lines"),
                  seed = 42
  ) 

p_0left




# ---------- subplot (b)   Kaplan--Meier survival probability of psychology entrants by gender.----------

lab_time  <- 8       # x-position where we place the "Male"/"Female" labels
palette_g <- c(female = "#ea801c", male = "#3594cc")

# ---- Data prep / sanity checks ----
req_cols <- c("career_length", "leaving", "gender", "cohort_1")
stopifnot(all(req_cols %in% names(df_KM)))

df_KM <- df_KM %>%
  mutate(
    gender   = factor(gender, levels = c("female", "male")),
    cohort_1 = factor(cohort_1, levels = c("2000", "2005", "2010")),
    leaving  = as.integer(leaving)  # ensure 0/1
  )

# ---- Fit models ----
km_fit         <- survfit(Surv(career_length, leaving == 1) ~ gender, data = df_KM)
km_fit_cohort  <- survfit(Surv(career_length, leaving == 1) ~ gender + cohort_1, data = df_KM)

# ---- Medians (overall by gender, and by gender×cohort if needed) ----
median_times <- survminer::surv_median(km_fit) %>%
  mutate(strata = sub("^gender=", "", strata))  # "female"/"male"

median_times_cohort <- survminer::surv_median(km_fit_cohort) %>%
  mutate(strata = sub("^gender=", "", strata))  # "female:cohort_1=..."

# ---- Base KM plot (no legend; p-value off; optional CI) ----
km_plot <- ggsurvplot(
  km_fit,
  data        = df_KM,
  conf.int    = TRUE,
  pval        = FALSE,
  risk.table  = FALSE,
  legend.title= "Gender",
  legend.labs = c("female", "male"),
  palette     = unname(palette_g[c("female", "male")]),
  xlab        = "Academic age",
  ylab        = "Survival probability",
  ylim        = c(0.3, 1.02),
  xlim        = c(1, 19),
  censor.size = 1.2
)

# ---- Axes & theme ----
km_plot$plot <- km_plot$plot +
  scale_x_continuous(breaks = c(2, 5, 10, 15), labels = c("2", "5", "10", "15")) +
  scale_y_continuous(limits = c(0.3, 1.02), expand = c(0, 0), breaks = seq(0.4, 1, by = 0.2)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 1.2) +
  theme(
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80"),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_blank(),
    axis.line          = element_line(color = "black"),  # left + bottom frame only
    panel.border       = element_blank(),
    legend.position    = "none",
    plot.margin        = margin(5.5, 30, 5.5, 5.5)
  )

# ---- Median markers: horizontal at 0.5 and vertical at each group's median ----
km_plot$plot <- km_plot$plot +
  geom_segment(
    data = median_times,
    aes(x = 0, y = 0.5, xend = median, yend = 0.5),
    inherit.aes = FALSE,
    linetype = "dashed", linewidth = 1, color = "grey40"
  ) +
  geom_segment(
    data = median_times,
    aes(x = median, y = 0.3, xend = median, yend = 0.5),
    inherit.aes = FALSE,
    linetype = "dashed", linewidth = 1, color = "grey40"
  )

# ---- Inline labels ("Male", "Female") with gentle separation ----
# Use extend=TRUE to guarantee values even if the exact time isn't in the step function.
summ <- summary(km_fit, times = lab_time, extend = TRUE)
label_df <- tibble(
  group_raw = as.character(summ$strata),          # e.g., "gender=female"
  group     = sub(".*=", "", group_raw),          # -> "female"/"male"
  x         = lab_time,
  y0        = summ$surv
) %>%
  mutate(
    label = ifelse(group == "female", "Female", "Male"),
    x     = ifelse(group == "female", x - 1, x + 1),         # small horizontal separation
    y     = ifelse(group == "female", y0 - 0.02, y0 + 0.02)  # small vertical separation
  )

km_plot$plot <- km_plot$plot +
  coord_cartesian(clip = "off") +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 7, hjust = -0.1
  )

# ---- Draw ----
print(km_plot)






# ---------- subplot (c) Annual transition rates by academic age and cohort group (2000–2004, 2005–2009, 2010–2014).----------

# ---- Checks ----
need <- c("next_age", "trans_female", "trans_male", "cohort_group")
stopifnot(all(need %in% names(df_career)))

# ---- Tidy to long, keep valid rows ----
df_long_yearly_cg <- df_career %>%
  pivot_longer(
    cols = c(trans_female, trans_male),
    names_to = "gender",
    values_to = "trans"
  ) %>%
  mutate(
    gender = ifelse(gender == "trans_female", "female", "male"),
    cohort_group = as.character(cohort_group)     # ensure text
  ) %>%
  filter(trans > 0.5)                              # your filter

# Keep only usable ages and create the grouping key used for color mapping
df_plot <- df_long_yearly_cg %>%
  filter(next_age >= 2) %>%
  mutate(
    # normalize cohort labels if needed; here we assume "2000","2005","2010"
    cohort_group = factor(cohort_group, levels = c("2000","2005","2010")),
    gender = factor(gender, levels = c("female","male")),
    grp = interaction(gender, cohort_group, drop = TRUE)
  )

# ---- Explicit names for end labels (legend-free labeling) ----
lab_names <- c(
  "female.2000" = "Female 2000–2004",
  "female.2005" = "Female 2005–2009",
  "female.2010" = "Female 2010–2014",
  "male.2000"   = "Male 2000–2004",
  "male.2005"   = "Male 2005–2009",
  "male.2010"   = "Male 2010–2014"
)

# ---- Colors (named to match grp levels) ----
cols_grp <- c(
  "female.2000" = "#a33d00", "female.2005" = "#ea801c", "female.2010" = "#ffc773",
  "male.2000"   = "#003366", "male.2005"   = "#3594cc", "male.2010"   = "#7bb8d9"
)

# ---- Plot ----
p1 <- ggplot(
  df_plot,
  aes(x = next_age, y = trans, color = grp)
) +
  # Highlight a window (4–6)
  annotate("rect", xmin = 4, xmax = 6, ymin = -Inf, ymax = Inf,
           alpha = 0.5, fill = "gray") +
  # Smooth lines (match loess settings with label computation below)
  geom_smooth(
    method = "loess", se = FALSE, linetype = "solid", linewidth = 1.5,
    method.args = list(span = 0.75, degree = 2)
  ) +
  # Axes
  scale_x_continuous(limits = c(2, 21), breaks = c(2, 5, 10, 15, 19)) +
  scale_y_continuous(limits = c(0.9, 1), breaks = seq(0.9, 1, by = 0.05)) +
  labs(x = "Academic age", y = "Translation rate", color = "Gender & Cohort") +
  scale_color_manual(values = cols_grp) +
  # Theme
  theme(
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_blank(),
    axis.line          = element_line(color = "black"),
    panel.border       = element_blank(),
    legend.position    = "none",                     # <- lowercase "none"
    legend.justification = c("right", "top"),
    legend.title       = element_text(size = 12, face = "bold"),
    legend.text        = element_text(size = 11),
    legend.key.height  = unit(1.5, "cm")
  )

# ---- Right-side labels: last x per group (+0.25), y from loess at last x ----
label_df <- df_plot %>%
  group_by(grp) %>%
  group_modify(~{
    x_last <- max(.x$next_age, na.rm = TRUE)
    fit    <- loess(trans ~ next_age, data = .x, span = 0.75, degree = 2)
    
    # Predict at x_last; if NA (edge cases), fall back to observed value at x_last
    y_pred <- as.numeric(predict(fit, newdata = data.frame(next_age = x_last)))
    y_obs  <- .x %>% filter(next_age == x_last) %>% summarise(y = mean(trans, na.rm = TRUE)) %>% pull(y)
    y_use  <- dplyr::coalesce(y_pred, y_obs)
    
    tibble(
      x     = x_last + 0.25,                          # your rule (kept inside axis limits)
      y     = y_use,
      label = lab_names[as.character(.y$grp)]
    )
  }) %>%
  ungroup()

# In case any label would exceed x-limit 21, keep it barely inside so it shows:
label_df <- label_df %>% mutate(x = pmin(x, 20.95))

p1 <- p1 +
  coord_cartesian(clip = "off") +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, nudge_x = 0,
    direction = "y", force = 5,
    box.padding = grid::unit(0.22, "lines"),
    point.padding = grid::unit(0.05, "lines"),
    segment.color = NA,
    size = 6,
    seed = 42
  ) +
  theme(plot.margin = margin(5.5, 70, 5.5, 5.5)) +
  guides(colour = "none")

p1






# ---------- subplot (d) Summary of gender gaps in transition rates across three time frames----------
# ---------- 1) Full captured career span (which varies by cohort)----------
# ---------- 2) Shared career window common to all cohort groups (academic age 2–9) ----------
# ---------- 3) Early-career period (academic age 4–6)  ----------

# ---- checks ----
need <- c("cohort_group","cohort","next_age","gender","trans")
stopifnot(all(need %in% names(df_long_yearly_cg)))

# ---- Ensure consistent types & cohort labels ----
df_long_yearly_cg <- df_long_yearly_cg %>%
  mutate(
    # cohort_group may be numeric or character; normalize to "2000","2005","2010"
    cohort_group = as.character(cohort_group),
    cohort_group = dplyr::recode(cohort_group,
                                 "2000"="2000–2004", "2005"="2005–2009", "2010"="2010–2014",
                                 .default = cohort_group),
    gender = factor(gender, levels = c("female","male"))
  )

# ---- Compute mean gap (male - female) + 95% CI for a filter/window ----
summarise_gap <- function(data, filter_expr, label_text) {
  data %>%
    filter({{ filter_expr }}) %>%
    select(cohort_group, next_age, gender, trans) %>%
    pivot_wider(names_from = gender, values_from = trans) %>%
    mutate(gender_gap = male - female) %>%
    group_by(cohort_group) %>%
    summarise(
      mean_gap = mean(gender_gap, na.rm = TRUE),
      se       = sd(gender_gap,   na.rm = TRUE) / sqrt(sum(!is.na(gender_gap))),
      lower    = mean_gap - 1.96 * se,
      upper    = mean_gap + 1.96 * se,
      .groups  = "drop"
    ) %>%
    mutate(label = label_text)
}

# ---- Build three panels: Overall career, Age 2–9, Age 4–6 ----
trans_overall <- summarise_gap(df_long_yearly_cg, TRUE,                 "Overall career")
trans_2_9     <- summarise_gap(df_long_yearly_cg, next_age <= 9,        "Age 2–9")
trans_4_6     <- summarise_gap(df_long_yearly_cg, next_age >= 4 & next_age <= 6, "Age 4–6")

trans_overall_origin <- bind_rows(trans_overall, trans_2_9, trans_4_6) %>%
  mutate(
    # Facet order on x-axis
    label = factor(label, levels = c("Overall career","Age 2–9","Age 4–6")),
    # Cohort legend order
    cohort_group = factor(cohort_group, levels = c("2000–2004","2005–2009","2010–2014"))
  ) %>%
  arrange(label)

# ---- Colors for cohorts ----
fill_cols <- c(
  "2000–2004" = "#004C6DFF",
  "2005–2009" = "#2E7F77FF",
  "2010–2014" = "#97D4BFFF"
)

# ---- Plot ----
p2 <- ggplot(trans_overall_origin,
             aes(x = label, y = mean_gap, fill = cohort_group)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.6) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.6),
    width = 0.2, linewidth = 1
  ) +
  labs(
    x = "", y = "Transition rate gap (male – female)",
    fill = NULL, title = ""
  ) +
  scale_fill_manual(values = fill_cols) +
  scale_y_continuous(limits = c(0, 0.03), breaks = seq(0, 0.03, by = 0.01)) +
  theme(
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_blank(),
    axis.line          = element_line(color = "black"),
    panel.border       = element_blank(),
    # Legend inside, clean background like your example
    legend.position    = c(0.82, 0.82),
    legend.justification = c(1, 1),
    legend.direction   = "vertical",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 17),
    legend.key.size    = unit(0.6, "cm"),
    legend.key         = element_rect(fill = "white", color = NA),
    legend.background  = element_rect(fill = scales::alpha("white", 0.6), color = NA)
  ) +
  guides(fill = guide_legend(override.aes = list(color = NA)))  # solid swatches

p2




# ---------- combine subplots for output----------

combined_plot_up <- p_0left+ plot_spacer() + km_plot$plot + 
  plot_layout(ncol = 3, widths = c(1,0.08,1)) +
  plot_annotation(
    tag_levels = "a",          # -> (a), (b), (c)...
    tag_prefix = "(", 
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.02, 0.98)  # top-left inside each panel
  )

combined_plot_down <-   (p1 + labs(tag = "(c)")) +
  (p2 + labs(tag = "(d)")) +
  plot_layout(widths = c(2, 1)) &
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

combined_plot <- combined_plot_up / combined_plot_down + plot_layout(heights = c(1, 1))  # 75% and 25%
combined_plot



ggsave(here::here("plots","fig_1.pdf"), combined_plot,
       width = 25, height = 14, dpi = 300,limitsize = FALSE)  # adjust size as needed





