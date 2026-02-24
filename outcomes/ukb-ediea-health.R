# UKBiobank health outcomes HR

## Environment setting
ukb_path <- "../UKB/dataset/"

library(tidyverse)
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(readr)
library(jsonlite)
library(mice)
library(survival)
library(ggplot2)
library(ggjournals)

load('./dat/ukb_data.rda')
disease_fieldids <- fromJSON("./UKB/disease_fieldids.json", simplifyVector = TRUE)


cox_covars <- c("age_rec", "sex", "edu_college", "eth_white", "bmi",
                "tdi", "smoking", "drinking", "actc")
cox_formula <- as.formula(
  paste0(
    "Surv(fu, event) ~ ediea",
    if (length(cox_covars) > 0) paste(" +", paste(cox_covars, collapse = " + ")) else ""
  )
)

# ------------------------ Outcome HR ------------------------
hr_list <- map_dfr(names(disease_fieldids) %>% unique(), function(g) {
  if(g!='Cancer_Melanoma') {
    cox_covars <- c("age_rec", "sex", "edu_college", "eth_white",
                    "tdi", "smoking", "drinking", "actc")
    
  } else {
    cox_covars <- c("age_rec", "sex", "edu_college", "eth_white", "bmi",
                    "tdi", "smoking", "drinking", "actc")
  }
  
  cox_formula <- as.formula(
    paste0(
      "Surv(fu, event) ~ ediea",
      if (length(cox_covars) > 0) paste(" +", paste(cox_covars, collapse = " + ")) else ""
    )
  )

  date_col <- paste0("date_", g)
  inc_col  <- paste0("inc_",  g)
  
  if (!date_col %in% names(ukb_data) || !inc_col %in% names(ukb_data)) {
    message("  >> no columns for ", g, ", skip")
    return(NULL)
  }
  
  eof <- as.Date(max(ukb_data[[date_col]], na.rm = TRUE))
  print(g)
  df_g <- ukb_data %>%
    filter(!is.infinite(.data$last_diet_date)) %>%
    mutate(
      baseline_date = as.Date(.data$last_diet_date),
      event_date = as.Date(.data[[date_col]]),
      t_end = if_else(!is.na(event_date), event_date,
                              if_else(!is.na(date_death), as.Date(date_death),
                                      eof)),
      
      event = if_else(!is.na(event_date), 1L, 0L),
      fu = as.numeric(t_end - baseline_date) / 365.25
    ) %>%
    filter(fu > 0)
  
  fit <- coxph(cox_formula, data = df_g)
  
  q_ediea   <- quantile(df_g$ediea, probs = c(0.1, 0.9), na.rm = TRUE)
  delta_90_10 <- unname(q_ediea[2] - q_ediea[1])
  
  tidy(fit) %>%
    filter(term == "ediea") %>%
    transmute(
      disease_group = g,
      term,
      HR   = exp(estimate * delta_90_10),
      LCL  = exp((estimate - 1.96 * std.error)*delta_90_10),
      UCL  = exp((estimate + 1.96 * std.error)*delta_90_10),
      p    = p.value,
      n    = nrow(df_g),
      events = sum(df_g$event)
    )
})

eof_mort <- as.Date(max(ukb_data$date_death, na.rm = TRUE))
df_mort <- ukb_data %>%
  filter(!is.infinite(last_diet_date)) %>%
  mutate(
    baseline_date = as.Date(last_diet_date),
    t_end = if_else(
      !is.na(date_death), as.Date(date_death), eof_mort
    ),
    event = if_else(!is.na(date_death), 1L, 0L),
    fu    = as.numeric(t_end - baseline_date) / 365.25
  ) %>%
  filter(fu > 0)

# delta = P90 - P10 for EDISEA
q <- quantile(df_mort$ediea, probs = c(0.10, 0.90), na.rm = TRUE, type = 2)
delta <- as.numeric(q[[2]] - q[[1]])

fit_mort <- coxph(cox_formula, data = df_mort)
hr_mort <- tidy(fit_mort) %>%
  filter(term == "ediea") %>%
  transmute(
    disease_group = "All-cause Mortality",
    term,
    HR    = exp(estimate*delta),
    LCL   = exp((estimate - 1.96 * std.error)*delta),
    UCL   = exp((estimate + 1.96 * std.error)*delta),
    p     = p.value,
    n     = nrow(df_mort),
    events = sum(df_mort$event)
  )

df_mort_extern <- ukb_data %>%
  filter(!is.infinite(last_diet_date)) %>%
  mutate(
    baseline_date = as.Date(last_diet_date),
    t_end = if_else(
      !is.na(date_traffic_death), as.Date(date_traffic_death), eof_mort
    ),
    event = if_else(!is.na(date_traffic_death), 1L, 0L),
    fu    = as.numeric(t_end - baseline_date) / 365.25
  ) %>%
  filter(fu > 0)

fit_mort_extern <- coxph(cox_formula, data = df_mort_extern)
hr_mort_extern <- tidy(fit_mort_extern) %>%
  filter(term == "ediea") %>%
  transmute(
    disease_group = "External-cause Mortality",
    term,
    HR    = exp(estimate*delta),
    LCL   = exp((estimate - 1.96 * std.error)*delta),
    UCL   = exp((estimate + 1.96 * std.error)*delta),
    p     = p.value,
    n     = nrow(df_mort),
    events = sum(df_mort$event)
  )

hr_results <- bind_rows(hr_list, hr_mort) %>%
  arrange(disease_group) %>%
  mutate(
    p_fdr = p.adjust(p, method = 'BH')
  )
hr_results

prefix_to_cat <- c(
  "Infection"     = "Diseases of infections (A01-B09)",
  "Cancer"     = "Cancer (C00-C97)",
  "BloodImmune"= "Blood and immune disorders (D50-D89)",
  "Endocrine"       = "Endocrine disorders (E00-E90)",
  "Mental"     = "Mental and behavioural disorders (F00-F99)",
  "Neuro"      = "Nervous system disorders (G00-G99)",
  "Eye"        = "Eye disorders (H00-H59)",
  "Ear"        = "Ear disorders (H60-H69)",
  "Circ"       = "Circulatory system disorders (I00-I99)",
  "Resp"       = "Respiratory system disorders (J00-J99)",
  "Digestive"     = "Digestive system disorders (K00-K93)",
  "Skin"       = "Skin disorders (L00-L99)",
  "MSK"        = "Musculoskeletal system disorders (M00-M99)",
  "GU"     = "Genitourinary system disorders (N00-N99)",
  "All-cause Mortality" = "All-cause Mortality"
  # "External-cause Mortality" = "External-cause Mortality (negative control)"
)
plot_dat <- hr_results %>%
  mutate(
    prefix = str_extract(disease_group, "^[^_]+"),
    category = prefix_to_cat[prefix],
    
    suffix = str_remove(disease_group, "^[^_]+_?"),
    
    text_color = ifelse(grepl("_any", disease_group, ignore.case = TRUE), "#71B5E6", "black"),
    outcome_label = if_else(str_detect(disease_group, "_any"), category, suffix),
    outcome_label = if_else(disease_group%in%c('All-cause Mortality', 'External-cause Mortality'), disease_group, outcome_label)
  ) %>%
  drop_na(category) %>%
  arrange(category, outcome_label) %>%
  mutate(id = row_number())

# ---------------------- FIGURE -------------------------
plot_dat_plot <- plot_dat %>%
  mutate(
    is_any = str_detect(disease_group, "_any")|str_detect(disease_group, "Mortality"),
    category = factor(
      category,
      levels = c(
        "All-cause Mortality",
        # "External-cause Mortality (negative control)",
        "Diseases of infections (A01-B09)",
        "Cancer (C00-C97)",
        "Endocrine disorders (E00-E90)",
        "Circulatory system disorders (I00-I99)",
        "Respiratory system disorders (J00-J99)",
        "Digestive system disorders (K00-K93)",
        "Genitourinary system disorders (N00-N99)",
        "Musculoskeletal system disorders (M00-M99)",
        "Nervous system disorders (G00-G99)",
        "Mental and behavioural disorders (F00-F99)",
        "Eye disorders (H00-H59)",
        "Ear disorders (H60-H69)",
        "Skin disorders (L00-L99)",
        "Blood and immune disorders (D50-D89)"
      )
    )
  ) %>%
  drop_na(category) %>%
  group_by(category) %>%
  arrange(desc(is_any), outcome_label, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(
    outcome_label = factor(outcome_label, levels = rev(unique(outcome_label))),
    sig = p_fdr < 0.05
  )

cat_bg <- plot_dat_plot %>%
  group_by(category) %>%
  mutate(
    y_pos = match(outcome_label, levels(droplevels(outcome_label)))
  ) %>%
  ungroup() %>%
  filter(is_any) %>%
  mutate(ymin = y_pos - 0.5, ymax = y_pos + 0.5)

x_min <- min(plot_dat_plot$LCL, na.rm = TRUE)
x_max <- 1.5
x_lim <- c(floor((x_min - 0.05) * 20) / 20, ceiling((x_max + 0.05) * 20) / 20)

p_ukb_hr <- ggplot(
  plot_dat_plot,
  aes(x = HR, y = outcome_label)
) +
  geom_rect(
    data = cat_bg,
    aes(xmin = -Inf, xmax = x_max, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey95",
    colour = NA
  ) +
  facet_grid(
    category ~ .,
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  ) +
  geom_vline(
    xintercept = 1,
    linetype   = "dashed",
    linewidth  = 0.4,
    colour     = "grey50"
  ) +
  geom_errorbarh(
    aes(xmin = LCL, xmax = if_else(UCL<1.5, UCL, 1.55)),
    height    = 0,
    linewidth = 0.5,
    colour    = "grey20"
  ) +
  geom_point(
    aes(fill = sig, colour = sig),
    shape  = 21,
    size   = 2.4,
    stroke = 0.4
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 21, `TRUE` = 21)
  ) +
  scale_fill_manual(
    values = c(
      `FALSE` = "white",
      `TRUE`  = "#71B5E6"
    )
  ) +
  geom_text(
    aes(
      x     = x_lim[2] - (x_lim[2] - x_lim[1]) * 0.02,
      label = sprintf("%.2f (%.2f-%.2f)", HR, LCL, UCL)
    ),
    hjust = 0,
    size  = 2.7
  ) +
  scale_x_continuous(
    limits = x_lim,
    breaks = pretty(x_lim, n = 4),
    expand = expansion(mult = c(0.01, 0.15))
  ) +
  
  labs(
    x = "HR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  theme(
    text        = element_text(size = 18),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey90"),

    strip.placement    = "outside",
    strip.text.y.left  = element_text(
      angle = 0,
      face  = "bold",
      size  = 9,
      hjust = 0
    ),
    strip.background   = element_rect(
      fill   = "grey95",
      colour = NA
    ),
    
    axis.text.y  = element_text(size = 9),
    axis.text.x  = element_text(size = 8),
    axis.title.x = element_text(size = 9, face = "bold", margin = margin(t = 6)),

    plot.title   = element_text(size = 11, face = "bold", hjust = 0),

    legend.position = "none",

    panel.spacing.y = unit(0.4, "lines"),
    panel.spacing.x = unit(0.8, "lines"),

  )

# ---FIGURE
library(patchwork)
library(cowplot)
library(grid)

x_min <- 0.6
x_max <- 1.2
x_text <- x_max + 0.02

plot_dat_plot2 <- plot_dat_plot %>%
  mutate(
    HR_clip  = pmin(pmax(HR,  x_min), x_max),
    LCL_clip = pmax(LCL, x_min),
    UCL_clip = pmin(UCL, x_max),
    left_trunc  = LCL < x_min,
    right_trunc = UCL > x_max
  )

dx <- 0.03 * (x_max - x_min)

left_arw  <- plot_dat_plot2 %>%
  filter(left_trunc) %>%
  transmute(category, outcome_label,
            x = x_min + dx, xend = x_min)

right_arw <- plot_dat_plot2 %>%
  filter(right_trunc) %>%
  transmute(category, outcome_label,
            x = x_max - dx, xend = x_max)

cat_bg2 <- plot_dat_plot2 %>%
  group_by(category) %>%
  mutate(y_pos = match(outcome_label, levels(droplevels(outcome_label)))) %>%
  ungroup() %>%
  filter(is_any) %>%
  mutate(ymin = y_pos - 0.5, ymax = y_pos + 0.5)

p_ukb_hr <- ggplot(
  plot_dat_plot2,
  aes(x = HR_clip, y = outcome_label)
) +
  geom_rect(
    data = cat_bg2,
    aes(xmin = x_min, xmax = x_max, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey95",
    colour = NA
  ) +
  facet_grid(
    category ~ .,
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  ) +
  geom_vline(
    xintercept = 1,
    linetype   = "dashed",
    linewidth  = 0.4,
    colour     = "grey50"
  ) +
  geom_errorbarh(
    aes(xmin = LCL_clip, xmax = UCL_clip),
    height    = 0,
    linewidth = 0.5,
    colour    = "grey20"
  ) +
  geom_segment(
    data = left_arw,
    aes(x = x, xend = xend, y = outcome_label, yend = outcome_label),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "grey20",
    arrow = arrow(type = "closed", length = unit(0.10, "in"))
  ) +
  geom_segment(
    data = right_arw,
    aes(x = x, xend = xend, y = outcome_label, yend = outcome_label),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "grey20",
    arrow = arrow(type = "closed", length = unit(0.10, "in"))
  ) +
  geom_point(
    aes(fill = sig, colour = sig),
    shape  = 21,
    size   = 2.4,
    stroke = 0.4
  ) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE`  = "#71B5E6")) +
  scale_colour_manual(values = c(`FALSE` = "grey20", `TRUE` = "#71B5E6")) +
  geom_text(
    aes(
      x     = x_text,
      label = sprintf("%.2f (%.2f-%.2f)", HR, LCL, UCL)
    ),
    hjust = 0,
    size  = 3
  ) +
  scale_x_continuous(
    breaks = seq(x_min, x_max, by = 0.2),
    expand = expansion(mult = c(0.01, 0.25))
  ) +
  coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
  labs(
    x = "HR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    text        = element_text(size = 18),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey90"),
    
    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 0, face = "bold", size = 9, hjust = 0),
    strip.background   = element_rect(fill = "grey95", colour = NA),
    
    axis.text.y  = element_text(size = 9),
    axis.text.x  = element_text(size = 8),
    axis.title.x = element_text(size = 9, face = "bold", margin = margin(t = 6)),
    
    legend.position = "none",
    panel.spacing.y = unit(0.4, "lines"),
    panel.spacing.x = unit(0.8, "lines"),

  )

# ---------（Chapter / Outcome / Hazard Ratio）---------
p_hdr <- cowplot::ggdraw() +
  cowplot::draw_label("Chapter", x = -0.89, y = 0.5, hjust = 0,
                      fontface = "bold", size = 10) +
  cowplot::draw_label("Outcome", x = -0.25, y = 0.5, hjust = 0,
                      fontface = "bold", size = 10) +
  cowplot::draw_label("Hazard Ratio", x = 0.83, y = 0.5, hjust = 0,
                      fontface = "bold", size = 10)

p_final <- p_hdr / p_ukb_hr + plot_layout(heights = c(0.05, 0.95))

ggsave(filename = "fig/figure_4.pdf", plot = p_final,
       width = 12, height = 12)

## Abstract
dat_any <- plot_dat_plot %>%
  filter(is_any) %>%
  mutate(
    disease = case_when(
      category == "All-cause Mortality" ~ "All-cause\nMortality",
      TRUE ~ str_trim(str_remove(category, "\\s*\\([^\\)]+\\)"))
    ),
    icd10 = case_when(
      category == "All-cause Mortality" ~ "Mortality",
      TRUE ~ str_match(category, "\\(([^\\)]+)\\)")[,2]
    )
  ) %>%
  arrange(HR) %>%
  mutate(
    disease = factor(disease, levels = disease),
    fill_rank = scales::rescale(row_number())
  )

dat_any <- dat_any %>%
  filter(disease_group!='External-cause Mortality') %>%
  mutate(icd10 = factor(icd10, levels = unique(icd10)))

pal15 <- c(
  "#0B3B3B", "#114444", "#164D4D", "#1C5656", "#265E5E",
  "#336567", "#416C70", "#4E7379", "#64828A", "#79909B",
  "#8F9FAC", "#A3AEBC", "#B4BDCA", "#C6CDD8", "#D7DCE6"
)

p_any_bar_v <- ggplot(dat_any, aes(x = disease, y = HR, fill = icd10)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  geom_col(width = 0.72, alpha = 0.95) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL),
                width = 0.18, linewidth = 0.45, colour = "grey25") +
  scale_y_log10(breaks = c(0.6,0.7,0.8,0.9,1),
                labels = c("0.6","0.7","0.8","0.9","1.0")) +
  scale_fill_manual(
    values = setNames(pal15, levels(dat_any$icd10)),
    name = "ICD-10"
  ) +
  labs(x = NULL, y = "HR") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.y = element_text(size = 9, face = "bold", margin = margin(r = 4)),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text  = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    legend.position = "right",
    plot.margin = margin(3, 3, 3, 3),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

p_any_bar_v