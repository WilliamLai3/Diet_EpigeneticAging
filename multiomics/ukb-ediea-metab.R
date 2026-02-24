# UKBiobank Proteomics

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
library(arrow)
library(ggrepel)
library(mediation)
library(STRINGdb)
library(igraph)
library(tidygraph)
library(ggraph)
library(lavaan)
library(glue)
library(ggforce)

parse_ukb_col <- function(x) {
  m <- str_match(x, "^p(\\d+)(?:_i(\\d+))?(?:_a(\\d+))?$")
  tibble(
    col = x,
    field_id = suppressWarnings(as.integer(m[, 2])),
    instance = suppressWarnings(as.integer(m[, 3])),
    array    = suppressWarnings(as.integer(m[, 4]))
  )
}

load("./dat/ukb_ediea_from_hrs.rda")
load("../UKB/data_gen/covars.rda")
load("../UKB/ukb_outcomes.rda") # Mort
nmr_inflam_data <- read_parquet("../UKB/dataset/participant_nmr_inflam.parquet")
fieldsum <- fread("../UKB/dataset/fieldsum.txt") %>%
  as_tibble() %>%
  mutate(
    field_id = as.integer(field_id),
    title = as.character(title),
    item_type = as.integer(item_type)
  )

ukb_data <- ukb_fg40 %>%
  left_join(data_cov, by = "eid") %>%
  left_join(disease_dates %>% 
              dplyr::select(eid, date_death, date_ad), 
            by = 'eid')
ukb_nmr_inflam <- ukb_data %>%
  inner_join(nmr_inflam_data %>%
               mutate(eid = as.numeric(eid)), by = "eid")

cox_covars <- c("age_rec", "sex", "edu_college", "eth_white", "bmi", 
                "tdi", "smoking", "drinking", "actc")
ukb_nmr_inflam[,cox_covars] <- ukb_nmr_inflam %>%
  dplyr::select(cox_covars) %>%
  mice(seed = 2025) %>%
  complete()

# ---------------- Multiple linear regression --------------
eof <- max(ukb_nmr_inflam$date_death, na.rm = TRUE)
ukb_nmr_inflam <- ukb_nmr_inflam %>%
  filter(!is.infinite(last_diet_date)) %>%
  mutate(
    baseline_date = as.Date(last_diet_date),
    event_date = as.Date(date_death),
    t_end = if_else(!is.na(event_date), event_date,
                    if_else(!is.na(date_death), as.Date(date_death),
                            as.Date(eof))),
    
    event = if_else(!is.na(event_date), 1L, 0L),
    fu = as.numeric(t_end - baseline_date) / 365.25
  ) %>%
  filter(fu > 0)

# ediea -> all-cause dementia
cox_formula_base <- as.formula(
  paste0(
    "Surv(fu, event) ~ ediea +",
    paste(cox_covars, collapse = " + ")
  )
)
fit_cox_base <- coxph(cox_formula_base, data = ukb_nmr_inflam)
beta_total <- coef(fit_cox_base)["ediea"]
tidy(fit_cox_base) %>%
  filter(str_detect(term, "ediea")) %>%
  transmute(
    term,
    HR   = exp(estimate),
    LCL  = exp(estimate - 1.96 * std.error),
    UCL  = exp(estimate + 1.96 * std.error),
    p    = p.value,
    n    = nrow(ukb_nmr_inflam),
    events = sum(ukb_nmr_inflam$event)
  )

# ============================================================
# 2) Define EDISEA, outcome for SEM, markers & metabolites
# ============================================================
category222 <- c(23652,23653,23654,23655,23651,20282,20283,23658,
                 23659,23649,23650,23660)
EDISEA_VAR <- "ediea"
inflam_markers <- c(
  "neutrophil", "lymphocyte", "monocyte", "platelet", "crp",
  "NLR", "PLR", "SII", "LMR"
)
inflam_markers <- inflam_markers[inflam_markers %in% names(ukb_nmr_inflam)]
# NMR metabolite columns
exclude_cols <- unique(c("eid", cox_covars, EDISEA_VAR,
                         inflam_markers, "date_death", "date_ad"))
num_cols <- names(ukb_nmr_inflam)[sapply(ukb_nmr_inflam, is.numeric)]
metab_candidates <- setdiff(num_cols, exclude_cols)
p_like <- metab_candidates[str_detect(metab_candidates, "^p\\d+(_i0)?(_a\\d+)?$")]
metab_cols <- parse_ukb_col(p_like) %>%
  left_join(fieldsum, by = "field_id") %>%
  filter(!str_detect(title, "QC Flag"), !str_detect(title, "count"), !str_detect(title, "C-reactive protein")) %>%
  filter(!field_id %in% category222)
metab_cols_qc <- parse_ukb_col(p_like) %>%
  left_join(fieldsum, by = "field_id") %>%
  filter(str_detect(title, "QC Flag"))

# ============================================================
# 3) Multiple linear regression
# ============================================================
# standardize
make_std_df <- function(dat, y, x, covars) {
  use_cols <- unique(c(y, x, covars))
  d <- dat %>% dplyr::select(all_of(use_cols))
  
  # drop NA rows for model variables
  d <- d %>% drop_na()
  
  # convert characters to factor
  for (nm in names(d)) {
    if (is.character(d[[nm]])) d[[nm]] <- as.factor(d[[nm]])
  }
  
  # standardize numeric columns among predictors/outcome
  # (so coefficient for x is standardized beta)
  std_cols <- c(y, x, covars)
  std_cols <- std_cols[std_cols %in% names(d)]
  for (nm in std_cols) {
    if (is.numeric(d[[nm]])) d[[nm]] <- as.numeric(scale(d[[nm]]))
  }
  
  d
}

fit_lm_std <- function(dat, y, x, covars) {
  d <- make_std_df(dat, y, x, covars)
  if (nrow(d) < 200) {
    return(tibble(y = y, n = nrow(d),
                  beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  fml <- as.formula(paste0(y, " ~ ", x, " + ", paste(covars, collapse = " + ")))
  fit <- lm(fml, data = d)
  
  sm <- summary(fit)$coefficients
  if (!x %in% rownames(sm)) {
    return(tibble(y = y, n = nrow(d),
                  beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  tibble(
    y = y, field_id = parse_ukb_col(y)$field_id[1], n = nrow(d),
    beta = unname(sm[x, "Estimate"]),
    se   = unname(sm[x, "Std. Error"]),
    p    = unname(sm[x, "Pr(>|t|)"])
  )
}

## ---- metabolites ------
message("Running LM for metabolites ...")
res_metab <- purrr::map_dfr(
  metab_cols$col,
  ~fit_lm_std(ukb_nmr_inflam, y = .x, x = EDISEA_VAR, covars = cox_covars)
) %>%
  mutate(fdr = p.adjust(p, method = "fdr"),
         fdr = pmax(fdr, 1e-300),
         sig = fdr<0.05,
         direction = if_else(beta<0, 'Inverse', 'Positive')) %>%
  arrange(fdr, p) %>% # 232
  left_join(fieldsum, by = 'field_id')

top10_metab_labels <- res_metab %>%
  filter(sig) %>%
  group_by(direction) %>%
  arrange(fdr, .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()

# FIGURE
plot_ediea_metab <- ggplot(
  res_metab,
  aes(x = beta, y = -log10(fdr))
) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(
    data = ~ filter(.x, !is.na(-log10(fdr))),
    color = "grey80",
    size  = 1.5
  ) +
  geom_point(
    data  = ~ filter(.x, sig),
    aes(colour = direction, alpha = sig),
    size = 2
  ) +
  theme_bmj() +
  scale_alpha_manual(
    values = c(`TRUE` = 1, `FALSE` = 0.3),
    guide = "none"
  ) +
  geom_text_repel(
    data = top10_metab_labels,
    aes(label = title),
    size          = 3,
    max.overlaps  = Inf,
    box.padding   = 0.4,
    point.padding = 0.2,
    segment.size  = 0.2,
    direction = "y",
    force = 20
  ) +
  scale_colour_manual(
    values = c(
      "Inverse" = "#2A6EBB",
      "Positive" = "#CD202C"
    ),
    name = "Association with EDISEA"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    x = "beta coefficient",
    y = "-log10(P_FDR value)",
    title = "a. Metabolites associated with EDISEA"
  )
plot_ediea_metab

# ---- inflammation markers ------
message("Running LM for inflammation markers ...")
res_inflam <- purrr::map_dfr(
  inflam_markers,
  ~fit_lm_std(ukb_nmr_inflam, y = .x, x = EDISEA_VAR, covars = cox_covars)
) %>%
  mutate(fdr = p.adjust(p, method = "fdr"), 
         sig = fdr<0.05,
         direction = if_else(beta<0, 'Inverse', 'Positive')
         ) %>%
  arrange(fdr, p)


# FIGURE
plot_ediea_inflam <- ggplot(
  res_inflam,
  aes(x = beta, y = -log10(fdr))
) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(
    data = ~ filter(.x, !is.na(-log10(fdr))),
    color = "grey80",
    size  = 1.5
  ) +
  geom_point(
    data  = ~ filter(.x, sig),
    aes(colour = direction, alpha = sig),
    size = 2
  ) +
  theme_bmj() +
  scale_alpha_manual(
    values = c(`TRUE` = 1, `FALSE` = 0.3),
    guide = "none"
  ) +
  geom_text_repel(
    data = ~ filter(.x, sig),
    aes(label = y),
    size          = 3,
    max.overlaps  = Inf,
    box.padding   = 0.4,
    point.padding = 0.2,
    segment.size  = 0.2
  ) +
  scale_colour_manual(
    values = c(
      "Inverse" = "#2A6EBB",
      "Positive" = "#CD202C"
    ),
    name = "Association with EDISEA"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    x = "beta coefficient",
    y = "-log10(P_FDR value)",
    title = "b. Inflammatory biomarkers associated with EDISEA"
  )
plot_ediea_inflam

# save
# save(res_metab,  file = "./dat/ukb_edisea_metabolites.rda")
# save(res_inflam, file = "./dat/ukb_edisea_inflammation.rda")

# ============================================================
# 4) SEM mediation (Metabolism PC + Inflammation PC as mediators)
# ============================================================
# choose metabolites for PCA
FDR_CUT <- 0.05
MAX_METAB_FOR_PCA <- 60

metab_pick <- res_metab %>%
  filter(!is.na(beta), !is.na(p)) %>%
  mutate(miss = sapply(y, \(nm) mean(is.na(ukb_nmr_inflam[[nm]])))) %>%
  arrange(fdr, miss) %>%
  filter(fdr < FDR_CUT)

metab_pca_cols <- metab_pick$y
message("Metabolites used for PCA: ", length(metab_pca_cols))

# inflammation PCA uses 9 markers
inflam_pca_cols <- inflam_markers
message("Inflammation used for PCA: ", length(inflam_pca_cols))

# ---- create SEM dataset
OUTCOME_VAR <- "event"

sem_cols <- unique(c("eid", EDISEA_VAR, OUTCOME_VAR, cox_covars,
                     metab_pca_cols, inflam_pca_cols))

sem_dat0 <- ukb_nmr_inflam %>%
  dplyr::select(all_of(sem_cols))

# keep complete cases for PCA variables + covars + x + outcome
sem_dat0 <- sem_dat0 %>% drop_na()

# ---- PCA for metabolism and inflammation
X_metab <- sem_dat0 %>% dplyr::select(all_of(metab_pca_cols)) %>% as.matrix()
X_infl  <- sem_dat0 %>% dplyr::select(all_of(inflam_pca_cols)) %>% as.matrix()

pc_metab <- prcomp(scale(X_metab), center = TRUE, scale. = TRUE)
pc_infl  <- prcomp(scale(X_infl),  center = TRUE, scale. = TRUE)

sem_dat <- sem_dat0 %>%
  mutate(
    Metabo_PC1 = as.numeric(pc_metab$x[, 1]),
    Inflam_PC1 = as.numeric(pc_infl$x[, 1])
  )

# standardize EDISEA & numeric covars & mediators
to_scale <- c(EDISEA_VAR, "Metabo_PC1", "Inflam_PC1", cox_covars)
for (nm in to_scale) {
  if (nm %in% names(sem_dat) && is.numeric(sem_dat[[nm]])) {
    sem_dat[[nm]] <- as.numeric(scale(sem_dat[[nm]]))
  }
}

# Decide if outcome is binary
is_binary <- is.numeric(sem_dat[[OUTCOME_VAR]]) &&
  all(na.omit(unique(sem_dat[[OUTCOME_VAR]])) %in% c(0, 1))

if (is_binary) {
  sem_dat[[OUTCOME_VAR]] <- ordered(sem_dat[[OUTCOME_VAR]])
}

# ---- SEM model: EDISEA -> (Metabo_PC1, Inflam_PC1) -> Outcome
multi_cat <- c("bmi","tdi","smoking","drinking","actc")
multi_cat <- intersect(multi_cat, names(sem_dat))

sem_dat2 <- sem_dat %>%
  mutate(across(all_of(multi_cat), ~ as.factor(.x)))

X_dummy <- model.matrix(~ . , data = sem_dat2[, multi_cat, drop = FALSE])
X_dummy <- X_dummy[, colnames(X_dummy) != "(Intercept)", drop = FALSE]

sem_dat2 <- sem_dat2 %>%
  select(-all_of(multi_cat)) %>%
  bind_cols(as.data.frame(X_dummy))

cox_covars2 <- setdiff(cox_covars, multi_cat)
cox_covars2 <- c(cox_covars2, colnames(X_dummy))
cov_str2 <- paste(cox_covars2, collapse = " + ")

model_sem <- paste0("
  Metabo_PC1 ~ a1*", EDISEA_VAR, " + ", cov_str2, "
  Inflam_PC1 ~ a2*", EDISEA_VAR, " + d*Metabo_PC1 + ", cov_str2, "
  ", OUTCOME_VAR, " ~ c*", EDISEA_VAR, " + b1*Metabo_PC1 + b2*Inflam_PC1 + ", cov_str2, "
  ind_met := a1*b1
  ind_inf := a2*b2
  ind_all := ind_met + ind_inf
  total   := c + ind_all
")

fit_sem <- sem(
  model_sem,
  data = sem_dat2,
  estimator = "WLSMV",
  ordered = OUTCOME_VAR
)

# saveRDS(fit_sem, "./dat/sem_fit.rds")

# ============================================================
# 5) Path plot with standardized coefficients + p-values (ggraph)
# ============================================================
pe <- parameterEstimates(fit_sem, standardized = TRUE) %>% as_tibble()

edge_df <- pe %>%
  filter(op == "~") %>%
  transmute(
    from = rhs,
    to   = lhs,
    est  = std.all,
    p    = pvalue
  ) %>%
  filter(
    from %in% c(EDISEA_VAR, "Metabo_PC1", "Inflam_PC1"),
    to   %in% c("Metabo_PC1", "Inflam_PC1", OUTCOME_VAR)
  )

fmt_p <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "(P < 0.001)", sprintf("(P = %.3f)", p)))
}
fmt_b <- function(b) sprintf("%.3f", b)

X_DIST <- 2.2
Y_DIST <- 1.8
node_df <- tibble(name = c(EDISEA_VAR, "Inflam_PC1", "Metabo_PC1", OUTCOME_VAR)) %>%
  mutate(
    plot_name = case_when(
      name == !!EDISEA_VAR ~ "EDISEA",
      name == "Inflam_PC1" ~ "Inflammation",
      name == "Metabo_PC1" ~ "Metabolism",
      name == !!OUTCOME_VAR ~ "Mortality",
      TRUE ~ name
    ),
    x = case_when(
      name == "Inflam_PC1" ~ -X_DIST,
      name == "Metabo_PC1" ~  X_DIST,
      name == !!EDISEA_VAR ~  0,
      name == !!OUTCOME_VAR ~  0
    ),
    y = case_when(
      name == !!EDISEA_VAR ~  Y_DIST,
      name == "Inflam_PC1" ~  0,
      name == "Metabo_PC1" ~  0,
      name == !!OUTCOME_VAR ~ -Y_DIST
    ),
    rx = case_when(name == !!EDISEA_VAR ~ 0.8, TRUE ~ 0.95),
    ry = 0.45,
    fill = case_when(
      name == !!EDISEA_VAR ~ "#B9D78A",
      name == "Inflam_PC1" ~ "#A9C7DB",
      name == "Metabo_PC1" ~ "#E6A3A0",
      name == !!OUTCOME_VAR ~ "#F0C680",
      TRUE ~ "grey90"
    )
  )

edge_plot <- edge_df %>%
  left_join(node_df %>% select(name, x, y, rx, ry), by = c("from" = "name")) %>%
  rename(x1 = x, y1 = y, rx1 = rx, ry1 = ry) %>%
  left_join(node_df %>% select(name, x, y, rx, ry), by = c("to" = "name")) %>%
  rename(x2 = x, y2 = y, rx2 = rx, ry2 = ry) %>%
  mutate(
    label_text = glue("{fmt_b(est)}\n{fmt_p(p)}"),
    
    angle = atan2(y2 - y1, x2 - x1),
    trim_r1 = (rx1 * ry1) / sqrt((ry1 * cos(angle))^2 + (rx1 * sin(angle))^2),
    x1_new = x1 + cos(angle) * trim_r1,
    y1_new = y1 + sin(angle) * trim_r1,
    
    trim_r2 = (rx2 * ry2) / sqrt((ry2 * cos(angle))^2 + (rx2 * sin(angle))^2),
    x2_new = x2 - cos(angle) * (trim_r2 + 0.05), 
    y2_new = y2 - sin(angle) * (trim_r2 + 0.05),
    
    type = case_when(
      from == !!EDISEA_VAR & to == "Inflam_PC1" ~ "top_left",
      from == !!EDISEA_VAR & to == "Metabo_PC1" ~ "top_right",
      from == !!EDISEA_VAR & to == !!OUTCOME_VAR ~ "vertical",
      from == "Inflam_PC1" & to == !!OUTCOME_VAR ~ "bottom_left",
      from == "Metabo_PC1" & to == !!OUTCOME_VAR ~ "bottom_right",
      (from == "Metabo_PC1" & to == "Inflam_PC1") | (from == "Inflam_PC1" & to == "Metabo_PC1") ~ "horizontal",
      TRUE ~ "other"
    ),
    
    curvature = case_when(
      type == "top_left"     ~ -0.3,
      type == "top_right"    ~  -0.3,
      type == "bottom_left"  ~ -0.3,
      type == "bottom_right" ~  0.3,
      type == "vertical"     ~  0.0,
      type == "horizontal"   ~  0.0,
      TRUE ~ 0
    ),

    mid_x = (x1 + x2) / 2,
    mid_y = (y1 + y2) / 2,
    
    lbl_x = case_when(
      type == "top_left"     ~ mid_x - 1.2,
      type == "top_right"    ~ mid_x + 1.2,
      type == "bottom_left"  ~ mid_x - 1.2,
      type == "bottom_right" ~ mid_x + 1.2,
      type == "vertical"     ~ mid_x - 0.7,
      type == "horizontal"   ~ mid_x + 0.8,
      TRUE ~ mid_x
    ),
    
    lbl_y = case_when(
      type == "top_left"     ~ mid_y + 0.5,
      type == "top_right"    ~ mid_y + 0.5,
      type == "bottom_left"  ~ mid_y - 0.2,
      type == "bottom_right" ~ mid_y - 0.2,
      type == "vertical"     ~ mid_y + 1.0,
      type == "horizontal"   ~ mid_y - 0.3,
      TRUE ~ mid_y
    )
  )

p_sem <- ggplot() +
  geom_curve(
    data = filter(edge_plot, curvature != 0, type %in% c('top_right', 'bottom_right')),
    aes(x = x1_new, y = y1_new, xend = x2_new, yend = y2_new),
    curvature = -0.4,
    colour = "grey50", linewidth = 0.8,
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    lineend = "round"
  ) +
  geom_curve(
    data = filter(edge_plot, curvature != 0, type %in% c('top_left', 'bottom_left')),
    aes(x = x1_new, y = y1_new, xend = x2_new, yend = y2_new),
    curvature = 0.4,
    colour = "grey50", linewidth = 0.8,
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    lineend = "round"
  ) +
  geom_segment(
    data = filter(edge_plot, curvature == 0),
    aes(x = x1_new, y = y1_new, xend = x2_new, yend = y2_new),
    colour = "grey50", linewidth = 0.8,
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    lineend = "round"
  ) +
  # 节点
  ggforce::geom_ellipse(
    data = node_df,
    aes(x0 = x, y0 = y, a = rx, b = ry, angle = 0, fill = fill),
    colour = "black", linewidth = 0.8
  ) +
  geom_text(
    data = node_df,
    aes(x = x, y = y, label = plot_name),
    size = 5, fontface = "plain"
  ) +
  geom_label(
    data = edge_plot,
    aes(x = lbl_x, y = lbl_y, label = label_text),
    fill = NA,
    label.size = 0,
    size = 4.5,
    lineheight = 0.9
  ) +
  scale_fill_identity() +
  coord_fixed(xlim = c(-3.5, 3.5), ylim = c(-2.5, 2.5), expand = FALSE) +
  theme_void() +
  labs(title = "c. Structural pathway analysis") + 
  theme(
    plot.title = element_text(face = "bold")
  )

p_sem

## gather plots
library(patchwork)
design <- "ABC"

plot_metab_inflam <- (plot_ediea_metab) +
  (plot_ediea_inflam) +
  (p_sem) +
  plot_layout(design = design, heights = c(1,1))
