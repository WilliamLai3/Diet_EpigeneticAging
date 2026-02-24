# UKBiobank Brain MRI

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
library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer('first', "dplyr")

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
load("../UKB/ukb_outcomes.rda")
brain_data <- read.csv("../UKB/dataset/participant_brain.csv")
confounders <- read_parquet("../UKB/dataset/participant_imagingConfounder.parquet")
names(brain_data) <- names(brain_data) |> str_replace("participant.", "")
names(confounders) <- names(confounders) |> str_replace("participant.", "")

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
ukb_brain <- ukb_data %>%
  inner_join(brain_data %>%
               mutate(eid = as.numeric(eid)), by = "eid") %>%
  inner_join(confounders %>%
               mutate(eid = as.numeric(eid)), by = "eid") %>%
  filter(!is.na(p26755_i2), !is.na(p25488_i2)) # thickness, ICV, FA, MD
  # filter(dieb_baseline=='no'&cvd_baseline=='no'&hbp_baseline=='no')

cox_covars <- c("age_rec", "sex", "edu_college", "eth_white", "bmi", 
                "tdi", "smoking", "drinking", "actc", "energy_kj",
                "dieb_baseline",'cvd_baseline','hbp_baseline',
                "p54_i2")
ukb_brain[,cox_covars] <- ukb_brain %>%
  dplyr::select(cox_covars) %>%
  mice(seed = 2025) %>%
  complete()

# ----------------------- MAIN ANALYSIS -------------------
idp_cols_all <- names(brain_data)
idp_cols_all <- setdiff(idp_cols_all, "eid")
idp_cols_all <- idp_cols_all[str_detect(idp_cols_all, "^p\\d+")]

idp_meta <- parse_ukb_col(idp_cols_all) %>%
  left_join(fieldsum %>% dplyr::select(field_id, title), by = "field_id") %>%
  filter(instance==2)

# ---- modality classification ----
THICK_IDS <- c(26755:26788,26856:26889)
AREA_IDS  <- c(26721:26754, 26822:26855)
CORTVOL_IDS <- c(26789:26821, 26890:26922)
CORTEX_VOL_EXTRA <- c(26552, 26583)   # cortex volume L/R
SUBCORT_IDS <- c(26558:26565, 26589:26596)
ICV_ID <- 26521
FA_IDS <- 25488:25514
MD_IDS <- 25515:25541

idp_meta <- idp_meta %>%
  mutate(
    modality = case_when(
      field_id %in% THICK_IDS ~ "Cortical thickness",
      field_id %in% AREA_IDS  ~ "Cortical area",
      field_id %in% CORTVOL_IDS ~ "Cortical volume",
      field_id %in% CORTEX_VOL_EXTRA ~ "Cortical volume",
      field_id %in% SUBCORT_IDS ~ "Subcortical volume",
      field_id %in% FA_IDS ~ "DTI FA",
      field_id %in% MD_IDS ~ "DTI MD",
      field_id == ICV_ID ~ "ICV",
      TRUE ~ "Other"
    )
  )

idp_meta <- idp_meta %>%
  filter(modality %in% c("Cortical thickness","Cortical area","Cortical volume",
                         "Subcortical volume","DTI FA","DTI MD"))

# ICV
icv_col <- parse_ukb_col(names(ukb_brain)) %>% filter(field_id == ICV_ID)

message("Brain IDPs kept: ", nrow(idp_meta))
message("ICV column: ", icv_col$col)

# ---------------------------
# 6) Run ROI-wise regressions (IDP ~ ediea + covars (+ ICV for area/volume))
# ---------------------------
base_covars <- c('ediea', cox_covars)

fit_lm_hc3_std <- function(d, y, x, covars) {
  cols <- unique(c(y, x, covars))
  dd <- d %>%
    dplyr::select(all_of(cols)) %>%
    mutate(across(all_of(c(y, x)), ~ suppressWarnings(as.numeric(.)))) %>%
    tidyr::drop_na(all_of(c(y, x)))
  
  n <- nrow(dd)
  if (n < 200) {
    return(tibble(y = y, n = n, beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  
  # z-standardize y and x (same as your original intent)
  dd[[y]] <- as.numeric(scale(dd[[y]]))
  dd[[x]] <- as.numeric(scale(dd[[x]]))
  
  # P90 - P10 on the SAME scale used in the model (x_z)
  delta_90_10 <- as.numeric(
    stats::quantile(dd[[x]], probs = 0.90, na.rm = TRUE) -
      stats::quantile(dd[[x]], probs = 0.10, na.rm = TRUE)
  )
  
  fml <- as.formula(paste0(y, " ~ ", x,
                           if (length(covars) > 0) paste0(" + ", paste(covars, collapse = " + ")) else ""))
  fit <- lm(fml, data = dd)
  
  V  <- sandwich::vcovHC(fit, type = "HC3")
  ct <- lmtest::coeftest(fit, vcov. = V)
  
  if (!x %in% rownames(ct)) {
    return(tibble(y = y, n = n, beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  
  beta_z <- unname(ct[x, "Estimate"])
  se_z   <- unname(ct[x, "Std. Error"])
  
  tibble(
    y = y, n = n,
    beta = beta_z * delta_90_10,   # 90th-10th effect
    se   = se_z   * delta_90_10,   # SE scales the same
    p    = unname(ct[x, "Pr(>|t|)"]),
    delta_90_10 = delta_90_10
  )
}

run_modality <- function(mod_name) {
  dat <- ukb_brain
  
  add_icv <- (mod_name %in% c("Cortical thickness", "Cortical area","Cortical volume","Subcortical volume"))
  covs <- c(cox_covars, if (add_icv) icv_col$col else NULL)
  
  ylist <- idp_meta %>% filter(modality == mod_name)
  
  res <- purrr::map_dfr(
    ylist$col,
    ~fit_lm_hc3_std(dat, y = .x, x = 'ediea', covars = covs)
  )
  
  res <- res %>%
    left_join(ylist %>% dplyr::select(col, field_id, title, modality),
              by = c("y" = "col")) %>%
    mutate(
      fdr = p.adjust(p, method = "fdr"),
      sig = !is.na(fdr) & fdr < 0.05
    )
  
  res
}

res_thick <- run_modality("Cortical thickness")
res_fa    <- run_modality("DTI FA") # NO significant assoc
res_md    <- run_modality("DTI MD")
res_subcortical_vol <- run_modality("Subcortical volume")

parse_fs_aparc <- function(title){
  x <- str_to_lower(title)

  hemi <- case_when(
    str_detect(x, "\\bleft\\b|\\blh\\b")  ~ "lh",
    str_detect(x, "\\bright\\b|\\brh\\b") ~ "rh",
    TRUE ~ NA_character_
  )

  region <- x %>%
    str_replace_all("\\(.*?\\)", " ") %>%
    str_replace_all("mean|average|thickness|cortical|mm", " ") %>%
    str_replace_all("\\bleft\\b|\\bright\\b|\\blh\\b|\\brh\\b", " ") %>%
    str_replace_all("of", " ") %>%
    str_replace_all("[^a-z\\s-]", " ") %>%
    str_squish()

  region_fs <- region %>% str_replace_all("[\\s-]", "")

  tibble(hemi = hemi, region_fs = region_fs)
}

plot_thick <- res_thick %>%
  mutate(tmp = map(title, parse_fs_aparc)) %>%
  unnest(tmp) %>%
  filter(!is.na(hemi), !is.na(region_fs), !is.na(beta))

# names=region label, value=beta
lh_vals <- plot_thick %>% filter(hemi=="lh") %>% select(region_fs, beta) %>% deframe()
rh_vals <- plot_thick %>% filter(hemi=="rh") %>% select(region_fs, beta) %>% deframe()
region_values <- list(lh = lh_vals, rh = rh_vals)


lim <- max(abs(plot_thick$beta), na.rm = TRUE)

## FIGURE
library(fsbrain)
subjects_dir <- "/Applications/freesurfer/8.0.0/subjects"
col_fn <- grDevices::colorRampPalette(c("#f7fbff", "#08306b"))

lh_list <- region_values$lh
rh_list <- region_values$rh

cm <- fsbrain::vis.region.values.on.subject(
  subjects_dir = subjects_dir,
  subject_id   = "fsaverage",
  atlas        = "aparc",
  lh_region_value_list = lh_list,
  rh_region_value_list = rh_list,
  views        = NULL,
  surface      = "pial",
  draw_colorbar = TRUE,
  makecmap_options = list(
    n = 200,
    colFn = col_fn,
    range = c(0, lim)
  ),
  silent = TRUE
)

fsbrain::export(
  cm,
  view_angles = fsbrain::get.view.angle.names(angle_set = "t4"),
  grid_like   = FALSE,
  output_img  = "fig/brain_thickness.png"
)

## subcortical
parse_aseg_subcort <- function(title, acc_region = NULL){
  x <- str_to_lower(title)
  
  hemi <- case_when(
    str_detect(x, "\\bleft\\b|\\blh\\b")  ~ "left",
    str_detect(x, "\\bright\\b|\\brh\\b") ~ "right",
    TRUE ~ NA_character_
  )
  
  region <- case_when(
    str_detect(x, "thalam")   ~ "thalamus proper",
    str_detect(x, "caudat")   ~ "caudate",
    str_detect(x, "putamen")  ~ "putamen",
    str_detect(x, "pallid")   ~ "pallidum",
    str_detect(x, "hippoc")   ~ "hippocampus",
    str_detect(x, "amygdal")  ~ "amygdala",
    str_detect(x, "ventral") | str_detect(x, "ventraldc") ~ "ventral DC",
    str_detect(x, "accumb")   ~ 'accumbens',
    TRUE ~ NA_character_
  )
  
  tibble(hemi = hemi, region = region)
}

subcort_plot <- res_subcortical_vol %>%
  mutate(tmp = purrr::map(title, parse_aseg_subcort)) %>%
  tidyr::unnest(tmp) %>%
  filter(!is.na(hemi), !is.na(region), !is.na(beta)) %>%
  mutate(
    side = "coronal",
    beta_sig = if_else(fdr < 0.05, beta, NA_real_),
    region_label = case_when(
      region == "thalamus proper" ~ "Thalamus",
      region == "accumbens area" ~ "Accumbens",
      region == "ventral dc"     ~ "Ventral DC",
      TRUE ~ str_to_title(region)
    ),
    hemi_label = if_else(hemi == "left", "L", "R"),
    roi_label  = paste0(hemi_label, " ", region_label),
    ci_low = beta - 1.96 * se,
    ci_high = beta + 1.96 * se
  )

library(ggseg)
data(aseg, package = "ggseg")
aseg_regions <- unique(na.omit(aseg$data$region))
subcort_plot <- subcort_plot %>%
  filter(region %in% aseg_regions)

lim_sub <- max(abs(subcort_plot$beta), na.rm = TRUE)
p_aseg <- ggplot() +
  geom_brain(
    atlas = aseg,
    data  = subcort_plot,
    aes(fill = beta_sig),
    side = "coronal",
    colour = "white",
    linewidth = 0.15
  ) +
  scale_fill_gradient2(
    midpoint = 0,
    limits = c(-lim_sub, lim_sub),
    na.value = "grey90",
    name = "β coefficient"
  ) +
  theme_brain()

p_aseg

ggsave("fig/brain_subcortical_aseg.pdf", p_aseg, width = 7.2, height = 3.2)