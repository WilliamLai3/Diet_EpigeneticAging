# EDIEA and mortality in NHANES
library(tidyverse)
library(survival)
library(mice)
library(survey)
library(linkET)
library(patchwork)
library(survminer)

load('./dat/nhanes_ediea_from_hrs.rda')
load("../NHANES/Diet&mor/MOR_A_J_19.rda")
epi_clock <- haven::read_sas("../NHANES/cov/dnmepi.sas7bdat")
data_eaa_nhanes_all <- data_eaa_nhanes_all %>%
  select(SEQN, names(data_eaa_nhanes_all)[2:39], ediea, CALOR_SUM)

load("../NHANES/cov/Cov_AJ.rda")

ALQ_AJ <- ALQ_AJ %>%
  mutate(
    SEQN = as.numeric(SEQN)
  ) %>%
  group_by(SEQN) %>%
  slice_tail(n = 1) %>%
  ungroup()

DEMO_select <- DEMO_AJ %>%
  select(SEQN, AGE, GENDER, RACE, EDU, PIR, SDMVPSU, SDMVSTRA, WTMEC2YR, WAVE) %>%
  left_join(
    SMQ_AJ,
    by = "SEQN"
  ) %>%
  left_join(
    ALQ_AJ %>% mutate(SEQN = as.numeric(SEQN)),
    by = "SEQN"
  ) %>%
  left_join(
    PA_AJ %>% mutate(SEQN = as.numeric(SEQN)),
    by = "SEQN"
  ) %>%
  left_join(
    BMI_AJ,
    by = "SEQN"
  ) %>%
  left_join(
    DIA_AJ,
    by = "SEQN"
  ) %>%
  left_join(
    BP_AJ,
    by = "SEQN"
  ) %>%
  left_join(
    CVD_AJ,
    by = "SEQN"
  )


health_status <- bind_rows(MOR_A, MOR_B, MOR_C, MOR_D, MOR_E,
                           MOR_F, MOR_G, MOR_H, MOR_I, MOR_J) %>%
  mutate(SEQN = as.numeric(seqn)) %>%
  select(-c('seqn'))

# merge datasets
nhanes <- data_eaa_nhanes_all %>% 
  # left_join(hrsfg_adj, by = c("SEQN")) %>%
  left_join(DEMO_select, by = c("SEQN")) %>%
  left_join(health_status, by = c("SEQN"))

nhanes <- nhanes %>% filter(!is.na(ediea)) # 90907 obs
nhanes <- nhanes %>% # 49838 obs
  filter(eligstat==1 & permth_int >=24 & AGE>=45) # 23830

nhanes <- nhanes %>% 
  mutate(
    cod = ifelse(!is.na(ucod_leading), ucod_leading, 0),
    fu_death = permth_int / 12,
    ad_death = ifelse(cod==6, 1, 0),
    ediea_t = cut_number(ediea, 3, labels = FALSE),
    # Grim_EAA_t = cut_number(Grim_EAA, 3, labels = FALSE),
    ediea_group = case_when(
      ediea_t==1 ~ 'low',
      ediea_t==2 ~ 'medium',
      ediea_t==3 ~ 'high',
      TRUE ~ NA_character_
    ),
    ediea_group = factor(ediea_group, levels = c("low","medium","high"))
    # ediea = scale(ediea)
  )

# Association with mortality
covs <- c("AGE", "GENDER", "RACE", "EDU", "PIR","METS",
          "SMO", "ALC", "BMIC", "CALOR_SUM", "CVD", "HBP", "DIA") # PIR, METS: >10% missing
nhanes[,covs] <- nhanes %>%
  select(covs) %>%
  mice(seed = 2025) %>%
  complete()


# Cox weighted
nhanesSvy <- svydesign(
  ids    = ~ SDMVPSU,
  strata = ~ SDMVSTRA,
  weights = ~ WTMEC2YR,
  nest   = TRUE,
  data   = nhanes
)


cph1_sum <- nhanes %>% 
  group_by(ediea_group) %>% 
  mutate(fol = permth_int/12) %>% 
  summarise(case = sum(mortstat),
            py = sum(fol))

cph_w_m1 <- svycoxph(
  Surv(fu_death, mortstat) ~ ediea_group + AGE + factor(GENDER) + CALOR_SUM + factor(WAVE),
  design = nhanesSvy
)
cph_w_m1_tidy <- tidy(cph_w_m1) %>% 
  filter(str_detect(term, "ediea_group")) %>%
  mutate(HR = sprintf("%.2f (%.2f-%.2f)",
                      exp(estimate),
                      exp(estimate - 1.96 * std.error),
                      exp(estimate + 1.96 * std.error)))

cph_w_m2 <- svycoxph(
  Surv(fu_death, mortstat) ~ ediea_group + AGE + factor(GENDER) + factor(RACE) + factor(EDU) +
    PIR + METS + factor(SMO) + factor(ALC) + factor(BMIC) + CALOR_SUM + factor(WAVE),
  design = nhanesSvy
)
cph_w_m2_tidy <- tidy(cph_w_m2) %>% 
  filter(str_detect(term, "ediea_group")) %>%
  mutate(HR = sprintf("%.2f (%.2f-%.2f)",
                      exp(estimate),
                      exp(estimate - 1.96 * std.error),
                      exp(estimate + 1.96 * std.error)))

cph_w_m3 <- svycoxph(
  Surv(fu_death, mortstat) ~ ediea_group + AGE + factor(GENDER) + factor(RACE) + factor(EDU) +
    PIR + METS + factor(SMO) + factor(ALC) + factor(BMIC) + CALOR_SUM + factor(WAVE) +DIA+HBP+CVD,
  design = nhanesSvy
)
nhanes_dth <- summary(cph_w_m3)$coefficients[1,]
cph_w_m3_tidy <- tidy(cph_w_m3) %>% 
  filter(str_detect(term, "ediea_group")) %>%
  mutate(HR = sprintf("%.2f, %.2f-%.2f",
                      exp(estimate),
                      exp(estimate - 1.96 * std.error),
                      exp(estimate + 1.96 * std.error)))

## K-M plots
nhanes$wts <- as.numeric(weights(nhanesSvy))
cox_w_m3_plot <- coxph(
  Surv(fu_death, mortstat) ~ ediea_group + AGE + factor(GENDER) +
    factor(RACE) + factor(EDU) + PIR + METS + factor(SMO) + factor(ALC) +
    factor(BMIC) + CALOR_SUM + factor(WAVE) +DIA+HBP+CVD,
  data    = nhanes,
  weights = wts
)

newdat <- nhanes |>
  summarise(
    AGE        = mean(AGE, na.rm = TRUE),
    PIR        = mean(PIR, na.rm = TRUE),
    METS       = mean(METS, na.rm = TRUE),
    CALOR_SUM  = mean(CALOR_SUM, na.rm = TRUE),
    GENDER     = first(GENDER),
    RACE       = first(RACE),
    EDU        = first(EDU),
    SMO        = first(SMO),
    ALC = first(ALC),
    BMIC       = first(BMIC),
    WAVE       = first(WAVE),
    HBP        = 0,
    DIA        = 0,
    CVD        = 0
  )

newdat <- bind_rows(
  mutate(newdat, ediea_group = factor("low",    levels = c("low","medium","high"))),
  mutate(newdat, ediea_group = factor("medium", levels = c("low","medium","high"))),
  mutate(newdat, ediea_group = factor("high",   levels = c("low","medium","high")))
)

fit_adj <- survfit(cox_w_m3_plot, newdata = newdat)
hr_med  <- cph_w_m3_tidy %>% filter(term == "ediea_groupmedium") %>% pull(HR)
hr_high <- cph_w_m3_tidy %>% filter(term == "ediea_grouphigh")   %>% pull(HR)
labs_edisea <- c(
  "Low \n(ref)",
  paste0("Medium \n(", hr_med,  ")"),
  paste0("High \n(",   hr_high, ")")
)
km_nhanes <- ggsurvplot(
  fit_adj,
  data        = newdat,
  # palette = 'lancet',
  conf.int    = FALSE,
  legend.labs  = labs_edisea,
  legend.title = "EDISEA (HR, 95% CI)",
  ylim = c(0.5,1.0),
  xlab        = "Follow-up time (years)",
  ylab        = "Adjusted survival probability",
  title       = "b. Mortality in NHANES",
  censor      = FALSE,
  censor.size = 1,
  censor.shape = 124
)
km_nhanes$plot <- km_nhanes$plot +
  theme(plot.title = element_text(face = "bold"))

km_nhanes$plot

fit_ediea_models <- function(data,
                             time,
                             event,
                             disease,
                             cov_m3_comorb = NULL
){
  dat <- data %>%
    filter(.data[[time]] > 0)
  
  ediea_vec <- dat[["ediea"]]
  ediea_min <- quantile(ediea_vec, 0.10, na.rm = TRUE)
  ediea_max <- quantile(ediea_vec, 0.90, na.rm = TRUE)
  ediea_step20 <- (ediea_max - ediea_min)

  cov_m2 <- c("AGE", "GENDER", "RACE", "EDU", "PIR","METS",
              "SMO", "ALC", "BMIC", "CALOR_SUM")

  f_m2 <- as.formula(
    paste0("Surv(", time, ", ", event, ") ~ ediea + ",
           paste(cov_m2, collapse = " + "))
  )
  
  models <- list(
    model2 = svycoxph(f_m2, design = nhanesSvy)
  )
  
  res_group <- map_dfr(names(models), function(mname) {
    models[[mname]] %>%
      tidy() %>%
      filter(str_detect(term, "ediea")) %>%
      mutate(
        outcome = disease,
        model   = mname,
        HR      = exp(estimate),
        LCI     = exp(estimate - 1.96 * std.error),
        UCI     = exp(estimate + 1.96 * std.error),
        HR_label = sprintf("%.2f (%.2f-%.2f)", HR, LCI, UCI),
        ediea_level = str_remove(term, "^ediea_group") %>% factor(levels=c('low','medium','high')),
        coef = estimate,
        se.coef. = std.error
      )
  })

  f_m2_cont <- as.formula(
    paste0("Surv(", time, ", ", event, ") ~ ediea + ",
           paste(cov_m2, collapse = " + "))
  )
  
  models_cont <- list(
    model2 = svycoxph(f_m2_cont, design = nhanesSvy)
  )

  res_per20 <- map_dfr(names(models_cont), function(mname) {
    tt <- models_cont[[mname]] %>%
      tidy() %>%
      filter(str_detect(term, 'ediea')) %>%
      clusterProfiler::slice(1)
    
    beta  <- tt$estimate
    se    <- tt$std.error

    logHR20  <- beta * ediea_step20
    se_log20 <- se   * ediea_step20
    
    HR  <- exp(logHR20)
    LCI <- exp(logHR20 - 1.96 * se_log20)
    UCI <- exp(logHR20 + 1.96 * se_log20)
    
    tibble(
      outcome = disease,
      model   = mname,
      term    = paste0("per20pct"),
      HR      = HR,
      LCI     = LCI,
      UCI     = UCI,
      HR_label = sprintf("%.2f (%.2f-%.2f)", HR, LCI, UCI),
      p.value = tt$p.value,
      ediea_level = factor("per20pct", levels = c("low","medium","high","per20pct")),
      coef = logHR20,
      se.coef. = se_log20
    )
  })

  bind_rows(res_group, res_per20)
}

res_nhanes <- bind_rows(
  # Mortality
  fit_ediea_models(nhanes, "fu_death", "mortstat",
                   disease = "Mortality", cov_m3_comorb = c("t2d", "hbp")),
)


cases_py_nhanes <- bind_rows(
  # Mortality
  nhanes %>%
    filter(fu_death > 0) %>%
    summarise(
      outcome = "Mortality",
      cases   = sum(mortstat == 1, na.rm = TRUE),
      py      = sum(fu_death, na.rm = TRUE)
    )
) %>%
  mutate(
    case_py_label = paste0(cases, " / ", round(py, 0))
  )

df_plot_hr_irr_nhanes <- df_plot_nhanes %>%
  filter(ediea_label == "90th-10th percentile difference") %>%
  left_join(
    cases_py_nhanes %>% select(outcome, case_py_label),
    by = "outcome"
  ) %>%
  mutate(outcome = factor(outcome, levels = unique(outcome)))