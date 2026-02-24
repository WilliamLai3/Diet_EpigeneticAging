# EDISEA and mortality in HRS
library(tidyverse)
library(survival)
library(mice)
library(ggjournals)
library(readr)
library(haven)
library(patchwork)
library(meta)
library(ggplot2)
library(scales)
library(grid)
source("./plot_rcs_linear.R")

load("dat/hrs_cov.rda")
load("dat/hrs_ediea.rda")
# load('dat/hrs_ediea_from_nhanes.rda'
load('dat/hrs_health.RData')
load('dat/hrsdem1422_dates.rda')
epi_clock <- haven::read_sas("dat/epiclocka_r.sas7bdat")
data_eaa <- h13_fg_adj %>% select(HHID, PN, ediea)
# data_eaa <- data_eaa_hrs_all %>% select(HHID, PN, ediea)

# merge datasets
hrs <- hrs_cov %>% 
  # left_join(h13_fg_adj, by = c("HHID", "PN")) %>%
  left_join(data_eaa, by = c("HHID", "PN")) %>%
  left_join(health_status, by = c("HHID", "PN"))

hrs <- hrs %>% filter(!is.na(ediea)) # 7410 obs
hrs <- hrs %>% filter(!is.na(death)) # 7398 obs

eof <- max(hrs$dod, na.rm = TRUE)
hrs <- hrs %>% 
  mutate(
    fu_death = ifelse(is.na(dod), eof-2013, dod-2013),
    ediea_t = cut_number(ediea, 3, labels = FALSE) %>% factor(),

    ediea_group = case_when(
      ediea_t==1 ~ 'low',
      ediea_t==2 ~ 'medium',
      ediea_t==3 ~ 'high',
      TRUE ~ NA_character_
    ),
    ediea_group = factor(ediea_group, levels = c("low","medium","high")),

    ageg = age %/% 5,
    ageg = ifelse(ageg<9, 9, ifelse(ageg>18, 18, ageg))
  )

# Association with mortality
covs <- c("age", "sex", "race", "educ", "hhinc", "smoke", 
          "drink", "actc", "bmi_cat", "toten", "cesd", "t2d", "hypten", "cvd")
hrs[,covs] <- hrs %>% 
  select(covs) %>% 
  mice(seed = 2025) %>%
  complete()

# Model 1/2/3 Mortality
cph1_sum <- hrs %>% 
  group_by(ediea_group) %>% 
  mutate(fol = fu_death) %>% 
  summarise(case = sum(death),
            py = sum(fol))

cph_m1 <- coxph(Surv(fu_death, death) ~ ediea_group + age + sex + toten,
                data = hrs)
cph_m1_tidy <- tidy(cph_m1) %>% 
  filter(str_detect(term, "ediea_group")) %>% 
  mutate(HR = sprintf("%.2f (%.2f-%.2f)", exp(estimate), 
                      exp(estimate - 1.96*std.error),
                      exp(estimate + 1.96*std.error)
  ))

cph_m2 <- coxph(Surv(fu_death, death) ~ ediea_group + age + sex + race + educ + hhinc + 
                  smoke + drink + actc + bmi_cat + toten, data = hrs)
cph_m2_tidy <- tidy(cph_m2) %>% 
  filter(str_detect(term, "ediea_group")) %>% 
  mutate(HR = sprintf("%.2f (%.2f-%.2f)", exp(estimate), 
                      exp(estimate - 1.96*std.error),
                      exp(estimate + 1.96*std.error)
  ))


cph_m3 <- coxph(Surv(fu_death, death) ~ ediea_group + age + sex + race + educ + hhinc + 
                 smoke + drink + actc + bmi_cat + toten + t2d + hypten + 
                 cvd, data = hrs)
cph_m3_tidy <- tidy(cph_m3) %>% 
  filter(str_detect(term, "ediea_group")) %>% 
  mutate(HR = sprintf("%.2f, %.2f-%.2f", exp(estimate), 
                      exp(estimate - 1.96*std.error),
                      exp(estimate + 1.96*std.error)
  ))

## K-M plots
library(survminer)

newdat <- hrs |>
  summarise(
    age    = mean(age, na.rm = TRUE),
    hhinc  = first(hhinc),
    toten  = mean(toten, na.rm = TRUE),
    cesd   = mean(cesd, na.rm = TRUE),
    sex    = first(sex),
    race   = first(race),
    educ   = first(educ),
    smoke  = first(smoke),
    drink  = first(drink),
    actc   = first(actc),
    bmi_cat = first(bmi_cat),
    t2d    = 0,
    hypten = 0,
    cvd    = 0
  )

newdat <- bind_rows(
  mutate(newdat, ediea_group = factor("low",    levels = c("low","medium","high"))),
  mutate(newdat, ediea_group = factor("medium", levels = c("low","medium","high"))),
  mutate(newdat, ediea_group = factor("high",   levels = c("low","medium","high")))
)

hr_med  <- cph_m3_tidy %>% filter(term == "ediea_groupmedium") %>% pull(HR)
hr_high <- cph_m3_tidy %>% filter(term == "ediea_grouphigh")   %>% pull(HR)
labs_edisea <- c(
  "Low \n(ref)",
  paste0("Medium \n(", hr_med,  ")"),
  paste0("High \n(",   hr_high, ")")
)
fit_adj <- survfit(cph_m3, newdata = newdat)
km_hrs <- ggsurvplot(
  fit_adj,
  data        = newdat,
  # palette     = "lancet",
  conf.int    = FALSE,
  legend.labs  = labs_edisea,
  legend.title = "EDISEA (HR, 95% CI)",
  ylim = c(0.8, 1),
  xlab        = "Follow-up time (years)",
  ylab        = "Adjusted survival probability",
  title       = "a. Mortality in HRS"
)
km_hrs$plot <- km_hrs$plot +
  theme(plot.title = element_text(face = "bold"))

# --- function: extract HR ---
extract_hr <- function(model_tidy){
  model_tidy %>%
    tidy() %>%
    filter(str_detect(term, "ediea_group")) %>%
    mutate(
      ediea_group = str_replace(term, "ediea_group", ""),
      ediea_group = recode(ediea_group,
                           "low"="low",
                           "medium"="medium",
                           "high"="high"),
      HR = sprintf("%.2f (%.2f-%.2f)",
                   exp(estimate),
                   exp(estimate - 1.96 * std.error),
                   exp(estimate + 1.96 * std.error))
    ) %>% 
    select(ediea_group, HR)
}
# --- extract HR tables ---
tbl_m1 <- extract_hr(cph_m1) %>% rename(Model1 = HR)
tbl_m2 <- extract_hr(cph_m2) %>% rename(Model2 = HR)
tbl_m3 <- extract_hr(cph_m3) %>% rename(Model3 = HR)
# --- merge into one wide table ---
final_tbl_hrs <- cph1_sum %>% 
  left_join(tbl_m1, by = "ediea_group") %>% 
  left_join(tbl_m2, by = "ediea_group") %>% 
  left_join(tbl_m3, by = "ediea_group") %>% 
  select(ediea_group, case, py, Model1, Model2, Model3)

final_tbl_hrs

# Mortality HR table
library(ggpubr)
tbl_for_plot <- final_tbl_hrs %>%
  mutate(
    ediea_group = factor(ediea_group,
                         levels = c("low", "medium", "high")),
    py = round(py, 0)
    # ediea_group = stringr::str_to_title(ediea_group)
  ) %>%
  arrange(ediea_group) %>%
  rename(
    Group  = ediea_group,
    Cases = case,
    `Person-yrs` = py
  )

table_plot_hrs <- ggtexttable(
  tbl_for_plot,
  rows  = NULL,
  theme = ttheme("classic")
)

# --------- Associations with mean scores, adjusted for baseline scores ----------
hrs_scores <- hrs %>%
  mutate(
    missing = is.na(date12) + is.na(date13) + is.na(date14) + is.na(date15),
    cesd_mean = rowMeans(select(hrs, c("cesd12", "cesd13",
                                       "cesd14", "cesd15")), na.rm = TRUE),
    cesd_chg13 = (cesd13 - cesd12)/2,
    cesd_chg14 = (cesd14 - cesd12)/4,
    cesd_chg15 = (cesd15 - cesd12)/6,

    tics_mean = rowMeans(select(hrs, c("tics27_12", "tics27_13",
                                       "tics27_14", "tics27_15")), na.rm = TRUE),
    tics_chg13 = (tics27_13 - tics27_12)/2,
    tics_chg14 = (tics27_14 - tics27_12)/4,
    tics_chg15 = (tics27_15 - tics27_12)/6,

    adl5_mean = rowMeans(select(hrs, c("adl5_12", "adl5_13",
                                       "adl5_14", "adl5_15")), na.rm = TRUE),
    adl5_chg13 = (adl5_13 - adl5_12)/2,
    adl5_chg14 = (adl5_14 - adl5_12)/4,
    adl5_chg15 = (adl5_15 - adl5_12)/6,

    iadl5_mean = rowMeans(select(hrs, c("iadl5_12", "iadl5_13",
                                        "iadl5_14", "iadl5_15")), na.rm = TRUE),
    iadl5_chg13 = (iadl5_13 - iadl5_12)/2,
    iadl5_chg14 = (iadl5_14 - iadl5_12)/4,
    iadl5_chg15 = (iadl5_15 - iadl5_12)/6,
    
    adl5_inc15   = pmax(6*adl5_chg15, 0),
    iadl5_inc15  = pmax(6*iadl5_chg15, 0),
    cesd_inc15 = pmax(6*cesd_chg15, 0)

  ) %>%
  filter(missing<4)

hrs_scores <- hrs_scores %>%
  mutate(
    cesd_chg = rowMeans(select(hrs_scores, c("cesd_chg13", "cesd_chg14", "cesd_chg15")), na.rm = TRUE),
    tics_chg = rowMeans(select(hrs_scores, c("tics_chg13", "tics_chg14", "tics_chg15")), na.rm = TRUE),
    adl5_chg = rowMeans(select(hrs_scores, c("adl5_chg13", "adl5_chg14", "adl5_chg15")), na.rm = TRUE),
    iadl5_chg = rowMeans(select(hrs_scores, c("iadl5_chg13", "iadl5_chg14", "iadl5_chg15")), na.rm = TRUE)
  )

# TICS
lm_tics <- lm(tics_chg15 ~ ediea + age +
                sex + race + educ + hhinc +
                smoke + drink + actc + bmi_cat + toten + tics27_12, data = hrs_scores)
summary(lm_tics)

plot_rcs_tics <- plot_rcs_linear(data = hrs_scores, x = "ediea", y = "tics_chg15",
                                 prob = 0.1,
                                 covs = c("age", "sex", "race", "educ", "hhinc",
                                          "smoke", "drink", "actc", "bmi_cat", "toten", "tics27_12")) +
  theme_bmj() +
  labs(title = "TICS-27 change 2014-2020 (HRS)",
       x = "EDISEA", y = "TICS-27 annual change")
plot_rcs_tics

# ADL
qpoi_adl <- lm(adl5_chg15 ~ediea + age +
                  sex + race + educ + hhinc +
                  smoke + drink + bmi_cat +
                  toten + adl5_12, data = hrs_scores %>%
                  filter(adl5_chg15>=0),
                family = quasipoisson(link = "log")
                )
summary(qpoi_adl)
hrs_adl <- summary(qpoi_adl)$coefficients[2,]

plot_rcs_adl <- plot_rcs_linear(data = hrs_scores, x = "ediea", y = "adl5_chg15",
                                 prob = 0.1,
                                 covs = c("age", "sex", "race", "educ", "hhinc",
                                          "smoke", "drink", "actc", "bmi_cat", "toten", "adl5_12")) +
  theme_bmj() +
  labs(title = "ADL change 2014-2020 (HRS)",
       x = "EDISEA", y = "ADL annual change")
plot_rcs_adl

# IADL
qpoi_iadl <- lm(iadl5_chg15 ~ediea + age +
                   sex + race + educ + hhinc +
                   smoke + drink + actc + bmi_cat +
                   toten + iadl5_12, data = hrs_scores %>%
                   filter(iadl5_chg15>=0),
                 family = quasipoisson(link = "log"))
summary(qpoi_iadl)

plot_rcs_iadl <- plot_rcs_linear(data = hrs_scores, x = "ediea", y = "iadl5_chg15",
                                prob = 0.1,
                                covs = c("age", "sex", "race", "educ", "hhinc",
                                         "smoke", "drink", "actc", "bmi_cat", "toten", "iadl5_12")) +
  theme_bmj() +
  labs(title = "IADL change 2014-2020 (HRS)",
       x = "EDISEA", y = "IADL annual change")
plot_rcs_iadl

# ----------------------- Associations with incident diseases ---------------------
eof_hbp <- max(hrs$date_hbp, na.rm = TRUE)
eof_diabetes <- max(hrs$date_diabetes, na.rm = TRUE)
eof_lung <- max(hrs$date_lung, na.rm = TRUE)
eof_heart <- max(hrs$date_heart, na.rm = TRUE)
eof_stroke <- max(hrs$date_stroke, na.rm = TRUE)
eof_cancer <- max(hrs$date_cancer, na.rm = TRUE)
eof_arthritis <- max(hrs$date_arthritis, na.rm = TRUE)
eof_dem <- max(hrs$date_dem, na.rm = TRUE)
dep_cut <- 4
baseline_date <- as.Date("2013-12-31")
hrs_diseases <- hrs %>%
  left_join(hrsdem_dates %>% select(HHID, PN, dem_new=dem, time2_dem), by=c('HHID', 'PN')) %>%
  mutate(
    # date_dem = if_else(is.na(date_dem), if_else(dem_new==1, time2_dem, NA_Date_), date_dem)
    date_dem = if_else(dem_new==1, time2_dem, date_dem),
  ) %>%
  mutate(
    missing = is.na(date12) + is.na(date13) + is.na(date14) + is.na(date15),
    last_date_dem = ifelse(!is.na(time2_dem), time2_dem, 
                       ifelse(!is.na(date15), date15,
                            ifelse(!is.na(date14), date14,
                                ifelse(!is.na(date13), date13, date12)))),
    last_date = ifelse(!is.na(date15), date15,
                    ifelse(!is.na(date14), date14,
                        ifelse(!is.na(date13), date13, date12))),
    e_hbp = ifelse(is.na(date_hbp), 0, 1),
    e_diab = ifelse(is.na(date_diabetes), 0, 1),
    e_lung = ifelse(is.na(date_lung), 0, 1),
    e_heart = ifelse(is.na(date_heart), 0, 1),
    e_stroke = ifelse(is.na(date_stroke), 0, 1),
    e_cancer = ifelse(is.na(date_cancer), 0, 1),
    e_arth = ifelse(is.na(date_arthritis), 0, 1),
    e_dem = ifelse(is.na(date_dem), 0, 1),
    e_ci = ifelse(is.na(date_ci), 0, 1),
    e_dem_ci = ifelse(!is.na(date_dem)|!is.na(date_ci), 1, 0),
    fu_hbp = ifelse(is.na(date_hbp),
                    as.Date(last_date)-as.Date("2013-12-31"),
                    date_hbp-as.Date("2013-12-31"))/365.25,
    fu_diab = ifelse(is.na(date_diabetes),
                     as.Date(last_date)-as.Date("2013-12-31"),
                     date_diabetes-as.Date("2013-12-31"))/365.25,
    fu_lung = ifelse(is.na(date_lung),
                     as.Date(last_date)-as.Date("2013-12-31"),
                     date_lung-as.Date("2013-12-31"))/365.25,
    fu_heart = ifelse(is.na(date_heart),
                      as.Date(last_date)-as.Date("2013-12-31"),
                      date_heart-as.Date("2013-12-31"))/365.25,
    fu_stroke = ifelse(is.na(date_stroke),
                       as.Date(last_date)-as.Date("2013-12-31"),
                       date_stroke-as.Date("2013-12-31"))/365.25,
    fu_cancer = ifelse(is.na(date_cancer),
                       as.Date(last_date)-as.Date("2013-12-31"),
                       date_cancer-as.Date("2013-12-31"))/365.25,
    fu_arth = ifelse(is.na(date_arthritis),
                     as.Date(last_date)-as.Date("2013-12-31"),
                     date_arthritis-as.Date("2013-12-31"))/365.25,
    fu_dem = ifelse(is.na(date_dem),
                    as.Date(last_date_dem)-as.Date("2013-12-31"),
                    date_dem-as.Date("2013-12-31"))/365.25,
    fu_ci = ifelse(is.na(date_ci),
                   as.Date(last_date)-as.Date("2013-12-31"),
                   as.Date(date_ci)-as.Date("2013-12-31"))/365.25,
    fu_dem_ci = pmin(fu_dem, fu_ci),
    
    dep_baseline = if_else(cesd12 >= dep_cut, 1L, 0L),
    date_dep = case_when(
      cesd12 >= dep_cut ~ date12,
      cesd13 >= dep_cut ~ date13,
      cesd14 >= dep_cut ~ date14,
      cesd15 >= dep_cut ~ date15,
      TRUE ~ NA_Date_
    ),
    e_dep = if_else(!is.na(date_dep), 1L, 0L),
    fu_dep = if_else(
      e_dep == 1L,
      as.numeric(as.Date(date_dep) - as.Date("2013-12-31")),
      as.numeric(as.Date(last_date) - as.Date("2013-12-31"))
    )/365.25
  ) %>%
  mutate(
    prev_hbp  = if_else(!is.na(date_hbp)      & date_hbp      <= baseline_date, 1L, 0L),
    prev_diab = if_else(!is.na(date_diabetes) & date_diabetes <= baseline_date, 1L, 0L),
    prev_lung = if_else(!is.na(date_lung)     & date_lung     <= baseline_date, 1L, 0L),
    prev_heart= if_else(!is.na(date_heart)    & date_heart    <= baseline_date, 1L, 0L),
    prev_stroke = if_else(!is.na(date_stroke) & date_stroke   <= baseline_date, 1L, 0L),
    prev_cancer = if_else(!is.na(date_cancer) & date_cancer   <= baseline_date, 1L, 0L),
    prev_arth   = if_else(!is.na(date_arthritis) & date_arthritis <= baseline_date, 1L, 0L),
    
    prev_chronic_n = prev_hbp + prev_diab + prev_lung + prev_heart + prev_stroke +
      prev_cancer + prev_arth,
    
    inc_hbp  = if_else(!is.na(date_hbp)      & date_hbp      > baseline_date, 1L, 0L),
    inc_diab = if_else(!is.na(date_diabetes) & date_diabetes > baseline_date, 1L, 0L),
    inc_lung = if_else(!is.na(date_lung)     & date_lung     > baseline_date, 1L, 0L),
    inc_heart= if_else(!is.na(date_heart)    & date_heart    > baseline_date, 1L, 0L),
    inc_stroke = if_else(!is.na(date_stroke) & date_stroke   > baseline_date, 1L, 0L),
    inc_cancer = if_else(!is.na(date_cancer) & date_cancer   > baseline_date, 1L, 0L),
    inc_arth   = if_else(!is.na(date_arthritis) & date_arthritis > baseline_date, 1L, 0L),
    
    inc_chronic_n = inc_hbp + inc_diab + inc_lung +
      inc_heart + inc_stroke + inc_cancer + inc_arth,
    
    fu_chronic_year = as.numeric(as.Date(last_date) - baseline_date) / 365.25
  )

# ----------------- health outcome HR ----------------
cov_m1_base   <- c("age", "sex", "toten")
cov_m2_extra  <- c("race", "educ", "hhinc", "smoke", "drink", "actc", "bmi_cat")

fit_ediea_models <- function(data,
                             time,   # 如 "fu_hbp"
                             event,  # 如 "e_hbp"
                             disease,               # 结局
                             cov_m3_comorb = NULL
){
  dat <- data %>%

  ediea_vec <- dat[["ediea"]]
  ediea_min <- quantile(ediea_vec, 0.10, na.rm = TRUE)
  ediea_max <- quantile(ediea_vec, 0.90, na.rm = TRUE)
  ediea_step20 <- (ediea_max - ediea_min)
  
  # Model 1
  cov_m1 <- cov_m1_base
  # Model 2
  cov_m2 <- c(cov_m1, cov_m2_extra)
  # Model 3
  cov_m3 <- c(cov_m2, cov_m3_comorb)
  
  f_m1 <- as.formula(
    paste0("Surv(", time, ", ", event, ") ~ ediea_group + ",
           paste(cov_m1, collapse = " + "))
  )
  f_m2 <- as.formula(
    paste0("Surv(", time, ", ", event, ") ~ ediea_group + ",
           paste(cov_m2, collapse = " + "))
  )
  f_m3 <- as.formula(
    paste0("Surv(", time, ", ", event, ") ~ ediea_group + ",
           paste(cov_m3, collapse = " + "))
  )
  
  models <- list(
    model1 = coxph(f_m1, data = dat),
    model2 = coxph(f_m2, data = dat),
    model3 = coxph(f_m3, data = dat)
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
  
  f_m1_cont <- as.formula(
    paste0("Surv(", time, ", ", event, ") ~ ediea + ",
           paste(cov_m1, collapse = " + "))
  )
  f_m2_cont <- as.formula(
    paste0("Surv(", time, ", ", event, ") ~ ediea + ",
           paste(cov_m2, collapse = " + "))
  )
  f_m3_cont <- as.formula(
    paste0("Surv(", time, ", ", event, ") ~ ediea + ",
           paste(cov_m3, collapse = " + "))
  )
  
  models_cont <- list(
    model1 = coxph(f_m1_cont, data = dat),
    model2 = coxph(f_m2_cont, data = dat),
    model3 = coxph(f_m3_cont, data = dat)
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

res_all <- bind_rows(
  # HBP
  fit_ediea_models(hrs_diseases, "fu_hbp", "e_hbp",
                   disease = "HBP", cov_m3_comorb = c("t2d", "cvd")),
  
  # Diabetes
  fit_ediea_models(hrs_diseases, "fu_diab", "e_diab",
                   disease = "Diabetes", cov_m3_comorb = c("hypten", "cvd")),
  
  # Lung disease
  fit_ediea_models(hrs_diseases, "fu_lung", "e_lung",
                   disease = "Lung disease", cov_m3_comorb = c("t2d", "hypten", "cvd")),
  
  # Heart disease
  fit_ediea_models(hrs_diseases, "fu_heart", "e_heart",
                   disease = "Heart disease", cov_m3_comorb = c("t2d", "hypten")),
  
  # Stroke
  fit_ediea_models(hrs_diseases, "fu_stroke", "e_stroke",
                   disease = "Stroke", cov_m3_comorb = c("t2d", "hypten", "cvd")),
  
  # Cancer
  fit_ediea_models(hrs_diseases, "fu_cancer", "e_cancer",
                   disease = "Cancer", cov_m3_comorb = c("t2d", "hypten", "cvd")),
  
  # Arthritis
  fit_ediea_models(hrs_diseases, "fu_arth", "e_arth",
                   disease = "Arthritis", cov_m3_comorb = c("t2d", "hypten", "cvd")),
  
  # Dementia
  fit_ediea_models(hrs_diseases, "fu_dem", "e_dem",
                   disease = "Dementia", cov_m3_comorb = c("t2d", "hypten", "cvd")),
  
  fit_ediea_models(hrs_diseases, "fu_death", "death",
                   disease = "Mortality", cov_m3_comorb = c("t2d", "hypten", "cvd")),
  
)

# ----------- ADL/IADL/CES-D/Chronic ----------
fit_irr <- function(df, outcome, offset_time = NULL, outcome_label, cov_m3 = NULL) {
  irr_base_cov_m1 <- c("age", "sex", "toten")
  irr_base_cov_m2 <- c("age", "sex", "toten", "race", "educ", "hhinc", "smoke", "drink", "actc", "bmi_cat")
  irr_base_cov_m3 <- c(irr_base_cov_m2, cov_m3)

  ediea_vec <- df[["ediea"]]
  ediea_min <- quantile(ediea_vec, 0.10, na.rm = TRUE)
  ediea_max <- quantile(ediea_vec, 0.90, na.rm = TRUE)
  ediea_step20 <- (ediea_max - ediea_min)
  
  form1 <- as.formula(
    paste0(outcome, " ~ ediea_group + ", paste(irr_base_cov_m1, collapse = " + "))
  )
  form2 <- as.formula(
    paste0(outcome, " ~ ediea_group + ", paste(irr_base_cov_m2, collapse = " + "))
  )
  form3 <- as.formula(
    paste0(outcome, " ~ ediea_group + ", paste(irr_base_cov_m3, collapse = " + "))
  )
  
  models <- list(
    model1 = glm(form1, family = poisson(link = "log"), offset = if (!is.null(offset_time)) log(df[[offset_time]]) else NULL, data = df),
    model2 = glm(form2, family = poisson(link = "log"), offset = if (!is.null(offset_time)) log(df[[offset_time]]) else NULL, data = df),
    model3 = glm(form3, family = poisson(link = "log"), offset = if (!is.null(offset_time)) log(df[[offset_time]]) else NULL, data = df)
  )
  
  res_group <- map_dfr(names(models), function(mname) {
    models[[mname]] %>%
      tidy() %>%
      filter(str_detect(term, "ediea")) %>%
      mutate(
        outcome = outcome_label,
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
  
  f_m1_cont <- as.formula(
    paste0(outcome, " ~ ediea + ", paste(irr_base_cov_m1, collapse = " + "))
  )
  f_m2_cont <- as.formula(
    paste0(outcome, " ~ ediea + ", paste(irr_base_cov_m2, collapse = " + "))
  )
  f_m3_cont <- as.formula(
    paste0(outcome, " ~ ediea + ", paste(irr_base_cov_m3, collapse = " + "))
  )
  
  models_cont <- list(
    model1 = glm(f_m1_cont, family = poisson(link = "log"), offset = if (!is.null(offset_time)) log(df[[offset_time]]) else NULL, data = df),
    model2 = glm(f_m2_cont, family = poisson(link = "log"), offset = if (!is.null(offset_time)) log(df[[offset_time]]) else NULL, data = df),
    model3 = glm(f_m3_cont, family = poisson(link = "log"), offset = if (!is.null(offset_time)) log(df[[offset_time]]) else NULL, data = df)
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
      outcome = outcome_label,
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

irr_adl <- fit_irr(hrs_scores, outcome = "adl5_inc15",
                   outcome_label = "ADL", cov_m3 = c("cesd12"))

irr_iadl <- fit_irr(hrs_scores, outcome = "iadl5_inc15",
                    outcome_label = "IADL", cov_m3 = c("iadl5_12"))

irr_cesd <- fit_irr(hrs_scores, outcome = "cesd_inc15",
                    outcome_label = "CES-D")

irr_chro <- fit_irr(hrs_diseases, outcome = 'inc_chronic_n', offset_time = "fu_chronic_year",
                    outcome_label = 'Chronic Diseases')

irr_all <- bind_rows(
  irr_adl,
  irr_iadl,
  # irr_cesd,
  irr_chro
)

## cases / person-years
cases_py_hrs <- bind_rows(
  # Heart disease
  hrs_diseases %>%
    filter(fu_heart > 0) %>%
    summarise(
      outcome = "Heart disease",
      cases   = sum(e_heart == 1, na.rm = TRUE),
      py      = sum(fu_heart, na.rm = TRUE)
    ),
  # Stroke
  hrs_diseases %>%
    filter(fu_stroke > 0) %>%
    summarise(
      outcome = "Stroke",
      cases   = sum(e_stroke == 1, na.rm = TRUE),
      py      = sum(fu_stroke, na.rm = TRUE)
    ),
  # Cancer
  hrs_diseases %>%
    filter(fu_cancer > 0) %>%
    summarise(
      outcome = "Cancer",
      cases   = sum(e_cancer == 1, na.rm = TRUE),
      py      = sum(fu_cancer, na.rm = TRUE)
    ),
  # Dementia
  hrs_diseases %>%
    filter(fu_dem > 0) %>%
    summarise(
      outcome = "Dementia",
      cases   = sum(e_dem == 1, na.rm = TRUE),
      py      = sum(fu_dem, na.rm = TRUE)
    ),
  # Mortality
  hrs_diseases %>%
    filter(fu_death > 0) %>%
    summarise(
      outcome = "Mortality",
      cases   = sum(death == 1, na.rm = TRUE),
      py      = sum(fu_death, na.rm = TRUE)
    )
) %>%
  mutate(
    case_py_label = paste0(cases, " / ", round(py, 0))
  )

df_plot_hr_irr_hrs <- bind_rows(df_plot, df_plot_irr) %>%
  filter(ediea_label == "90th-10th percentile difference") %>%
  left_join(
    cases_py_hrs %>% select(outcome, case_py_label),
    by = "outcome"
  ) %>%
  mutate(outcome = factor(outcome, levels = unique(outcome)))

# ------------------------------ META, HRS & ELSA --------------------------------
meta_hrs_elsa <- bind_rows(
  df_plot_hr_irr_hrs %>% 
    filter(ediea_level == "per20pct") %>% 
    mutate(Cohort = "HRS"),
  df_plot_hr_irr_elsa %>% 
    filter(ediea_level == "per20pct") %>% 
    mutate(Cohort = "ELSA"),
  df_plot_hr_irr_nhanes %>% 
    filter(ediea_level == "per20pct") %>% 
    mutate(Cohort = "NHANES")
) %>% 
  select(outcome, Cohort, coef, se.coef., ediea_level, case_py_label)

meta_pooled_all <- meta_hrs_elsa %>% 
  group_by(outcome) %>% 
  group_modify(~{
    m <- metagen(TE    = .x$coef, seTE  = .x$se.coef.,
      studlab = .x$Cohort, common  = FALSE, random  = TRUE
    )
    bind_rows(
      .x,
      tibble(
        outcome     = .x$outcome[1],
        Cohort      = "Pooled",
        coef        = m$TE.random,
        se.coef.    = m$seTE.random,
        ediea_level = "per20pct"
      )
    )
  }) %>% 
  ungroup() %>% 
  mutate(
    HR      = exp(coef),
    LCI  = exp(coef - 1.96 * se.coef.),
    UCI = exp(coef + 1.96 * se.coef.),
    HR_text = sprintf("%.2f (%.2f-%.2f)", HR, LCI, UCI),
    p       = 2 * pnorm(-abs(coef / se.coef.)),
    is_sig = p<0.05,
    
    Cohort = factor(Cohort, levels = c("Pooled", "NHANES", "ELSA", "HRS"))
  )


library(tidytext)

x_cap_min  <- 0.1
x_cap_max  <- 1.6
x_axis_min <- 0.1
x_axis_max <- 1.7
cohort_order <- c("HRS","NHANES","ELSA","Pooled")
df <- meta_pooled_all %>%
  filter(!is.na(HR), !is.na(LCI), !is.na(UCI)) %>%
  mutate(
    Cohort  = trimws(as.character(Cohort)),
    outcome = forcats::fct_inorder(outcome)
  ) %>%
  group_by(outcome) %>%
  mutate(
    rank_index = match(Cohort, rev(cohort_order)),
    Cohort_plot = tidytext::reorder_within(Cohort, rank_index, outcome),
    y = dense_rank(rank_index),
    
    LCI_cap = pmax(LCI, x_cap_min),
    UCI_cap = pmin(UCI, x_cap_max),
    trunc_l = LCI < x_cap_min,
    trunc_u = UCI > x_cap_max,
    is_sig  = as.logical(is_sig)
  ) %>%
  ungroup()


bg <- df %>% distinct(outcome, y) %>% mutate(is_even = y %% 2 == 0, .groups = "drop")
hdr <- df %>%
  group_by(outcome) %>%
  summarise(y_hdr = max(y) + 0.65, .groups = "drop") %>%
  slice_head(n=2)

p_forest_meta <- ggplot(df) +
  facet_wrap(~ outcome, nrow = 4, ncol = 2, scales = "free_y") + 
  geom_rect(
    data = bg %>% filter(is_even),
    aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
    inherit.aes = FALSE, fill = "grey97", colour = NA
  ) +
  geom_vline(xintercept = seq(x_cap_min, x_cap_max, by = 0.2), linewidth = 0.35, colour = "grey90") +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.7, colour = "grey35") +
  geom_segment(
    aes(x = LCI_cap, xend = UCI_cap, y = Cohort_plot, yend = Cohort_plot),
    linewidth = 0.75, colour = "grey25"
  ) +
  geom_segment(
    data = df %>% filter(trunc_u),
    aes(x = x_cap_max - 0.02, xend = x_cap_max, y = Cohort_plot, yend = Cohort_plot),
    inherit.aes = FALSE, linewidth = 0.75, colour = "grey25",
    arrow = arrow(type = "closed", length = unit(0.10, "in"))
  ) +
  geom_segment(
    data = df %>% filter(trunc_l),
    aes(x = x_cap_min + 0.02, xend = x_cap_min, y = Cohort_plot, yend = Cohort_plot),
    inherit.aes = FALSE, linewidth = 0.75, colour = "grey25",
    arrow = arrow(type = "closed", length = unit(0.10, "in"))
  ) +
  geom_point(
    data = df %>% filter(Cohort != "Pooled"),
    aes(x = HR, y = Cohort_plot, fill = is_sig),
    shape = 21, size = 2.6, stroke = 0.75, colour = "black"
  ) +
  geom_point(
    data = df %>% filter(Cohort == "Pooled"),
    aes(x = HR, y = Cohort_plot, fill = is_sig),
    shape = 23, size = 3.2, stroke = 0.85, colour = "black"
  ) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#2F6FAE")) +
  geom_text(
    aes(x = x_axis_min + 0.01, y = Cohort_plot, label = case_py_label),
    hjust = 0, size = 3.0, colour = "grey10", vjust = -1
  ) +
  geom_text(
    aes(x = x_axis_max - 0.01, y = Cohort_plot, label = HR_text),
    hjust = 1, size = 3.0, colour = "grey10", vjust = -1
  ) +
  geom_text(
    data = hdr,
    aes(x = x_axis_min + 0.01, y = y_hdr, label = "Cases / person-years"),
    inherit.aes = FALSE, hjust = 0, size = 3.1, fontface = "bold", colour = "grey10",
    vjust = -1.0
  ) +
  geom_text(
    data = hdr,
    aes(x = x_axis_max - 0.01, y = y_hdr, label = "HR (95% CI)"),
    inherit.aes = FALSE, hjust = 1, size = 3.1, fontface = "bold", colour = "grey10",
    vjust = -1.0
  ) +
  scale_x_continuous(limits = c(x_axis_min, x_axis_max)) +
  tidytext::scale_y_reordered() + 
  labs(
    x = "HR/IRR (90th vs 10th percentile, 95% CI)",
    y = NULL, title = "c. Incidence of age−related health outcomes"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    strip.background = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    axis.ticks.length = unit(2.5, "mm"),
    panel.spacing = unit(1.0, "lines"),
    plot.margin = margin(8, 18, 8, 8)
  )

p_forest_meta

# Abstract Graph
df_pooled <- meta_pooled_all %>%
  filter(Cohort == "Pooled") %>%
  filter(!is.na(HR), !is.na(LCI), !is.na(UCI)) %>%
  mutate(
    outcome = fct_rev(fct_inorder(outcome)),
    LCI_cap = pmax(LCI, x_cap_min),
    UCI_cap = pmin(UCI, x_cap_max),
    trunc_l = LCI < x_cap_min,
    trunc_u = UCI > x_cap_max
  )

p_pooled_only <- ggplot(df_pooled, aes(y = outcome)) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.6, colour = "grey35") +
  geom_segment(aes(x = LCI_cap, xend = UCI_cap, yend = outcome),
               linewidth = 0.7, colour = "grey25") +
  geom_segment(
    data = df_pooled %>% filter(trunc_u),
    aes(x = x_cap_max - 0.02, xend = x_cap_max, yend = outcome),
    inherit.aes = FALSE, linewidth = 0.7, colour = "grey25",
    arrow = arrow(type = "closed", length = unit(0.08, "in"))
  ) +
  geom_segment(
    data = df_pooled %>% filter(trunc_l),
    aes(x = x_cap_min + 0.02, xend = x_cap_min, yend = outcome),
    inherit.aes = FALSE, linewidth = 0.7, colour = "grey25",
    arrow = arrow(type = "closed", length = unit(0.08, "in"))
  ) +
  geom_point(aes(x = HR, fill = is_sig),
             shape = 23, size = 2.6, stroke = 0.7, colour = "black") +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#2F6FAE")) +
  geom_text(aes(x = x_axis_max, label = HR_text),
            hjust = 1, size = 2.6, colour = "grey10", nudge_y = 0.18) +
  scale_x_continuous(limits = c(x_axis_min, x_axis_max)) +
  labs(x = "HR (95% CI)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(face = "bold", size = 9),
    axis.text.y  = element_text(size = 8, colour = "black"),
    axis.text.x  = element_text(size = 8, colour = "black"),
    axis.ticks.length = unit(2.0, "mm"),
    plot.margin = margin(4, 28, 4, 4)
  )
