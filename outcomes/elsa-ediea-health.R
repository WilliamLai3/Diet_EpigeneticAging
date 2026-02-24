# EDISEA and outcomes in ELSA
library(tidyverse)
library(survival)
library(mice)
library(ggjournals)
library(readr)
library(haven)
library(patchwork)
library(haven)
library(mice)
library(survminer)
source("./plot_rcs_linear.R")
source("./plot_rcs_pr_rms.R")

DEM_THRESHOLD <- 1.5

yn01 <- function(x) {
  x <- as.numeric(x)
  ifelse(x == 1, 1,
         ifelse(x %in% c(0, 2, 3, 4, 5, 7, 8, 9), 0, NA_real_))
}

rightwrong01 <- function(x) {
  ans <- ifelse(x %in% c(-1,-8,-9), NA_integer_, as.numeric(x))
  return (ans)
}

year_month2date <- function(year, month) {
  date = if_else(
    !is.na(year) & !is.na(month),
    as.Date(
      sprintf("%04d-%02d-15", year, month),
      format = "%Y-%m-%d"
    ),
    as.Date(NA_character_)
  )
  return (date)
}

death <- haven::read_sav('../ELSA-5050/spss/spss28/index_file_wave_0-wave_5_v2.sav')

load('./dat/elsa_ediea_from_hrs.rda')
h_elsa <- haven::read_sav('../ELSA-5050/spss/spss28/h_elsa_g3.sav')
# h_elsa <- haven::read_dta('../ELSA-2/elsa_main/stata/stata13_se/h_elsa_g2.dta')

wave_10 <- haven::read_sav('../ELSA-5050/spss/spss28/wave_10_elsa_data_eul_v4.sav')
wave_11 <- haven::read_sav('../ELSA-5050/spss/spss28/wave_11_elsa_data_eul_v1.sav')

# -------------------------- Population Cognitive SD ------------------------
waves <- 1:9
get_var <- function(stub, w, default = NA_real_) {
  v <- sprintf("r%d%s", w, stub)
  n <- nrow(h_elsa)
  if (v %in% names(h_elsa)) {
    h_elsa[[v]]
  } else {
    rep(default, n)
  }
}
## long table
make_wave <- function(w) {
  data.frame(
    IDAUNIQ  = h_elsa$idauniq,
    wave     = w,
    age      = h_elsa[[sprintf("r%dagey",   w)]],
    edu      = h_elsa$raeduc_e,  # 1–3：教育层
    orient   = h_elsa[[sprintf("r%dorient", w)]],   # 定向力, 0-4
    memory   = h_elsa[[sprintf("r%dtr20",   w)]],   # 总记忆（立即+延迟）, 0-20
    fluency  = get_var("verbf", w),   # 语义流畅, 0-56
    numeracy = get_var("numer_e", w),   # 计算/数理, 0-6
    adl      = get_var("adltot6", w),    # 6 项 ADL 总分, 0-6
    iadl     = get_var("iadltot1_e", w)
  )
}

cog_long <- map_dfr(waves, make_wave) %>%
  mutate(
    edu = as.numeric(edu),
    age_band = cut(
      age, breaks = c(50, 60, 70, 81),
      include.lowest = TRUE, right = FALSE,
      labels = c("50-59", "60-69", "70-80")
    )
  )

norms <- cog_long %>%
  filter(age >= 50, age <= 80, !is.na(edu)) %>%
  group_by(edu) %>%
  summarise(
    mean_orient   = mean(orient,   na.rm = TRUE),
    sd_orient     = sd(orient,     na.rm = TRUE),
    mean_memory   = mean(memory,   na.rm = TRUE),
    sd_memory     = sd(memory,     na.rm = TRUE),
    mean_fluency  = mean(fluency,  na.rm = TRUE),
    sd_fluency    = sd(fluency,    na.rm = TRUE),
    mean_numeracy = mean(numeracy, na.rm = TRUE),
    sd_numeracy   = sd(numeracy,   na.rm = TRUE),
    .groups = "drop"
  )

cog_long <- cog_long %>%
  left_join(norms, by = c("edu")) %>%
  mutate(
    z_orient   = (orient   - mean_orient)   / sd_orient,
    z_memory   = (memory   - mean_memory)   / sd_memory,
    z_fluency  = (fluency  - mean_fluency)  / sd_fluency,
    z_numeracy = (numeracy - mean_numeracy) / sd_numeracy,
    
    imp_orient   = !is.na(z_orient)   & z_orient   <= -DEM_THRESHOLD,
    imp_memory   = !is.na(z_memory)   & z_memory   <= -DEM_THRESHOLD,
    imp_fluency  = !is.na(z_fluency)  & z_fluency  <= -DEM_THRESHOLD,
    imp_numeracy = !is.na(z_numeracy) & z_numeracy <= -DEM_THRESHOLD,

    n_imp   = imp_orient + imp_memory + imp_fluency + imp_numeracy,
    cog_imp = n_imp >= 2,

    adl_imp = !is.na(adl) & adl >= 1,
    iadl_imp = !is.na(iadl) & iadl >= 1,
    
    alg_dem_wave = as.integer(cog_imp & (adl_imp|iadl_imp))
  )

alg_w9 <- cog_long %>%
  filter(wave == 9) %>%
  transmute(
    idauniq = IDAUNIQ,
    alg_dem9 = alg_dem_wave
  )

# -------------------------- Readin Wave10/11 -------------------------------
wave_10 <- wave_10 %>%
  transmute(
    idauniq,
    age10 = indager,
    age10_cat = cut(
      age10,
      breaks = c(50, 60, 70, 81),
      right  = FALSE,
      include.lowest = TRUE,
      labels = c("50-59", "60-69", "70-80")
    ),

    r10iwindy_num = as.numeric(iintdaty),
    r10iwindm_num = as.numeric(iintdatm),
    date10 = year_month2date(r10iwindy_num, r10iwindm_num),
    
    hbp10 = yn01(heeverbp),
    t2d10    = yn01(heacd),
    stroke10 = yn01(heeverst),
    cancer10 = yn01(heeverca),
    lung10   = yn01(heevercl),
    arth10   = yn01(heeverar),
    dem10    = yn01(heeverdm),
    
    heart10_raw = pmax(
      yn01(heeveran),  # angina
      yn01(heevermi),  # heart attack
      yn01(heeverhf),  # heart failure
      yn01(heeveroh),  # other heart condition
      na.rm = TRUE
    ),
    heart10 = ifelse(is.infinite(heart10_raw), NA_integer_, heart10_raw),
    
    ## Algorithm Dementia
    # Memory
    memory_score10 = rightwrong01(cflisen)+rightwrong01(cflisd), # 0-20
    # Fluency
    fluency_score10 = rightwrong01(cfani),
    # Executive
    serial7_correct10 = 
      as.numeric((rightwrong01(cfsva) == 93) +
      (rightwrong01(cfsvb) == 86) +
      (rightwrong01(cfsvc) == 79) +
      (rightwrong01(cfsvd) == 72) +
      (rightwrong01(cfsve) == 65)),
    countdown_correct10 = case_when(
      rightwrong01(cfc20frst) == 1 ~ 1,
      rightwrong01(cfc20fscnd) == 1 ~ 1,
      TRUE ~ 0
    ),
    exec_score10 = serial7_correct10 + countdown_correct10,
    # Orientation
    orient_score10 = yn01(cfdatd) + yn01(cfdatm) + 
      yn01(cfdaty) + yn01(cfday),
    # ADL
    adl_score10 = rightwrong01(headldr) + rightwrong01(headlwa) + 
      rightwrong01(headlba) + rightwrong01(headlea) +
      rightwrong01(headlbe) + rightwrong01(headlwc),
    # IADL
    iadl_score10 = rightwrong01(headlma) + rightwrong01(headlda) +
      rightwrong01(headlpr) + rightwrong01(headlsh) +
      rightwrong01(headlph) + rightwrong01(headlsp) +
      rightwrong01(headlme) + rightwrong01(headlho) +
      rightwrong01(headlmo)
    
  ) %>%
  select(-heart10_raw)

wave_11 <- wave_11 %>%
  transmute(
    idauniq,
    age11 = indager,
    age11_cat = cut(
      age11,
      breaks = c(50, 60, 70, 81),
      right  = FALSE,
      include.lowest = TRUE,
      labels = c("50-59", "60-69", "70-80")
    ),
    
    r11iwindy_num = as.numeric(iintdaty),
    r11iwindm_num = as.numeric(iintdatm),
    date11 = year_month2date(r11iwindy_num, r11iwindm_num),
    
    hbp11 = yn01(hediia1),
    t2d11    = yn01(hedids12),
    stroke11 = yn01(hediia7),
    cancer11 = yn01(hedids5),
    lung11   = yn01(hedids1),
    arth11   = yn01(hedids3),
    dem11    = (yn01(hedids9)|yn01(hedids8)) %>% as.numeric(),
    
    ## angina / heart attack / heart failure / other heart condition
    heart11_raw = pmax(
      yn01(hediia2),  # angina
      yn01(hediia3),  # heart attack
      yn01(hediia4),  # heart failure
      yn01(hediia5),  # heart murmur
      yn01(hediia6),  # abnormal heart rhythm
      yn01(hediia95), # other heart condition
      na.rm = TRUE
    ),
    heart11 = ifelse(is.infinite(heart11_raw), NA_integer_, heart11_raw),
    
    ## Algorithm Dementia
    # Memory
    memory_score11 = rightwrong01(cflisen)+rightwrong01(cflisd),
    # Fluency
    fluency_score11 = rightwrong01(cfani),
    # Executive
    serial7_correct11 = 
      as.numeric((rightwrong01(cfsva) == 93) +
                   (rightwrong01(cfsvb) == 86) +
                   (rightwrong01(cfsvc) == 79) +
                   (rightwrong01(cfsvd) == 72) +
                   (rightwrong01(cfsve) == 65)),
    countdown_correct11 = case_when(
      rightwrong01(cfc20frst) == 1 ~ 1,
      rightwrong01(cfc20fscnd) == 1 ~ 1,
      TRUE ~ 0
    ),
    exec_score11 = serial7_correct11 + countdown_correct11,
    # Orientation
    orient_score11 = yn01(cfdatd) + yn01(cfdatm) + 
      yn01(cfdaty) + yn01(cfday),
    # ADL
    adl_score11 = rightwrong01(headldr) + rightwrong01(headlwa) + 
      rightwrong01(headlba) + rightwrong01(headlea) +
      rightwrong01(headlbe) + rightwrong01(headlwc),
    # IADL
    iadl_score11 = rightwrong01(headlma) + rightwrong01(headlda) +
      rightwrong01(headlpr) + rightwrong01(headlsh) +
      rightwrong01(headlph) + rightwrong01(headlsp) +
      rightwrong01(headlme) + rightwrong01(headlho) +
      rightwrong01(headlmo)
  ) %>%
  select(-heart11_raw)

# ---------------------------- Define incident outcome --------------------------
elsa_cov <- h_elsa %>%
  transmute(

    idauniq,

    r9iwindy_num = as.numeric(r9iwindy),
    r9iwindm_num = as.numeric(r9iwindm),
    date9 = year_month2date(r9iwindy_num, r9iwindm_num),

    age    = as.numeric(r9agey),
    age9_cat = cut(
      age,
      breaks = c(50, 60, 70, 81),
      right  = FALSE,
      include.lowest = TRUE,
      labels = c("50-59", "60-69", "70-80")
    ),
    sex    = as.numeric(ragender),
    race   = as.numeric(raracem),
    educ   = as.numeric(raeduc_e),
    hhinc  = case_when(
      !is.na(h9itot) & h9itot < 50000 ~ '<50000',
      !is.na(h9itot) & h9itot >= 50000 ~ '>=50000',
      TRUE ~ NA_character_
      ) %>% factor(),

    smoke = case_when(
      r9smoken == 1                      ~ 2L,        # current smoker
      r9smokev == 1 & r9smoken != 1      ~ 1L,        # former
      r9smokev == 0                      ~ 0L,        # never
      TRUE                               ~ NA_integer_
    ) %>% factor(),
    drink = as.numeric(r9drink),
    actc  = case_when(
      r9vgactx_e %in% c(4,5) ~ 'Low',
      r9vgactx_e %in% c(3) ~ 'Medium',
      r9vgactx_e %in% c(2) ~ 'High',
      TRUE ~ NA_character_
    ) %>% factor(),
    bmi_cat  = as.numeric(r8mbmicat),
    bmi_cat = case_when(
      bmi_cat %in% c(1,2) ~ "Normal",
      bmi_cat %in% c(3) ~ "Overweight",
      bmi_cat %in% c(4,5,6) ~ "Obesity",
      TRUE ~ NA_character_
    ) %>% factor(),

    cesd = r9cesd,

    heart     = yn01(r9hearte),
    stroke    = yn01(r9stroke),
    cancer    = yn01(r9cancre),
    hbp = yn01(r9hibpe),
    t2d = yn01(r9diabe),
    lung = yn01(r9lunge),
    arth = yn01(r9arthre),
    dementia = yn01(r9demene),
    cvd = case_when(
      heart == 1L | stroke == 1L ~ 1L,
      heart == 0L & stroke == 0L ~ 0L,
      TRUE                       ~ NA_integer_
    ),
    
    adl_score9 = r9adltot6,
    iadl_score9 = r9iadltot1_e
  ) %>%
  left_join(wave_10, by = 'idauniq') %>%
  left_join(wave_11, by = 'idauniq') %>%
  mutate(
    inc_hbp = ifelse(hbp==0&(hbp10==1|hbp11==1), 1, 0),
    inc_heart = ifelse(heart==0&(heart10==1|heart11==1), 1, 0),
    inc_stroke = ifelse(stroke==0&(stroke10==1|stroke11==1), 1, 0),
    inc_cancer = ifelse(cancer==0&(cancer10==1|cancer11==1), 1, 0),
    inc_t2d = ifelse(t2d==0&(t2d10==1|t2d11==1), 1, 0),
    inc_lung = ifelse(lung==0&(lung10==1|lung11==1), 1, 0),
    inc_arth = ifelse(arth==0&(arth10==1|arth11==1), 1, 0),
    inc_dem = ifelse(dementia==0&(dem10==1|dem11==1), 1, 0),
  ) %>%
  mutate(
    last_date = dplyr::coalesce(date11, date10),

    hbp_event_date = case_when(
      hbp == 0 & hbp10 == 1                 ~ date9 + (date10 - date9) / 2,
      hbp == 0 & (is.na(hbp10) | hbp10 == 0) & hbp11 == 1 ~ date10 + (date11 - date10) / 2,
      hbp == 0 & inc_hbp == 0              ~ last_date,   # 截尾
      TRUE                                 ~ as.Date(NA_character_)
    ),
    fu_hbp = as.numeric(hbp_event_date - date9) / 365.25,

    heart_event_date = case_when(
      heart == 0 & heart10 == 1 ~ date9 + (date10 - date9) / 2,
      heart == 0 & (is.na(heart10) | heart10 == 0) & heart11 == 1 ~ date10 + (date11 - date10) / 2,
      heart == 0 & inc_heart == 0 ~ last_date,
      TRUE ~ as.Date(NA_character_)
    ),
    fu_heart = as.numeric(heart_event_date - date9) / 365.25,

    stroke_event_date = case_when(
      stroke == 0 & stroke10 == 1 ~ date9 + (date10 - date9) / 2,
      stroke == 0 & (is.na(stroke10) | stroke10 == 0) & stroke11 == 1 ~ date10 + (date11 - date10) / 2,
      stroke == 0 & inc_stroke == 0 ~ last_date,
      TRUE ~ as.Date(NA_character_)
    ),
    fu_stroke = as.numeric(stroke_event_date - date9) / 365.25,
    
    cancer_event_date = case_when(
      cancer == 0 & cancer10 == 1 ~ date9 + (date10 - date9) / 2,
      cancer == 0 & (is.na(cancer10) | cancer10 == 0) & cancer11 == 1 ~ date10 + (date11 - date10) / 2,
      cancer == 0 & inc_cancer == 0 ~ last_date,
      TRUE ~ as.Date(NA_character_)
    ),
    fu_cancer = as.numeric(cancer_event_date - date9) / 365.25,

    t2d_event_date = case_when(
      t2d == 0 & t2d10 == 1 ~ date9 + (date10 - date9) / 2,
      t2d == 0 & (is.na(t2d10) | t2d10 == 0) & t2d11 == 1 ~ date10 + (date11 - date10) / 2,
      t2d == 0 & inc_t2d == 0 ~ last_date,
      TRUE ~ as.Date(NA_character_)
    ),
    fu_t2d = as.numeric(t2d_event_date - date9) / 365.25,

    lung_event_date = case_when(
      lung == 0 & lung10 == 1 ~ date9 + (date10 - date9) / 2,
      lung == 0 & (is.na(lung10) | lung10 == 0) & lung11 == 1 ~ date10 + (date11 - date10) / 2,
      lung == 0 & inc_lung == 0 ~ last_date,
      TRUE ~ as.Date(NA_character_)
    ),
    fu_lung = as.numeric(lung_event_date - date9) / 365.25,

    arth_event_date = case_when(
      arth == 0 & arth10 == 1 ~ date9 + (date10 - date9) / 2,
      arth == 0 & (is.na(arth10) | arth10 == 0) & arth11 == 1 ~ date10 + (date11 - date10) / 2,
      arth == 0 & inc_arth == 0 ~ last_date,
      TRUE ~ as.Date(NA_character_)
    ),
    fu_arth = as.numeric(arth_event_date - date9) / 365.25,
    
    dem_event_date = case_when(
      dementia == 0 & dem10 == 1 ~ date9 + (date10 - date9) / 2,
      dementia == 0 & (is.na(dem10) | dem10 == 0) & dem11 == 1 ~ date10 + (date11 - date10) / 2,
      dementia == 0 & inc_dem == 0 ~ last_date,
      TRUE ~ as.Date(NA_character_)
    ),
    fu_dem = as.numeric(dem_event_date - date9) / 365.25
    
  ) %>%
  left_join(alg_w9, by = 'idauniq') %>%
  mutate(educ_num = as.numeric(educ)) %>%
  left_join(norms, by = c("educ_num" = "edu")) %>%
  mutate(
    educ = case_when(
      educ %in% c(1, 3, 4) ~ "Below College",
      educ %in% c(5) ~ "College and above",
      TRUE ~ NA_character_
    ) %>% factor(),

    z_mem10 = (memory_score10  - mean_memory)  / sd_memory,
    z_flu10 = (fluency_score10 - mean_fluency) / sd_fluency,
    z_ori10 = (orient_score10  - mean_orient)  / sd_orient,
    z_numeracy10 = (exec_score10 - mean_numeracy) / sd_numeracy,
    
    imp_mem10 = !is.na(z_mem10) & z_mem10 <= -DEM_THRESHOLD,
    imp_flu10 = !is.na(z_flu10) & z_flu10 <= -DEM_THRESHOLD,
    imp_ori10 = !is.na(z_ori10) & z_ori10 <= -DEM_THRESHOLD,
    imp_num10 = !is.na(z_numeracy10) & z_numeracy10 <= -DEM_THRESHOLD,
    
    n_imp10   = imp_mem10 + imp_flu10 + imp_ori10 + imp_num10,
    cog_imp10 = n_imp10 >= 2,
    adl_imp10 = !is.na(adl_score10) & adl_score10 >= 1,
    iadl_imp10 = !is.na(iadl_score10) & iadl_score10 >= 1,
    
    alg_dem10 = as.integer(cog_imp10 & (adl_imp10|iadl_imp10)),

    z_mem11 = (memory_score11  - mean_memory)  / sd_memory,
    z_flu11 = (fluency_score11 - mean_fluency) / sd_fluency,
    z_ori11 = (orient_score11  - mean_orient)  / sd_orient,
    z_numeracy11 = (exec_score11 - mean_numeracy) / sd_numeracy,
    
    imp_mem11 = !is.na(z_mem11) & z_mem11 <= -DEM_THRESHOLD,
    imp_flu11 = !is.na(z_flu11) & z_flu11 <= -DEM_THRESHOLD,
    imp_ori11 = !is.na(z_ori11) & z_ori11 <= -DEM_THRESHOLD,
    imp_num11 = !is.na(z_numeracy11) & z_numeracy11 <= -DEM_THRESHOLD,
    
    n_imp11   = imp_mem11 + imp_flu11 + imp_ori11 + imp_num11,
    cog_imp11 = n_imp11 >= 2,
    adl_imp11 = !is.na(adl_score11) & adl_score11 >= 1,
    iadl_imp11 = !is.na(iadl_score11) & iadl_score11 >= 1,
    
    alg_dem11 = as.integer(cog_imp11 & (adl_imp11|iadl_imp11)),
    

    inc_alg_dem = ifelse(alg_dem9 == 0L & (alg_dem10 == 1L | alg_dem11 == 1L), 1L, 0L),
    
    alg_dem_event_date = case_when(
      alg_dem9 == 0L & alg_dem10 == 1L ~ date9 + (date10 - date9) / 2,
      alg_dem9 == 0L & (is.na(alg_dem10) | alg_dem10 == 0L) & alg_dem11 == 1L ~ date10 + (date11 - date10) / 2,
      alg_dem9 == 0L & inc_alg_dem == 0L ~ last_date,
      TRUE ~ as.Date(NA_character_)
    ),
    fu_alg_dem = as.numeric(alg_dem_event_date - date9) / 365.25,
    
    inc_all_dem = case_when(
      inc_dem == 1L | inc_alg_dem == 1L ~ 1L,
      inc_dem == 0L & inc_alg_dem == 0L ~ 0L,
      TRUE                                ~ NA_integer_
    ),
    all_dem_event_date = case_when(
      inc_all_dem == 1L ~ pmin(
        if_else(inc_dem   == 1L, dem_event_date,   as.Date(NA_character_)),
        if_else(inc_alg_dem == 1L, alg_dem_event_date, as.Date(NA_character_)),
        na.rm = TRUE
      ),
      inc_all_dem == 0L ~ last_date,
      TRUE              ~ as.Date(NA_character_)
    ),
    fu_all_dem = as.numeric(all_dem_event_date - date9) / 365.25
    
  ) %>%
  mutate(
    adl_chg = ifelse(!is.na(adl_score11), adl_score11-adl_score9, 
                     ifelse(!is.na(adl_score10), adl_score10-adl_score9, 0)),
    adl_chg = pmax(adl_chg, 0),
    iadl_chg = ifelse(!is.na(iadl_score11), iadl_score11-iadl_score9, 
                     ifelse(!is.na(iadl_score10), iadl_score10-iadl_score9, 0)),
    iadl_chg = pmax(iadl_chg, 0),
    
    inc_chronic = inc_hbp + inc_heart + inc_stroke + inc_cancer + inc_t2d + inc_lung +
      inc_arth,
    fu_chronic = as.numeric(as.Date(last_date) - date9) / 365.25
  )

covs = c("age", "sex", "race", "educ", "hhinc",
         "smoke", "drink", "actc", "bmi_cat", "energy_kcal", "t2d", "hbp", "cvd",
         "adl_score9", "iadl_score9")
elsa_dat <- elsa_fg40 %>%
  left_join(elsa_cov, by = 'idauniq') %>%
  mutate(
    ediea_t = cut_number(ediea, 3, labels = FALSE) %>% factor(),

    ediea_group = case_when(
      ediea_t==1 ~ 'low',
      ediea_t==2 ~ 'medium',
      ediea_t==3 ~ 'high',
      TRUE ~ NA_character_
    ),
    ediea_group = factor(ediea_group, levels = c("low","medium","high")),
  )
covs_for_mice <- elsa_dat %>%
  dplyr::select(all_of(covs)) %>%
  mutate(
    across(
      where(haven::is.labelled),
      ~ as.numeric(.)
    )
  )
imp <- mice(covs_for_mice, seed = 2025)
elsa_dat[, covs] <- complete(imp)

# --------------------------------- HR -----------------------------------------
cov_m1_base   <- c("age", "sex", "energy_kcal")
cov_m2_extra  <- c("race", "educ", "hhinc", "smoke", "drink", "actc", "bmi_cat")

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

cph_dem <- elsa_dat %>%
  filter(fu_alg_dem>0) %>%
  group_by(ediea_group) %>% 
  mutate(fol = fu_alg_dem) %>% 
  summarise(case = sum(inc_alg_dem),
            py = sum(fol))

agebase <- 60
cov_set60 <- elsa_dat %>%
  mutate(
    age_start_60 = ifelse(age<agebase, agebase, age),
    surv_time_60 = ifelse(age<agebase, fu_alg_dem + age - agebase, fu_alg_dem),
    age_end_60 = age_start_60 + surv_time_60,
  ) %>%
  filter(age_end_60>=agebase)

survf <- survfit(Surv(age, age+fu_alg_dem, inc_alg_dem) ~ ediea_group, data = elsa_dat %>% filter(fu_alg_dem>0))
ggsurvplot(survf, fun = "cumhaz",
           xlim = c(45,100),
           ylim = c(0.0,0.5))

res_all <- bind_rows(
  # Heart disease
  fit_ediea_models(elsa_dat, "fu_heart", "inc_heart",
                   disease = "Heart disease", cov_m3_comorb = c("t2d", "hbp")),
  
  # Stroke
  fit_ediea_models(elsa_dat, "fu_stroke", "inc_stroke",
                   disease = "Stroke", cov_m3_comorb = c("t2d", "hbp")),
  
  # Cancer
  fit_ediea_models(elsa_dat, "fu_cancer", "inc_cancer",
                   disease = "Cancer", cov_m3_comorb = c("t2d", "hbp", "cvd")),

  # Dementia
  fit_ediea_models(cov_set60, "surv_time_60", "inc_alg_dem",
                   disease = "Dementia", cov_m3_comorb = c("t2d", "hbp", "cvd"))
)

hr_dem <- coxph(Surv(age_start_60, age_end_60, inc_alg_dem) ~ ediea_group + sex + energy_kcal + race +
                  educ + hhinc + smoke + drink + actc + bmi_cat,
                data = cov_set60 %>% filter(surv_time_60>0))
summary(hr_dem)

# ----------- ADL/IADL/CES-D/Chronic ----------
fit_irr <- function(df, outcome, offset_time = NULL, outcome_label, cov_m3 = NULL) {
  irr_base_cov_m1 <- c("age", "sex", "energy_kcal")
  irr_base_cov_m2 <- c("age", "sex", "energy_kcal", "race", "educ", "hhinc", "smoke", "drink", "actc", "bmi_cat")
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

irr_adl <- fit_irr(elsa_dat, outcome = "adl_chg",
                   outcome_label = "ADL", cov_m3 = c("adl_score9"))

irr_iadl <- fit_irr(elsa_dat, outcome = "iadl_chg",
                    outcome_label = "IADL", cov_m3 = c("iadl_score9"))

irr_chro <- fit_irr(elsa_dat, outcome = 'inc_chronic', offset_time = "fu_chronic",
                    outcome_label = 'Chronic Diseases')

irr_all <- bind_rows(
  irr_adl,
  irr_iadl,
  irr_chro
)

df_plot_irr <- irr_all %>%
  filter(model == "model3") %>%
  mutate(
    ediea_label = case_when(
      ediea_level == "medium"   ~ "T2",
      ediea_level == "high"     ~ "T3",
      ediea_level == "per20pct" ~ "90th-10th percentile difference",
      TRUE                      ~ NA_character_
    )
  ) %>%
  select(outcome, ediea_level, ediea_label,
         HR, LCI, UCI, HR_label, p.value, coef, se.coef.)

df_plot_irr <- df_plot_irr %>%
  mutate(ediea_label = factor(ediea_label,
                              levels = c("T1 (Ref.)",
                                         "T2", "T3",
                                         "90th-10th percentile difference"))) %>%
  complete(outcome,
           ediea_label,
           fill = list(HR = NA_real_, LCI = NA_real_, UCI = NA_real_,
                       HR_label = NA_character_, p.value = NA_real_)) %>%
  mutate(
    is_ref = ediea_label == "T1 (Ref.)",
    HR_text = case_when(
      is_ref              ~ "Ref.",
      !is.na(HR_label)    ~ HR_label,
      TRUE                ~ ""
    ),
    P_text = case_when(
      is_ref              ~ "",
      is.na(p.value)      ~ "",
      p.value < 0.001     ~ "<0.001",
      TRUE                ~ sprintf("%.3f", p.value)
    ),
    ediea_label = fct_rev(ediea_label)
  ) %>%
  mutate(
    is_sig = p.value < 0.05
  )

## cases / person-years
cases_py_elsa <- bind_rows(
  # Heart disease
  elsa_dat %>%
    filter(fu_heart > 0) %>%
    summarise(
      outcome = "Heart disease",
      cases   = sum(inc_heart == 1, na.rm = TRUE),
      py      = sum(fu_heart, na.rm = TRUE)
    ),
  # Stroke
  elsa_dat %>%
    filter(fu_stroke > 0) %>%
    summarise(
      outcome = "Stroke",
      cases   = sum(inc_stroke == 1, na.rm = TRUE),
      py      = sum(fu_stroke, na.rm = TRUE)
    ),
  # Cancer
  elsa_dat %>%
    filter(fu_cancer > 0) %>%
    summarise(
      outcome = "Cancer",
      cases   = sum(inc_cancer == 1, na.rm = TRUE),
      py      = sum(fu_cancer, na.rm = TRUE)
    ),
  # Dementia
  cov_set60 %>%
    filter(surv_time_60 > 0) %>%
    summarise(
      outcome = "Dementia",
      cases   = sum(inc_alg_dem == 1, na.rm = TRUE),
      py      = sum(surv_time_60, na.rm = TRUE)
    )
) %>%
  mutate(
    case_py_label = paste0(cases, " / ", round(py, 0))
  )

df_plot_hr_irr_elsa <- bind_rows(df_plot, df_plot_irr) %>%
  filter(ediea_label == "90th-10th percentile difference") %>%
  left_join(
    cases_py_elsa %>% select(outcome, case_py_label),
    by = "outcome"
  ) %>%
  mutate(outcome = factor(outcome, levels = unique(outcome)))
  
p_forest_hr_irr_elsa <- ggplot() +
  facet_grid(outcome ~ ., scales = "free_y", space = "free_y") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # CI
  geom_errorbarh(
    data = df_plot_hr_irr_elsa %>% filter(!is.na(HR)),
    aes(y = ediea_label, xmin = LCI, xmax = UCI, color = is_sig),
    height = 0.15, size = 0.5
  ) +
  geom_point(
    data = df_plot_hr_irr_elsa %>% filter(!is.na(HR)),
    aes(y = ediea_label, x = HR, fill = is_sig),
    shape = 22, size = 2.5
  ) +
  # HR(95%CI)
  geom_text(
    data = df_plot_hr_irr_elsa,
    aes(y = ediea_label, x = 1.35, label = HR_text),
    hjust = 0, size = 3
  ) +
  # cases / person-yrs
  geom_text(
    data = df_plot_hr_irr_elsa,
    aes(y = ediea_label, x = 1.35, vjust=-2, label = case_py_label),
    hjust = 0, size = 3
  ) +
  scale_x_continuous(
    limits = c(0.2, 1.8),
    breaks = c(0.5, 1.0, 1.5)
  ) +
  labs(x = "HR/IRR (90th vs 10th percentile, 95% CI)", y = NULL, title = 'e. in ELSA') +
  theme_bmj() +
  theme(legend.position = "none", axis.text.y  = element_blank(),
        plot.title   = element_text(hjust = 0.5))

p_forest_hr_irr_elsa

