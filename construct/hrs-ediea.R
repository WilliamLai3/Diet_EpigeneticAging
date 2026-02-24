# Generating epirical dietary index for epigenetic aging

# Load libraries
library(tidyverse)
library(ggjournals)
library(broom)
library(cowplot)
library(patchwork)
library(mice)
library(conflicted)
library(rms)
library(survival)
source('./plot_rcs_linear.R')

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer('first', "dplyr")
conflict_prefer('unname', "base")


# Load data
load("dat/hrs_cov.rda")
load("dat/h13_fg.RData")
load('dat/hrs_health.RData')
epi_clock <- haven::read_sas("dat/epiclocka_r.sas7bdat")
epiclock_list <- names(epi_clock)[-c(1:2,16:17)]

gen_eaaz <- function(chron_age, epi_age){
  lmodel <- lm(epi_age ~ chron_age)
  z_res <- scale(lmodel$residuals)
  return(as.numeric(z_res))
}

gen_ageaccel <- function(chron_age, epi_age){
  fit <- lm(epi_age ~ chron_age)
  as.numeric(residuals(fit))
}

data_eaa <- left_join(hrs_cov, h13_fg_adj, by = c("HHID", "PN")) %>% 
  left_join(epi_clock, by = c("HHID", "PN"))

data_eaa <- data_eaa %>% filter(ffq_val) # 7410 in FFQ
data_eaa <- data_eaa %>% filter(!is.na(DNAMGRIMAGE)) # 1711 with epigenetic age

data_eaa <- data_eaa %>% 
  mutate(Horvath_EAA = gen_eaaz(age, HORVATH_DNAMAGE),
         Hannum_EAA = gen_eaaz(age, HANNUM_DNAMAGE),
         Levine_EAA = gen_eaaz(age, LEVINE_DNAMAGE),
         HorvathSkin_EAA = gen_eaaz(age, HORVATHSKIN_DNAMAGE),
         Lin_EAA = gen_eaaz(age, LIN_DNAMAGE),
         Weidner_EAA = gen_eaaz(age, WEIDNER_DNAMAGE),
         Vidalbralo_EAA = gen_eaaz(age, VIDALBRALO_DNAMAGE),
         Yang_EAA = gen_eaaz(age, YANG_DNAMAGE),
         Zhang_EAA = gen_eaaz(age, ZHANG_DNAMAGE),
         Bocklandt_EAA = gen_eaaz(age, BOCKLANDT_DNAMAGE),
         Garagnani_EAA = gen_eaaz(age, GARAGNANI_DNAMAGE),
         Grim_EAA = gen_eaaz(age, DNAMGRIMAGE),
         mPOA_EAA = scale(MPOA-1),
         
         Grim_EAA_r = gen_ageaccel(age, DNAMGRIMAGE),
         mPOA_EAA_r = MPOA-1,
         Levine_EAA_r = gen_ageaccel(age, LEVINE_DNAMAGE)
  )

covs <- c("age", "sex", "race", "educ", "hhinc",
          "smoke", "actc", "bmi_cat", "toten")
data_eaa[,covs] <- data_eaa %>%
  select(covs) %>%
  mice(seed = 2025) %>%
  complete()

# Food group and GrimAgeAccel
fg_var <- names(h13_fg_adj)[8:45]
fg_map <- data.frame(foodgroup = fg_var, 
                     fg_name =c("Low-fat dairy", "High-fat dairy", "Butter", "Margarine", "Fruits", 
                                "Fruit juice", "Tomato", "Green leafy vegetables", "Cruciferous vegetables", "Dark yellow vegetables", 
                                "Legumes", "Other vegetables", "Eggs", "Red meat", "Processed meat", 
                                "Poultry", "Organ meat", "Fish and other seafoods", "Whole grains", "Refined grains", 
                                "Pizza", "Potato", "Fries", "Snack", "Artificial sweetened beverages", 
                                "Sugar-sweetened beverages", "Beer", "Wine", "Liquor", "Tea", 
                                "Coffee", "Sweets", "Nuts", "Soup", "Condiments", "Garlic", 
                                "Mayonnaise", "Oils"))

# ------------- 1) linear for each foodgroup -----------
lm_fg <- function(fg_var, covs) {
  x <- data_eaa[[fg_var]]
  q <- quantile(x, probs = c(0.10, 0.90), na.rm = TRUE)
  range_90_10 <- diff(q)
  
  df <- data_eaa %>%
    mutate(x_90_10 = .data[[fg_var]] / range_90_10) %>%
    select(Grim_EAA, x_90_10, all_of(covs))
  
  lm_sum <- lm(Grim_EAA ~ ., data = df) %>% summary()
  
  out <- c(foodgroup = fg_var,
           lm_sum[["coefficients"]][2, ])  # 第一个自变量就是 x_90_10
  return(out)
}

table_fg <- map_df(fg_var, 
             ~lm_fg(fg_var = .x, covs = covs)) %>% 
  mutate(std.error = `Std. Error` %>% as.numeric(),
         Beta = Estimate %>% as.numeric(),
         LCI = Beta - 1.96 * std.error,
         UCI = Beta + 1.96 * std.error,
         P = `Pr(>|t|)` %>% as.numeric(),
         `P-value` = ifelse(P < 0.05, "P < 0.05", "P >= 0.05")) %>% 
  left_join(fg_map, by = c("foodgroup"))


# nonlinear test
library(splines)

lm_fg_nonlinear <- function(fg_var, covs, spline_df = 4) {
  x <- data_eaa[[fg_var]]
  q <- quantile(x, probs = c(0.10, 0.90), na.rm = TRUE)
  range_90_10 <- diff(q)

  if (!is.finite(range_90_10) || range_90_10 == 0) {
    return(tibble(
      foodgroup = fg_var,
      Estimate = NA_real_, `Std. Error` = NA_real_, `Pr(>|t|)` = NA_real_,
      P_nonlinear = NA_real_
    ))
  }
  
  df <- data_eaa %>%
    mutate(x_90_10 = .data[[fg_var]] / range_90_10) %>%
    select(Grim_EAA, x_90_10, all_of(covs)) %>%
    filter(is.finite(x_90_10), is.finite(Grim_EAA))
  
  fit_lin <- lm(Grim_EAA ~ x_90_10 + ., data = df)
  fit_spl <- lm(Grim_EAA ~ ns(x_90_10, df = spline_df) + ., data = df)
  
  p_nl <- anova(fit_lin, fit_spl)$`Pr(>F)`[2]
  
  co <- summary(fit_lin)$coefficients["x_90_10", , drop = FALSE]
  
  tibble(
    foodgroup = fg_var,
    Estimate = unname(co[1, "Estimate"]),
    `Std. Error` = unname(co[1, "Std. Error"]),
    `Pr(>|t|)` = unname(co[1, "Pr(>|t|)"]),
    P_nonlinear = as.numeric(p_nl)
  )
}

table_fg_nonlinear <- map_df(fg_var, ~ lm_fg_nonlinear(fg_var = .x, covs = covs, spline_df = 4)) %>%
  mutate(
    std.error = as.numeric(`Std. Error`),
    Beta = as.numeric(Estimate),
    LCI = Beta - 1.96 * std.error,
    UCI = Beta + 1.96 * std.error,
    P = as.numeric(`Pr(>|t|)`),
    `P-value` = ifelse(P < 0.05, "P < 0.05", "P >= 0.05"),
    P_nonlinear_fdr = p.adjust(P_nonlinear, method = "fdr"),
    `Non-linearity (FDR<0.05)` = ifelse(!is.na(P_nonlinear_fdr) & P_nonlinear_fdr < 0.05, "Yes", "No")
  ) %>%
  left_join(fg_map, by = "foodgroup")


table_fg_plot <- table_fg %>%
  left_join(
    table_fg_nonlinear %>% dplyr::select(foodgroup, P_nonlinear, P_nonlinear_fdr),
    by = "foodgroup"
  ) %>%
  mutate(
    # choose which p to display:
    # p_disp = P_nonlinear,          # raw P-nonlinear
    p_disp = P_nonlinear_fdr,        # <- use FDR if you prefer
    p_nl_label = dplyr::case_when(
      is.na(p_disp) ~ "",
      p_disp < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_disp)
    )
  )

# x-position of the p column (a bit to the right of the CI)
x_pcol <- with(table_fg_plot, max(UCI, na.rm = TRUE) + 0.12 * diff(range(c(LCI, UCI), na.rm = TRUE)))
n_fg   <- nrow(table_fg_plot)

plot_fg_eaa <- ggplot(
  table_fg_plot,
  aes(y = fct_reorder(fg_name, Beta), x = Beta, color = `P-value`)
) +
  geom_point() +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # P-nonlinear column
  geom_text(aes(x = x_pcol, label = p_nl_label),
          inherit.aes = TRUE, color = "black", size = 5, hjust = 0) +
  annotate("text", x = x_pcol, y = n_fg + 1,
           label = "p-nonlinear", hjust = 0, fontface = "bold", size = 6) +
  
  theme_bmj() +
  scale_color_manual(values = c("P < 0.05" = "#2A6EBB",
                                "P >= 0.05" = "#747678")) +
  labs(y = NULL,
       x = "GrimAgeAccel",
       title = "a. Food groups and GrimAge Acceleration") +
  theme(
    legend.position = "none",
    plot.margin = margin(6, 60, 6, 6)  # give room for the right column
  ) +
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = expansion(add = c(0.4, 1.2)))  # room for the header row

plot_fg_eaa

# --------- 2) elastic-net model to predict epigenetic age acceleration --------
library(glmnet)
library(stringr)
set.seed(2025)

n <- nrow(data_eaa)
train_id <- sample(1:n, size = round(0.7 * n))
test_id  <- base::setdiff(1:n, train_id)
data_train <- data_eaa[train_id, ]
data_test  <- data_eaa[test_id, ]

fg_cols <- fg_var

## z-score on Training set ------
X_train_raw <- as.matrix(data_train[, fg_cols])
X_train_z   <- scale(X_train_raw, center = TRUE, scale = TRUE)
mu_fg <- attr(X_train_z, "scaled:center")
sd_fg <- attr(X_train_z, "scaled:scale")

# zscore test/whole set
X_test_z <- scale(as.matrix(data_test[, fg_cols]),
                  center = mu_fg, scale = sd_fg)
X_all_z <- scale(as.matrix(data_eaa[, fg_cols]),
                 center = mu_fg, scale = sd_fg)

y_train <- data_train$Grim_EAA
y_all   <- data_eaa$Grim_EAA

alpha_grid <- seq(0.2, 1, by = 0.2)
cv_list <- map(alpha_grid, ~ cv.glmnet(
  X_train_z, y_train,
  alpha = .x,
  family = "gaussian",
  nfolds = 10,
  standardize = FALSE,
  dfmax = 11
))
cv_errors <- sapply(cv_list, function(cvfit) min(cvfit$cvm))
best_idx  <- which.min(cv_errors)
best_alpha <- alpha_grid[best_idx]
best_lambda <- cv_list[[best_idx]]$lambda.min
best_alpha; best_lambda

cvfit <- glmnet(
  X_train_z, y_train,
  alpha = best_alpha,
  lambda = best_lambda,
  family = "gaussian",
  standardize = FALSE
  # type.measure = "mse"
)
lambda_opt <- cvfit$lambda.min

ediea_model <- glmnet(
  X_all_z, y_all,
  alpha = best_alpha,
  lambda = best_lambda,
  family = "gaussian",
  standardize = FALSE
)
coef_all <- as.matrix(coef(ediea_model))
# leave intercept
beta_all <- coef_all[-1, 1]
names(beta_all) <- fg_cols
# nonzero
beta_sel <- beta_all[beta_all != 0]
sel_fg   <- names(beta_sel)

save(ediea_model, file = "./dat/EDISEA_model_HRS_train.rda")
ediea_coef <- -coef(ediea_model) %>%
  as.matrix() %>%
  data.frame() %>%
  rownames_to_column("foodgroup")
names(ediea_coef)[2] <- "Estimate"
ediea_coef <- ediea_coef %>%
  inner_join(fg_map, by = c('foodgroup')) %>%
  filter(abs(Estimate) > 0.0)

X_all_sel <- X_all_z[, sel_fg, drop = TRUE]
beta_sel_vec <- beta_sel[sel_fg]
EDISEA_cont <- -as.numeric(X_all_sel %*% beta_sel_vec)
data_eaa$ediea <- EDISEA_cont

# Correlation on Test set
EDISEA_test  <- EDISEA_cont[test_id]
Grim_test   <- data_eaa$Grim_EAA[test_id]
cor_spearman <- cor(EDISEA_test, Grim_test, method = "spearman", use = "complete.obs")
cor_spearman

plot_fg_ediea <- ggplot(ediea_coef, aes(y = fct_reorder(fg_name, Estimate), x = Estimate)) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.3f", Estimate), x = -0.12), size = 5) +
  # lolipop plot
  geom_segment(aes(xend = 0, yend = fg_name), size = 1, color = "#2A6EBB") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bmj() + 
  coord_cartesian(xlim = c(-0.13, 0.05)) +
  labs(y = NULL, x = "Weight",
       title = "b. Weights of food groups in EDISEA")
plot_fg_ediea

# --------------------- 3) Train/Test quantile -------------------
hrsfg_adj_z <- scale(as.matrix(h13_fg_adj[, fg_cols]),
                     center = mu_fg, scale = sd_fg)
hrsfg_adj_z <- hrsfg_adj_z[, sel_fg, drop = FALSE]
h13_fg_adj$ediea <- -as.numeric(hrsfg_adj_z %*% beta_sel)
save(h13_fg_adj, ediea_model, data_eaa, file = "./dat/hrs_ediea.rda")
hist(data_eaa$ediea)
cor(data_eaa$Grim_EAA, data_eaa$ediea, method = "spearman")

# ------------------------ Correlations with other epigenetic clocks --------------------
library(linkET)
eaa_list <- c("Grim_EAA_r",
              # "Zhang_EAA",
              "mPOA_EAA_r",
              "Levine_EAA_r", # PhenoAge
              "Horvath_EAA",
              "Hannum_EAA",

              # "Lin_EAA",
              # "Weidner_EAA",
              "Vidalbralo_EAA",

              "HorvathSkin_EAA",
              "Yang_EAA",
              "Bocklandt_EAA",
              "Garagnani_EAA"
)
corr_mat <- data_eaa %>% select(eaa_list)
corr_ediea <- cor(corr_mat, method = "spearman")
mantel <- data_eaa %>% select("ediea", eaa_list) %>%
  cor(method = "spearman") %>%
  as.data.frame() %>%
  mutate(spec = "EDISEA",
         env = rownames(.),
         `Correlation with EDISEA` = ediea %>% as.numeric()
  ) %>%
  dplyr::select(spec, env, `Correlation with EDISEA`)
plot_cor <- qcorrplot(corr_ediea, type = "lower", diag = FALSE) +
  geom_square() +
  scale_fill_gradient2(name = "rho between AgeAccels",
                       low = "#2A6EBB", mid = "white",
                       high = "#CD202C", midpoint = 0) +
  geom_couple(aes(colour = `Correlation with EDISEA`),
              data = mantel, size = 3,
              curvature = nice_curvature()) +
  scale_colour_gradient2(name = "rho with EDISEA", 
                         low = "#7D5CC6", mid = "white",
                         high = "#E37222", midpoint = 0) +
  theme(legend.position = "right", plot.title = element_text(face = "bold")) +
  labs(title = "c. Correlations between AgeAccels in HRS")

data_eaa_long <- data_eaa %>%
  select(HHID, PN, ediea, eaa_list[1:2]) %>%
  pivot_longer(cols = eaa_list[1:2], names_to = "eaa", values_to = "eaa_val")
data_eaa_long <- data_eaa_long %>%
  mutate(
    eaa_name = case_when(
      eaa == "Grim_EAA_r" ~ paste0("GrimAgeAccel\n r=",
                                 cor(data_eaa$Grim_EAA, data_eaa$ediea,
                                     method = "spearman") %>% round(3)),
      eaa == "Horvath_EAA" ~ paste0("HorvathAgeAccel, r=",
                                    cor(data_eaa$Horvath_EAA, data_eaa$ediea,
                                        method = "spearman") %>% round(3)),
      eaa == "Hannum_EAA" ~ paste0("HannumAgeAccel, r=",
                                   cor(data_eaa$Hannum_EAA, data_eaa$ediea,
                                       method = "spearman") %>% round(3)),
      eaa == "Levine_EAA_r" ~ paste0("PhenoAgeAccel\n r=",
                                   cor(data_eaa$Levine_EAA, data_eaa$ediea,
                                       method = "spearman") %>% round(3)),
      eaa == "Lin_EAA" ~ paste0("LinAgeAccel, r=",
                                cor(data_eaa$Lin_EAA, data_eaa$ediea,
                                    method = "spearman") %>% round(3)),
      eaa == "Zhang_EAA" ~ paste0("ZhangAgeAccel, r=",
                                  cor(data_eaa$Zhang_EAA, data_eaa$ediea,
                                      method = "spearman") %>% round(3)),
      eaa == "Weidner_EAA" ~ paste0("WeidnerAgeAccel, r=",
                                    cor(data_eaa$Weidner_EAA, data_eaa$ediea,
                                        method = "spearman") %>% round(3)),
      eaa == "Vidalbralo_EAA" ~ paste0("VidalBraloAgeAccel, r=",
                                       cor(data_eaa$Vidalbralo_EAA, data_eaa$ediea,
                                           method = "spearman") %>% round(3)),
      eaa == "mPOA_EAA_r" ~ paste0("DunedinPoAm\n r=",
                                 cor(data_eaa$mPOA_EAA, data_eaa$ediea,
                                     method = "spearman") %>% round(3))
    ),
    ediea_q = cut_number(ediea, 5, labels = FALSE) %>% factor(),
    eaa_name = factor(
      eaa_name,
      levels = c(
        eaa_name[eaa == "Grim_EAA"][1],
        base::setdiff(unique(eaa_name), eaa_name[eaa == "Grim_EAA"][1])
      )
    )
  )

df_poam <- data_eaa_long %>%
  filter(eaa == "mPOA_EAA_r") %>%
  mutate(ediea_q = cut_number(ediea, 5, labels = FALSE) %>% factor())
df_other <- data_eaa_long %>%
  filter(eaa != "mPOA_EAA_r") %>%
  mutate(ediea_q = cut_number(ediea, 5, labels = FALSE) %>% factor())

df_other %>%
  group_by(ediea_q) %>%
  filter(eaa=='Grim_EAA_r') %>%
  summarise(    n = sum(is.finite(eaa_val)),
                mean = mean(eaa_val, na.rm = TRUE),
                sd   = sd(eaa_val, na.rm = TRUE),
                med  = median(eaa_val, na.rm = TRUE),
                q50  = quantile(eaa_val, 0.5, na.rm = TRUE),
                q75  = quantile(eaa_val, 0.75, na.rm = TRUE),
                .groups = "drop")

ymin <- -15
ymax <-  20

poam_min <- quantile(df_poam$eaa_val, 0.01, na.rm = TRUE)
poam_max <- quantile(df_poam$eaa_val, 0.99, na.rm = TRUE)

poam_to_primary <- function(x) (x - poam_min) / (poam_max - poam_min) * (ymax - ymin) + ymin
primary_to_poam <- function(x) (x - ymin) / (ymax - ymin) * (poam_max - poam_min) + poam_min

df_poam <- df_poam %>%
  mutate(eaa_val_primary = poam_to_primary(eaa_val))

p_other <- ggplot(df_other, aes(x = ediea_q, y = eaa_val)) +
  geom_boxplot() +
  facet_wrap(~ eaa_name, nrow = 1) +
  theme_bmj() +
  labs(x = "EDISEA quintiles", y = "AgeAccel",
       title = "c. EDISEA and 2nd-generation epigenetic aging clocks in HRS") +
  coord_cartesian(ylim = c(ymin, ymax))

poam_breaks <- pretty(c(poam_min, poam_max), n = 5)
p_poam <- ggplot(df_poam, aes(x = ediea_q, y = eaa_val_primary)) +
  geom_boxplot(fill = "#5F7FA2") +
  facet_wrap(~ eaa_name, nrow = 1) +
  theme_bmj() +
  labs(x = "EDISEA quintiles", y = "AgeAccel") +
  scale_y_continuous(
    limits = c(ymin, ymax),
    breaks = poam_to_primary(poam_breaks),
    
    sec.axis = sec_axis(
      ~ primary_to_poam(.), 
      name = "DunedinPoAm",
      breaks = poam_breaks 
    )
  ) +
  theme(
    axis.text.y.left  = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.line.y.left  = element_blank(),
    axis.title.y = element_blank(), 
    axis.title.x = element_blank()
  )

plot_ediea_gen2 <- (p_other | p_poam) + plot_layout(widths = c(1, 1))
plot_ediea_gen2

# --------------------------------- External Validation -------------------------
# Load data
load("../NHANES/cov/Cov_AJ.rda")
DEMO_select <- DEMO_AJ %>%
  select(SEQN, AGE, GENDER, RACE, EDU, PIR) %>%
  left_join(
    SMQ_AJ,
    by = "SEQN"
  ) %>%
  left_join(
    PA_AJ %>% mutate(SEQN = as.numeric(SEQN)),
    by = "SEQN"
  ) %>%
  left_join(
    BMI_AJ,
    by = "SEQN"
  )

load("../NHANES/Diet&mor/fped_eq/hrsfg_adj.rda")
epi_clock <- haven::read_sas("../NHANES/cov/dnmepi.sas7bdat")
epiclock_list <- c("HorvathAge", "HannumAge", "SkinBloodAge", "PhenoAge",
                   "ZhangAge", "LinAge", "WeinerAge", "VidalBraloAge", "GrimAgeMort", "GrimAge2Mort")

data_eaa_nhanes <- left_join(hrsfg_adj, DEMO_select, by = c("SEQN")) %>%
  left_join(epi_clock, by = c("SEQN")) # 4449

data_eaa_nhanes_all <- data_eaa_nhanes
data_eaa_nhanes <- data_eaa_nhanes %>% filter(!is.na(GrimAgeMort) & !is.na(AGE)) # 2311 with epigenetic age

data_eaa_nhanes <- data_eaa_nhanes %>%
  mutate(Grim_EAA = gen_eaaz(AGE, GrimAgeMort),
         Horvath_EAA = gen_eaaz(AGE, HorvathAge),
         Hannum_EAA = gen_eaaz(AGE, HannumAge),
         SkinBlood_EAA = gen_eaaz(AGE, SkinBloodAge),
         Pheno_EAA = gen_eaaz(AGE, PhenoAge),
         Zhang_EAA = gen_eaaz(AGE, ZhangAge),
         Lin_EAA = gen_eaaz(AGE, LinAge),
         Weidner_EAA = gen_eaaz(AGE, WeidnerAge),
         VidalBralo_EAA = gen_eaaz(AGE, VidalBraloAge),
         DunedinPoAm_EAA = scale(DunedinPoAm-1),
         HorvathTelo_EAA = gen_eaaz(AGE, HorvathTelo),
         
         Grim_EAA_r = gen_ageaccel(AGE, GrimAgeMort),
         Pheno_EAA_r = gen_ageaccel(AGE, PhenoAge),
         DunedinPoAm_EAA_r = DunedinPoAm-1
  )

# Food group and GrimAgeAccel
fg_var <- names(hrsfg_adj)[2:39]
covs <- c("AGE", "GENDER", "RACE", "EDU", "PIR", "METS",
          "SMO", "BMIC", "CALOR_SUM") # PIR, METS: >10% missing
data_eaa_nhanes[,covs] <- data_eaa_nhanes %>%
  select(covs) %>%
  mice(seed = 2025) %>%
  complete()

fg_map <- data.frame(foodgroup = fg_var,
                     fg_name =c("lfdai","hfdai","buttr","marga","fruit","fruju","tomat","gleaf","cruci","dkyel",
                                "legum","othvg","eggs","rmeat","promt","poult","organ","fishs","whgrn","rfgrn",
                                "pizza","potat","fries","snack","artsb","sugsb","beer","wine","liqur","tea",
                                "coffe","sweet","nuts","soup","condi","garlc","mayon","oils"))
name_lookup <- setNames(fg_map$fg_name, fg_map$foodgroup)
newdata_fg <- data_eaa_nhanes %>%
  select(fg_var)
names(newdata_fg) <- name_lookup[names(newdata_fg)]

newdata_fg_all <- data_eaa_nhanes_all %>%
  select(fg_var)
names(newdata_fg_all) <- name_lookup[names(newdata_fg_all)]
# --------- 2) validation --------
library(ggpubr)
# my_comparisons <- combn(paste0(1:5), 2, simplify = FALSE)
my_comparisons <- list(
  c("1", "2"),
  c("1", "3"),
  c("1", "4"),
  c("1", "5")
)

y_max <- 2.2

X_ext_raw <- as.matrix(newdata_fg[, fg_cols])
X_ext_z <- scale(X_ext_raw)
X_ext_z <- X_ext_z[, sel_fg, drop = FALSE]
# X_ext_z <- as.matrix(newdata_fg[, sel_fg])
data_eaa_nhanes$ediea <- -as.numeric(X_ext_z %*% beta_sel)

X_ext_raw_all <- as.matrix(newdata_fg_all[, fg_cols])
X_ext_z_all <- scale(X_ext_raw_all)
X_ext_z_all <- X_ext_z_all[, sel_fg, drop = FALSE]
data_eaa_nhanes_all$ediea <- -as.numeric(X_ext_z_all %*% beta_sel)
save(data_eaa_nhanes_all, file='./dat/nhanes_ediea_from_hrs.rda')

eaa_list <- c("Grim_EAA_r",
              # "Zhang_EAA",
              "DunedinPoAm_EAA_r",
              "Pheno_EAA_r",
              "HorvathTelo_EAA",
              "Horvath_EAA",
              "Hannum_EAA",
              "Lin_EAA",
              "Weidner_EAA",
              "VidalBralo_EAA",

              "SkinBlood_EAA"
)

data_eaa_long_nhanes <- data_eaa_nhanes %>%
  select(SEQN, ediea, eaa_list[1:2]) %>%
  pivot_longer(cols = eaa_list[1:2], names_to = "eaa", values_to = "eaa_val")
data_eaa_long_nhanes <- data_eaa_long_nhanes %>%
  mutate(
    eaa_name = case_when(
      eaa == "Grim_EAA_r" ~ paste0("GrimAgeAccel\n r=",
                                 cor(data_eaa_nhanes$Grim_EAA, data_eaa_nhanes$ediea,
                                     method = "spearman") %>% round(3)),
      eaa == "Zhang_EAA" ~ paste0("ZhangAgeAccel\n r=",
                                  cor(data_eaa_nhanes$Zhang_EAA, data_eaa_nhanes$ediea,
                                      method = "spearman") %>% round(3)),
      eaa == "Horvath_EAA" ~ paste0("HorvathAgeAccel\n r=",
                                    cor(data_eaa_nhanes$Horvath_EAA, data_eaa_nhanes$ediea,
                                        method = "spearman") %>% round(3)),
      eaa == "Pheno_EAA_r" ~ paste0("PhenoAgeAccel\n r=",
                                  cor(data_eaa_nhanes$Pheno_EAA, data_eaa_nhanes$ediea,
                                      method = "spearman") %>% round(3)),
      eaa == "DunedinPoAm_EAA_r" ~ paste0("DunedinPoAm\n r=",
                                        cor(data_eaa_nhanes$DunedinPoAm_EAA, data_eaa_nhanes$ediea,
                                            method = "spearman") %>% round(3)),
      eaa == "HorvathTelo_EAA" ~ paste0("HorvathTeloAge\n r=",
                                        cor(data_eaa_nhanes$HorvathTelo_EAA, data_eaa_nhanes$ediea,
                                            method = "spearman") %>% round(3)),
    ),
    ediea_q = cut_number(ediea, 5, labels = FALSE) %>% factor(),
    eaa_name = factor(
      eaa_name,
      levels = c(
        eaa_name[eaa == "Grim_EAA"][1],
        base::setdiff(unique(eaa_name), eaa_name[eaa == "Grim_EAA"][1])
      )
    )
  )

df_poam_nhanes <- data_eaa_long_nhanes %>%
  filter(eaa == "DunedinPoAm_EAA_r") %>%
  mutate(ediea_q = cut_number(ediea, 5, labels = FALSE) %>% factor())
df_other_nhanes <- data_eaa_long_nhanes %>%
  filter(eaa != "DunedinPoAm_EAA_r") %>%
  mutate(ediea_q = cut_number(ediea, 5, labels = FALSE) %>% factor())

df_other_nhanes %>%
  group_by(ediea_q) %>%
  filter(eaa=='Grim_EAA_r') %>%
  summarise(    n = sum(is.finite(eaa_val)),
                mean = mean(eaa_val, na.rm = TRUE),
                sd   = sd(eaa_val, na.rm = TRUE),
                med  = median(eaa_val, na.rm = TRUE),
                q50  = quantile(eaa_val, 0.5, na.rm = TRUE),
                q75  = quantile(eaa_val, 0.75, na.rm = TRUE),
                .groups = "drop")

poam_min_nhanes <- quantile(df_poam_nhanes$eaa_val, 0.01, na.rm = TRUE)
poam_max_nhanes <- quantile(df_poam_nhanes$eaa_val, 0.99, na.rm = TRUE)

poam_to_primary_nhanes <- function(x) (x - poam_min_nhanes) / (poam_max_nhanes - poam_min_nhanes) * (ymax - ymin) + ymin
primary_to_poam_nhanes <- function(x) (x - ymin) / (ymax - ymin) * (poam_max_nhanes - poam_min_nhanes) + poam_min_nhanes

df_poam_nhanes <- df_poam_nhanes %>%
  mutate(eaa_val_primary = poam_to_primary_nhanes(eaa_val))

p_other_nhanes <- ggplot(df_other_nhanes, aes(x = ediea_q, y = eaa_val)) +
  geom_boxplot() +
  facet_wrap(~ eaa_name, nrow = 1) +
  theme_bmj() +
  labs(x = "EDISEA quintiles", y = "AgeAccel",
       title = "d. EDISEA and 2nd-generation epigenetic clocks in NHANES") +
  coord_cartesian(ylim = c(ymin, ymax))

poam_breaks_nhanes <- pretty(c(poam_min, poam_max), n = 5)
p_poam_nhanes <- ggplot(df_poam_nhanes, aes(x = ediea_q, y = eaa_val_primary)) +
  geom_boxplot(fill = "#5F7FA2") +
  facet_wrap(~ eaa_name, nrow = 1) +
  theme_bmj() +
  labs(x = "EDISEA quintiles", y = "AgeAccel") +
  scale_y_continuous(
    limits = c(ymin, ymax),
    breaks = poam_to_primary_nhanes(poam_breaks_nhanes),
    
    sec.axis = sec_axis(
      ~ primary_to_poam_nhanes(.), 
      name = "DunedinPoAm",
      breaks = poam_breaks_nhanes 
    )
  ) +
  theme(
    axis.text.y.left  = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.line.y.left  = element_blank(),
    axis.title.y = element_blank(), 
    axis.title.x = element_blank()
  )

plot_ediea_extern <- (p_other_nhanes | p_poam_nhanes) + plot_layout(widths = c(1, 1))
plot_ediea_extern

## gather plots
design2 <- "
AB
AC
AD
EF
"
base_size <- 18

big_theme <- theme(
  text        = element_text(size = base_size),
  plot.title  = element_text(size = base_size + 2, face = "bold", hjust = 0),
  axis.title  = element_text(size = base_size),
  axis.text   = element_text(size = base_size - 2),
  legend.title= element_text(size = base_size - 1),
  legend.text = element_text(size = base_size - 2),
  strip.text  = element_text(size = base_size - 1)
)

plot_eaa_val <-
  free(plot_fg_eaa) +
  free(plot_fg_ediea) +
  free(plot_ediea_gen2) +
  free(plot_ediea_extern) +
  free(plot_ediea_eaa_bar_hrs2) + free(plot_ediea_eaa_bar_nhanes2) +
  plot_layout(design = design2, widths = c(1,1))

plot_eaa_val <- plot_eaa_val & big_theme

ggsave(plot_eaa_val,
       filename = "fig/figure_2.pdf",
       width = 18, height = 18)

