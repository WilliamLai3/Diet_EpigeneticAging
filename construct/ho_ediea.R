library(tidyverse)

load('../UKB/ukb_fg40.rda')
load('./dat/EDISEA_model_HRS_train.rda')


fg_var <- names(ukb_fg40)[5:42]
fg_map <- data.frame(foodgroup = fg_var, 
                     fg_name =c("Low-fat dairy", "High-fat dairy", "Butter", "Margarine", "Fruits", 
                                "Fruit juice", "Tomato", "Green leafy vegetables", "Cruciferous vegetables", "Dark yellow vegetables", 
                                "Legumes", "Other vegetables", "Eggs", "Red meat", "Processed meat", 
                                "Poultry", "Organ meat", "Fish and other seafoods", "Whole grains", "Refined grains", 
                                "Pizza", "Potato", "Fries", "Snack", "Artificial sweetened beverages", 
                                "Sugar-sweetened beverages", "Beer", "Wine", "Liquor", "Tea", 
                                "Coffee", "Sweets",  "Nuts", "Soup", "Condiments", "Garlic", 
                                "Mayonnaise", "Oils"))

coef_all <- as.matrix(coef(ediea_model))
# leave intercept
beta_all <- coef_all[-1, 1]
names(beta_all) <- fg_var
# nonzero
beta_sel <- beta_all[beta_all != 0]
sel_fg   <- names(beta_sel)

ediea_coef <- -coef(ediea_model) %>%
  as.matrix() %>%
  data.frame() %>%
  rownames_to_column("foodgroup")
names(ediea_coef)[2] <- "Estimate"
ediea_coef <- ediea_coef %>%
  inner_join(fg_map, by = c('foodgroup')) %>%
  filter(abs(Estimate) > 0.0)

## 连续型 EDISEA 指数：Z-score 加权求和 -----
X_all = as.matrix(ukb_fg40[, fg_var])
X_all_z <- scale(X_all)
X_all_z <- X_all_z[, sel_fg, drop = FALSE]

EDISEA_cont <- -as.numeric(X_all_z %*% beta_sel)
ukb_fg40$ediea <- EDISEA_cont

save(ukb_fg40, file = './dat/ukb_ediea_from_hrs.rda')
