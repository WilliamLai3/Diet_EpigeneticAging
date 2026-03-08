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
library(scales)
library(tidytext)

load("./dat/ukb_ediea_from_hrs.rda")
load("../UKB/data_gen/covars.rda")
load("../UKB/ukb_outcomes.rda") # inc_Mental_dementia
olink_data <- read_parquet("../UKB/dataset/participant_olink.parquet")
names(olink_data) <- names(olink_data) |> str_replace("olink_instance_0.", "")

ukb_data <- ukb_fg40 %>%
  left_join(data_cov, by = "eid") %>%
  left_join(disease_dates %>% 
              dplyr::select(eid, date_death, date_ad), 
            by = 'eid')
ukb_olink <- ukb_data %>%
  inner_join(olink_data %>%
               mutate(eid = as.numeric(eid)), by = "eid")

cox_covars <- c("age_rec", "sex", "edu_college", "eth_white", "bmi", 
                "tdi", "smoking", "drinking", "actc", "p30901","p30902",
                "p74","p54","p3166","p21842")
ukb_olink[,cox_covars] <- ukb_olink %>%
  dplyr::select(cox_covars) %>%
  mice(seed = 2025) %>%
  complete()

# ---------------- Multiple linear regression, EDISEA ~ protein NPX --------------
protein_cols <- setdiff(names(olink_data), "eid")

ukb_olink_scaled <- ukb_olink %>%
  mutate(ediea_group = as.numeric(cut_number(ediea, 3)) %>% factor()) %>%
  mutate(across(all_of(protein_cols), ~ as.numeric(scale(.x))))

eof <- max(ukb_olink_scaled$date_death, na.rm = TRUE)
ukb_olink_scaled <- ukb_olink_scaled %>%
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
fit_cox_base <- coxph(cox_formula_base, data = ukb_olink_scaled)
beta_total <- coef(fit_cox_base)["ediea"]
tidy(fit_cox_base) %>%
  filter(str_detect(term, "ediea")) %>%
  transmute(
    term,
    HR   = exp(estimate),
    LCL  = exp(estimate - 1.96 * std.error),
    UCL  = exp(estimate + 1.96 * std.error),
    p    = p.value,
    n    = nrow(ukb_olink_scaled),
    events = sum(ukb_olink_scaled$event)
  )

run_mediation_for_protein <- function(prot_name, data, sims = 3000) {
  message("Running mediation for: ", prot_name)

  dat <- data %>%
    dplyr::filter(!is.na(.data[[prot_name]])) %>%
    dplyr::filter(
      !is.na(ediea),
      !is.na(fu),
      !is.na(event),
      dplyr::if_all(dplyr::all_of(cox_covars), ~ !is.na(.x))
    )
  
  if (nrow(dat) < 100) {
    return(tibble(
      protein          = prot_name,
      n                = nrow(dat),
      beta_ediea_prot  = NA_real_,
      se_ediea_prot    = NA_real_,
      p_ediea_prot     = NA_real_,
      beta_prot_dem    = NA_real_,
      se_prot_dem      = NA_real_,
      p_prot_dem       = NA_real_,
      beta_ediea_direct= NA_real_,
      HR_ediea_direct  = NA_real_,
      prop_med         = NA_real_,
      acme             = NA_real_,
      acme_p           = NA_real_,
      ade              = NA_real_,
      ade_p            = NA_real_,
      total_effect     = NA_real_,
      total_p          = NA_real_,
      prop_med_qb      = NA_real_,
      prop_med_qb_p    = NA_real_
    ))
  }
  
  lm_formula <- as.formula(
    paste0(
      prot_name, " ~ ediea + ",
      paste(cox_covars, collapse = " + ")
    )
  )
  
  fit_lm <- stats::lm(lm_formula, data = dat)
  fit_lm$call$formula <- stats::formula(fit_lm)
  
  sm_lm  <- summary(fit_lm)
  
  beta_ediea_prot <- sm_lm$coefficients["ediea", "Estimate"]
  se_ediea_prot   <- sm_lm$coefficients["ediea", "Std. Error"]
  p_ediea_prot    <- sm_lm$coefficients["ediea", "Pr(>|t|)"]

  cox_formula <- as.formula(
    paste0(
      "Surv(fu, event) ~ ediea + ",
      prot_name, " + ",
      paste(cox_covars, collapse = " + ")
    )
  )
  
  fit_cox <- survreg(cox_formula, data = dat)
  fit_cox$call$formula <- stats::formula(fit_cox)
  
  sm_cox  <- summary(fit_cox)
  print(sm_cox$coefficients[prot_name])
  
  beta_ediea_direct <- sm_cox$coefficients["ediea"]
  # HR_ediea_direct   <- sm_cox$coefficients["ediea", "exp(coef)"]
  tab <- sm_cox$table
  beta_prot_dem     <- sm_cox$coefficients[prot_name]
  se_prot_dem       <- tab[prot_name, "Std. Error"]
  p_prot_dem        <- tab[prot_name, "p"]
  
  med_out <- mediation::mediate(
        model.m   = fit_lm,
        model.y   = fit_cox,
        treat     = "ediea",
        mediator  = prot_name,
        # sims      = sims,
        boot = FALSE
      )
  acme         <- med_out$d0
  acme_p       <- med_out$d0.p
  ade          <- med_out$z0
  ade_p        <- med_out$z0.p
  total_effect <- med_out$tau.coef
  total_p      <- med_out$tau.p
  prop_med_qb  <- med_out$n0
  prop_med_qb_p<- med_out$n0.p
  prop_med     <- prop_med_qb

  tibble(
    protein           = prot_name,
    n                 = nrow(dat),
    beta_ediea_prot   = beta_ediea_prot,
    se_ediea_prot     = se_ediea_prot,
    p_ediea_prot      = p_ediea_prot,
    beta_prot_dem     = beta_prot_dem,
    se_prot_dem       = se_prot_dem,
    p_prot_dem        = p_prot_dem,
    beta_ediea_direct = beta_ediea_direct,
    prop_med          = prop_med,
    acme              = acme,
    acme_p            = acme_p,
    ade               = ade,
    ade_p             = ade_p,
    total_effect      = total_effect,
    total_p           = total_p,
    prop_med_qb       = prop_med_qb,
    prop_med_qb_p     = prop_med_qb_p
  )
}

# prot_results <-   map_dfr(
#   protein_cols,
#   ~ run_mediation_for_protein(.x, data = ukb_olink_scaled)
# )
# save(prot_results, file = './dat/ukb_prot.rda')

## ---------------------- PROTEOMICS ASSOC -------------------------
load('./dat/ukb_prot.rda')
prot_results <- prot_results %>%
  mutate(
    protein = protein %>% toupper(),
    fdr_ediea_prot = p.adjust(p_ediea_prot, method = "fdr"),
    fdr_prot_dem   = p.adjust(p_prot_dem,   method = "fdr"),
    fdr_acme_p = p.adjust(acme_p, method = "fdr"),
    HR_prot_dem    = exp(beta_prot_dem),
    fdr_ediea_prot2 = pmax(fdr_ediea_prot, 1e-300),
    neglog10_fdr    = -log10(fdr_ediea_prot2),
    direction       = if_else(beta_ediea_prot < 0, "Inverse", "Positive"),
    sig             = fdr_ediea_prot < 0.05
  ) # 810 significant proteins

top10_labels <- prot_results %>%
  filter(sig) %>%
  group_by(direction) %>%
  arrange(fdr_ediea_prot, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

## FIGURE: EDISEA ~ proteomics
plot_ediea_proteomics <- ggplot(
  prot_results,
  aes(x = beta_ediea_prot, y = neglog10_fdr)
) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(
    data = ~ filter(.x, !is.na(neglog10_fdr)),
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
    data = top10_labels,
    aes(label = protein),
    size          = 5,
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
    legend.position = "none",
    axis.title.x = element_text(face = "plain"),
    axis.title.y = element_text(face = "plain")
  ) +
  labs(
    x = "beta coefficient",
    y = "-log10(P_FDR value)",
    title = "a. Plasma proteomic profiling of the EDISEA"
  )

plot_ediea_proteomics

# --------------------- PROTEOMICS Mediator ------------------------
mediators <- prot_results %>%
  filter(
    !is.na(beta_ediea_prot),
    !is.na(beta_prot_dem),
    !is.na(prop_med),
    fdr_ediea_prot < 0.05,
    fdr_prot_dem < 0.05,
    fdr_acme_p < 0.05,
    prop_med > 0.1,
    # direction=="Positive"
  ) %>%
  mutate(
    direction = if_else(beta_ediea_prot < 0, "Inverse", "Positive")
  )


top_mediators <- mediators %>%
  arrange(desc(prop_med)) %>%
  slice_head(n = 10)
top_mediators %>%
  dplyr::select(protein, beta_ediea_prot, HR_prot_dem, prop_med,
         fdr_ediea_prot, fdr_prot_dem)

## --------------------------- ENRICHMENT -------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_cran <- c("patchwork","stringr","forcats")
pkgs_bioc <- c("clusterProfiler","org.Hs.eg.db","enrichplot","DOSE","msigdbr")

invisible(lapply(pkgs_cran, function(p) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)))
invisible(lapply(pkgs_bioc, function(p) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, update=FALSE, ask=FALSE)))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(forcats)
ratio_to_num <- function(x){
  sapply(strsplit(as.character(x), "/"), \(z) as.numeric(z[1]) / as.numeric(z[2]))
}

sig_prots <- prot_results %>%
  filter(
    !is.na(beta_ediea_prot),
    !is.na(beta_prot_dem),
    !is.na(prop_med),
    fdr_ediea_prot < 0.05,
    fdr_prot_dem < 0.05,
    fdr_acme_p < 0.05
    # prop_med > 0.1,
    # direction=="Positive"
  ) %>%
  mutate(
    direction = if_else(beta_ediea_prot < 0, "Inverse", "Positive")
  ) %>%
  dplyr::pull(protein) %>% unique() %>% toupper()

bg_prots <- unique(protein_cols) %>% toupper()

map_to_entrez <- function(x) {
  x <- unique(na.omit(x))
  
  m1 <- suppressMessages(
    bitr(x, fromType="SYMBOL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)
  )
  
  miss <- setdiff(x, m1$SYMBOL)
  m2 <- tibble::tibble()
  if (length(miss) > 0) {
    m2 <- suppressMessages(
      bitr(miss, fromType="ALIAS", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)
    )
  }
  
  dplyr::bind_rows(m1, m2) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
}

map_sig <- map_to_entrez(sig_prots)
map_bg  <- map_to_entrez(bg_prots)

gene_sig <- unique(map_sig$ENTREZID)
gene_bg  <- unique(map_bg$ENTREZID)

message("Mapped mediators: ", length(gene_sig), " / ", length(sig_prots))
message("Mapped background: ", length(gene_bg), " / ", length(bg_prots))

# GO
run_go <- function(ont) {
  ego <- enrichGO(
    gene          = gene_sig,
    universe      = gene_bg,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ont,          # "BP" "CC" "MF"
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  }
  ego
}

ego_bp <- run_go("BP")
ego_cc <- run_go("CC")
ego_mf <- run_go("MF")

# KEGG
ekegg <- tryCatch({
  enrichKEGG(
    gene          = gene_sig,
    universe      = gene_bg,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  ) |> setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}, error = function(e) NULL)

go_bp_df <- as.data.frame(ego_bp)
go_cc_df <- as.data.frame(ego_cc)
go_mf_df <- as.data.frame(ego_mf)
kegg_df <- as.data.frame(ekegg)

go_all <- dplyr::bind_rows(
  go_bp_df %>% dplyr::mutate(ont = "BP"),
  go_cc_df %>% dplyr::mutate(ont = "CC"),
  go_mf_df %>% dplyr::mutate(ont = "MF"),
  kegg_df %>% dplyr::mutate(ont = "KEGG"),
)

go_top5 <- go_all %>%
  group_by(ont) %>%
  slice_min(order_by = p.adjust, n=5) %>%
  ungroup() %>%
  mutate(
    ont = factor(ont, levels = c("BP", "CC", "MF", "KEGG")),
    term = reorder_within(Description, ratio_to_num(GeneRatio), ont)
  )

count_lim <- range(go_top5$Count, na.rm = TRUE)
padj_lim  <- range(go_top5$p.adjust, na.rm = TRUE)

p_ef <- ggplot(go_top5, aes(x = ratio_to_num(GeneRatio), y = term)) +
  geom_point(aes(size = Count, colour = p.adjust)) +
  facet_grid(ont ~ ., scales = "free_y") +
  scale_y_reordered() +
  scale_size_continuous(limits = count_lim, range = c(2.0, 7.0)) +
  scale_color_gradientn(
    colours = c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C"),
    trans = "log10",
    breaks = scales::log_breaks(n = 5),
    labels = scales::label_scientific(),
    name = "p.adjust"
  ) +
  guides(
    size = guide_legend(override.aes = list(shape = 21, fill = NA, colour = "grey20", alpha = 1)),
    fill = guide_colorbar(barheight = unit(35, "mm"))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.10))) +
  theme_bw() +
  theme(
    strip.placement = "outside",
    strip.text.y.right = element_text(face = "bold", size = 10),
    strip.text.y.left = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(6, 10, 6, 6)
  ) +
  labs(
    title = "c. Pathways enrichment analysis",
    x = "GeneRatio"
  )

p_ef

## --------------------------- PROTEOMICS PPI ---------------------------------
sig_prots <- mediators %>% dplyr::pull(protein) %>% unique()

genes <- tibble(gene = sig_prots) %>%
  distinct()

# 2) INIT STRING（9606）
string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = 0,
  input_directory = ""
)

genes_df <- data.frame(
  gene = as.character(genes$gene), stringsAsFactors = FALSE
)

genes_df$gene <- toupper(trimws(genes_df$gene))
genes_df <- genes_df[genes_df$gene != "" & !is.na(genes_df$gene), , drop = FALSE]
genes_df <- unique(genes_df)

mapped <- string_db$map(genes_df, "gene", removeUnmappedRows = TRUE)

ppi_raw <- string_db$get_interactions(mapped$STRING_id)

ppi <- ppi_raw %>%
  as_tibble() %>%
  dplyr::rename(from = from, to = to, score = combined_score) %>%
  filter(from %in% mapped$STRING_id, to %in% mapped$STRING_id) %>%
  filter(score >= 700)

id2gene <- mapped %>% dplyr::select(STRING_id, gene)
ppi2 <- ppi %>%
  left_join(id2gene, by = c("from" = "STRING_id")) %>% dplyr::rename(from_gene = gene) %>%
  left_join(id2gene, by = c("to"   = "STRING_id")) %>% dplyr::rename(to_gene   = gene) %>%
  dplyr::select(from_gene, to_gene, score)

ppi2 %>% arrange(desc(score)) %>% head(10)


key_nodes <- mediators %>%
  arrange(desc(prop_med)) %>%
  slice_head(n = 20) %>%
  pull(protein) %>% unique() %>% toupper()

node_df <- tibble(name = sort(unique(c(ppi2$from_gene, ppi2$to_gene)))) %>%
  left_join(mediators %>% dplyr::select(protein, prop_med, beta_ediea_prot, fdr_acme_p, direction),
            by = c("name" = "protein")) %>%
  mutate(
    is_key = name %in% key_nodes,

    node_group = if_else(is_key, "Key mediators", "Other mediators")
  )

g <- graph_from_data_frame(ppi2, directed = FALSE, vertices = node_df)

comp <- components(g)
g <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))

set.seed(2025)
xy <- tibble(name = V(g)$name,
             node_group = if_else(name %in% key_nodes, "Key mediators", "Other mediators")) %>%
  group_by(node_group) %>%
  mutate(theta = seq(0, 2*pi, length.out = n() + 1)[-1],
         r     = if_else(node_group == "Key mediators", 0.1, 0.2),
         x     = r * cos(theta),
         y     = r * sin(theta)) %>%
  ungroup()

V(g)$prop_med  <- ifelse(is.na(V(g)$prop_med), 0, V(g)$prop_med)
V(g)$direction <- ifelse(is.na(V(g)$direction), "Other", V(g)$direction)

dx <- 0.02 * diff(range(xy$x, na.rm = TRUE))
dy <- 0.02 * diff(range(xy$y, na.rm = TRUE))
p_ppi <- ggraph(g, layout = "manual", x = xy$x, y = xy$y) +
  geom_edge_bend(aes(alpha = score), colour = "grey60", show.legend = FALSE) +
  geom_node_point(
    aes(size=prop_med, fill = node_group),
    shape=21, colour = 'black', stroke=0.6, alpha=0.95
  ) +
  scale_size_continuous(
    range = c(2.5, 14),
    trans = "sqrt",
    breaks = pretty_breaks(4),
    labels = percent_format(accuracy = 1),
    name = "Proportion"
  ) +
  geom_node_text(
    aes(label = name, vjust = -2.0),
    colour = "black",
    size = 4
  ) +
  scale_fill_manual(values = c("Key mediators" = "#D1A679", "Other mediators" = "#AFC8E8"),
                    guide = 'none') +
  scale_edge_alpha(range = c(0.2, 0.9)) +
  theme_void() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold")) +
  labs(title = "b. Protein-protein interactions of key mediators")

p_ppi
ggsave(p_ppi, file = 'fig/figure_1/ppi.pdf', width = 8, height = 6)


## gather plots
library(patchwork)
design <- "ABBB
CCDD"

plot_prot <- free(plot_ediea_proteomics) +
 (p_ppi) +
  (p_ef) +
  free(p_corr) +
  plot_layout(design = design)