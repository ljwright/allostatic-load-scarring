library(tidyverse)
library(glue)
library(mice)
library(Hmisc)
library(officer)
library(flextable)
library(magrittr)
library(haven)
library(gallimaufr)

rm(list = ls())

## Weights don't seem to be right
## CC for GHQ gives large value for allostatic mediation
## DECIDE ON FINAL COVARIATES (DROP THOSE NOT IN USE)
## ADD SCA PLOTS
## RUN SCA AND SAVE REGRESSION OBJECTS FROM 3_REGRESSIONS TO ENSURE CONSISTENCY IN 4_SCA


# 1. Load Data ----
load("Data/df_raw.Rdata")
load("Data/imp.Rdata")

bio_pretty <- function(var){
  raw_vars <- str_subset(names(all_vars), var)
  clean_name <- pretty_bio[[var]]
  
  case_when(str_detect(raw_vars, "Quartile") ~ glue("{clean_name} (High-Risk)"),
            str_detect(raw_vars, "Index") ~ glue("{clean_name} (Index)"),
            TRUE ~ glue("{clean_name} (Z-Score)")) %>%
    as.character() %>%
    set_names(raw_vars)
}

pretty_bio <- c(Allostatic = "Allostatic Load", HbA1c = "HbA1c", Insulin = "Insulin-like Growth Factor 1", 
                CReactive = "C-reactive Protein", ClaussFib = "Fibrinogen", DHEAS = "DHEA-S",
                Creatinine = "Creatinine Clearence Rate",  Trig = "Triglycerids",
                HDL_Ratio = "Total-to-HDL Cholesterol Ratio", Pulse = "Pulse", 
                Systolic = "Systolic Blood Pressure", Diastolic = "Diastolic Blood Pressure",
                WtH_Ratio = "Waist-to-Height Ratio")

pretty_lbls <- c(
  n = "n",
  Unem_Age = "Youth Unemployment",
  map(names(pretty_bio), bio_pretty) %>% unlist(),
  Age = "Age",
  Gender = "Gender",
  GHQ_Likert = "GHQ-12 Likert",
  Income = "(Log) Household Income",
  NSSEC5 = "Current NS-SEC",
  Tenure = "Housing Tenure",
  Smoke_W2 = "Smoking Status",
  Alcohol_W2 = "Alcohol Use",
  Education = "Education",
  FatherNSSEC5_14 = "Father's NS-SEC",
  ParentEdu = "Parental Education",
  ParentsHH_14 = "Family Composition",
  Ethnicity = "Ethnicity",
  Foreign = "Birth Country",
  Generation = "Immigrant Generation",
  Wave = "Wave",
  `(Intercept)` = "Constant",
  gender = "Genders",
  allostatic = "Allostatic Load",
  behav = "Health Behaviours", 
  mental = "Mental Health",
  sep = "Current SEP",
  obs = "Observations",
  imps = "Imputations"
) %>%
  make_lbls(df_raw) %>%
  mutate(cat_clean = case_when(cat_clean == "0" ~ "Low Risk",
                               cat_clean == "1" ~ "High Risk",
                               TRUE ~ cat_clean))
save(pretty_lbls, file = "Data/pretty_lbls.Rdata")


# 2. Sample ----
df_s <- read_dta("Data/df_sample.dta") %>%
  add_count(name = "n_nurse") %>%
  filter(Blood_Sample == 1) %>%
  add_count(name = "n_blood") %>%
  filter(between(Age, 25, 64)) %>%
  add_count(name = "n_age") %>%
  filter(!is.na(Blood_Weight_X),
         Blood_Weight_X > 0) %>%
  add_count(name = "n_weight") %>%
  filter(is.na(CReactive) | CReactive <= 10) %>%
  add_count(name = "n_crp") %>%
  left_join(select(df_raw, pidp, GHQ_Likert), by = "pidp") %>%
  mutate(n_ghq = sum(!is.na(GHQ_Likert))) %>%
  select(-c(pidp, Blood_Weight_X, Blood_Sample, Age, GHQ_Likert)) %>%
  drop_na() %>%
  add_count(name = "n_outcome") %>%
  select(matches("^n_")) %>%
  distinct() %>%
  pivot_longer(everything()) %>%
  mutate(p = case_when(row_number() == 1 ~ 100,
                       row_number() < n() ~ value*100/lag(value),
                       row_number() == n() ~ value*100/lag(value, 2)),
         p = round(p, 2),
         value = format(value, big.mark = ",") %>% trimws(),
         string = glue("{value} ({p}%)")) %>%
  select(name, string)

df_s

# 2. Descriptives ----
desc <- list()

# Observed
desc$miss <- df_raw %>% 
  summarise(across(everything(), ~ 100*sum(is.na(.x))/n())) %>%
  pivot_longer(everything(), names_to = "var", values_to = "miss") %>%
  mutate(miss = glue("{round(miss, 2)}%"),
         var = ifelse(var == "Unem_Age", "n", var))

desc$cc_n <- count(df_raw, Unem_Age) %>%
  drop_na() %>%
  mutate(p = 100*n/sum(n),
         n = format(n, big.mark = ",") %>% trimws(),
         p = round(p, 2),
         string = glue("{n} ({p}%)"),
         var = "n", cat = "n") %>%
  select(group_var = Unem_Age, var, cat, string)

desc$cc <- get_desc(df_raw, id_var = "pidp", group_var = "Unem_Age") %>%
  drop_na(group_var) %>%
  select(-miss) %>%
  filter(var != "n") %>%
  bind_rows(desc$cc_n) %>%
  left_join(desc$miss, by = "var")

# Imputed
desc$mi_n <- make_long("Allostatic_Index") %>%
  filter(imp > 0, 
         pidp %in% obs_pidp$Allostatic_Index) %>%
  mutate(Blood_Weight_X = Blood_Weight_X*n()/sum(Blood_Weight_X)) %>%
  count(Unem_Age, wt = Blood_Weight_X/max(imp)) %>%
  mutate(p = 100*n/sum(n),
         n = round(n, 2) %>% format(big.mark = ",") %>% trimws(),
         p = round(p, 2),
         string = glue("{n} ({p}%)"),
         var = "n", cat = "n") %>%
  select(group_var = Unem_Age, var, cat, string)

desc$mi_covar <- make_long("Allostatic_Index") %>%
  filter(imp > 0, 
         pidp %in% obs_pidp$Allostatic_Index) %>%
  select(-dep_var) %>%
  get_desc(id_var = "pidp", imp_var = "imp", weight_var = "Blood_Weight_X",
           group_var = "Unem_Age") %>%
  select(-miss) %>%
  filter(var != "n") %>%
  bind_rows(desc$mi_n)

desc$get_dep <- function(var){
  make_long(var) %>%
    filter(imp > 0,
           pidp %in% obs_pidp[[!!var]]) %>%
    select(pidp, imp, Blood_Weight_X, Unem_Age, Income, FatherEdu, dep_var) %>%
    get_desc(id_var = "pidp", imp_var = "imp", weight_var = "Blood_Weight_X",
             group_var = "Unem_Age") %>%
    filter(var == "dep_var") %>%
    mutate(var = !!var,
           cat = ifelse(cat == "dep_var", !!var, cat))
}  

desc$mi_dep <- str_subset(names(obs_pidp), "GHQ_Likert", TRUE) %>%
  map_dfr(desc$get_dep) %>%
  select(-miss)

# LR Tests
lr_test <- function(var, imp_var = "Allostatic_Index"){
  weights <- make_long(imp_var) %>%
    filter(imp == 1,
           pidp %in% obs_pidp$Allostatic_Index) %>%
    pull(Blood_Weight_X)
  
  form <- paste("Unem_Age ~", var) %>%
    as.formula()
  
  fit1 <- make_long(imp_var) %>%
    filter(imp > 0,
           pidp %in% obs_pidp$Allostatic_Index) %>%
    nest(data = -imp) %>%
    mutate(model = map(data,
                       function(x)   glm(form, family = "binomial",
                                         data = x, weights = weights))) %>%
    pull(model) %>%
    as.mira() %>%
    D3()
  
  p <- summary(fit1)$comparison[["p.value"]]
  
  ifelse(p <= 0.05, "*", "")
}

lr_tests <- tibble(var = names(df_raw)) %>%
  filter(!(var %in% c("Unem_Age", "pidp", "PSU", "Blood_Weight_X"))) %>%
  mutate(imp_var = ifelse(var %in% all_vars, var, "Allostatic_Index"),
         run_var = ifelse(var %in% all_vars, "dep_var", var),
         significant = map2_chr(run_var, imp_var, lr_test)) %>%
  select(-imp_var, -run_var)

save(lr_tests, file = "Data/lr_tests.Rdata")
load("Data/lr_tests.Rdata")
desc$lr_tests <- lr_tests

# Flextable
desc$df <- bind_rows(desc$mi_dep, desc$mi_covar) %>%
  rename(mi = string) %>%
  full_join(rename(desc$cc, cc = string), by = c("group_var", "var", "cat")) %>%
  filter(!(var %in% c("Blood_Weight_X", "PSU", "FatherEdu"))) %>%
  mutate(group_var = ifelse(group_var == levels(df_raw$Unem_Age)[1], "0", "1")) %>%
  left_join(pretty_lbls, by = c("var", "cat" = "desc_cat")) %>%
  left_join(desc$lr_tests, by = "var") %>%
  mutate(mi = ifelse((level == 1 | cat == "1") & group_var == "1" & var != "n",
                     glue("{mi}{significant}"), mi),
         miss = as.character(miss),
         miss = case_when(str_detect(var_clean, "Z-Score") ~ "",
                          str_detect(var_clean, "High-Risk") ~ miss,
                          level > 1 ~ "", 
                          TRUE ~ miss)) %>%
  arrange(index) %>%
  select(group_var, var_clean, cat_clean, cc, mi, miss)

desc$headers <- list(var_clean = "", cat_clean = "Variable",
                     cc_0 = "<6 Months Unemployment", 
                     cc_1 = "6+ Months Unemployment",
                     miss = "% Missing",
                     mi_0 = "<6 Months Unemployment", 
                     mi_1 = "6+ Months Unemployment")
desc$span <- list(var_clean = "", cat_clean = "",
                  cc_0 = "Unweighted Observed Data", 
                  cc_1 = "Unweighted Observed Data",
                  miss = "",
                  mi_0 = "Weighted Imputed Data", 
                  mi_1 = "Weighted Imputed Data")

desc$flx <- desc$df %>%
  pivot_wider(names_from = group_var, values_from = c(mi, cc)) %>%
  select(var_clean, cat_clean, cc_0, cc_1, miss, mi_0, mi_1) %>%
  filter(cat_clean != "Low Risk") %>%
  mutate(cat_clean = case_when(str_detect(cat_clean, "Index") ~ "Index",
                               str_detect(cat_clean, "Z-Score") ~ "Z-Score",
                               TRUE ~ cat_clean),
         var_clean = str_replace(var_clean, " \\((Index|Z\\-Score|High\\-Risk)\\)", ""),
         var_clean = ifelse(cat_clean == var_clean, "", var_clean)) %>%
  make_flx(desc$headers, desc$span)
desc$flx
save_as_docx(desc$flx, path = "Tables/descriptive_statistics.docx")


# 3. Format Regression Results ----
load("Data/main_results.Rdata")
load("Data/mod_specs.Rdata")

format_res <- function(res){
  map(res, list) %>%
    as_tibble_row()
}

df_outcome <- mod_specs %>%
  distinct(outcome) %>%
  left_join(pretty_lbls %>%
              rename_with(~ str_replace(.x, "var", "outcome")) %>%
              select(outcome, outcome_clean, index),
            by = "outcome") %>%
  mutate(type = case_when(str_detect(outcome_clean, "Index") ~ "quartile",
                          str_detect(outcome_clean, "Z-Score") ~ "z",
                          str_detect(outcome_clean, "High-Risk") ~ "quartile"),
         type_clean = case_when(str_detect(outcome_clean, "Index") ~ "Index",
                                str_detect(outcome_clean, "Z-Score") ~ "Z-Score",
                                str_detect(outcome_clean, "High-Risk") ~ "High-Risk"),
         type_axis = case_when(type_clean == "Index" ~ "IRR",
                               type_clean == "Z-Score" ~ "Effect Size",
                               type_clean == "High-Risk" ~ "Odds Ratio"),
         outcome_clean = str_replace(outcome_clean, " \\((Index|Z\\-Score|High\\-Risk)\\)", ""),
         outcome_fct = fct_reorder(outcome_clean, index)) %>%
  select(-index)

mod_cleaner <- c(basic = "Basic Model",
               allostatic = "+ Allostatic Load",
               behav = "+ Health Behaviours",
               sep = "+ Socio-Economic Position",
               all = "+ All Mediators") %>%
  fct_inorder()

main_res <- mod_specs %>%
  mutate(res = map(main_res, format_res)) %>%
  left_join(df_outcome, by = "outcome") %>%
  mutate(mod_clean = mod_cleaner[mod],
         sex_clean = str_to_title(sex))
         
# 4. Figures ----
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

get_col <- function(col){
  main_res %>%
    mutate(res = map(res, ~ .x %>% select(any_of(!!col)))) %>%
    unnest(res) %>%
    unnest(!!col)
}

# Figure 1
fig_1 <- function(type, sample){
  df <- get_col("pool") %>%
    filter(unem_age == FALSE, mod == "basic", term == "a2_0") %>%
    mutate(type_axis = fct_rev(type_axis),
           outcome_fct = fct_rev(outcome_fct)) %>%
    filter(type == !!type, sample == !!sample)
  
  yint <- ifelse(type == "quartile", 1, 0)
  if (type == "quartile"){
    df <- mutate(df, across(c(beta, lci, uci), exp)) 
  }
  
  p <- ggplot(df) +
    aes(x = outcome_fct, y = beta, ymin = lci, ymax = uci) +
    facet_grid(type_axis ~ sex_clean, scales = "free_y", switch = "y", space = "free_y") +
    geom_hline(yintercept = yint) +
    geom_pointrange(color = cbbPalette[4]) +
    coord_flip() +
    labs(x = NULL, y = "Marginal Effect") +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0))
  
  if (type == "quartile"){
    p <- p + scale_y_log10()
  }
  
  return(p)
}

main_res %>%
  drop_na() %>%
  distinct(type, sample) %$%
  map2(type, sample, fig_1)

# Figure 2
fig_2 <- function(type, sample){
  df <- get_col("pool") %>%
    filter(str_detect(outcome, "Allostatic"), term == "a2_0", unem_age == FALSE) %>%
    filter(type == !!type, sample == !!sample)
  
  yint <- ifelse(type == "quartile", 1, 0)
  if (type == "quartile"){
    df <- mutate(df, across(c(beta, lci, uci), exp)) 
  }
  
  p <- ggplot(df) +
    aes(x = sex_clean, y = beta, ymin = lci, ymax = uci, 
        color = mod_clean, shape = mod_clean) +
    geom_hline(yintercept = yint) +
    geom_pointrange(position = position_dodge(0.6)) +
    labs(x = NULL, y = "Marginal Effect", color = NULL, shape = NULL) +
    scale_color_brewer(palette = "Set1") +
    scale_shape_manual(values = 15:19) +
    theme_minimal() +
    theme(legend.position = c(.85, .85),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0))
  
  if (type == "quartile"){
    p <- p + scale_y_log10()
  }
  
  return(p)
}

main_res %>%
  drop_na() %>%
  distinct(type, sample) %$%
  map2(type, sample, fig_2)

# Figure 3
fig_3 <- function(type, sample){
  df <- get_col("margins") %>%
    filter(str_detect(outcome, "Allostatic"), unem_age == TRUE, mod == "basic") %>%
    filter(type == !!type, sample == !!sample)
  
  yint <- ifelse(type == "quartile", 1, 0)
  if (type == "quartile"){
    df <- mutate(df, across(c(beta, lci, uci), exp)) 
  }
  
  p <- ggplot(df) +
    aes(x = Age, y = beta, ymin = lci, ymax = uci) +
    facet_grid(~ sex_clean, scales = "free_y", switch = "y", space = "free_y") +
    geom_hline(yintercept = yint) +
    geom_ribbon(color = NA, alpha = 0.2, fill = cbbPalette[7]) +
    geom_line(color = cbbPalette[7]) +
    labs(x = NULL, y = "Marginal Effect") +
    theme_minimal()
  
  if (type == "quartile"){
    p <- p + scale_y_log10()
  }
  
  return(p)
}

main_res %>%
  drop_na() %>%
  distinct(type, sample) %$%
  map2(type, sample, fig_3)


# Figure 4
fig_4 <- function(type, sample){
  df <- get_col("predict") %>%
    filter(str_detect(outcome, "Allostatic"), unem_age == TRUE, mod == "basic") %>%
    filter(type == !!type, sample == !!sample)
  
  yint <- ifelse(type == "quartile", 1, 0)
  if (type == "quartile"){
    df <- mutate(df, across(c(beta, lci, uci), exp)) 
  }
  
  p <- ggplot(df) +
    aes(x = Age, y = beta, ymin = lci, ymax = uci, 
        color = Unem_Age, fill = Unem_Age) +
    facet_grid(~ sex_clean, scales = "free_y", switch = "y", space = "free_y") +
    geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    scale_color_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = "Predicted Allostatic Load",
         color = NULL, fill = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

main_res %>%
  drop_na() %>%
  distinct(type, sample) %$%
  map2(type, sample, fig_4)

type <- "z"
sample <- "mi"
yint <- 0
imp_var <- "Allostatic_Index"
rm(type, sample, yint)

# Figure 5
fig_5 <- function(imp_var, sample){
  df <- get_col("pool") %>%
    filter(outcome == "GHQ_Likert", term == "a2_0", unem_age == FALSE,
           (mediate == TRUE & mod %in% c("all", "basic")) |
             mediate == FALSE & mod != "all") %>%
    mutate(mod = ifelse(mediate == FALSE | mod == "all", mod, "allostatic"),
           mod_clean = mod_cleaner[mod]) %>%
  filter(imp_var == !!imp_var, sample == !!sample)
  
  p <- ggplot(df) +
    aes(x = sex_clean, y = beta, ymin = lci, ymax = uci, 
        color = mod_clean, shape = mod_clean) +
    geom_hline(yintercept = yint) +
    geom_pointrange(position = position_dodge(0.6)) +
    labs(x = NULL, y = "Marginal Effect", color = NULL, shape = NULL) +
    scale_color_brewer(palette = "Set1") +
    scale_shape_manual(values = 15:19) +
    theme_minimal() +
    theme(legend.position = c(.85, .85),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0))
  
  return(p)
}

main_res %>%
  filter(outcome == "GHQ_Likert") %>%
  distinct(imp_var, sample) %$%
  map2(imp_var, sample, fig_5)

# 5. Specification Curve Analysis ----
# Format Results
load("Data/sca_results.Rdata")

# 4. Plots ----
ggplot(sca) + 
  aes(x = rank, y = beta, color = signif) +
  facet_wrap(~ sex) +
  geom_hline(yintercept = 0) +
  geom_point()

sca %>%
  unnest(vars) %>%
  ggplot() +
  aes(x = rank_std, y = vars, color = signif) +
  facet_wrap(~ sex) +
  geom_jitter(alpha = 0.05)

sca %>%
  unnest(vars) %>%
  ggplot() +
  aes(x = rank_std, y = factor(n_vars), color = signif) +
  facet_wrap(~ sex) +
  geom_jitter(alpha = 0.05)

# Figure 6
# Figure 7 
# Figure 8
   

# FIGURE 1: MAIN EFFECT ESTIMATES; BY Z-SCORE & SAMPLE
# FIGURE 2: MEDIATION OF MAIN EFFECT; BY Z-SCORE
# FIGURE 3: AGE EFFECTS (MARGINS); BY Z-SCORE
# FIGURE 4: AGE EFFECTS (PREDICTION); BY Z-SCORE
# FIGURE 5: MEDIATION
# FIGURE 6: SCA CURVE
# FIGURE 7: SCA BY NUMBER OF VARIABLES
# FIGURE 8: SCA BY PRESENCE OF VARIABLE



