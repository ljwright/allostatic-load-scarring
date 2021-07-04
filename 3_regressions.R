library(tidyverse)
library(magrittr)
library(broom)
library(glue)
library(miceadds)
library(Hmisc)
library(tictoc)
library(furrr)
library(splines)
library(fastglm)
library(estimatr)

rm(list = ls())

# 1. Load Data ----
load("Data/df_raw.Rdata")
load("Data/imp.Rdata")

# Age Splines
make_ns <- function(x, deg_free){
  ns(x, deg_free) %>%
    as_tibble() %>%
    rename_with(~ glue("a{deg_free}_{.x}")) %>%
    mutate(across(everything(), as.numeric))
}

age_long <- tibble(Age = seq(min(df_raw$Age), max(df_raw$Age))) %>%
  mutate(make_ns(Age, 3),
         make_ns(Age, 2)) %>%
  expand_grid(Unem_Age = unique(df_raw$Unem_Age)) %>%
  drop_na() %>%
  mutate(a2_0 = as.numeric(Unem_Age) - 1,
         a2_1 = a2_0*a2_1,
         a2_2 = a2_0*a2_2) %>%
  arrange(Age, Unem_Age) %>%
  relocate(Unem_Age, .after = 1)

age_specs <- age_long %>%
  pivot_longer(-c(Age, Unem_Age), names_to = "term")

join_splines <- function(df){
  df %>%
    left_join(age_long, by = c("Age", "Unem_Age"))
}

# 2. Model Arguments ----
# Sex and Model Forms
get_func <- function(outcome, mediate, unem_age, controls){
  base_form <- glue("dep_var ~ a3_1 + a3_2 + a3_3 + a2_0")
  if (outcome == "GHQ_Likert"){
    base_form <- glue("GHQ_Likert ~ a3_1 + a3_2 + a3_3 + a2_0")
    if (mediate == TRUE){
      base_form <- glue("{base_form} + dep_var")
    }
  }
  if (unem_age ==  TRUE){
    base_form <-  glue("{base_form} + a2_1 + a2_2")
  }
  
  mod_form <- glue_collapse(controls, " + ") %>%
    glue(base_form, " + ", .)
}


covariates <- c("Ethnicity", "Foreign", "Education",
                "ParentEdu", "ParentsHH_14", "Wave",
                "FatherNSSEC5_14")
behav <- c("Smoke_W2", "Alcohol_W2")
sep <- c("NSSEC5", "Tenure", "Income")

mod_forms <- list(basic = covariates,
                  behav = c(covariates, behav),
                  sep = c(covariates, sep), 
                  all = c(covariates, sep, behav))
sexes <- list(male = "Male", female = "Female",
              all = c("Male", "Female"))

# Weighted Means and Modes
wtd_mean_or_mode <- function(var, weight){
  if (is.factor(var)){
    wtd.table(var, weight, type = "table") %>%
      enframe() %>%
      arrange(desc(value)) %>%
      slice(1) %>%
      pull(name) %>%
      factor(levels = levels(var))
  }
  else{
    wtd.mean(var, weight)
  }
}

convert_factor <- function(var, name){
  if (is.factor(var)){
    name <- paste0(name, as.character(var))
    var <- 1
  }
  set_names(var, name)
}

df_m <- make_long("Allostatic_Z") %>%
  filter(imp > 0)

v_means <- map(df_m, wtd_mean_or_mode,
               df_m$Blood_Weight_X) %>%
  imap(convert_factor) %>%
  set_names(NULL) %>%
  do.call("c", .)

rm(df_m)

# Model Specificatons
biomarkers <- str_subset(names(df_raw), "_Quartile") %>%
  c(., str_replace(., "_Quartile", ""))

all_vars <- str_subset(names(df_raw), "Allostatic") %>%
  c(biomarkers) %>%
  set_names(., .)

mod_specs <- expand_grid(outcome = c(all_vars, "GHQ_Likert"),
                         mediate = c(FALSE, TRUE),
                         unem_age = c(FALSE, TRUE),
                         mod = names(mod_forms),
                         sex = names(sexes),
                         sample = c("cc", "mi")) %>%
  uncount(ifelse(outcome == "GHQ_Likert", 2, 1), .id = "n") %>%
  mutate(family = case_when(str_detect(outcome, "Quartile") ~ "binomial",
                            outcome == "Allostatic_Index" ~ "poisson",
                            TRUE ~ "gaussian"),
         imp_var = case_when(outcome != "GHQ_Likert" ~ outcome,
                             n == 1 ~ "Allostatic_Index", 
                             n == 2 ~ "Allostatic_Z")) %>%
  filter(outcome == "GHQ_Likert" | mediate == FALSE,
         !(outcome %in% biomarkers & mod != "basic"),
         !(outcome %in% c(biomarkers, "GHQ_Likert") & unem_age == TRUE)) %>%
  mutate(id = row_number(), .before = 1) %>%
  select(-n)

save(all_vars, mod_specs, file = "Data/mod_specs.Rdata")


# 3. Model Functions ----
# Single Model
get_fast <- function(df, mod_form, family){
  mod_form <- mod_form %>% as.formula()
  
  if (family == "gaussian"){
    
    mod <- lm_robust(mod_form, df, Blood_Weight_X)
    
  } else{
    
    mod_f <- model.frame(mod_form, data = df, weights = Blood_Weight_X)
    x <- model.matrix(mod_form, mod_f)
    y <- model.response(mod_f)
    if (is.factor(y)){
      y <- as.numeric(y) - 1
    }
    weights <- model.weights(mod_f)
    
    mod <- fastglm(x, y, family = family, weights = weights, method = 3)
    
  }
  
  return(mod)
}

get_coefs <- function(mod){
  enframe(coef(mod), name = "term", value = "estimate")
}

# MARGINS
get_margins <- function(boots){
  age_specs %>%
    filter(Unem_Age == "6+ Months Unemployment") %>%
    filter(str_detect(term, "^a2")) %>%
    inner_join(boots, by = "term") %>%
    group_by(Age, Unem_Age, boot) %>%
    summarise(estimate = sum(value*estimate),
              .groups = "drop") 
}

# PREDICT
get_predict <- function(boots){
  age_specs %>%
    distinct(Age, Unem_Age) %>%
    expand_grid(boots) %>%
    left_join(age_specs, by = c("Age", "Unem_Age", "term")) %>%
    mutate(value = case_when(!is.na(value) ~ value,
                             !is.na(v_means[term]) ~ v_means[term],
                             TRUE ~ 0)) %>%
    group_by(Age, Unem_Age, boot) %>%
    summarise(estimate = sum(value*estimate),
              .groups = "drop") 
}

# Transformed Margins
get_transformed <- function(boots, family){
  get_predict(boots) %>%
    mutate(family = !!family,
           estimate = exp(estimate),
           estimate = ifelse(family == "binomial",
                             estimate/(1+estimate),
                             estimate)) %>%
    select(-family) %>%
    arrange(Unem_Age) %>%
    group_by(Age, boot) %>%
    mutate(estimate = estimate - lag(estimate)) %>%
    slice(2) %>%
    ungroup()
}


# Bootstraps
get_result <- function(df, mod_form, family, unem_age){
  res <- list()
  
  # 1. Model
  res$coefs <- get_fast(df, mod_form, family) %>%
    get_coefs()
  
  # 2. Bootstraps
  res$boots <- map_dfr(1:500, 
                       ~ df %>%
                         sample_frac(replace = TRUE) %>%
                         get_fast(mod_form, family) %>%
                         get_coefs(),
                       .id = "boot") %>%
    mutate(boot = as.integer(boot))
  
  # 3. Sample Variance
  res$rubin <- res$boots %>%
    group_by(term) %>%
    summarise(se = sd(estimate)) %>%
    full_join(res$coefs, ., by = "term")
  res$coefs <- NULL
  
  # 4. Margins & Predict
  if (unem_age == TRUE){
    res$margins <- get_margins(res$boots)
    res$predict <- get_predict(res$boots)
  }
  if (family %in% c("binomial", "poisson")){
    res$transform <- get_transformed(res$boots, family)
  }
  
  # 4. Return
  return(res)
}

# Clean Results
# GET CI
get_ci <- function(x){
  quantile(x, c(.5, .025, .975)) %>%
    as_tibble_row() %>%
    rename(beta = 1, lci = 2, uci = 3)
}

get_pool <- function(boots){
  boots %>%
    group_by(Age, Unem_Age) %>%
    summarise(get_ci(estimate),
              .groups = "drop")
}

# POOL USING RUBIN'S RULES
tidy_pool <- function(object, alpha = 0.05){
  crit <- stats::qt(alpha/2, object$df, lower.tail = FALSE)
  
  tibble(term = names(object$qbar),
         beta = object$qbar, 
         se = sqrt(diag(object$t)),
         lci = beta - crit * se,
         uci = beta + crit * se) %>%
    select(-se)
}


pool_rubin <- function(rubin){
  rubin %>%
    arrange(imp, term) %>%
    nest(res = -imp) %>%
    mutate(res = map(res, 
                     ~  .x %$%
                       tibble(beta = set_names(estimate, term) %>% list(),
                              se = list(se)))) %>%
    unnest(res) %$%
    pool_mi(beta, se = se) %>%
    tidy_pool()
}

# Get Results
get_results <- function(spec_id){
  print(spec_id)
  
  spec <- mod_specs %>%
    slice(spec_id) %>%
    as.list()
  
  v_imps <- 1:imp[[1]]$m
  
  mod_form <- get_func(spec$outcome, spec$mediate, 
                       spec$unem_age, mod_forms[[spec$mod]])
  if (spec$sex == "all") mod_form <- glue("{mod_form} + Gender")
  if (spec$sample == "cc"){
    mod_form <- str_replace(mod_form, "Wave", "1")
    v_imps <- 0
  }
  
  get_df <- function(imp, imp_var, dep_var){
    make_long(imp_var) %>%
      filter(imp == !!imp,
             Gender %in% sexes[[spec$sex]]) %>%
      join_splines() %>%
      filter(pidp %in% obs_pidp[[!!dep_var]])
  }
  
  get_res <- function(imp){
    get_df(imp, spec$imp_var, spec$outcome) %>%
      get_result(mod_form, spec$family, spec$unem_age)
  }
  
  res <- map(v_imps, get_res) %>%
    transpose() %>%
    map(bind_rows, .id = "imp")
  
  # Add Pooled
  res$pool <- res$boots %>%
    group_by(term) %>%
    summarise(get_ci(estimate),
              .groups = "drop")
  
  if (spec$sample == "cc"){
    res$rubin <- res$rubin %>%
      uncount(2, .id = "imp")
  }
  res$rubin <- pool_rubin(res$rubin)
  
  # Add Margins, Predict if sensible
  if (spec$unem_age){
    res$predict <- get_pool(res$predict)
    res$margins <- get_pool(res$margins) %>%
      select(-Unem_Age)
  }
  if (spec$family %in% c("binomial", "poisson")){
    res$transform <- get_pool(res$transform) %>%
      select(-Unem_Age)
  }
  
  # Observations
  res$obs <- get_df(min(v_imps), spec$imp_var, spec$outcome) %>%
    get_fast(mod_form, spec$family)
  res$obs <- ifelse(spec$family == "gaussian", res$obs$nobs, length(res$obs$n))
  
  # Return Object
  res$boots <- NULL
  gc()
  
  return(res)
}


# 4. Run Models ----
plan(multisession, workers = 2)
Sys.time()
tic()
set.seed(1)
main_res <- mod_specs %>%
  sample_frac() %>%
  mutate(res = future_map(id, get_results,
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)))
toc()
future:::ClusterRegistry("stop")

save(main_res, file = "Data/main_results.Rdata")

Sys.time()
tic()
main_res <- map(mod_specs$id, get_results)
toc()

save(main_res, file = "Data/main_results.Rdata")
