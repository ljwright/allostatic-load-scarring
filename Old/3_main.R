library(tidyverse)
library(magrittr)
library(broom)
library(glue)
library(survey)
library(margins)
library(prediction)
library(miceadds)
library(Hmisc)
library(tictoc)
library(furrr)

rm(list = ls())

# 1. Load Data ----
load("Data/df_analysis.Rdata")
load("Data/imp_long.Rdata")

# imp_long <- imp_long %>%
#   filter(imputation == "mice_all")

# 2. Model Arguments ----
# Sex and Model Forms
get_func <- function(outcome, mediate, unem_age, controls){
  base_form <- glue("{outcome} ~ Unem_Age + Age_C + I(Age_C^2) + I(Age_C^3)")
  if (outcome == "GHQ_Likert" & mediate == TRUE){
    base_form <- glue("{base_form} + Allostatic")
  }
  if (unem_age ==  TRUE){
    base_form <-  glue("{base_form} + Unem_Age*Age_C + Unem_Age*I(Age_C^2)")
  }
  
  mod_form <- glue_collapse(controls, " + ") %>%
    glue(base_form, " + ", .)
}


covariates <- c("NonWhite", "Foreign", "Education",
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

# Weighted Means
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

df_m <- imp_long %>%
  filter(imputation == "index",
         imp > 0) %>%
  select(-c(imputation, imp, pidp, PSU))

df_means <- map(df_m, wtd_mean_or_mode,
                df_m$Blood_Weight_X) %>%
  as_tibble()
rm(df_m)

# 3. Model Functions ----
get_result <- function(df, mod_form, family){
  dsgn <- svydesign(id = ~ PSU,
                    weights = ~ Blood_Weight_X,
                    data = df) 
  
  mod <- mod_form %>%
    as.formula() %>%
    svyglm(family = family, design = dsgn)
  
  chop_res <- function(res){
    res %$%
      tibble(beta = set_names(beta, term) %>% list(),
             se = list(se), lci = list(lci), uci = list(uci))
  }
  
  # Coefficients
  coefs <- broom::tidy(mod, conf.int = TRUE) %>%
    select(term, beta = estimate, se = std.error,
           lci = conf.low, uci = conf.high) %>%
    chop_res()
  
  # Margins
  mrgn <- margins(mod, data = df_means,
                  variables = "Unem_Age",
                  at = list(Age_C = 0:39), type = "link",
  ) %>%
    summary() %>%
    as_tibble() %>%
    select(term = Age_C, beta = AME, se = SE,
           lci = lower, uci = upper) %>%
    chop_res()
  
  # Predict
  prdct <- prediction(mod, data = df_means, type = "link",
                      at = list(Age_C = 0:39, Unem_Age = levels(df$Unem_Age))) %>%
    as_tibble() %>%
    mutate(term = glue("{Age_C}_{Unem_Age}"),
           lci = NA) %>%
    select(term, beta = fitted, se = se.fitted) %>%
    mutate(lci = qnorm(0.025, beta, se),
           uci = qnorm(0.975, beta, se)) %>%
    chop_res()
  
  tibble(coefs = list(coefs), 
         margins = list(mrgn),
         predict = list(prdct),
         obs = nobs(mod))
}

get_results <- function(imputation, imp, outcome, mediate, unem_age, 
                        mod, sex, family){
  
  mod_form <- get_func(outcome, mediate, unem_age, mod_forms[[mod]])
  if (sex == "all") mod_form <- glue("{mod_form} + Female")
  if (imp == 0) mod_form <- str_replace(mod_form, "Wave", "1")
  
  df_obs <- df %>%
    rename(Allostatic = Allostatic_Index) %>%
    select(pidp, x = all_of(outcome)) %>%
    filter(!is.na(x))
  
  imp_long %>%
    filter(imputation == !!imputation,
           imp == !!imp,
           Female %in% sexes[[sex]]) %>%
    semi_join(df_obs, by = "pidp") %>%
    get_result(mod_form, family)
}

# 4. Run Models ----
plan(multisession)
tic()
main_res <- imp_long %>%
  select(imputation, imp) %>%
  distinct() %>%
  expand_grid(outcome = c("Allostatic", "GHQ_Likert"),
              mediate = c(FALSE, TRUE),
              unem_age = c(FALSE, TRUE),
              mod = names(mod_forms),
              sex = names(sexes)) %>%
  mutate(family = ifelse(outcome == "GHQ_Likert" | imputation == "z",
                         "gaussian", "poisson")) %>%
  filter(outcome == "GHQ_Likert" | mediate == FALSE) %>%
  mutate(res = future_pmap(list(imputation = imputation,
                                imp = imp, outcome = outcome,
                                mediate = mediate, unem_age = unem_age,
                                mod = mod, sex = sex, family = family),
                           get_results, .progress = TRUE)) %>%
  unnest(res)
toc()
future:::ClusterRegistry("stop")

save(main_res, file = "Data/main_results.Rdata")


# 5. Pool Results ----
pool_results <- function(res, family){
  
  make_exp <- function(df){
    if (family[1] == "poisson"){
      df <- df %>%
        mutate(across(c(beta, lci, uci), exp))
    }
    return(df)
  }
  
  if (length(res) > 1){
    pooled <- res %>%
      bind_rows() %$%
      pool_mi(beta, se = se) %>%
      summary() %>%
      as_tibble(rownames = "term") %>%
      select(term, beta = results, p, 
             lci = `(lower`, uci = `upper)`) %>%
      make_exp() %>%
      list()
    
  } else{
    pooled <- res %>%
      bind_rows() %>%
      rowwise() %>%
      mutate(beta = list(enframe(beta, name = "term", value = "beta")),
             df = tibble(term = beta$term, beta = beta$beta,
                         se = se, lci = lci, uci = uci) %>%
               make_exp() %>%
               list()) %>%
      pull(df)
  }
  return(pooled)
}

main_pool <- main_res  %>%
  mutate(sample = ifelse(imp == 0, "cc", "mi")) %>%
  group_by(imputation, sample, outcome, mediate, unem_age, mod, sex) %>%
  summarise(across(c(coefs, margins, predict), ~ pool_results(.x, family)),
            obs = max(obs), imps = n(), .groups = "drop")
save(main_pool, file = "Data/main_pooled.Rdata")
