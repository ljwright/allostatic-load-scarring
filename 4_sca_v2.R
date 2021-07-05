library(tidyverse)
library(mice)
library(glue)
library(survey)
library(furrr)
library(tictoc)
library(Hmisc)
library(gallimaufr)

rm(list = ls())


# 1. Load Data ----
load("Data/df_raw.Rdata")
load("Data/mod_funcs.Rdata")

obs_pidp <- df_raw %>%
  drop_na(Allostatic_Index) %>%
  pull(pidp)

biomarkers <- str_subset(names(df_raw), "_Quartile") %>%
  str_replace("_Quartile", "")

keep_vars <- c("pidp", "Blood_Weight_X", "PSU",
               str_subset(names(age_long), "^a3_"), "a2_0",
               biomarkers, glue("{biomarkers}_Quartile"),
               covariates, "Gender")

df_raw <- df_raw %>%
  mutate(across(matches("Quartile"), ~ as.numeric(.x) - 1)) %>%
  join_splines() %>%
  select(all_of(keep_vars))

rm(keep_vars, get_func, join_splines, age_long)

# 2. Model Functions
get_imp <- function(quartile, vars){
  if (quartile == TRUE) vars <- paste0(vars, "_Quartile")
  
  dep_var <- df_raw %>%
    select(all_of(vars)) %>%
    rowSums()
  
  df_mice <- df_raw %>%
    mutate(dep_var = dep_var) %>%
    select(-all_of(c(biomarkers, glue("{biomarkers}_Quartile"))))
  
  pred <- make.predictorMatrix(df_mice)
  pred[, c("pidp", "PSU", "Gender")] <- 0
  
  imp <- vector(mode = "list", length = 2) %>%
    set_names(c("Male", "Female"))
  
  for (gender in c("Male", "Female")){
    df_g <- df_mice %>%
      filter(Gender == gender)
    
    imp[[gender]] <- mice(df_g, predictorMatrix = pred,
                          maxit = 10, m = 1, seed = 1, 
                          printFlag = FALSE)
  }
  
  df <- map_dfr(imp, complete) %>%
    as_tibble() %>%
    filter(pidp %in% obs_pidp) %>%
    group_by(Gender) %>%
    mutate(dep_var = wtd_scale(dep_var, Blood_Weight_X)) %>%
    ungroup()
  
  return(df)
}

get_lm <- function(sex, df){
  df_mod <- df %>%
    filter(Gender %in% sexes[[sex]])
  
  dsgn <- svydesign(id = ~ PSU,
                    weights = ~ Blood_Weight_X,
                    data = df_mod)
  
  mod_form <- glue_collapse(covariates, ' + ')
  mod_form <- glue("dep_var ~ a2_0 + a3_1 + a3_2 + a3_3 + {mod_form}")
  if (sex == "all") mod_form <- glue("{mod_form} + Gender")
  
  mod <- mod_form %>%
    as.formula() %>%
    svyglm(family = "gaussian", design = dsgn) %>%
    broom::tidy(conf.int = TRUE) %>%
    filter(term == "a2_0") %>%
    mutate(sex = !!sex) %>%
    select(sex, beta = estimate, p = p.value,
           lci = conf.low, uci = conf.high) 
  
  return(mod)
}

# 3. Run Models
get_estimates <- function(spec_id){
  quartile <- mod_specs$quartile[spec_id]
  vars <- mod_specs$vars[[spec_id]]
  
  df <- get_imp(quartile, vars)
  res <- map_dfr(names(sexes), get_lm, df) %>%
    mutate(id = spec_id)
  
  return(res)
}

mod_specs <- map(6:length(biomarkers),
                 ~ combn(biomarkers, .x, simplify = FALSE)) %>%
  flatten() %>%
  expand_grid(quartile = c(TRUE, FALSE), vars = .) %>%
  mutate(id = row_number())

# Run
set.seed(1)
tic()
plan(multisession, workers = 4)
sca_res <- future_map_dfr(sample(mod_specs$id[1:4]), get_estimates, .progress = TRUE)
future:::ClusterRegistry("stop")
toc()

save(sca_res, file = "Data/sca_raw.Rdata")


# Clean
sca <- left_join(sca_res, mod_specs, by = "id") %>%
  mutate(n_vars = map_int(vars, length)) %>%
  group_by(sex) %>%
  mutate(rank = dense_rank(beta),
         signif = ifelse(p < 0.05, "p < 0.05", "p >= 0.05")) %>%
  ungroup()

save(sca, file = "Data/sca_results.Rdata")