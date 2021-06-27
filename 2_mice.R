library(tidyverse)
library(haven)
library(mice)
library(glue)
library(tictoc)
library(furrr)
library(Hmisc)

rm(list = ls())

# 1. Load data ----
df_raw <- read_dta("Data/df_analysis.dta") %>%
  as_factor() %>%
  zap_formats() %>%
  zap_label() %>%
  mutate(across(matches("_Quartile"), as.factor))
save(df_raw, file = "Data/df_raw.Rdata")


# 2. mice Arguments ----
all_vars <- str_subset(names(df_raw), "_Quartile") %>%
  c(., str_replace(., "_Quartile", "")) %>%
  c(., str_subset(names(df_raw), "Allostatic")) %>%
  set_names(., .)

map_mice <- function(var){
  df_mice <- df_raw %>%
    rename(dep_var = all_of(!!var)) %>%
    select(-any_of(!!all_vars))
  
  pred <- make.predictorMatrix(df_mice)
  pred[, c("pidp", "PSU")] <- 0
  
  
  imp <- mice(df_mice, method = "rf",
              predictorMatrix = pred,
              maxit = 10, m = 60,
              seed = 1, printFlag = FALSE)
  return(imp)
}

get_obs <- function(var){
  df_raw %>%
    filter(!is.na({{ var }})) %>%
    pull(pidp)
}

obs_pidp <- map(all_vars, ~ get_obs(Allostatic_Z)) %>%
  c(., list(GHQ_Likert = get_obs(GHQ_Likert)))

make_long <- function(var){
  imp[[var]] %>%
    complete("long", TRUE) %>%
    as_tibble() %>%
    rename(imp = .imp) %>%
    select(-.id)
}


# 3. Individual Imputations ----
tic()
plan(multisession, workers = 4)
imp <- sample(all_vars) %>%
  future_map(map_mice, .progress = TRUE)
future:::ClusterRegistry("stop")
toc()

save(imp, make_long, obs_pidp,
     file = "Data/imp.Rdata")
