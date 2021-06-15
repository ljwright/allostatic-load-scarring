library(tidyverse)
library(haven)
library(mice)
library(glue)
library(tictoc)
library(parallel)
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

biomarkers <- str_subset(names(df_raw), "_Quartile") %>%
  str_replace("_Quartile", "") %>%
  glue_collapse("|")

df <- df_raw %>%
  select(-matches(biomarkers))
save(df, file = "Data/df_analysis.Rdata")

df_bio <- df_raw %>%
  select(-matches("Allostatic"))
save(df_bio, biomarkers, file = "Data/df_bio.Rdata")


# 2. mice Random Forest ----
m <- 64
cores <- 4
m_core <- ceiling(m/cores)

map_mice <- function(df_mice){
  pred <- make.predictorMatrix(df_mice)
  pred[, c("pidp", "PSU")] <- 0
  
  plan(multisession)
  imp <- future_map(1:cores, 
                    ~ mice(df_mice, method = "rf",
                           predictorMatrix = pred,
                           maxit = 10, m = m_core,
                           seed = .x, printFlag = FALSE))
  future:::ClusterRegistry("stop")
  return(imp)
}

make_long <- function(imps){
  map_dfr(imps, ~ complete(.x, "long", TRUE), .id = "iter") %>%
    as_tibble() %>%
    rename(imp = .imp) %>%
    select(-.id) %>%
    mutate(iter = as.numeric(iter),
           imp = ifelse(imp == 0, 0, 
                        max(imp)*(iter-1)+imp)) %>%
    select(-iter) %>%
    distinct()
}

get_obs <- function(var){
  df %>%
    filter(!is.na({{ var }})) %>%
    pull(pidp)
}

obs_pidp <- biomarkers %>%
  str_split("\\|") %>%
  pluck(1) %>%
  c("Allostatic", .) %>%
  set_names(., .) %>%
  map(~ get_obs(Allostatic_Z)) %>%
  c(., list(GHQ_Likert = get_obs(GHQ_Likert)))

save(m, cores, m_core, map_mice, 
     make_long, obs_pidp, 
     file = "Data/mice_args.Rdata")


# 3. Allostatic Imputations ----
load("Data/mice_args.Rdata")
load("Data/df_analysis.Rdata")
load("Data/df_bio.Rdata")

imp <- list()

tic()
imp$index <-  df %>%
  select(- Allostatic_Z) %>%
  rename(Allostatic = Allostatic_Index) %>%
  map_mice()
toc()

tic()
imp$z <-  df %>%
  select(- Allostatic_Index) %>%
  rename(Allostatic = Allostatic_Z) %>%
  map_mice()
toc()

tic()
imp$bio_q <- df_bio %>%
  select(-matches(glue("^({biomarkers})$"))) %>%
  rename_with(~ str_replace(.x, "_Quartile", "")) %>%
  mutate(across(matches(biomarkers), as.factor)) %>%
  map_mice()
toc()

tic()
imp$bio_z <- df_bio %>%
  select(-matches("_Quartile")) %>%
  map_mice()
toc()

save(imp, make_long, obs_pidp, 
     file = "Data/imp.Rdata")


# 4. Invidual Imputations ----
load("Data/mice_args.Rdata")
load("Data/df_raw.Rdata")

all_vars <- str_subset(names(df_raw), "_Quartile") %>%
  c(., str_replace(., "_Quartile", "")) %>%
  c(., str_subset(names(df_raw), "Allostatic")) %>%
  set_names(., .)

get_obs <- function(var){
  df_raw %>%
    filter(!is.na({{ var }})) %>%
    pull(pidp)
}

obs_pidp <- map(all_vars, ~ get_obs(Allostatic_Z)) %>%
  c(., list(GHQ_Likert = get_obs(GHQ_Likert)))

make_long <- function(var){
  imp_ind[[var]] %>%
    complete("long", TRUE) %>%
    as_tibble() %>%
    rename(imp = .imp) %>%
    select(-.id)
}

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

plan(multisession, workers = 6)
imp_ind <- sample(all_vars) %>%
  future_map(map_mice, .progress = TRUE)
future:::ClusterRegistry("stop")

save(imp_ind, make_long, obs_pidp, 
     file = "Data/imp_ind.Rdata")

# 5. Format Results ----
load("Data/imp.Rdata")

imp_long <- map_dfr(imp, make_long, .id = "imputation")

df_std <- df %>%
  filter(!is.na(Allostatic_Index)) %>%
  select(pidp) %>%
  left_join(imp_long, by = "pidp") %>%
  filter(imputation == "z") %>%
  summarise(mean = wtd.mean(Allostatic, Blood_Weight_X),
            sd = wtd.var(Allostatic, Blood_Weight_X) %>% sd())

imp_long %>%
  mutate(Allostatic = ifelse(imputation == "z",
                             (Allostatic - df_std$mean)/df_std$sd, 
                             Allostatic)) %>%
  group_by(imputation) %>%
  summarise(mean = wtd.mean(Allostatic, Blood_Weight_X),
            sd = wtd.var(Allostatic, Blood_Weight_X) %>% sd())

save(imp_long, file = "Data/imp_long.Rdata")
