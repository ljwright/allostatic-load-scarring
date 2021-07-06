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

df_raw <- df_raw %>%
  mutate(across(matches("Quartile"), ~ as.numeric(.x) - 1)) %>%
  filter(is.na(Allostatic_Index))

biomarkers <- str_subset(names(df_raw), "_Quartile") %>%
  str_replace("_Quartile", "")


# 2. SCA Arguments ---- # CHANGE MOD FORM
mod_form <- "
Allostatic ~ Unem_Age + Age_C + I(Age_C^2) + I(Age_C^3) +
Foreign + FatherNSSEC5_14 + NonWhite +
Female + ParentEdu + Education + Wave + ParentsHH_14"

sexes <- list(all = c("Male", "Female"),
              male = "Male", female = "Female")

get_lm <- function(vars, sex, quartile){
  if (quartile == TRUE) vars <- glue("{var}_Quartile")
  
  allostatic <- df_raw %>% 
    select(all_of(biomarkers)) %>%
    rowSums() 
  
  df_mod <- df_raw %>%
    mutate(Allostatic = !!allostatic) %>%
    filter(Female %in% sexes[[sex]]) %>%
    mutate(Allostatic = wtd_scale(Allostatic, Blood_Weight_X))
  
  dsgn <- svydesign(id = ~ PSU,
                    weights = ~ Blood_Weight_X,
                    data = df_mod)
  
  if (sex != "all") mod_form <- str_replace(mod_form, "Female", "1")
  
  mod <- mod_form %>%
    as.formula() %>%
    svyglm(family = "gaussian", design = dsgn) %>%
    broom::tidy(conf.int = TRUE) %>%
    filter(str_detect(term, "Unem_Age")) %>%
    select(beta = estimate, p = p.value,
           lci = conf.low, uci = conf.high)
}

# 3. Run Models ----
set.seed(1)
tic()
plan(multisession)
sca <- map(6:length(biomarkers),
           ~ combn(biomarkers, .x, simplify = FALSE)) %>%
  flatten() %>%
  expand_grid(sex = c("all", "male", "female"), quartile = c(TRUE, FALSE), vars = .) %>%
  # group_by(sex) %>% sample_n(100) %>% ungroup() %>%
  mutate(res = future_pmap(list(vars, sex, quartile), get_lm, .progress = TRUE)) %>%
  mutate(n_vars = map_int(vars, length)) %>%
  unnest(res) %>%
  group_by(sex) %>%
  mutate(id = row_number(),
         rank = dense_rank(beta),
         signif = ifelse(p < 0.05, "p < 0.05", "p >= 0.05")) %>%
  ungroup()
future:::ClusterRegistry("stop")
toc()

save(sca, file = "Data/sca_results.Rdata")
