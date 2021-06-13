library(tidyverse)
library(mice)
library(glue)
library(survey)
library(furrr)
library(tictoc)
library(Hmisc)

rm(list = ls())


# 1. Load Data ----
load("Data/df_raw.Rdata")

# MICE RF ADDING ALL BIOMARKERS TOGETHER?
biomarkers <- str_subset(names(df_raw), "_Quartile") %>%
  str_replace("_Quartile", "")

df <- df_raw %>%
  filter(!is.na(Allostatic_Index))


# 2. SCA Arguments ----
mod_form <- "
Allostatic ~ Unem_Age + Age_C + I(Age_C^2) + I(Age_C^3) +
Foreign + FatherNSSEC5_14 + NonWhite +
Female + ParentEdu + Education + Wave + ParentsHH_14"

sexes <- list(all = c("Male", "Female"),
              male = "Male", female = "Female")

get_lm <- function(vars, sex, quartile){
  if (quartile == TRUE) vars <- glue("{var}_Quartile")
  
  df_mod <- df %>%
    filter(Female %in% sexes[[sex]]) %>%
    rowwise() %>%
    mutate(Allostatic = sum(c_across(all_of(vars)))) %>%
    ungroup()
  
  wtd_sd <- wtd.var(df_mod$Allostatic, 
                    df_mod$Blood_Weight_X) %>% sqrt()
  
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
           lci = conf.low, uci = conf.high) %>%
    mutate(sd = wtd_sd)
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
  mutate(res = future_pmap(list(vars = vars, sex = sex, quartile = quartile), get_lm, .progress = TRUE)) %>%
  mutate(n_vars = map_int(vars, length)) %>%
  unnest(res) %>%
  group_by(sex) %>%
  mutate(id = row_number(),
         rank = dense_rank(beta),
         rank_std = dense_rank(beta/sd),
         rank_range = dense_rank(beta/n_vars),
         signif = ifelse(p < 0.05, "p < 0.05", "p >= 0.05")) %>%
  ungroup()
future:::ClusterRegistry("stop")
toc()

save(sca, file = "Data/sca_results.Rdata")

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
