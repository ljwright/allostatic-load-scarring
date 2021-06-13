library(quantreg)
library(tidyverse)
library(glue)
library(Hmisc)
library(tictoc)
library(mice)
library(magrittr)
library(furrr)
library(twilio)

rm(list = ls())

# 1. Load Data ----
send_msg <- function(message = "AWS computations complete"){
  tw_send_message(to = "+447956701440",
                  from = "+12054966817",
                  body = message)
}

load("Data/imp_long.Rdata")
load("Data/df_analysis.Rdata")
imp_long <- imp_long %>%
  filter(between(imp, 1, 30)) %>%
  semi_join(filter(df, !is.na(Allostatic_Index)), by = "pidp")
rm(df)

controls <- c("NonWhite", "Foreign", "Education",
              "ParentEdu", "ParentsHH_14", "Wave", 
              "FatherNSSEC5_14")

# 2. Model Functions ----
get_rq <- function(mod_form, df_reg, tau){
  df_reg <- df_reg %>%
    sample_frac(1, TRUE)
  
  coefs <- rq(mod_form, tau = tau, 
              data = df_reg, weights = Blood_Weight_X) %>%
    coef() %>%
    enframe(name = "cat")
}

boot_rq <- function(imputation, imp, tau, sex){
  sex <- as.character(sex)
  if (sex == "all") controls <- c(controls, "Female")
  
  mod_form <- get_func(controls) %>%
    as.formula()
  
  df_reg <- imp_long %>%
    filter(imputation == !!imputation,
           imp == !!imp,
           Female %in% sexes[[sex]])
  
  map_dfr(1:500, ~ get_rq(mod_form, df_reg, tau), .id = "b") %>%
    mutate(across(where(is.character), factor))
}

get_func <- function(controls){
  base_form <- "Allostatic ~ Unem_Age + Age_C + I(Age_C^2) + I(Age_C^3)" 
  
  mod_form <- glue_collapse(controls, " + ") %>%
    glue(base_form, " + ", .)
}

sexes <- list(male = "Male", female = "Female",
              all = c("Male", "Female"))


# 3. Run Regressions ----
tic()
plan(multisession, workers = 8)
boot_res <- expand_grid(imputation = unique(imp_long$imputation),
                        imp = unique(imp_long$imp),
                        tau = 1:9/10,
                        sex = names(sexes) %>% factor()) %>%
  mutate(coef = future_pmap(list(imputation = imputation, imp = imp,
                                 tau = tau, sex = sex),
                            boot_rq, .progress = TRUE)) %>%
  unnest(coef)
future:::ClusterRegistry("stop")
toc()
send_msg("Bootstraps complete")


# 4. Format Results ----
quantile_coef <- boot_res %>%
  group_by(imputation, sex, cat, tau) %>%
  summarise(value = quantile(value, c(.5, .025, .975)),
            .groups = "drop") %>%
  mutate(stat = rep(c("beta", "lci", "uci"), length.out = n())) %>%
  pivot_wider(names_from = "stat", values_from = "value") %>%
  get_ci()

save(quantile_coef, file = "Data/quantile_results.Rdata")
rm(boot_res)