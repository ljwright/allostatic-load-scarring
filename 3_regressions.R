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

rm(list = ls())

# 1. Load Data ----
load("Data/df_raw.Rdata")
load("Data/imp_ind.Rdata")

# Age Splines
make_ns <- function(x, deg_free){
  ns(x, deg_free) %>%
    as_tibble() %>%
    rename_with(~ glue("a{deg_free}_{.x}")) %>%
    mutate(across(everything(), as.numeric))
}

age_long <- tibble(Age_C = seq(min(df_raw$Age_C), max(df_raw$Age_C))) %>%
  mutate(make_ns(Age_C, 3),
         make_ns(Age_C, 2)) %>%
  expand_grid(Unem_Age = unique(df_raw$Unem_Age)) %>%
  drop_na() %>%
  mutate(a2_0 = as.numeric(Unem_Age) - 1,
         a2_1 = a2_0*a2_1,
         a2_2 = a2_0*a2_2) %>%
  arrange(Age_C, Unem_Age) %>%
  mutate(Age = Age_C + 25, .after = 1) %>%
  relocate(Unem_Age, .after = 2)

age_specs <- age_long %>%
  select(-Age_C) %>%
  pivot_longer(-c(Age, Unem_Age), names_to = "term")

join_splines <- function(df){
  df %>%
    left_join(age_long, by = c("Age_C", "Unem_Age")) %>%
    select(-Age_C)
}

# 2. Model Arguments ----
# Sex and Model Forms
get_func <- function(outcome, mediate, unem_age, controls){
  base_form <- glue("dep_var ~ a3_1 + a3_2 + a3_3 + a2_0")
  if (outcome == "GHQ_Likert" & mediate == TRUE){
    base_form <- glue("{base_form} + Allostatic")
  }
  if (unem_age ==  TRUE){
    base_form <-  glue("{base_form} + a2_1 + a2_2")
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
  mutate(family = case_when(str_detect(outcome, "Quartile") ~ "binomial",
                            outcome == "Allostatic_Index" ~ "poisson",
                            TRUE ~ "gaussian")) %>%
  filter(outcome == "GHQ_Likert" | mediate == FALSE,
         !(outcome %in% biomarkers & (mod != "basic" | unem_age != FALSE))) %>%
  mutate(id = row_number(), .before = 1)


# 3. Model Functions ----
# Single Model
get_fast <- function(df, mod_form, family){
  mod_form <- mod_form %>% as.formula()
  
  mod_f <- model.frame(mod_form, data = df, weights = Blood_Weight_X)
  x <- model.matrix(mod_form, mod_f)
  y <- model.response(mod_f)
  weights <- model.weights(mod_f)
  
  fastglm(x, y, family = family, weights = weights, method = 3)
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


# Bootstraps
get_result <- function(df, mod_form, family, unem_age){
  res <- list()
  
  # 1. Model
  res$coefs <- get_fast(df, mod_form, family) %>%
    get_coefs()
  
  # 2. Bootstraps
  res$boots <- map_dfr(1:10, 
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
  
  # 4. Margins & Predict
  if (unem_age == TRUE){
    res$margins <- get_margins(res$boots)
    res$predict <- get_predict(res$boots)
  }
  if (family == "binomial"){
    res$predict <- get_predict(res$boots)
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

get_pdiff <- function(boots){
  convert_odds <- function(log_odds){
    odds <- exp(log_odds)
    p <- odds/(1 + odds)
  }
  
  boots %>%
    mutate(estimate = convert_odds(estimate)) %>%
    group_by(Age, boot) %>%
    mutate(estimate = estimate - lag(estimate)) %>%
    slice(2) %>%
    ungroup() %>%
    get_pool()
}

# POOL USING RUBIN'S RULES
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
    summary() %>%
    as_tibble(rownames = "term") %>%
    select(term, beta = 2, lci = 6, uci = 7)
}




# Get Results
get_results <- function(spec_id){
  spec <- mod_specs %>%
    slice(spec_id) %>%
    as.list()
  
  v_imps <- 1:3 # 1:max(imp_long$imp)
  
  mod_form <- get_func(spec$outcome, spec$mediate, 
                       spec$unem_age, mod_forms[[spec$mod]])
  if (spec$sex == "all") mod_form <- glue("{mod_form} + Female")
  if (spec$sample == "cc"){
    mod_form <- str_replace(mod_form, "Wave", "1")
    v_imps <- 0
  }
  
  get_res <- function(imp, var){
    make_long(var) %>%
      filter(imp == !!imp,
             Female %in% sexes[[spec$sex]]) %>%
      join_splines() %>%
      filter(pidp %in% obs_pidp[[var]]) %>%
      get_result(mod_form, spec$family, spec$unem_age)
  }
  
  res <- map(v_imps, get_res, spec$outcome) %>%
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
    res$margins <- get_pool(res$margins)
  }
  if (family == "binomial"){
    res$predict <- get_pdiff(res$predict)
  }
  
  # Observations
  res$obs <- imp_long %>%
    filter(imputation == !!spec$imputation,
           imp == min(v_imps),
           Female %in% sexes[[spec$sex]]) %>%
    semi_join(df_obs, by = "pidp") %>%
    get_fast(mod_form, spec$family) %>%
    pluck("n") %>%
    length()
  
  # Return Object
  res$boots <- NULL
  return(res)
}

# 4. Run Models ----
plan(multisession)
tic()
set.seed(1)
main_res <- mod_specs %>%
  filter(imputation == "index", sex == "all",
         mod == "all", outcome == "Allostatic",
         imp <= 3) %>%
  mutate(res = future_pmap(list(imputation = imputation,
                                imp = imp, outcome = outcome,
                                mediate = mediate, unem_age = unem_age,
                                mod = mod, sex = sex, family = family),
                           get_results, 
                           .progress = TRUE,
                           .options = furrr_options(seed = TRUE))) %>%
  select(id, unem_age, imp, res) %>%
  unnest(res)
toc()
future:::ClusterRegistry("stop")

save(main_res, file = "Data/main_results.Rdata")


chop_res <- function(res){
  
}







# 5. Pool Results ----
mod_ids <- mod_specs %>%
  select(-imp) %>%
  distinct()


# Cofficients
main_res %>%
  select(id, coefs) %>%
  unnest(coefs) %>%
  group_by(id, term) %>%
  summarise(beta = mean(coef),
            .groups = "drop") %>%
  ungroup()

main_res %>%
  select(id, imp, boots) %>%
  unnest(boots) %>%
  group_by(id, imp, term) %>%
  summarise(beta = mean(coef),
            var = var(coef),
            .groups = "drop")

main_res %>%
  select(id, boots) %>%
  unnest(boots) %>%
  group_by(id, term) %>%
  summarise(get_ci(coef),
            .groups = "drop") %>%
  ungroup()

main_res %>%
  select(id, imp, boots) %>%
  unnest(boots)



# GET PREDICT
age_spec <- age_long %>%
  pivot_longer(-c(Age, Unem_Age), names_to = "term") %>%
  nest(spec = -c(Age, Unem_Age))

get_pred <- function(spec, res){
  left_join(res, spec, by = "term") %>%
    mutate(value = case_when(!is.na(value) ~ value,
                             !is.na(v_means[term]) ~ v_means[term],
                             TRUE ~ 0)) %$%
    sum(estimate*value)
}

get_preds <- function(res){
  age_spec %>%
    mutate(estimate = map_dbl(spec, get_pred, res)) %>%
    select(-spec)
}

res_pred <- main_res %>%
  filter(unem_age == TRUE) %>% # SAVE TIME BY FILTERING MORE EFFECTIVELY
  select(id, imp, boots) %>%
  unnest(boots) %>%
  mutate(res = map(res, get_preds)) %>%
  unnest(res) %>%
  group_by(id, Age, Unem_Age) %>%
  summarise(get_ci(estimate),
            .groups = "drop")


# GET MARGINS
mrg_long <- age_long %>%
  filter(Unem_Age == "6+ Months Unemployment") %>%
  select(Age, matches("^a2")) %>%
  pivot_longer(-Age, names_to = "term")

get_ci <- function(x){
  quantile(x, c(.5, .025, .975)) %>%
    as_tibble_row() %>%
    rename(beta = 1, lci = 2, uci = 3)
}

main_res %>%
  filter(unem_age == TRUE) %>%
  select(id, imp, boots)
unnest(boots) %>%
  filter(str_detect(term, "^a2")) %>%
  left_join(mrg_long, by = "term") %>%
  group_by(id, imp, boot, Age) %>%
  summarise(estimate = sum(coef*value),
            .groups = "drop") %>%
  group_by(id, Age) %>%
  summarise(get_ci(estimate),
            .groups = "drop")





main_pool <- main_res  %>%
  mutate(sample = ifelse(imp == 0, "cc", "mi")) %>%
  group_by(imputation, sample, outcome, mediate, unem_age, mod, sex) %>%
  summarise(across(c(coefs, margins, predict), ~ pool_results(.x, family)),
            obs = max(obs), imps = n(), .groups = "drop")
save(main_pool, file = "Data/main_pooled.Rdata")

main_pool %>% 
  unnest(margins) %>%
  select(unem_age, sex, term, beta, lci, uci, sample) %>%
  # separate(term, c("age", "unem"), remove = FALSE, sep = "_") %>%
  mutate(age = as.numeric(term)) %>%
  filter(unem_age == FALSE, sample == "mi") %>%
  ggplot() +
  aes(x = age, y = beta, ymin = lci, ymax = uci) +
  # aes(color = unem, fill = unem) +
  facet_wrap(~ sex) +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line()

library(msm)
?deltamethod

deltamethod(~ exp(x2), coef(mod), vcov(mod))
exp(coef(mod)[2])


mod <- lm(GHQ_Likert ~ Allostatic_Index + Female, df)

qnorm(c(.025, .5, .975), 
      exp(coef(mod)[2]),
      deltamethod(~ exp(x2), coef(mod), vcov(mod)))

mod %>%
  tidy(conf.int = TRUE) %>%
  slice(2) %>%
  mutate(across(where(is.numeric), exp))

x <- map_dbl(1:10000,
             ~ df %>%
               sample_frac(replace = TRUE) %>%
               lm(GHQ_Likert ~ Allostatic_Index + Female, .) %>%
               tidy() %>%
               slice(2) %>%
               pull(estimate) %>%
               exp())
sd(x)
mean(x)
quantile(x, c(.025, .5, .975))
deltamethod(~ exp(x2), coef(mod), vcov(mod))


dsgn <- svydesign(id = ~ PSU,
                  weights = ~ Blood_Weight_X,
                  data = df) 

mod <- svyglm(Allostatic_Index ~ Unem_Age + Age + I(Age^2) + Female,
              family = "poisson", design = dsgn)
tidy(mod, conf.int = TRUE)

mod <- glm(Allostatic_Index ~ Unem_Age + Age + I(Age^2) + Female,
           family = "poisson", df, weights = Blood_Weight_X)
tidy(mod, conf.int = TRUE)


get_glm <- function(x){
  df_x <- sample_frac(df, replace = TRUE)
  
  mod <- glm(Allostatic_Index ~ Unem_Age + Age + I(Age^2) + Female,
             family = "poisson", df_x, weights = df_x$Blood_Weight_X)
  
  coef(mod)[2]
}
tic()
x <- map_dbl(1:5000, get_glm)
toc()

mean(x)
sd(x)
quantile(x, c(.025, .5, .975))


get_fast <- function(index){
  df_x <- sample_frac(df, replace = TRUE)
  
  mod_form <- Allostatic_Index ~ Unem_Age + Age + I(Age^2) + Female
  mod_f <- model.frame(Allostatic_Index ~ Unem_Age + Age + I(Age^2) + Female, 
                       data = df_x, weights = Blood_Weight_X)
  x <- model.matrix(mod_form, mod_f)
  y <- model.response(mod_f)
  weights <- model.weights(mod_f)
  
  mod <- fastglm(x, y, family = "poisson", weights = weights)
  
  coef(mod)[2]
}
tic()
x <- map_dbl(1:5000, get_fast)
toc()

mean(x)
sd(x)
quantile(x, c(.025, .5, .975))