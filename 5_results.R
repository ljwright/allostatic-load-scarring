library(tidyverse)
library(glue)
library(mice)
library(Hmisc)
library(officer)
library(flextable)
library(magrittr)
library(haven)

rm(list = ls())


# 1. Load Data ----
load("Data/df_analysis.Rdata")
load("Data/imp_long.Rdata")

df <- df %>%
  mutate(Age = Age_C + 25,
         Allostatic = Allostatic_Index)

imp_long <- imp_long %>%
  filter(imputation == "index",
         imp > 0) %>%
  semi_join(filter(df, !is.na(Allostatic_Index)), by = "pidp") %>%
  mutate(Age = Age_C + 25) %>%
  arrange(imputation, imp, pidp)

pretty_lbls <- c(
  n = "n",
  Unem_Age = "Youth Unemployment",
  "Unem_Age6+ Months Unemployment:Age_C" = "Youth Unemployment x Age",
  "Unem_Age6+ Months Unemployment:I(Age_C^2)" = "Youth Unemployment x Age^2",
  Allostatic = "Allostatic Load (Index)",
  Age = "Age",
  Age_C = "Age",
  "I(Age_C^2)" = "Age^2",
  "I(Age_C^3)" = "Age^3",
  Female = "Gender",
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
  NonWhite = "Ethnicity",
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
  enframe(name = "var", value = "var_clean") %>%
  mutate(cat = map(var, ~ paste0(.x,levels(df[[.x]])))) %>%
  unnest(cols = cat) %>%
  mutate(
    cat_clean = ifelse(var==cat, 
                       var_clean,
                       str_replace(cat, var, "")),
    cat_clean = case_when( var == "NonWhite" & cat_clean == "No" ~ "White",
                           var == "NonWhite" & cat_clean == "Yes" ~ "Non-White",
                           var == "Foreign" & cat_clean == "No" ~ "UK-born",
                           var == "Foreign" & cat_clean == "Yes" ~ "Foreign-born",
                           var == "Wave" & cat_clean == "2" ~ "Wave 2",
                           var == "Wave" & cat_clean == "3" ~ "Wave 3",
                           TRUE ~ cat_clean),
    index = row_number()) %>%
  group_by(var) %>%
  mutate(level = row_number(), levels = n()) %>%
  ungroup()
save(pretty_lbls, file = "Data/pretty_lbls.Rdata")

# 2. Sample ----
df_s <- read_dta("Data/df_sample.dta")

df_s %>%
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
  left_join(select(df, pidp, GHQ_Likert), by = "pidp") %>%
  mutate(n_ghq = sum(!is.na(GHQ_Likert))) %>%
  select(-c(pidp, Blood_Weight_X, Blood_Sample, Age, GHQ_Likert)) %>%
  drop_na() %>%
  add_count(name = "n_outcome") %>%
  select(matches("^n_")) %>%
  distinct()

bio_clean <- c(
  Pulse = "Pulse rate", Systolic = "Sytolic blood pressure",
  Diastolic = "Diastolic blood pressure", WtH_Ratio = "Waist-to-Height Ratio",
  Insulin = "IGF-1", Creatinine = "Creatinine Clearance Rate",
  DHEAS = "DHEA-S", ClaussFib = "Fibrinogen", CReactive = "C-reactive protein",
  HDL_Ratio = "Cholesterol-to-HDL ratio", Trig = "Triglycerides",
  HbA1c = "HbA1c")
save(bio_clean, file = "Data/bio_clean.Rdata")

df_sample <- df_s %>%
  filter(Blood_Sample == 1,
         between(Age, 25, 64),
         Blood_Weight_X > 0,
         is.na(CReactive) | CReactive <= 10) %>%
  select(-c(pidp, Blood_Weight_X, Blood_Sample, Age))
  
df_pat <- md.pattern(df_sample, plot = FALSE) %>%
  as_tibble(rownames = "n") %>%
  mutate(n = as.numeric(n)) %>%
  arrange(desc(n)) %>%
  slice(1:10) %>%
  mutate(id = row_number()) %>%
  select(-V13) %>%
  pivot_longer(-c(n, id)) %>%
  mutate(prop = 100*n/nrow(df_sample),
         n = format(n, big.mark = ","),
         string = glue("{n}\n({round(prop, 2)}%)"),
         name = bio_clean[name],
         value = ifelse(value == 0, "Missing", "Observed") %>%
           fct_rev())

df_pat %>%
  ggplot() +
  aes(x = name, y = id, fill = value) +
  geom_tile(color = "grey50") +
  coord_flip() +
  scale_y_continuous(breaks = 1:10, 
                     labels = df_pat %>% 
                       select(id, string) %>% 
                       distinct() %>% 
                       pull(string)) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  theme_minimal() +
  labs(x = NULL, y = NULL, fill = NULL)
ggsave("Images/biomarker_missingness.png",
       height = 9.9, width = 21, units = "cm")  
  

sample_tbl <- df_sample %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(missing = 100*mean(is.na(value)),
            .groups = "drop") %>%
  mutate(name = bio_clean[name],
         missing = glue("{round(missing, 2)}%")) %>%
  flextable() %>%
  set_header_labels(name = "Biomarker", missing = "Missing") %>%
  fontsize(size = 11, part = "all") %>%
  align(j = 1, align = "right", part = "all") %>%
  align(j = 2, align = "center", part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit()
sample_tbl
save_as_docx(sample_tbl, path = "Tables/missing_biomarkers.docx")

rm(sample_tbl, df_s, df_sample, df_pat)


# 3. Descriptive Statistics ----
get_descriptives <- function(df, norm = 1, weight = FALSE){
  if (weight == FALSE) df <- mutate(df, Blood_Weight_X = 1)
  df <- df %>%
    filter(Blood_Weight_X > 0) %>%
    mutate(Blood_Weight_X = Blood_Weight_X*n()/sum(Blood_Weight_X)) %>%
    rename(unem = Unem_Age)
  
  df_cat <- df %>%
    select(-matches("^Unem_")) %>%
    select(unem, where(is.factor), Blood_Weight_X) %>%
    pivot_longer(-c(unem, Blood_Weight_X)) %>%
    filter(!is.na(value)) %>%
    count(unem, name, value, wt = Blood_Weight_X) %>%
    group_by(unem, name) %>%
    mutate(prop = 100*n/sum(n),
           n = n/norm,
           across(c(n, prop), round, 2),
           n = format(n, big.mark = ",") %>%
             trimws(),
           string = glue("{n} ({prop}%)"),
           cat = glue("{name}{value}")) %>%
    ungroup() %>%
    select(unem, cat, string)
  
  df_cont <- df %>%
    select(unem, where(is.numeric), Blood_Weight_X) %>%
    pivot_longer(-c(unem, Blood_Weight_X)) %>%
    filter(!is.na(value)) %>%
    group_by(unem, name) %>%
    summarise(mean = wtd.mean(value, Blood_Weight_X),
              sd = wtd.var(value, Blood_Weight_X) %>% sqrt(),
              .groups = "drop") %>%
    mutate(across(c(mean, sd), round, 2),
           string =glue("{mean} ({sd})"),
           value = name) %>%
    select(unem, cat = name, string)
  
  df_n <- df %>%
    count(unem, wt = Blood_Weight_X) %>%
    filter(!is.na(unem)) %>%
    mutate(prop = 100*n/sum(n),
           n = n/norm,
           across(c(n, prop), round, 2),
           n = format(n, big.mark = ",") %>%
             trimws(),
           string = glue("{n} ({prop}%)"),
           cat = "n") %>%
    select(unem, cat, string)
  
  bind_rows(df_cat, df_cont, df_n) %>%
    ungroup() %>%
    filter(!is.na(unem))
}

desc <- list()
desc$cc <- get_descriptives(df)
desc$mi <- get_descriptives(imp_long, max(imp_long$imp), TRUE)

lr_test <- function(var){
  weights <- imp_long %>%
    filter(imp == 1) %>%
    pull(Blood_Weight_X)
  
  form <- paste("Unem_Age ~", var) %>%
    as.formula()
  
  fit1 <- imp_long %>%
    nest(data = -imp) %>%
    mutate(model = map(data,
                       function(x)   glm(form, family = "binomial",
                                         data = x, weights = weights))) %>%
    pull(model) %>%
    as.mira() %>%
    D3()
}
lr_tests  <- imp_long %>%
  select(-c(imp, Unem_Age, Blood_Weight_X, imputation, PSU, pidp)) %>%
  names() %>% structure(., names = .) %>%
  map(lr_test) %>%
  map_dfr(~ summary(.x)$comparison, .id = "var") %>%
  mutate(significant = ifelse(p.value<=0.05, "*", "")) %>%
  select(var, significant)
save(lr_tests, file = "Data/lr_tests.Rdata")
load("Data/lr_tests.Rdata")
desc$lr_tests <- lr_tests

desc$missing <- df %>%
  summarise(across(everything(), 
                   ~ 100*sum(is.na(.x)/n()))) %>%
  pivot_longer(everything(), names_to = "var", values_to = "prop") %>%
  mutate(prop_missing = glue("{round(prop, 2)}%")) %>%
  select(var, prop_missing)

desc$table <- bind_rows(cc = desc$cc, mi = desc$mi, .id = "data") %>%
  filter(cat != "Age_C") %>%
  inner_join(pretty_lbls, by = "cat") %>%
  arrange(index, level) %>%
  left_join(desc$lr_tests, by = "var") %>%
  left_join(desc$missing, by = "var") %>%
  mutate(unem = ifelse(unem == levels(df$Unem_Age)[1],
                       glue("unem1_{data}"), glue("unem2_{data}")),
         prop_missing = ifelse(level>1, "", prop_missing),
         significant = ifelse(level>1, "", significant)) %>%
  select(-data) %>%
  pivot_wider(names_from = "unem", values_from = "string") %>%
  mutate(unem2_mi = ifelse(level == 1, glue("{unem2_mi}{significant}", .na = ""), unem2_mi),
         var_clean = ifelse(var_clean == cat_clean, "", var_clean)) %>%
  select(var_clean, cat_clean, unem1_cc, unem2_cc, prop_missing, unem1_mi, unem2_mi)


desc$flextable <- flextable(desc$table) %>%
  merge_v(j = ~ var_clean) %>%
  set_header_labels(var_clean = "", cat_clean = "Variable", unem1_cc = "<6 + Months\nUnemployment",
                    unem2_cc = "6+ Months\nUnemployment", prop_missing = "% Missing",
                    unem1_mi = "<6 + Months\nUnemployment", unem2_mi = "6+ Months\nUnemployment") %>%
  border_remove() %>%
  add_header(var_clean = "", cat_clean = "", unem1_cc = "Unweighted Observed Data",
             unem2_cc = "Unweighted Observed Data",
             unem1_mi = "Weighted Imputed Data",  unem2_mi = "Weighted Imputed Data") %>%
  fontsize(size = 11, part = "all") %>%
  merge_h(part = "header") %>%
  border_inner_h(border = fp_border(color="gray30", style = "dashed")) %>%
  hline_top(border = fp_border(color="black", width = 2), part = "all") %>%
  hline_bottom(border = fp_border(color="black", width = 2), part = "all") %>%
  fix_border_issues(part = "all") %>%
  align(j=3:7, align="center", part = "all") %>%
  valign(j = 1, valign = "top") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit()
desc$flextable
save_as_docx(desc$flextable, path = "Tables/descriptive_statistics.docx")


# 3. Main Results ----
load("Data/main_pooled.Rdata")

mod_clean <- c(basic = "Basic Model",
               behav = "+ Health Behaviours",
               sep = "+ Socio-Economic Position",
               all = "+ All Mediators")

clean_coefs <- function(coefs, family){
  clean_coef <- function(x, cat, family){
    case_when(str_detect(cat, "Age_C") & family == "poisson" ~ x^39,
              str_detect(cat, "Age_C") & family != "poisson" ~ x*39,
              TRUE ~ x) %>%
      round(2)
  }
  
  coefs %>%
    rename(cat = term) %>%
    mutate(across(c(beta, lci, uci), 
                  ~ clean_coef(.x, cat, family)))
}

main_pool <- main_pool %>%
  mutate(mod_clean = factor(mod_clean[mod], mod_clean),
         allostatic = ifelse(imputation == "index", "Index", "Z-Score"),
         Sex = str_to_title(sex),
         imps = ifelse(imps == 1, "-", imps),
         margins = map(margins,
                       ~ mutate(.x, Age = as.numeric(term) + 25)),
         coefs = map2(coefs, family, clean_coefs),
         predict = map(predict, 
                       ~ separate(.x, term, c("Age", "unem"), "_") %>%
                         mutate(Age = as.numeric(Age) + 25))
  )

# Plots
allostatic_plot <- list()

allostatic_plot$make <- function(func){
  main_pool %>%
    select(sample, imputation) %>%
    distinct() %$%
    map2(sample, imputation, 
         allostatic_plot[[func]])
}

allostatic_plot$main <- function(sample, imputation){
  if (imputation == "index"){
    yintercept <- 1
    ylab <- "Marginal Effect (IRR)"
  } else{
    yintercept <- 0
    ylab <- "Marginal Effect (SD)"
  }
  
  p <- main_pool %>%
    filter(outcome == "Allostatic", unem_age == FALSE,
           sample == !!sample, imputation == !!imputation) %>%
    unnest(coefs) %>%
    filter(str_detect(cat, "^Unem")) %>%
    ggplot() +
    aes(x = Sex, y = beta, ymin = lci, ymax = uci,
        color = mod_clean, shape = mod_clean) +
    geom_hline(yintercept = yintercept) + 
    geom_pointrange(position = position_dodge(0.8)) +
    scale_color_brewer(palette = "Set1") +
    scale_shape_manual(values = 15:19) +
    labs(x = NULL, color = NULL, shape = NULL,
         y = ylab) +
    theme_minimal() +
    theme(legend.position = c(.85,.85))
  
  glue("Images/allostatic_main_{sample}_{imputation}.png") %>%
    ggsave(plot = p, dpi = 300,  width = 21, 
           height = 9.9, units = "cm")  
  
  return(p)
}
allostatic_plot$make("main")


allostatic_plot$age <- function(df, col_var, palette = "Set1",
                      h_line = 1,
                      y_lab = "Marginal Effect (IRR)",
                      col_lab = NULL,
                      leg_pos = c(.85, .85), nrow = 1){
  ggplot(df) +
    aes(x = Age, y = beta, ymin = lci, ymax = uci,
        color = {{ col_var }}, fill = {{ col_var }},
        linetype = {{ col_var }}) +
    facet_wrap( ~ Sex, nrow = nrow) +
    geom_hline(yintercept = h_line) +
    geom_line() +
    geom_ribbon(color = NA, alpha = 0.2) +
    scale_x_continuous(breaks = 2.5:6.5*10) +
    scale_color_brewer(palette = palette) +
    scale_fill_brewer(palette = palette) +
    labs(x = "Age", y = y_lab,
         color = col_lab, fill = col_lab,
         linetype = col_lab) +
    theme_minimal() +
    theme(legend.position = leg_pos)
}


allostatic_plot$age_basic <- function(sample, imputation){
  if (imputation == "index"){
    h_line <- 1
    y_lab <- "Marginal Effect (IRR)"
  } else{
    h_line <- 0
    y_lab <- "Marginal Effect (SD)"
  }
  
  p <- main_pool %>%
    unnest(margins) %>%
    filter(outcome == "Allostatic", imputation == !!imputation,
           sample == !!sample, mod == "basic", unem_age == TRUE) %>%
    allostatic_plot$age(Sex, h_line = h_line, y_lab = y_lab) +
    aes(linetype = NULL) +
    guides(color = FALSE, fill = FALSE)
  
  glue("Images/allostatic_age_{sample}_{imputation}.png") %>%
    ggsave(plot = p, dpi = 300,  width = 21, 
           height = 9.9, units = "cm")  
  
  return(p)
}
allostatic_plot$make("age_basic")

allostatic_plot$age_all <- function(sample, imputation){
  if (imputation == "index"){
    h_line <- 1
    y_lab <- "Marginal Effect (IRR)"
  } else{
    h_line <- 0
    y_lab <- "Marginal Effect (SD)"
  }
  
  p <- main_pool %>%
    unnest(margins) %>%
    filter(sample == !!sample, imputation == !!imputation,
           unem_age == TRUE, outcome == "Allostatic") %>%
    allostatic_plot$age(mod_clean, leg_pos = "bottom", 
                        nrow = 3, col_lab = "Model",
                        h_line = h_line, y_lab = y_lab)
  
  glue("Images/allostatic_age_all_{sample}_{imputation}.png") %>%
    ggsave(plot = p, dpi = 300,  width = 21, 
           height = 29.7, units = "cm")  
  
  return(p)
}
allostatic_plot$make("age_all")

allostatic_plot$age_predict <- function(sample, imputation){
  p <- main_pool %>%
    unnest(predict) %>%
    filter(sample == !!sample, imputation == !!imputation,
           unem_age == TRUE, outcome == "Allostatic", mod == "basic") %>%
    allostatic_plot$age(unem, "Dark2", NULL, "Predicted Allostatic Load")
  
  glue("Images/allostatic_predicted_{sample}_{imputation}.png") %>%
    ggsave(plot = p, dpi = 300,  width = 21, 
           height = 9.9, units = "cm")  
  
  return(p)
}
allostatic_plot$make("age_predict")

allostatic_plot$age_compare <- function(sample, imputation){
  h_line <- as.numeric(imputation == "index")
  
  main_basic <- main_pool %>%
    unnest(margins) %>%
    filter(sample == !!sample, imputation == !!imputation, 
           mod == "basic", unem_age == TRUE, outcome == "Allostatic") %>%
    select(-mod_clean)
  
  p <- main_pool %>%
    unnest(margins) %>%
    filter(sample == !!sample, imputation == !!imputation, 
           mod != "basic", unem_age == TRUE, outcome == "Allostatic") %>%
    ggplot() +
    aes(x = Age, y = beta, ymin = lci, ymax = uci, group = mod) +
    facet_grid(mod_clean ~ Sex, switch = "y") +
    geom_hline(yintercept = h_line) +
    geom_line(data = main_basic, color = "grey50", linetype = "dashed") +
    geom_line(aes(color = mod_clean)) +
    geom_ribbon(aes(color = mod_clean, fill = mod_clean),
                color = NA, alpha = 0.2) +
    scale_x_continuous(breaks = 2.5:6.5*10) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Age", y = NULL,
         color = NULL, fill = NULL) +
    guides(color = FALSE, fill = FALSE) +
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.spacing.y = unit(1, "lines"),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0))
  
  glue("Images/allostatic_age_compare_{sample}_{imputation}.png") %>%
    ggsave(plot = p, dpi = 300,  width = 21, 
           height = 21, units = "cm")  
  
  return(p)
}
allostatic_plot$make("age_compare")


# Mediation Plots
med_clean <- c(basic = "Basic model", allostatic = "+ Allostatic",
             behav = "+ Health Behaviours", sep = "+ Socio-Economic Position", 
             all = "+ All Mediators")

med_res <- main_pool %>%
  filter(outcome == "GHQ_Likert") %>%
  mutate(med_mod = case_when(mediate == TRUE & mod == "basic" ~ "allostatic",
                             mediate == TRUE & mod == "all" ~ "all",
                             mediate == FALSE & mod != "all" ~ mod),
         med_clean = factor(med_clean[med_mod], med_clean)) %>%
  filter(!is.na(med_mod))

med_plot <- list()

med_plot$make <- function(func){
  med_res %>%
    select(sample, imputation) %>%
    distinct() %$%
    map2(sample, imputation, 
         med_plot[[func]])
}

med_plot$main <- function(sample, imputation){
  p <- med_res %>%
    unnest(margins) %>%
    filter(sample == !!sample,
           imputation == !!imputation, unem_age == FALSE) %>%
    filter(Age == 25) %>%
    ggplot() +
    aes(x = Sex, y = beta, ymin = lci, ymax = uci, 
        color = med_clean, shape = med_clean) +
    geom_hline(yintercept = 0) + 
    geom_pointrange(position = position_dodge(0.8)) +
    scale_color_brewer(palette = "Set1") +
    scale_shape_manual(values = 15:19) +
    labs(x = NULL, color = NULL, shape = NULL,
         y = "Marginal Effect") +
    theme_minimal() +
    theme(legend.position = c(.85,.85))
  
  glue("Images/ghq_main_{sample}_{imputation}.png") %>%
    ggsave(plot = p, dpi = 300,  width = 21, 
           height = 9.9, units = "cm") 
  
  return(p)
}
med_plot$make("main")

med_plot$age <- function(sample, imputation){
  med_basic <- med_res %>%
    unnest(margins) %>%
    filter(sample == !!sample, imputation == !!imputation, 
           med_mod == "basic", unem_age == TRUE) %>%
    select(-med_clean)
  
  p <- med_res %>%
    unnest(margins) %>%
    filter(sample == !!sample, imputation == !!imputation, 
           med_mod != "basic", unem_age == TRUE) %>%
    ggplot() +
    aes(x = Age, y = beta, ymin = lci, ymax = uci, group = med_mod) +
    facet_grid(med_clean ~ Sex, switch = "y") +
    geom_hline(yintercept = 0) +
    geom_line(data = med_basic, color = "grey50", linetype = "dashed") +
    geom_line(aes(color = med_clean)) +
    geom_ribbon(aes(color = med_clean, fill = med_clean),
                color = NA, alpha = 0.2) +
    scale_x_continuous(breaks = 2.5:6.5*10) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Age", y = NULL,
         color = NULL, fill = NULL) +
    guides(color = FALSE, fill = FALSE) +
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.spacing.y = unit(1, "lines"),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0))
  
  glue("Images/ghq_age_{sample}_{imputation}.png") %>%
    ggsave(plot = p, dpi = 300,  width = 21,
           height = 29.7, units = "cm")
  
  return(p)
}
med_plot$make("age")




# Allostatic Tables
mod_args <- tribble(
  ~mod, ~behav, ~sep,
  "basic", "NO", "NO",
  "behav", "YES", "NO",
  "sep", "NO", "YES",
  "all", "YES", "YES"
) %>%
  pivot_longer(-mod, names_to = "cat", values_to = "string")

clean_mod <- function(df){
  df_c <- slice(df, 1)
  
  args <- filter(mod_args, mod == df_c$mod)
  
  tibble(cat = c("obs", "imps", args$cat),
         string = c(format(df_c$obs, big.mark = ",") %>% trimws(),
                    df_c$imps, args$string)) %>%
    expand_grid(select(df_c, imputation, sample, mod, unem_age, sex)) %>%
    bind_rows(df)
}

allostatic_coefs <- main_pool %>%
  filter(outcome == "Allostatic") %>%
  unnest(coefs) %>%
  mutate(string = glue("{beta}\n({lci}, {uci})")) %>%
  group_split(imputation, sample, mod, unem_age, sex) %>%
  map_dfr(clean_mod) %>%
  left_join(pretty_lbls, by = "cat") %>%
  mutate(var_clean = ifelse(levels <= 2, "", var_clean)) %>%
  select(var_clean, cat_clean, cat, imputation, sample, 
         mod, unem_age, sex, string, index) %>%
  arrange(imputation, sample, mod, sex, index) %>%
  pivot_wider(names_from = "mod", values_from = "string", values_fill = "") %>%
  select(var_clean, cat_clean, cat, imputation, sample, unem_age, sex, all_of(names(mod_clean)))

allostatic_tbl <- list()
allostatic_tbl$short_func <- function(sample, imputation, unem_age){
  tbl <- allostatic_coefs %>%
    filter(str_detect(cat, "(Unem_Age|behav|sep|obs|imps)"),
           imputation == !!imputation, sample == !!sample, unem_age == !!unem_age) %>%
    select(sex, everything(), -c(imputation, var_clean, unem_age, sample, cat)) %>%
    mutate(sex = str_to_title(sex)) %>%
    flextable() %>%
    set_header_labels(sex = "Gender", cat_clean = "Variable", basic = "(1)",
                      behav = "(2)", sep = "(3)", all = "(4)") %>%
    merge_v(1) %>%
    align(j=3:6, align="center", part = "all") %>%
    valign(j = 1, valign = "top") %>%
    border(i = ~ cat_clean == "Imputations",
           border.bottom = fp_border(color="black", style = "solid")) %>%
    border(i = ~ cat_clean == "Health Behaviours",
           border.top = fp_border(color="gray70", style = "dashed")) %>%
    fix_border_issues(part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    autofit()
  
  if (unem_age == TRUE){
    file_name <- glue("Tables/allostatic_short_{sample}_{imputation}_age.docx")
  } else{
    file_name <- glue("Tables/allostatic_short_{sample}_{imputation}_main.docx")
  }
  
  save_as_docx(tbl, path = file_name)
  return(tbl)
}
allostatic_coefs %>%
  select(sample, imputation, unem_age) %>%
  distinct() %$%
  pmap(list(sample, imputation, unem_age), allostatic_tbl$short_func)
allostatic_tbl$short_func("mi", "z", TRUE)

allostatic_tbl$long_func <- function(sample, sex, unem_age, imputation){
  tbl <- allostatic_coefs  %>%
    mutate(var_clean = ifelse(str_detect(cat, "(Unem_Age|gender|behav|mental|sep|obs|imps)"),
                              " ", var_clean)) %>%
    filter(imputation == !!imputation, sample == !!sample, 
           sex == !!sex, unem_age == !!unem_age) %>%
    select(-c(imputation, sex, sample, cat, unem_age)) %>%
    flextable() %>%
    set_header_labels(var_clean = "", cat_clean = "Variable", basic = "(1)",
                      behav = "(2)", mental = "(3)",
                      sep = "(4)", all = "(5)") %>%
    merge_v(1) %>%
    align(j=3:6, align="center", part = "all") %>%
    valign(j = 1, valign = "top") %>%
    border_inner_h(border = fp_border(color="gray70", style = "dashed")) %>%
    border(i = ~ cat_clean == "Constant", 
           border.bottom = fp_border(color="gray30", style = "solid")) %>%
    fix_border_issues(part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    autofit()
  
  if (unem_age == TRUE){
    file_name <- glue("Tables/allostatic_long_{sample}_{imputation}_{sex}_age.docx")
  } else{
    file_name <- glue("Tables/allostatic_long_{sample}_{imputation}_{sex}_main.docx")
  }
  
  save_as_docx(tbl, path = file_name)
  tbl
}

allostatic_coefs %>%
  select(sample, sex, unem_age, imputation) %>%
  distinct() %$%
  pmap(list(sample, sex, unem_age, imputation), allostatic_tbl$long_func)

# Mediation Table
med_args <- tribble(
  ~mod, ~allostatic, ~behav, ~sep,
  "basic", "NO", "NO", "NO",
  "allostatic", "YES", "NO", "NO",
  "behav", "NO", "YES","NO", 
  "sep", "NO", "NO", "YES",
  "all", "YES", "YES", "YES"
) %>%
  pivot_longer(-mod, names_to = "cat", values_to = "string")

clean_med <- function(df){
  df_c <- slice(df, 1)
  
  args <- filter(med_args, mod == df_c$med_mod)
  
  tibble(cat = c("gender", "obs", "imps", args$cat),
         string = c(str_to_title(df_c$sex),
                    format(df_c$obs, big.mark = ",") %>% trimws(),
                    df_c$imps, args$string)) %>%
    expand_grid(select(df_c, imputation, sample, med_mod, unem_age, sex)) %>%
    bind_rows(df)
}

med_coefs <- main_pool %>%
  filter(outcome == "GHQ_Likert") %>%
  mutate(med_mod = case_when(mediate == TRUE & mod == "basic" ~ "allostatic",
                             mediate == TRUE & mod == "all" ~ "all",
                             mediate == FALSE & mod != "all" ~ mod),
         med_clean = factor(med_clean[med_mod], med_clean)) %>%
  filter(!is.na(med_mod)) %>%
  unnest(coefs) %>%
  mutate(string = glue("{beta}\n({lci}, {uci})")) %>%
  group_split(imputation, sample, med_mod, unem_age, sex) %>%
  map_dfr(clean_med) %>%
  left_join(pretty_lbls, by = "cat") %>%
  mutate(var_clean = ifelse(levels <= 2, "", var_clean)) %>%
  select(var_clean, cat_clean, cat, imputation, sample, 
         med_mod, unem_age, sex, string, index) %>%
  arrange(imputation, sample, med_mod, sex, index) %>%
  pivot_wider(names_from = "med_mod", values_from = "string", values_fill = "") %>%
  select(var_clean, cat_clean, cat, imputation, sample, unem_age, sex, all_of(names(med_clean)))

med_tbl <- list()
med_tbl$short_func <- function(sample, imputation, unem_age){
  tbl <- med_coefs %>%
    filter(str_detect(cat, "(Unem_Age|allostatic|behav|sep|obs|imps)"),
           imputation == !!imputation, sample == !!sample, unem_age == !!unem_age) %>%
    select(sex, everything(), -c(imputation, var_clean, unem_age, sample, cat)) %>%
    mutate(sex = str_to_title(sex)) %>%
    flextable() %>%
    set_header_labels(sex = "Gender", cat_clean = "Variable", basic = "(1)",
                      allostatic = "(2)", behav = "(3)", sep = "(4)", all = "(5)") %>%
    merge_v(1) %>%
    align(j=3:7, align="center", part = "all") %>%
    valign(j = 1, valign = "top") %>%
    border(i = ~ cat_clean == "Imputations",
           border.bottom = fp_border(color="black", style = "solid")) %>%
    border(i = ~ cat_clean == "Allostatic Load",
           border.top = fp_border(color="gray70", style = "dashed")) %>%
    fix_border_issues(part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    autofit()
  
  if (unem_age == TRUE){
    file_name <- glue("Tables/med_short_{sample}_{imputation}_age.docx")
  } else{
    file_name <- glue("Tables/med_short_{sample}_{imputation}_main.docx")
  }
  
  save_as_docx(tbl, path = file_name)
  return(tbl)
}
med_coefs %>%
  select(sample, imputation, unem_age) %>%
  distinct() %$%
  pmap(list(sample, imputation, unem_age), med_tbl$short_func)

med_tbl$long_func <- function(sample, sex, unem_age, imputation){
  tbl <- med_coefs  %>%
    mutate(var_clean = ifelse(str_detect(cat, "(Unem_Age|gender|allostatic|behav|mental|sep|obs|imps)"),
                              " ", var_clean)) %>%
    filter(imputation == !!imputation, sample == !!sample, 
           sex == !!sex, unem_age == !!unem_age) %>%
    select(-c(imputation, sex, sample, cat, unem_age)) %>%
    flextable() %>%
    set_header_labels(var_clean = "", cat_clean = "Variable", basic = "(1)",
                      allostatic = "(2)", behav = "(3)", sep = "(4)", all = "(5)") %>%
    merge_v(1) %>%
    align(j=3:7, align="center", part = "all") %>%
    valign(j = 1, valign = "top") %>%
    border_inner_h(border = fp_border(color="gray70", style = "dashed")) %>%
    border(i = ~ cat_clean == "Constant", 
           border.bottom = fp_border(color="gray30", style = "solid")) %>%
    fix_border_issues(part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    autofit()
  
  if (unem_age == TRUE){
    file_name <- glue("Tables/allostatic_long_{sample}_{imputation}_{sex}_age.docx")
  } else{
    file_name <- glue("Tables/allostatic_long_{sample}_{imputation}_{sex}_main.docx")
  }
  
  save_as_docx(tbl, path = file_name)
  tbl
}
med_coefs %>%
  select(sample, sex, unem_age, imputation) %>%
  distinct() %$%
  pmap(list(sample, sex, unem_age, imputation), med_tbl$long_func)


# 4. Quantile Results ----
load("Data/quantile_results.Rdata")

quantile_coef <- quantile_coef %>%
  left_join(pretty_lbls, by = "cat") %>%
  arrange(sex, tau, index) %>%
  mutate(across(c(beta, lci, uci),
                ~ ifelse(str_detect(cat, "Age_C"), .x*39, .x) %>%
                  round(2)),
         string = glue("{beta}\n({lci}, {uci})"), 
         q = glue("Q{tau*100}"),
         sex_i = glue("{sex}_{index}"),
         var_clean = ifelse(levels <= 2, "", var_clean)) 

# Plots
quantile_coef %>%
  filter(str_detect(cat, "^Unem")) %>%
  mutate(sex_f = str_to_title(sex),
         imp = ifelse(imputation == "z", "Z-Score", "Index")) %>%
  ggplot() +
  aes(x = tau, y = beta, ymin = lci, ymax = uci,
      color = sex_f, fill = sex_f) +
  facet_grid(imp ~ sex_f, scales = "free_y") +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.2) +
  scale_x_continuous(breaks = 1:9/10) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(panel.spacing.y = unit(1, "lines"), 
        strip.text.y = element_text(angle=0)) +
  guides(color = FALSE, fill = FALSE) +
  labs(x = "Quantile", y = "Marginal Effect")
ggsave("Images/quantile.png", dpi = 300,
       width = 21, height = 16, units = "cm")

# Tables
# Main Effects
quantile_tbl <- list()

quantile_tbl$short <- quantile_coef %>%
    filter(str_detect(cat, "Unem_Age"), 
           unem_age == FALSE) %>%
    mutate(sex_imp = glue("{sex}_{imputation}")) %>%
    select(q, string, sex_imp) %>%
    pivot_wider(names_from = sex_imp, values_from = string) %>%
    select(q, matches("index"), matches("z")) %>%
    flextable() %>%
    border_remove() %>%
    set_header_labels(q = "Quantile", all_index = "All", female_index = "Female", male_index = "Male",
                      all_z = "All", female_z = "Female", male_z = "Male") %>%
    add_header(values = list(q = "", all_index = "Index", female_index = "Index", male_index = "Index",
                             all_z = "Z-Score", female_z = "Z-Score", male_z = "Z-Score")) %>%
    merge_h(1, part = "header") %>%
    align(j=2:7, align="center", part = "all") %>%
    hline_top(border = fp_border(color="black", width = 2), part = "all") %>%
    hline_bottom(border = fp_border(color="black", width = 2), part = "all") %>%
    border(j = c(1, 4), border.right = fp_border(color="grey70", style = "dashed")) %>%
    fix_border_issues(part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    autofit()
quantile_tbl$short 
save_as_docx(quantile_tbl$short, path = glue("Tables/quantile_short.docx"))

# Age Effects
quantile_names <- map(c("all", "female", "male"), ~ glue("{.x}_{3:5}")) %>% 
  unlist() %>%
  c("q", .)
quantile_header_1 <- c("Quantile", rep(c("Main Effect", "x Age", "x Age^2"), 3)) %>%
  set_names(quantile_names) %>% as.list()
quantile_header_2 <- c("", rep(c("All", "Female", "Male"), each = 3)) %>%
  set_names(quantile_names) %>% as.list()
quantile_tbl$short_age <- function(imputation){
  tbl <- quantile_coef %>%
    filter(str_detect(cat, "Unem_Age"), 
           imputation == !!imputation, 
           unem_age == TRUE) %>%
    select(q, string, sex_i) %>%
    pivot_wider(names_from = sex_i, values_from = string) %>%
    flextable() %>%
    border_remove() %>%
    set_header_labels(values = quantile_header_1) %>%
    add_header_row(values = quantile_header_2) %>%
    merge_h(1, part = "header") %>%
    align(j=2:10, align="center", part = "all") %>%
    hline_top(border = fp_border(color="black", width = 2), part = "all") %>%
    hline_bottom(border = fp_border(color="black", width = 2), part = "all") %>%
    border(j = c(1, 4, 7),
           border.right = fp_border(color="grey70", style = "dashed")) %>%
    fix_border_issues(part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    autofit()
  
  save_as_docx(tbl, path = glue("Tables/quantile_short_age_{imputation}.docx"))
  
  return(tbl)
}
quantile_tbl$short_age("index")
quantile_tbl$short_age("z")

quantile_tbl$long <- function(sex, imputation, unem_age){
  tbl <- quantile_coef %>%
    filter(sex == !!sex,
           imputation == !!imputation,
           unem_age == !!unem_age) %>%
    select(var_clean, cat_clean, string, q) %>%
    pivot_wider(names_from = q, values_from = string) %>%
    flextable() %>%
    set_header_labels(var_clean = "", cat_clean = "Variable") %>%
    merge_v(1, part = "body") %>%
    align(j=3:11, align="center", part = "all") %>%
    border_inner_h(border = fp_border(color="gray70", style = "dashed")) %>%
    fix_border_issues(part = "all") %>%
    fontsize(size = 11, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    autofit()
  
  save_as_docx(tbl, path = glue("Tables/quantile_long_{sex}.docx"))
  
  tbl
}
quantile_coef %>%
  select(sex, imputation, unem_age) %>%
  distinct() %$%
  pmap(list(sex, imputation, unem_age), quantile_tbl$long)


# 5. Specification Curve Analysis ----
load("Data/sca_results.Rdata")

sca <- sca %>%
  mutate(type = ifelse(quartile, "Index", "Z-Score"),
         sex_u = str_to_title(sex),
         across(c(beta, lci, uci),
                list(std = ~ .x/sd,
                     range = ~ .x/n_vars)))

sca %>%
  group_by(sex, signif) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n*100/sum(n))

sca %>%
  group_by(type, sex, signif) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n*100/sum(n))

sca %>%
  group_by(sex) %>%
  summarise(lci = quantile(beta_std, .025),
            uci = quantile(beta_std, .975))

sca %>%
  group_by(type, sex) %>%
  summarise(lci = quantile(beta_std, .025),
            uci = quantile(beta_std, .975))
  
load("Data/bio_clean.Rdata")

sca_main <- sca %>% filter(n_vars == max(n_vars)) %>%
  select(-signif, -type)

# Plots
sca_plot <- function(p, y_lab = NULL){
  p +
    scale_x_continuous(labels = scales::comma) +
    scale_color_manual(values = c("red", "black")) +
    scale_fill_manual(values = c("red", "black")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.spacing.x = unit(1.25, "lines")) +
    labs(x = "Rank", y = y_lab, color = NULL) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) 
}


p <- ggplot(sca) + 
  aes(x = rank, y = beta_std, ymin = lci_std,
      ymax = uci_std, color = signif) +
  facet_wrap(~ sex_u) +
  geom_hline(yintercept = 0) +
  geom_pointrange(size = 0.05, alpha = 0.15) +
  geom_point(size = 0.05) +
  geom_point(data = sca_main, shape = 23, 
             color = "black", fill = "white",
             size = 3, stroke = 1.5)
sca_plot(p)
ggsave("Images/sca_main.png", dpi = 300,
       width = 21, height = 9.9, units = "cm")


p <- ggplot(sca) + 
  aes(x = rank_q, y = beta_std, ymin = lci_std,
      ymax = uci_std, color = type) +
  facet_wrap(~ sex_u) +
  geom_hline(yintercept = 0) +
  geom_point(size = 0.05) +
  geom_point(data = sca_main, shape = 23, 
             color = "black", fill = "white",
             size = 3, stroke = 1.5)
sca_plot(p) +
  scale_color_brewer(palette = "Accent") +
  guides(color = guide_legend(override.aes = list(size = 1, alpha = 1))) +
  labs(color = "Allostatic Load")
ggsave("Images/sca_comparison.png", dpi = 300,
       width = 21, height = 9.9, units = "cm")

p <- sca %>%
  unnest(vars) %>%
  mutate(vars = str_replace(vars, "_Quartile", ""),
         var_clean = bio_clean[vars]) %>%
  ggplot() +
  aes(x = rank, y = var_clean) +
  facet_wrap( ~ sex_u) +
  geom_boxplot(width = 0.6) +
  geom_vline(xintercept = 2510, linetype = "dashed", color = "grey70") 
sca_plot(p) +
  theme(panel.grid = element_blank())
ggsave("Images/sca_vars.png", dpi = 300,
       width = 21, height = 9.9, units = "cm")

p <- ggplot(sca) +
  aes(x = rank, y = factor(n_vars), color = signif) +
  facet_wrap(~ sex_u) +
  geom_jitter(alpha = 0.3)
sca_plot(p, "Variables in Index")
ggsave("Images/sca_index.png", dpi = 300,
       width = 21, height = 9.9, units = "cm")
