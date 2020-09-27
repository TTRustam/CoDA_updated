# !!! Check main results and comments after line 220. 
# -------------------------------------------------------------------------------------------------- #  
library("here")
library("tidyverse")
options(dplyr.summarise.inform = FALSE)
# -------------------------------------------------------------------------------------------------- #  
# !!! The most important parameters. Training window and forecast horizon.
# 1990 is more or less fine, but still it seems that vecm model cannot handle forecast if
# training period sis too short
# 1992 is the shortest perspective we can take to fit the 10 years in sample data, hence 16 years is min.
years <-  1975:2008
ih    <-  10
# -------------------------------------------------------------------------------------------------- #  
# Other parameters. Age structure, which ages to keep, which causes exist, etc.
# !!! The order of causes in important.
ages      <- c(0, 2.5, seq(7, 97, 5), 100)
ages_kept <- c(as.character(ages[-c(1:6, 19:22)]), "85")
m         <- length(ages)
n         <- length(years)
years_for <- c(years, (max(years) + 1):(max(years) + ih))
nmCoD_m <- c("Other cancer",
             "Other disease",
             "Prostate",
             "Stomach",
             "Trachea, bronchus and lungs")

nmCoD_f <- c("Cervix uteri",
             "Female breast",
             "Other cancer",
             "Other disease",
             "Stomach", 
             "Trachea, bronchus and lungs")

# Here many things should be checked
k_f         <- length(nmCoD_f)
k_m         <- length(nmCoD_m)
min_age     <- 25 - 1 
o           <- length(ages_kept)
var_names_m <- str_c(rep(nmCoD_m, each = 13), ages_kept)
var_names_f <- str_c(rep(nmCoD_f, each = 13), ages_kept)
sp_ac_m     <- split(var_names_m, ceiling(seq_along(var_names_m) / o))
sp_ac_f     <- split(var_names_f, ceiling(seq_along(var_names_f) / o))
oldest_ages <- "85"
# -------------------------------------------------------------------------------------------------- #  
# Data, in order: Total Mx, cause-specific Mx, population, life-table and cause-specific cancer mortality
# Total Mx
tot_death <- readRDS("mortality/tot_death.rds") %>% 
  dplyr::select(sex, age, num_range(range = years, prefix = ""))

# cause-specific Mx
df_mort   <- readRDS("mortality/mortality.rds") %>%
  filter(year %in% years)

# population
df_exp    <- readRDS("exposure/exposure.rds") %>%
  filter(year %in% years)

# life-table (Preston method)
lt        <- readRDS("life_table/life_table.rds") %>%
  mutate(data = map(data, ~ .x %>% 
                      filter(year %in% years)))

# cause-specific cancer mortality
m_cancer  <- readRDS("mortality/cancer_mortality.rds") %>% 
  mutate(data = map(data, ~ .x %>%
                      dplyr::select(age, num_range(range = years, prefix = ""))))
# -------------------------------------------------------------------------------------------------- #  
# functions. Inner function consist of different types of SVD and Life-table constructor.
# Main functions are coda forecast
# Auxiliary functions are small and big helpers. Wherever I could, I made a function (if needed > 2 times)
# lot_data is for preparing the plots and ggplot function
source("coda_functions/for_russia_auxillary_functions.R")
source("coda_functions/for_russia_inner_functions.R")
source("coda_functions/for_russia_main_functions.R")
source("coda_functions/plot_data.R")
# -------------------------------------------------------------------------------------------------- #
# call a function from plot_data
data <- prepare_data()
# -------------------------------------------------------------------------------------------------- #
# This should be changed before final iteration. SInce causes are different for M and F,
# We should specify the exact set of causes to be used for naming variables DF in other function
# I know how to change, just takes some time, to wrap the head around it.
nm_1 <-  var_names_f
nm_2 <-  nmCoD_f
# -------------------------------------------------------------------------------------------------- #
# !!! These cooments could be neglected. SImple descriotion of 2s model.
# model example
# This model is being done in 2 steps (hence is the name)
# First it calculates and applies weights to the rows of data (years). 
# Thus ages with more occurrences take more wights and recent years are give more weight. 0.9 to 0.0003
# Then takes acomp of lt deaths subtracts geom. mean takes clr of this 
# Then calculates the SVD of this function with weights for rows and columns 
# row weights by reverse of w_scaled, col_weights by age weight. Takes first approximation. 1.
# bx = first right singular vector
# kt = first singular value of x * first left singular vector
# squares of first 4 singular values divided by sum of squares of these 4 values
# then forecast kt with 1 1 1 arima with drift and Max. Likelihood + forecast ih years ahead
# -------------------------------------------------------------------------------------------------- #
# Models for females
model_fit_for_2s_f   <- CoDa_2step(dx         = data$dx_com_f,
                                   w_base     = 0.1,
                                   age_weight = data$f_col,
                                   ih         = ih,
                                   k          = k_f,
                                   years      = years,
                                   ages       = ages_kept)

model_fit_for_ct_f   <- CoDa_CT(dx    = data$dx_com_f,
                                ih    = ih,
                                k     = k_f,
                                years = years,
                                ages  = ages_kept)

# The warning here (in case of short fitting period) indicates that 
# Johanson trace test produces some NaN results. 
# Best be avoided But I have no Idea how to deal with it.
# This depends on a chosen fitting period, if it is small, this problem pops out.
# If it is too small, the function does not work a t all since test results return complex numbers 
model_fit_for_vecm_f <- CoDa_vecm(dx    = data$dx_com_f,
                                  ir    = 3,
                                  ih    = ih,
                                  k     = k_f,
                                  years = years,
                                  ages  = ages_kept)

model_fit_for_lc_f   <- COD_LC(mx    = data$f_cause,
                               t     = ih,
                               k     = k_f,
                               ages  = ages_kept,
                               years = years)
# -------------------------------------------------------------------------------------------------- #
# same for males. Note that k parameter (number of causes) also varies
nm_1 <-  var_names_m
nm_2 <-  nmCoD_m  

model_fit_for_2s_m   <- CoDa_2step(dx         = data$dx_com_m,
                                   w_base     = 0.1,
                                   age_weight = data$m_col,
                                   ih         = ih,
                                   k          = k_m,
                                   years      = years,
                                   ages       = ages_kept)

model_fit_for_ct_m   <- CoDa_CT(dx    = data$dx_com_m,
                                ih    = ih,
                                k     = k_m,
                                years = years,
                                ages  = ages_kept)

model_fit_for_vecm_m <- CoDa_vecm(dx    = data$dx_com_m,
                                  ir    = 3,
                                  ih    = ih,
                                  k     = k_m,
                                  years = years,
                                  ages  = ages_kept)

model_fit_for_lc_m   <- COD_LC(mx    = data$m_cause,
                               t     = ih,
                               k     = k_m,
                               ages  = ages_kept,
                               years = years)

# -------------------------------------------------------------------------------------------------- #
# Prepare the data for plotting. Namely cause-specific dx.
plot_data_f <- prepare_to_plot(mod_2s = model_fit_for_2s_f,
                               LC     = model_fit_for_lc_f,
                               vecm   = model_fit_for_vecm_f,
                               ct     = model_fit_for_ct_f,
                               obs    = data$dx_com_f,
                               sex    = "Females") 
plot_data_m <- prepare_to_plot(mod_2s = model_fit_for_2s_m,
                               LC     = model_fit_for_lc_m,
                               vecm   = model_fit_for_vecm_m,
                               ct     = model_fit_for_ct_m,
                               obs    = data$dx_com_m,
                               sex    = "Males")  
# -------------------------------------------------------------------------------------------------- #
# Plots
females_plot <- plot_function(plot_data_f)
males_plot   <- plot_function(plot_data_m)
# If plots need to be saved + 
# ggsave(filename = "20_forecast_fem.jpeg", width = 6, height = 4, units = "in", scale = 2, 
#        path = here("figures"))
# -------------------------------------------------------------------------------------------------- #
# plot of compositional data
comp_data <- rel_and_neo_distr(dat = data)
comp_plots <- plot_comp_data()
# to save
# ggsave(filename = "distribution.jpeg", width = 6, height = 4, units = "in", scale = 2, 
#        path = here("figures"))
# -------------------------------------------------------------------------------------------------- #
# kt and bx plots
kt_bx <- kt_and_bx_plot()
# -------------------------------------------------------------------------------------------------- #
# in sample comparison
in_sample <- rmse_in_sample(start = max(years)) 
# -------------------------------------------------------------------------------------------------- #
# RMSE
# females
f_vecm <- rmse_models(model = "VECM CoDA", data = "RMSE/females_observed", sx = "Females")
f_ct   <- rmse_models(model = "CT CoDA",   data = "RMSE/females_observed", sx = "Females")
f_2s   <- rmse_models(model = "2S CoDA",   data = "RMSE/females_observed", sx = "Females")
f_LC   <- rmse_models(model = "LC",        data = "RMSE/females_observed", sx = "Females")

# males
m_vecm <- rmse_models(model = "VECM CoDA", data = "RMSE/males_observed", sx = "Males")
m_ct   <- rmse_models(model = "CT CoDA",   data = "RMSE/males_observed", sx = "Males")
m_2s   <- rmse_models(model = "2S CoDA",   data = "RMSE/males_observed", sx = "Males")
m_LC   <- rmse_models(model = "LC",        data = "RMSE/males_observed", sx = "Males")

# RMSE resulting table
rmse <- tibble(rmse  = c(f_2s, f_LC, f_vecm, f_ct, m_2s, m_LC, m_vecm, m_ct),
               model = c("2S CoDA", "LC", "VECM CoDA", "CT CoDA", "2S CoDA", "LC", "VECM CoDA", "CT CoDA"),
               sex   = c(rep("Females", 4), rep("Males", 4))) %>%
  pivot_wider(names_from  = "sex", 
              values_from = "rmse")                           %>% 
  mutate(model = factor(model, levels = c("2S CoDA", "VECM CoDA", "CT CoDA", "LC")))

# -------------------------------------------------------------------------------------------------- #
# Results printed 
# 10 years ahead forecast. This plots only make sence if we forecast beyonf 2018
females_plot
males_plot

# 10 years in sample comparisons. Only makes sence if we forecast in sample e.g. 2008:2018
in_sample

# compositional plots 
comp_plots

# kt and bx
kt_bx

# RMSE
rmse
# END.
# -------------------------------------------------------------------------------------------------- #