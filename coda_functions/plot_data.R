# -------------------------------------------------------------------------------------------------- #  
# This function prepares all the data we need to make a CoDA forecast.
prepare_data <- function() { 
  
  # Remove ages below 24 and years before 1965 from life tables
  f_dx <- lt_create_RT(gender = "Female", variable = "dx")
  f_Lx <- lt_create_RT(gender = "Female", variable = "Lx")
  m_dx <- lt_create_RT(gender = "Male",   variable = "dx")
  m_Lx <- lt_create_RT(gender = "Male",   variable = "Lx")
  
  # divide each cancer localization deaths by total deaths, then multiply by total dx
  m_rel <- rel_dat_func(dx = m_dx, gender = "Male")
  f_rel <- rel_dat_func(dx = f_dx, gender = "Female")
  
  # this one to CoDA the model. just reshaping data
  dx_com_m <- com_func(m_rel)
  dx_com_f <- com_func(f_rel)
  
  # total sum of relative deaths by sex and years 
  dx_tot_m    <- tot_func(dat = m_rel)
  dx_tot_f    <- tot_func(dat = f_rel)
  
  # mean relative deaths at each age
  col_means_m <- colMeans(dx_tot_m)
  col_means_f <- colMeans(dx_tot_f)
  
  # divide mean relative (cancer by causes) deaths by age BY sum of the total means by age (colmeans) (0 to 1)
  m_col <- col_func(dat = m_rel, means = col_means_m)
  f_col <- col_func(dat = f_rel, means = col_means_f)
  
  # divide cancer by causes relative data by Lx // For LC CoDA model
  m_cause <- cause_func(dat = m_rel, Lx = m_Lx)
  f_cause <- cause_func(dat = f_rel, Lx = f_Lx)
  
  return(tibble::lst(f_dx, f_Lx, m_dx, m_Lx, m_rel, f_rel, dx_com_m, dx_com_f, 
                     col_means_m, col_means_f, m_col, f_col, m_cause, f_cause))
}
# -------------------------------------------------------------------------------------------------- #  
library("scales")
library("ggthemes")
# -------------------------------------------------------------------------------------------------- #
# Prepare data for plotting
prepare_to_plot <- function(mod_2s, LC, vecm, ct, obs, sex) { 
  
  
  perc_d_obs <- obs                                   %>% 
    pivot_longer(everything(),
                 names_to  = "var",
                 values_to = "val")                   %>% 
    separate(var, c("cause", "age"), "_", remove = F) %>% 
    dplyr::select(-age)                               %>% 
    group_by(cause)                                   %>% 
    nest()                                            %>%
    mutate(data = map(data, ~ .x                         %>% 
                        group_by(var)                    %>% 
                        mutate(year = years)             %>% 
                        pivot_wider(names_from  = "var",
                                    values_from = "val") %>% 
                        dplyr::select(-year)             %>% 
                        ungroup))                     %>% 
    perc_d_func()                                     %>% 
    mutate(year = years)                              %>% 
    pivot_longer(-year, 
                 names_to  = "cause", 
                 values_to = "value")                 %>% 
    filter(cause != "Other disease")                  %>% 
    mutate(model  = "Observed")
  
  perc_d_2s   <- mod_2s$dx_forcast    %>%
    perc_d_func()                     %>% 
    mutate(year = years_for)          %>% 
    pivot_longer(-year, 
                 names_to  = "cause", 
                 values_to = "value") %>% 
    filter(cause != "Other disease")  %>% 
    mutate(model  = "2S CoDA")
  
  
  perc_d_ct   <- ct$dx_forcast[,-3]   %>%
    perc_d_func()                     %>% 
    mutate(year = years_for)          %>% 
    pivot_longer(-year,
                 names_to  = "cause",
                 values_to = "value") %>% 
    filter(cause != "Other disease")  %>% 
    mutate(model  = "CT CoDA")
  
  perc_d_vecm <- vecm$dx_forcast      %>% 
    perc_d_func()                     %>% 
    mutate(year = years_for)          %>% 
    pivot_longer(-year,
                 names_to  = "cause",
                 values_to = "value") %>% 
    filter(cause != "Other disease")  %>% 
    mutate(model  = "VECM CoDA")
  
  perc_d_lc   <- LC$dx_forcast        %>%
    perc_d_func()                     %>% 
    mutate(year = years_for)          %>% 
    pivot_longer(-year, 
                 names_to  = "cause", 
                 values_to = "value") %>% 
    filter(cause != "Other disease")  %>% 
    mutate(model  = "LC")
  
  
  to_plot <- perc_d_obs    %>% 
    full_join(perc_d_2s)   %>% 
    full_join(perc_d_ct)   %>% 
    full_join(perc_d_vecm) %>% 
    full_join(perc_d_lc)   %>%
    mutate(sex = sex) %>% 
    mutate(model = factor(model, 
                          levels = c("Observed", "2S CoDA", 
                                     "VECM CoDA", "CT CoDA", "LC")))
  return(to_plot)
  
  }
# -------------------------------------------------------------------------------------------------- #  
# The ggplot function itself
plot_function <- function(data) {
  data %>% 
    ggplot(aes(x = year, y = value, group = interaction(cause, model)))    +
    geom_line(aes(linetype = model, color = cause), size = 1, alpha = 0.7) +    
    # facet_grid(ind ~ sex, scales = "free_y", switch = "y")               +
    scale_x_continuous(breaks = seq(min(years), max(years), 5))            +
    scale_color_tableau()                                                  +
    expand_limits(y = 0)                                                   +
    # scale_size_continuous(range = c(1, 1.5))                             +
    scale_y_continuous(breaks = pretty_breaks())                           +
    theme_bw()                                                             +
    theme(
      legend.position   = c(.85, .62),
      legend.direction  = "vertical",
      axis.text.x       = element_text(colour = "black",       size = 12), 
      legend.key        = element_rect(colour = "transparent"),
      legend.background = element_rect(colour = "transparent", fill = alpha('white', 0)),
      axis.text.y       = element_text(colour = "black",       size = 12), 
      legend.title      = element_blank(),
      legend.text       = element_text(size = 12, face = "bold", color = "black"),
      strip.text.x      = element_text(size = 12, face = "bold", color = "black"),
      strip.text.y      = element_text(size = 12, face = "bold", color = "black"),
      axis.title        = element_blank(), 
      strip.background  = element_blank(),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width  = unit(0.8, "cm")) + 
    guides(size = FALSE) +
    geom_vline(xintercept = max(years))
  
  }
# -------------------------------------------------------------------------------------------------- #  
# Here we prepare compositional data fro plotting.
rel_and_neo_distr <- function(dat) {
  
  # relative distances for males and females
  d_rel_dist_f <- dat$col_means_f / sum(dat$col_means_f)
  d_rel_dist_m <- dat$col_means_m / sum(dat$col_means_m)
  rel_dist_f   <- rep(d_rel_dist_f, times = length(nmCoD_f))
  rel_dist_m   <- rep(d_rel_dist_m, times = length(nmCoD_m))
  
# -------------------------------------------------------------------------------------------------- #
# ------------ Aitchinson composition sum age == 1 ------------- #
  neo_close_m <- neo_func(dat$m_rel)
  neo_close_f <- neo_func(dat$f_rel)
  
  geo_mean_neo_f  <- geometricmeanCol(neo_close_f)
  neo_cent_f      <- neo_close_f - geo_mean_neo_f
  mean_neo_cent_f <- geometricmeanCol(neo_cent_f)
  close_dx_p_f    <- acomp(dat$f_dx)
  geo_mean_neo_m  <- geometricmeanCol(neo_close_m)
  neo_cent_m      <- neo_close_m - geo_mean_neo_m
  mean_neo_cent_m <- geometricmeanCol(neo_cent_m)
  close_dx_p_m    <- acomp(dat$m_dx)
  
  ax_p_f                 <- geometricmeanCol(close_dx_p_f)
  dx_cent_p_f            <- close_dx_p_f - ax_p_f
  clr_cent_p_f           <- clr(dx_cent_p_f)
  rownames(clr_cent_p_f) <- ages_kept
  colnames(clr_cent_p_f) <- as.character(years)  
  ax_p_m                 <- geometricmeanCol(close_dx_p_m)
  dx_cent_p_m            <- close_dx_p_m - ax_p_m
  clr_cent_p_m           <- clr(dx_cent_p_m)
  rownames(clr_cent_p_m) <- ages_kept
  colnames(clr_cent_p_m) <- as.character(years)
  
  ############## ALSO PLOT NEOCENT f
  distr_m  <- plot_prep(dat = neo_close_m, dat2 = geo_mean_neo_m)
  distr_f  <- plot_prep(dat = neo_close_f, dat2 = geo_mean_neo_f)
  neo_df_m <- neo_plot(neo_cent_m)
  neo_df_f <- neo_plot(neo_cent_f)
  
  distribution <- distr_m$distr %>% 
    mutate(sex = "Male") %>% 
    full_join(distr_f$distr) %>% 
    mutate(sex = ifelse(is.na(sex), "Female", sex), 
           sex = factor(sex, levels = c("Male", "Female")))
  
  geometr_lines <- distr_m$geom_line %>% 
    mutate(sex = "Male") %>% 
    full_join(distr_f$geom_line) %>% 
    mutate(sex = ifelse(is.na(sex), "Female", sex)) %>% 
    mutate(sex = factor(sex, levels = c("Male", "Female")))
  
  neo_final <- neo_df_m %>% 
    mutate(sex = "Male") %>% 
    full_join(neo_df_f) %>% 
    mutate(sex = ifelse(is.na(sex), "Female", sex), 
           sex = factor(sex, levels = c("Male", "Female")))
  
  
  return(tibble::lst(distr_m, distr_f, neo_df_m, neo_df_f, 
                     distribution, geometr_lines, neo_final))
  
}

# -------------------------------------------------------------------------------------------------- #
# function to prepare and plot compositional data
plot_comp_data <- function() {

one <- ggplot()                                                                            + 
  geom_line(data  = comp_data$distribution, 
            aes(x = age, y = val,   group = year, color = year), size = 1)                 +
  geom_line(data  = comp_data$geometr_lines, 
            aes(x = age, y = value, group = meh), color = "black", size = 1, linetype = 2) +
  facet_grid(. ~ sex, scales = "free_y", switch = "y")                                     +
  theme_bw()                                                                               +
  scale_y_continuous(breaks = pretty_breaks())                                             +
  scale_color_tableau()                                                                    +
  theme(legend.position = "bottom",
        strip.placement = "horizontal",
        axis.text.x  = element_text(colour = "black", size = 12), 
        axis.text.y  = element_text(colour = "black", size = 12), 
        legend.title = element_blank(),
        legend.text  = element_text(size = 12, face = "bold", color = "black"),
        strip.text.x = element_text(size = 12, face = "bold", color = "black"),
        strip.text.y = element_text(size = 12, face = "bold", color = "black"),
        axis.title   = element_blank(), 
        strip.background  = element_blank(),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width  = unit(0.8, "cm"))                                               + 
  guides(size = FALSE) 

two <- ggplot()                                                          + 
  geom_line(data  = comp_data$neo_final, 
            aes(x = age, y = val, group = year, color = year), size = 1) + 
  facet_grid(. ~ sex, scales = "free_y", switch = "y")                   +
  theme_bw()                                                             +
  scale_y_continuous(breaks = pretty_breaks())                           +
  scale_color_tableau()                                                  +
  theme(legend.position = "bottom",
        strip.placement = "horizontal",
        axis.text.x  = element_text(colour = "black", size = 12), 
        axis.text.y  = element_text(colour = "black", size = 12), 
        legend.title = element_blank(),
        legend.text  = element_text(size = 12, face = "bold", color = "black"),
        strip.text.x = element_text(size = 12, face = "bold", color = "black"),
        strip.text.y = element_text(size = 12, face = "bold", color = "black"),
        axis.title   = element_blank(), 
        strip.background  = element_blank(),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width  = unit(0.8, "cm"))                              + 
  guides(size = FALSE)

return(list(one, two))
}

# -------------------------------------------------------------------------------------------------- #
# kt and bx plot function. Also should generalize this later.
kt_and_bx_plot <- function() {
  # bx
  
  bx_func <- function(data, name, cod) {
    
    dat <- data                         %>% 
      as_tibble()                       %>% 
      set_names(cod)                    %>% 
      mutate(age = ages_kept)           %>% 
      pivot_longer(-age, 
                   names_to  = "cause", 
                   values_to = "value") %>% 
      mutate(model = name)
    
    return(dat)  
  }
  
# kt
  kt_func <- function(data, name, cod) {
    
    dat <- data                         %>% 
      as_tibble()                       %>% 
      set_names(cod)                    %>% 
      mutate(year = years_for)          %>% 
      pivot_longer(-year, 
                   names_to  = "cause", 
                   values_to = "value") %>% 
      mutate(model = name)
    
    return(dat)  
  }
  
# bx
  bx_step_m <- bx_func(model_fit_for_2s_m$bxs,    "2S-CoDA",   nmCoD_m)
  bx_vecm_m <- bx_func(model_fit_for_vecm_m$bxs,  "VECM-CoDA", nmCoD_m)  
  bx_lc_m   <- bx_func(model_fit_for_lc_m$bx_all, "LC",        nmCoD_m)  
  bx_ct_m   <- bx_func(model_fit_for_ct_m$bxs,    "CT-CoDA",   nmCoD_m)  
  bx_step_f <- bx_func(model_fit_for_2s_f$bxs,    "2S-CoDA",   nmCoD_f)
  bx_vecm_f <- bx_func(model_fit_for_vecm_f$bxs,  "VECM-CoDA", nmCoD_f)  
  bx_lc_f   <- bx_func(model_fit_for_lc_f$bx_all, "LC",        nmCoD_f)  
  bx_ct_f   <- bx_func(model_fit_for_ct_f$bxs,    "CT-CoDA",   nmCoD_f)  
  
# kt
  kt_step_m <- kt_func(model_fit_for_2s_m$kt_in_all,    "2S-CoDA",   nmCoD_m)
  kt_vecm_m <- kt_func(model_fit_for_vecm_m$kt_all_for, "VECM-CoDA", nmCoD_m)  
  kt_lc_m   <- kt_func(model_fit_for_lc_m$kt_all,       "LC",        nmCoD_m)  
  kt_step_f <- kt_func(model_fit_for_2s_f$kt_in_all,    "2S-CoDA",   nmCoD_f)
  kt_vecm_f <- kt_func(model_fit_for_vecm_f$kt_all_for, "VECM-CoDA", nmCoD_f)  
  kt_lc_f   <- kt_func(model_fit_for_lc_f$kt_all,       "LC",        nmCoD_f)  
  
  kt_ct_f   <- as_tibble(c(model_fit_for_ct_f$kt,
                           model_fit_for_ct_f$kt_for)) %>% 
    mutate(year  = years_for, 
           model = "CT-CoDA", 
           cause = "general")
  kt_ct_m   <- as_tibble(c(model_fit_for_ct_m$kt,
                           model_fit_for_ct_m$kt_for)) %>% 
    mutate(year  = years_for, 
           model = "CT-CoDA", 
           cause = "general")
  
  bx_females <- bx_step_f %>% 
    full_join(bx_vecm_f)  %>% 
    full_join(bx_lc_f)    %>% 
    full_join(bx_ct_f)    %>% 
    mutate(sex = "Female")
  
  bx_males <- bx_step_m  %>% 
    full_join(bx_ct_m)   %>% 
    full_join(bx_lc_m)   %>% 
    full_join(bx_vecm_m) %>% 
    mutate(sex = "Male")
  
  kt_females <- kt_step_f %>% 
    full_join(kt_vecm_f)  %>% 
    full_join(kt_lc_f)    %>% 
    full_join(kt_ct_f)    %>% 
    mutate(sex = "Female")
  
  kt_males <- kt_step_m  %>% 
    full_join(kt_ct_m)   %>% 
    full_join(kt_lc_m)   %>% 
    full_join(kt_vecm_m) %>% 
    mutate(sex = "Male")

  bx <- bx_females                                          %>% 
    full_join(bx_males)                                     %>% 
    filter(cause != "Other disease")                        %>% 
    mutate(cause = as.factor(cause))                        %>% 
    mutate(sex = factor(sex, levels = c("Male", "Female"))) %>% 
    ggplot(aes(x = age, y = value, group = cause, color = cause)) + 
    geom_line(size = 1)                                           + 
    facet_grid(model ~ sex, switch = "y", scales = "free")        +
    theme_bw()                                                    +
    scale_y_continuous(breaks = pretty_breaks())                  +
    scale_color_tableau()                                         +
    theme(legend.position   = "bottom",
          strip.placement   = "horizontal",
          axis.text.x       = element_text(colour = "black", size = 12), 
          axis.text.y       = element_text(colour = "black", size = 12), 
          legend.title      = element_blank(),
          legend.text       = element_text(size = 12, face = "bold", color = "black"),
          strip.text.x      = element_text(size = 12, face = "bold", color = "black"),
          strip.text.y      = element_text(size = 12, face = "bold", color = "black"),
          axis.title        = element_blank(), 
          strip.background  = element_blank(),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width  = unit(0.8, "cm"))                    + 
    guides(size = FALSE)                                          + 
    geom_hline(yintercept = 0)  

    kt <- kt_females                                             %>% 
    full_join(kt_males)                                          %>% 
    filter(cause != "Other disease")                             %>% 
    mutate(cause = ifelse(cause == "general", "General", cause)) %>% 
    mutate(cause = as.factor(cause), 
           cause = fct_relevel(cause, "General"))                %>% 
    mutate(sex = factor(sex, levels = c("Male", "Female")))      %>% 
    ggplot(aes(x = year, y = value, group = cause, color = cause)) + 
      geom_line(size = 1)                                          + 
      facet_grid(model ~ sex, switch = "y", scales = "free")       +
      theme_bw()                                                   +
      scale_y_continuous(breaks = pretty_breaks())                 +
      scale_color_tableau()                                        +
      theme(legend.position   = "bottom",
            strip.placement   = "horizontal",
            axis.text.x       = element_text(colour = "black", size = 12), 
            axis.text.y       = element_text(colour = "black", size = 12), 
            legend.title      = element_blank(),
            legend.text       = element_text(size = 12, face = "bold", color = "black"),
            strip.text.x      = element_text(size = 12, face = "bold", color = "black"),
            strip.text.y      = element_text(size = 12, face = "bold", color = "black"),
            axis.title        = element_blank(), 
            strip.background  = element_blank(),
            legend.key.height = unit(0.5, "cm"),
            legend.key.width  = unit(0.8, "cm"))                   + 
      guides(size = FALSE)                                         + 
      geom_hline(yintercept = 0)
  
  return(list(kt = kt, 
              bx = bx))
}

# -------------------------------------------------------------------------------------------------- #
# Plot is sample comparisons + failed RMSE
# RMSE and In sample? By hands and seems to be wrong
rmse_in_sample <- function(start) {

perc_d_obs_f <- readRDS("RMSE/females_observed") %>% 
  mutate(sex = "Females")

perc_d_obs_m <- readRDS("RMSE/males_observed") %>% 
  mutate(sex = "Males")

perc_d_obs <- perc_d_obs_f %>% 
  full_join(perc_d_obs_m)

in_sample_10_years <- plot_data_f                              %>% 
  full_join(plot_data_m)                                       %>% 
  filter(model != "Observed")                                  %>% 
  full_join(perc_d_obs)                                        %>% 
  group_by(year, model, sex)                                   %>% 
  summarise(value = sum(value))                                %>% 
  ungroup()                                                    %>% 
  pivot_wider(names_from  = model, 
              values_from = value)                             %>%
  filter(year > start)                                         %>%
  group_by(sex)                                                %>% 
  summarise(`2s`   = sqrt(mean((`2S CoDA`   - Observed) ^ 2)),
            `vecm` = sqrt(mean((`VECM CoDA` - Observed) ^ 2)),
            `ct`   = sqrt(mean((`CT CoDA`   - Observed) ^ 2)),
            `lc`   = sqrt(mean((LC          - Observed) ^ 2))) %>%
  pivot_longer(-sex,
               names_to  = "model",
               values_to = "value")                            %>%
  mutate(value = value * 100000)                               %>% 
  pivot_wider(names_from  = "sex", 
              values_from = "value")

plots <- plot_data_f                                                      %>% 
  full_join(plot_data_m)                                                  %>% 
  full_join(perc_d_obs)                                                   %>% 
  mutate(sex   = factor(sex, levels = c("Males", "Females")),
         model = factor(model, levels = c("Observed", "2S CoDA", 
                                          "VECM CoDA", "CT CoDA", "LC"))) %>% 
  ggplot(aes(x = year, y = value, group = interaction(cause, model)))    +
  geom_line(aes(linetype = model, color = cause), size = 1, alpha = 0.7) +
  scale_x_continuous(breaks = seq(1965, 2018, 5))                        +
  scale_color_tableau(guide = guide_legend(ncol = 2))                    +
  expand_limits(y = 0)                                                   +
  facet_wrap(sex ~ ., scales = "free", ncol = 1)                         +
  scale_y_continuous(breaks = pretty_breaks())                           +
  theme_bw()                                                             +
  theme(legend.position   = c(0.9, 0.24),
        legend.direction  = "vertical",
        axis.text.x       = element_text(colour = "black", size = 12), 
        legend.key        = element_rect(colour = "transparent"),
        legend.background = element_rect(colour = "transparent", fill = alpha('white', 0)),
        axis.text.y       = element_text(colour = "black", size = 12),
        legend.title      = element_blank(),
        legend.text       = element_text(size = 12, face = "bold", color = "black"),
        strip.text.x      = element_text(size = 12, face = "bold", color = "black"),
        strip.text.y      = element_text(size = 12, face = "bold", color = "black"),
        axis.title        = element_blank(), 
        strip.background  = element_blank(),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width  = unit(0.8, "cm"))                             + 
  guides(size     = FALSE, 
         linetype = guide_legend(ncol = 2))                              +
  geom_vline(xintercept = start)

return(plot = plots)
}

# -------------------------------------------------------------------------------------------------- #
# Correct RMSE, on the rolling bases ts
rmse_models <- function(model, sx, data) {
  
  perc_d_obs <- readRDS(data)     %>% 
    group_by(year)                %>% 
    summarise(value = sum(value)) %>% 
    ungroup()                     %>%
    dplyr::select(value)          %>% 
    ts(start     = min(years), 
       end       = max(years) + 10, 
       frequency = 1)
  
  in_sample_10_years <- plot_data_f   %>%
    full_join(plot_data_m)            %>% 
    filter(model != "Observed")       %>% 
    group_by(year, model, sex)        %>% 
    summarise(value = sum(value))     %>% 
    ungroup()                         %>% 
    pivot_wider(names_from  = model, 
                values_from = value)
  
  fit <- in_sample_10_years      %>% 
    filter(sex == sx)            %>% 
    dplyr::select(!! sym(model)) %>% 
    ts(start     = min(years), 
       end       = max(years) + 10, 
       frequency = 1)
  
# choose only RMSE statistics
  acc <- accuracy(fit, perc_d_obs, (length(fit) - 10):length(fit))[2] 
  
  return(acc)
}
# -------------------------------------------------------------------------------------------------- #
