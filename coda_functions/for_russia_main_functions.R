library("compositions")
library("FactoMineR")
library("demography")
library("urca") # For coda vector-error-correction model (VECM)
library("vars")

# -------------------------------------------------------------------------------------------------- #
CoDa_2step <- function(dx, w_base, age_weight, ih, k, years, ages) { 

  m           <- length(ages)
  s           <- k # Why. Remove later !!!!!!!!!!!!!!!! was ses in f args here and other f. No need.
  w_base1     <- w_base
# calculated age weights from the master file
  age_weight1 <- age_weight
# x = 1:number of rows in dx file
# nrow = n years
  x           <- (1:nrow(dx))
# weighted rows (rows = years)
  w           <- w_base * (1 - w_base) ^ x
  w_scaled    <- w / sum(w)
  
# -------------------------------------------------------------------------------------------------- #
# START FROM HERE
# -------------------------------------------------------------------------------------------------- #
#--------------- fit step 1 ---------------#
# returns list with parameters ax, bx, kt and ALSO dx and par(all param) and clr_cent data
  fit_one <- CoDa_multi_para_wei_both(dx           = dx,
                                      nr_rank      = 1,
                                      w_base       = w_base1,
                                      age_weight_2 = age_weight1)
  
  R2_one <- to_bind(data = fit_one,
                    indx = 1:4, 
                    var  = "vs") %>% 
    t()
  
# fit Arima model to kt
  fcst_one <- forecast(Arima(fit_one$kt,
                             order         = c(1, 1, 1) ,
                             include.drift = TRUE, 
                             method        = "ML"), h = ih)
  
  kt_fit_one <- fcst_one$fitted
  kt_for_one <- fcst_one$mean
  ax_all     <- matrix(fit_one$ax, length(ages_kept), k)
  bx_all     <- matrix(fit_one$bx, length(ages_kept), k)
  
# -------------------------------------------------------------------------------------------------- #
#------------------- step 2 -------------------# 
# initial ?
  clr_proj_fit_one <- as_tibble(matrix(fit_one$kt, 
                                       length(years), 1) %*% 
                                  t(fit_one$bx))
  
# forecast ? 
  clr_proj_for_one <- as_tibble(matrix(c(fit_one$kt, kt_for_one), 
                                       length(years) + ih, 1) %*% 
                                  t(fit_one$bx)) 
  
# years, ages, causes
  data_clr <- as.data.frame.rmult(fit_one$clr_cent)   %>%
    mutate(year  = years)                             %>% 
    pivot_longer(-year, 
                 names_to  = "inx", 
                 values_to = "value")                 %>%
    mutate(age   = str_extract(inx, "[[:digit:]]+"), 
           cause = str_remove(inx,  "_[[:digit:]]+")) %>% 
    dplyr::select(-inx)                               %>% 
    pivot_wider(names_from  = "age", 
                values_from = "value")                %>% 
    dplyr::select(-year)                              %>% 
    group_by(cause)                                   %>% 
    nest()                                            %>% 
    ungroup()
  
# -------------------------------------------------------------------------------------------------- #
# ESSENTIALLY HERE WE SHOULD CHANGE
  clr_proj_fit_one_ar <- clr_proj_fit_one %>% 
    cst_next(period = years)
  
# changes to transmute
  res_1_in <- data_clr                             %>%
    left_join(clr_proj_fit_one_ar, by = "cause")   %>%
    transmute(data_new = map2(data.x, 
                              data.y, ~  .x - .y)) %>%
    mutate(par_1 = map(data_new, ~ svd(.x)),
           d     = map(par_1, "d")) 
  
# -------------------------------------------------------------------------------------------------- #
  res_1_in <- res_1_in %>% 
    mutate(d1            = map(d,                     ~ .x[1:4]),
           R2_one        = map2(d1, d,                ~ .x ^ 2 / sum(.y ^ 2)),
           R2            = map(R2_one,                ~ .x                  %>%
                                 enframe()                                  %>%
                                 pivot_wider(names_from  = "name",
                                             values_from = "value")))       %>% 
    mutate(U_in          = map(par_1, "u"), 
           V_in          = map(par_1, "v"),
           S_in          = map(d,                     ~ diag(.x)), 
           bx_in         = map(V_in,                  ~ .x[ ,1]), 
           kt_in         = map2(S_in, 
                                U_in,                 ~ .x[1,1] * .y[ ,1]),
           bx2_in        = map(V_in,                  ~ .x[ ,2]),
           kt_2_in       = map2(S_in, 
                                U_in,                 ~ .x[2,2] * .y[ ,2])) %>% 
    mutate(kt_in_for_1   = map(kt_in,                 ~ forecast(auto.arima(.x, max.d = 1), h = ih)$mean), 
           kt_in_for     = map2(kt_in, 
                                kt_in_for_1,          ~ .x %>% c(.y)), 
           kt_2_in_for_1 = map(kt_2_in,               ~ forecast(auto.arima(.x, max.d = 1), h = ih)$mean),
           kt_2_in_for   = map2(kt_2_in, 
                                kt_2_in_for_1,        ~ .x %>% c(.y)))      %>% 
    dplyr::select(-c(kt_in_for_1, kt_2_in_for_1))                           %>% 
    mutate(bx_in         = map(bx_in,                 ~ .x %>% t()),
           bx2_in        = map(bx2_in,                ~ .x %>% t()),
           clr_proj_for_in    = pmap(list(a = kt_in_for,
                                          b = bx_in,
                                          c = kt_2_in_for,
                                          d = bx2_in),
                                     function(a, b, c, d) matrix(a, (n + ih), 1) %*% 
                                       b + matrix(c, (n + ih), 1) %*% d),
           clr_proj_fit_in    = pmap(list(a = kt_in,
                                          b = bx_in,
                                          c = kt_2_in,
                                          d = bx2_in),
                                     function(a, b, c, d) matrix(a, n, 1) %*% 
                                       b + matrix(c, n, 1) %*% d))
  
  kt_in_all <- res_1_in$kt_in_for %>% 
    set_names(nm_2)              %>% 
    bind_rows()
  
  bx_in_all <-  res_1_in$bx_in    %>% 
    set_names(nm_2)              %>% 
    bind_rows()
  
  R2_in <- res_1_in$R2 %>% 
    bind_rows()
  
# -------------------------------------------------------------------------------------------------- #  
#--- calculate forecast ---------------#
  clr_proj_for_in_all  <- res_1_in %>% 
    res_func(var = "clr_proj_for_in")
  clr_proj_fit_in_all  <- res_1_in %>% 
    res_func(var = "clr_proj_fit_in")  
  
# projections
  clr_proj_for_one_all <- clr_proj_for_one + clr_proj_for_in_all
  clr_proj_fit_one_all <- clr_proj_fit_one + clr_proj_fit_in_all
  
# Inv clr
  BK_proj_fit_one      <- clrInv(clr_proj_fit_one_all)
  BK_proj_for_one      <- clrInv(clr_proj_for_one_all)
  
# Add geometric mean
  proj_fit_one         <- BK_proj_fit_one + fit_one$ax
  proj_for_one         <- BK_proj_for_one + fit_one$ax
  
# All functions are roughly the same, will combine into one or two
  proj_for_one_ar <- as.data.frame.acomp(proj_for_one) %>%
    `*` (100000)                                       %>% 
    cst_next(period = years_for)
  
  proj_fit_one_ar <- as.data.frame.acomp(proj_fit_one) %>%
    `*` (100000)                                       %>% 
    cst_next(period = years)
  
  dx_ar  <- dx %>%
    cst_next(period = years)
  
  dx_ar1 <- as.data.frame.acomp(acomp(dx)) %>%
    cst_next(period = years)

# ----- calculate explained varation -------------#
  var_joint <- vars_func(data  = clr_proj_fit_one_ar, 
                         data2 = data_clr,
                         var   = "data")
  
  var_int   <- vars_func(data  = dplyr::select(res_1_in, clr_proj_fit_in) %>% 
                           mutate(cause = nm_2),
                         data2 = data_clr,
                         var   = "clr_proj_fit_in")
  
  var_res   <- 1 - var_joint - var_int
  exp_var   <- cbind(var_joint, var_int, var_res) %>% 
    set_names(c("var_joint", "var_int", "var_res"))
  
# Life table setup
  a_one     <- c(rep(2.5, length(ages))) 
  n_one     <- c(rep(5, length(ages)))
  radix_one <- 100000
  
  proj_for_one_tot <- proj_for_one_ar %>% 
    for_proj_lt()
  
  for_Lx_one <- LT_get_Lx(data = proj_for_one_tot, 
                          radix = radix_one, 
                          a = a_one, 
                          n = n_one)$Lx
  
  for_mx_one <- proj_for_one_ar %>%
    ungroup() %>% 
    transmute(cause      = nm_2,
              for_mx_one = map(data, ~ .x / for_Lx_one))
  
  return(list(dx         = dx,
              R2         = R2_one,
              R2_in      = R2_in,
              axs        = ax_all,
              bxs        = bx_all,
              kt         = fit_one$kt,
              kt_for     = kt_for_one,
              kt_in_all  = kt_in_all,
              bx_in_all  = bx_in_all, 
              dx_forcast = proj_for_one_ar,
              BK_fit     = clr_proj_fit_one_all,
              dx_fit     = proj_fit_one_ar,
              dx_obs     = dx_ar,
              exp_var    = exp_var,
              for_mx_one = for_mx_one))
}

# -------------------------------------------------------------------------------------------------- #
# CT-CoDa function 
CoDa_CT <- function(dx, ih, k, years, ages, radix_one) {
# dx life table deaths stacked horizontaly with dimensions time x age for each stackked matrix 
# years = The years included in fitting period
# ages = single ages    
  m       <- length(ages)
  fit_one <- CoDa_multi_para(dx      = dx, 
                             nr_rank = 3) 
  
  R2_one  <- to_bind(data = fit_one, 
                     indx = 1:4, 
                     var  = "vs") %>% 
    t()
  
# fit Arima model to kt
  fcst_one   <- arim_func(data = fit_one, 
                          var = "kt", 
                          ord = c(1, 1, 1))
  
  kt_fit_one <- fcst_one$fitted
  kt_for_one <- fcst_one$mean
  
  fcst_two   <- arim_func(data = fit_one,
                          var = "kt2", 
                          ord = c(1, 0, 1))
  kt_fit_two <- fcst_two$fitted
  kt_for_two <- fcst_two$mean
  
  fcst_3     <- arim_func(data = fit_one,
                          var = "kt3", 
                          ord = c(1, 0, 1))
  kt_fit_3   <- fcst_3$fitted
  kt_for_3   <- fcst_3$mean
  
  ax_all  <- matrix(fit_one$ax,  length(ages), k)
  bx_all  <- matrix(fit_one$bx,  length(ages), k)
  bx2_all <- matrix(fit_one$bx2, length(ages), k)
  
# projections
  clr_proj_fit_one <- matrix(fit_one$kt, length(years), 1) %*% 
    t(fit_one$bx)                                           + 
    matrix(fit_one$kt2, length(years), 1)                  %*% 
    t(fit_one$bx2)                                          + 
    matrix(fit_one$kt3,length(years), 1)                   %*% 
    t(fit_one$bx3)
  
  clr_proj_for_one <- matrix(c(fit_one$kt, kt_for_one), length(years) + ih, 1) %*% 
    t(fit_one$bx)                                                               + 
    matrix(c(fit_one$kt2, kt_for_two), length(years) + ih, 1)                  %*% 
    t(fit_one$bx2) + matrix(c(fit_one$kt3, kt_for_3),
                            length(years) + ih, 1)                             %*% 
    t(fit_one$bx3)
  
# Inv clr
  BK_proj_fit_one <- clrInv(clr_proj_fit_one)
  BK_proj_for_one <- clrInv(clr_proj_for_one)
  
# Add geometric mean
  proj_fit_one    <- BK_proj_fit_one + fit_one$ax
  proj_for_one    <- BK_proj_for_one + fit_one$ax
  
  proj_for_one_ar <- as.data.frame.acomp(proj_for_one) %>% 
    `*` (100000)                                       %>%
    cst_next(period = years_for)
  proj_fit_one_ar <- as.data.frame.acomp(proj_fit_one) %>% 
    `*` (100000)                                       %>%
    cst_next(period = years)
  
  dx_ar  <- as_tibble(dx)  %>%
    cst_next(period = years)
  dx_ar1 <- as.data.frame.acomp(acomp(dx)) %>%
    cst_next(period = years)
  
# Life table setup
# should make causes arbitrary !!!!!!
  proj_for_one_tot <- proj_for_one_ar %>% 
    for_proj_lt()
  
  a_one      <- c(rep(2.5, length(ages))) 
  n_one      <- c(rep(5, length(ages)))
  radix_one  <- 100000
  for_Lx_one <- LT_get_Lx(proj_for_one_tot,
                          radix_one,
                          a_one,
                          n_one)$Lx
  
  proj_for_one_ar <- proj_for_one_ar %>%
    mutate(for_mx_one = map(data, ~ .x / for_Lx_one))
  
  return(list(dx               = dx,
              R2               = R2_one,
              axs              = ax_all,
              bxs              = bx_all,
              kt               = fit_one$kt,
              kt_for           = kt_for_one,
              kt2              = fit_one$kt2,
              kt_for2          = kt_for_two, 
              dx_forcast       = proj_for_one_ar,
              bx2_all          = bx2_all,
              dx_fit           = proj_fit_one_ar,
              dx_obs           = dx_ar,
              sv               = fit_one$par$d,
              BK_obs           = fit_one$clr_cent,
              BK_fit           = clr_proj_fit_one,
              for_mx_one       = proj_for_one_ar,
              proj_for_one_tot = proj_for_one_tot))
}

# -------------------------------------------------------------------------------------------------- #
CoDa_vecm <- function(dx, ih, k, years, ages, ir) { 
# - kt_for
# dx life table deaths stackked horizontial with dimensions time x age for each stackked matrix 
# years = The years included in fitting period
# ages = single ages    
  m             <- length(ages)
  close_dx      <- acomp(dx)
  ax            <- geometricmeanCol(close_dx)
  dx_cent       <- close_dx - ax
  clr_cent      <- as.data.frame.rmult(clr(dx_cent))
  
  clr_cent_ar_2 <- clr_cent                                 %>%
    cst_next(period = years)                                %>% 
    mutate(par = map(data,  ~ svd(.x, nu = 2, nv = 2)), 
           U   = map(par,   ~ .x$u), 
           V   = map(par,   ~ .x$v),
           S   = map(par,   ~ diag(.x$d)), 
           bx  = map(V,     ~ .x[ ,1]),
           kt  = map2(S, U, ~ .x[1,1] * .y[ ,1]), 
           bx2 = map(V,     ~ .x[ ,2]),
           kt2 = map2(S, U, ~ .x[2,2] * .y[ ,2]), 
           d   = map(par,   ~ .x$d), 
           R2  = map(d,     ~ .x[1:3] ^ 2 / sum(.x ^ 2) %>%
                       t()                              %>%
                       as_tibble()))
  
  kt_all   <- all_func(data     = clr_cent_ar_2, 
                       variable = "kt")
  bx_all   <- all_func(data     = clr_cent_ar_2, 
                       variable = "bx")
  kt_all_2 <- all_func(data     = clr_cent_ar_2, 
                       variable = "kt2")
  bx_all_2 <- all_func(data     = clr_cent_ar_2, 
                       variable = "bx2")
  
  R2       <- clr_cent_ar_2 %>%
    ungroup() %>% 
    dplyr::select(R2)       %>%
    unnest(cols = R2)
  
  ax_all   <- matrix(ax, length(ages), k)
  jo_test  <- jo_func(kt_all)
  
# projections
# HERE sometimes function vec2var cannot coerce NaN and complex numbers
  vecvar     <- vec2var(jo_test, r = ir) 
  
  for_kt_all <- predict(vecvar, n.ahead = ih)$fcst %>%
    map(~ as_tibble(.x))
  
# definitely check
  kt_all_for <- kt_all                %>% 
    full_join(map(for_kt_all, "fcst") %>% 
                bind_cols()           %>% 
                set_names(nm_2), 
              by = nm_2) 
  
# unknown attribute demands suppress warning ?????? bug
  kt_all_for_2 <- kt_all_2                                                                %>% 
    pivot_longer(everything(), 
                 names_to  = "cause", 
                 values_to = "val")                                                       %>% 
    nest(data = val)                                                                      %>% 
    mutate(data2 = map(data, ~ as_tibble(as.numeric(forecast(auto.arima(.x, max.d = 0), h = ih)$mean))),
           data3 = map2(data, 
                        data2, ~ suppressWarnings(.x %>%
                                                    full_join(.y, by = c("val" = "value"))))) %>% 
    dplyr::select(cause, data3)                                                           %>% 
    unnest(cols = data3)                                                                  %>% 
    group_by(cause)                                                                       %>% 
    mutate(row = row_number())                                                            %>% 
    pivot_wider(names_from = "cause", 
                values_from = "val")                                                      %>% 
    dplyr::select(-row)                                                                   %>% 
    ungroup()
  
  clr_data <- tibble(cause            = nm_2, 
                     clr_proj_for_one = NA, 
                     clr_proj_fit_one = NA) %>%
    mutate(clr_proj_for_one  = map(cause, ~ b_create(data1 = kt_all_for,
                                                     data2 = kt_all_for_2,
                                                     data3 = bx_all, 
                                                     data4 = bx_all_2,
                                                     indx = .x)), 
           clr_proj_fit_one = map(cause,  ~ b_create(data1 = kt_all,
                                                     data2 = kt_all_2,
                                                     data3 = bx_all, 
                                                     data4 = bx_all_2,
                                                     indx = .x)))
  
  clr_proj_for_one_1 <- for_clr(data     = clr_data, 
                                variable = "clr_proj_for_one")
  clr_proj_fit_one_1 <- for_clr(data     = clr_data, 
                                variable = "clr_proj_fit_one")
  
# Inv clr
  BK_proj_for_one <- clrInv(clr_proj_for_one_1)
  BK_proj_fit_one <- clrInv(clr_proj_fit_one_1)
  
# Add geometric mean
  proj_for_one <- BK_proj_for_one + ax
  proj_fit_one <- BK_proj_fit_one + ax  

  proj_for_one_ar <- as.data.frame.acomp(proj_for_one) %>%
    `*` (100000)                                       %>% 
    cst_next(period = years_for)
  
  proj_fit_one_ar <- as.data.frame.acomp(proj_fit_one)  %>%
    `*` (100000)                                        %>% 
    cst_next(period = years)
  
  dx_ar  <- as_tibble(dx) %>%
    cst_next(period = years)
  
  dx_ar1 <- as.data.frame.acomp(acomp(dx)) %>%
    cst_next(period = years)
  
# Life table setup
  a_one     <- c(rep(2.5, length(ages))) 
  n_one     <- c(rep(5, length(ages)))
  radix_one <- 100000
  
  proj_for_one_tot <- proj_for_one_ar %>% 
    for_proj_lt()
  
  for_Lx_one <- LT_get_Lx(proj_for_one_tot, 
                          radix_one, 
                          a_one, 
                          n_one)$Lx
  
  for_mx_one <- proj_for_one_ar %>%
    mutate(for_mx_one = map(data, ~ .x / for_Lx_one))
  
  return(list(dx           = dx,
              axs          = ax_all,
              bxs          = bx_all,
              kt           = kt_all,
              kt_all_for   = kt_all_for,
              R2           = R2,
              kt_all_for_2 = kt_all_for_2,
              bx_all_2     = bx_all_2,
              dx_fit       = proj_fit_one_ar,
              BK_fit       = clr_data[ ,-2],
              dx_forcast   = proj_for_one_ar,
              dx_obs       = dx_ar,
              jo_test      = jo_test,
              for_mx_one   = for_mx_one))
}

# -------------------------------------------------------------------------------------------------- #
COD_LC <- function(mx, t, k, ages, years) {
  
  m <- length(ages)
  
  lc_temp <- mx                           %>% 
    mutate(data         = map(data,  ~ .x %>%
                                t()))        %>% 
    mutate(data2        = map(data,  ~ .x %>%
                                LC_no_jump(t = t)), 
           LC_mx_for_ar = map(data2, ~ .x$forefit_mx), 
           kt           = map(data2, ~ .x$kt), 
           bx           = map(data2, ~ .x$bx), 
           R2           = map(data2, ~ .x$R2)) %>% 
    ungroup()

  for_all <- function(data, variable) {
    
    quo_variable <-  quo(variable)
    
    dat          <- lc_temp                  %>%
      dplyr::select(cause, all_of(variable)) %>%
      unnest(cols = !! quo_variable)         %>%
      group_by(cause)                        %>%
      mutate(row  = row_number())            %>%
      pivot_wider(names_from  = "cause",
                  values_from = variable)    %>%
      dplyr::select(-row) %>%
      ungroup()
    
    return(dat)
  }
  
  kt_all <- for_all(data     = lc_temp, 
                    variable = "kt")
  bx_all <- for_all(data     = lc_temp, 
                    variable = "bx")
  R2     <- for_all(data     = lc_temp, 
                    variable = "R2")
  
  LC_mx_for_all <- lc_temp      %>% 
    ungroup()                   %>% 
    dplyr::select(LC_mx_for_ar) %>% 
    t()                         %>% 
    reduce(`+`)
  
  lc_temp <- lc_temp %>% 
    mutate(q_i =  map(LC_mx_for_ar, ~ (.x * 5) / (1 + (5 - 2.5) * LC_mx_for_all)))
  
  qi_all <- lt_create(data = LC_mx_for_all) %>% 
    map("qx")                               %>% 
    set_names(years_for)                    %>% 
    bind_rows(.id = NULL)                   %>% 
    t()
  
  Lx_all <- lt_create(data = LC_mx_for_all) %>% 
    map("nLx")                              %>% 
    set_names(years_for)                    %>% 
    bind_cols(.id = NULL)                   %>%
    t()
  
  lc_temp    <- lc_temp %>%
    mutate(d_i_for_lc = map(LC_mx_for_ar, ~ .x * Lx_all))
  
  dx_forcast <- lc_temp %>%
    transmute(cause = cause, 
              data  = map(d_i_for_lc,     ~ .x %>%
                            as_tibble()))
  
  return(list(dx_forcast = dx_forcast,
              kt_all     = kt_all,
              bx_all     = bx_all,
              for_mx_one = lc_temp$LC_mx_for_ar,
              R2         = R2))   
}
# -------------------------------------------------------------------------------------------------- #
