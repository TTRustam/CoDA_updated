
lt_create_RT <- function(gender, variable, data) {
  f_dx <- lt                                        %>% 
    filter(sex == gender)                           %>% 
    unnest(data)                                    %>% 
    ungroup()                                       %>% 
    dplyr::select(year, age, all_of(variable))      %>% 
    filter(year > 1964, 
           !age %in% c( "0", "01-04", "05-09",
                        "10-14", "15-19", "20-24")) %>% 
    pivot_wider(names_from  = "year", 
                values_from = variable)             %>% 
    dplyr::select(-age)                             %>% 
    dplyr::select(all_of(as.character(years)))
  
  return(f_dx)
}


rel_dat_func <- function(dx, gender) {
  
  deaths <- tot_death     %>% 
    filter(sex == gender) %>% 
    dplyr::select(-c(sex, age))
  
  data <- m_cancer        %>%
    filter(sex == gender) %>%
    mutate(data = map(data, ~ .x            %>% 
                        dplyr::select(-age) %>% 
                        `/` (deaths)        %>%
                        `*` (dx)            %>% 
                        t()))
  
  return(data)  
}

com_func <- function(dat) {
  
  dat <- dat                        %>%
    mutate(data = map(data, ~ .x    %>%
                        t()         %>% 
                        as_tibble() %>% 
                        set_names(years))) %>% 
    unnest(cols = data)                    %>%
    group_by(sex, cause)                   %>% 
    mutate(age = ages_kept)                %>% 
    unite("ind", c("cause", "age"))        %>%
    ungroup()                              %>% 
    pivot_longer(-c(sex, ind), 
                 names_to   = "year",
                 values_to  = "value")     %>% 
    pivot_wider(names_from  = "ind", 
                values_from = "value")     %>% 
    dplyr::select(-c(sex, year))
  
  return(dat)
}

tot_func <- function(dat) {
  
  dat <- dat                                  %>%
    mutate(data = map(data, ~ .x    %>%
                        as_tibble() %>% 
                        set_names(ages_kept))) %>% 
    unnest(cols = data)                        %>% 
    group_by(sex, cause)                       %>% 
    mutate(year = years)                       %>% 
    group_by(sex, year)                        %>% 
    summarise_at(vars(-cause), sum)            %>% 
    ungroup()                                  %>% 
    dplyr::select(-c(sex, year))
  
  return(dat)
}

col_func <- function(dat, means) {
  
  
  dat   <- dat                     %>%
    mutate(data = map(data, ~ .x   %>% 
                        colMeans() %>% 
                        `/` (sum(means)))) %>% 
    unnest(cols = data)                          %>% 
    pull(data)
  
  return(dat)
}

cause_func <- function(dat, Lx){
  
  
  dat <- dat %>%
    mutate(data = map(data, ~ t(.x) / Lx)) %>% 
    ungroup() %>% 
    dplyr::select(-sex)
  
  return(dat)
}

neo_func <- function(dat) {
  
  dat <- dat                                     %>% 
    mutate(data = map(data, ~ as_tibble(.x)))    %>% 
    unnest(cols = data)                          %>%
    filter(cause != "Other disease")             %>% 
    group_by(cause)                              %>% 
    mutate(year = years)                         %>% 
    group_by(year)                               %>% 
    dplyr::select(-c(cause, sex))                %>%
    summarise(across(c(starts_with("V")), sum))  %>% 
    ungroup()                                    %>% 
    dplyr::select(-year)                         %>%
    acomp()
  
  return(dat)
}

plot_prep <- function(dat, dat2) {
  
  dat <- as.data.frame.acomp(dat)   %>% 
    set_names(ages_kept)            %>%
    mutate(year = years)            %>% 
    pivot_longer(-year, 
                 names_to  = "age", 
                 values_to = "val") %>% 
    filter(year %in% c(1965, 1980, 
                       2000, 2018)) %>% 
    mutate(year = as.factor(year))
  
  geom_line <- as_tibble(dat2) %>% 
    mutate(age = ages_kept, 
           meh = as.factor(111))
  
  return(list(distr     = dat, 
              geom_line = geom_line))
}

neo_plot <- function(dat) {
  
  dat <- as.data.frame.acomp(dat)   %>% 
    set_names(ages_kept)            %>%
    mutate(year = years)            %>% 
    pivot_longer(-year, 
                 names_to  = "age", 
                 values_to = "val") %>% 
    filter(year %in% c(1965, 1980, 
                       2000, 2018)) %>% 
    mutate(year = as.factor(year))
  
  return(dat)
}

to_bind <- function(data, indx, var) {
  
  data <- data$par[[eval(var)]][indx] ^ 2 / 
    sum(  data$par[[eval(var)]]       ^ 2)
  
  return(data)
}

cst_next <- function(data, period) {
  data <- data                                  %>%
    set_names(nm_1)                             %>%
    mutate(year = period)                       %>%
    pivot_longer(-year,
                 names_to  = "cause",
                 values_to = "val")             %>%
    mutate(age   = cause,
           cause = str_remove(cause, "[0-9]+")) %>%
    group_by(cause)                             %>%
    nest()                                      %>%
    mutate(data = map(data, ~ .x                         %>%
                        pivot_wider(names_from  = "age",
                                    values_from = "val") %>%
                        dplyr::select(-year)))           %>%
    mutate(data = map(data, ~ .x %>%
                        select_if(~ !(all(is.na(.))))))

  return(data)
}

res_func <- function(data, var) {
  
  quo_var <- quo(var)
  enq_var <- enquo(var)
  map_var <- rlang::sym(var)
  
  data <- data                       %>%
    dplyr::select(!! quo_var)        %>%
    mutate(!! var := map(!! map_var, ~ .x  %>%
                           as_tibble()     %>% 
                           set_names(ages_kept)),
           cause   = nm_2)           %>%
    unnest(cols = !! quo(var))       %>%
    pivot_longer(-cause, 
                 names_to  = "age",
                 values_to = "val")  %>%
    group_by(cause, age)             %>%
    mutate(n = row_number())         %>%
    unite(indx, c("cause", "age"))   %>%
    pivot_wider(names_from  = "indx",
                values_from = "val") %>%
    dplyr::select(-n)                %>%
    ungroup()
  
  return(data)
}

vars_func <- function(data,
                      data2, 
                      var) {
  
  map_var <- rlang::sym(var)
  
  data    <- data                          %>% 
    full_join(rename(data2, data2 = data)) %>% 
    mutate(data  = map(!! map_var, ~ .x    %>% 
                        as.matrix()        %>% 
                        norm(type = "F")   %>% 
                        `^` (2)),
           data2 = map(data2, ~ .x         %>%
                         norm(type = "F")  %>%
                         `^` (2)),
           data3 = map2(data,
                        data2, ~ .x / .y)) %>% 
    ungroup()                              %>% 
  dplyr::select(data3)                     %>% 
  unnest(cols = data3)
  
  return(data)
}

for_proj_lt <- function(data) {
  
  data <- data                                 %>% 
    pivot_wider(., 
                names_from  = "cause", 
                values_from = "data")          %>% 
    unnest(everything())                       %>% 
    mutate(year = years_for)                   %>% 
    pivot_longer(-year, 
                 names_to  = "cause", 
                 values_to = "val")            %>% 
    mutate(age = str_extract(cause, "[0-9]+")) %>% 
    group_by(year, age)                        %>% 
    summarise(x = sum(val))                    %>% 
    pivot_wider(names_from  = "age", 
                values_from = "x")             %>% 
    ungroup()                                  %>% 
    dplyr::select(-year)
  
  return(data)
}

LT_get_Lx <- function(data, radix, a, n) {
  
  lx <- apply(data[, -ncol(data)], 1, cumsum) %>% 
    t()                                       %>% 
    as_tibble()                               %>% 
    mutate_all(~ (radix - .x))                %>%
    add_column(x       = radix, 
               .before = 1)                   %>%
    set_names(ages_kept)
  
  Lx <- (lx * n) - (data * (n - a))
  
  Tx <- t(Lx)                         %>% 
    as_tibble()                       %>% 
    mutate_all(~ rev(cumsum(rev(.)))) %>% 
    t()                               %>% 
    as_tibble()                       %>% 
    replace(. == 0, NA)               %>% 
    set_names(ages_kept)
  
  return(list(dx = data, 
              lx = lx, 
              Lx = Lx))
}

perc_d_func <- function(data) {
  
  data <- data                                %>%
    mutate(data = map(data, ~ .x %>%
                        apply(1, sum)))       %>%
    mutate(data = map(data, ~ as_tibble(.x))) %>%
    unnest(cols = data)                       %>%
    group_by(cause)                           %>%
    mutate(value = value / 100000,
           n     = row_number())              %>%
    pivot_wider(names_from  = "cause",
                values_from = "value")        %>%
    dplyr::select(-n)                         %>% 
    
    ungroup()
  
  return(data)
}

arim_func <- function(data, var, ord) {
  
  order_one <- ord
  fcs       <- forecast(auto.arima(data[[eval(var)]], max.d = 1), h = ih)
  
  return(fcs)
}

# This thing causes problems
jo_func <- partial(ca.jo, type  = "trace",
                   ecdet        = "trend",
                   K            = 2,
                   spec         = "transitory")

all_func <- function(data, variable) {
  
  quo_variable <- quo(variable)
  
  result       <- data                    %>%
    dplyr::select(cause, !! quo_variable) %>%
    unnest(cols =  !! quo_variable)       %>%
    group_by(cause) %>%
    mutate(row = row_number())            %>%
    pivot_wider(names_from  = "cause",
                values_from = variable)   %>%
    dplyr::select(-row)                   %>% 
    set_names(nm_2)                       %>% 
    ungroup()
  
  return(result)
}

b_create <- function(data1, data2, 
                     data3, data4, indx) {
  
  data <- pmap_dfc(list(a = data1[ ,indx],
                        b = t(data3[ ,indx]),
                        c = data2[ ,indx],
                        d = t(data4[ ,indx])),
                   function(a, b, c, d) a * b + c * d)
  return(data)
}

for_clr <- function(data, variable) {
  
  quo_var <- quo(variable)
  
  data <- data                          %>%
    dplyr::select(cause, !! quo_var)    %>% 
    unnest(cols = !! quo_var)           %>% 
    pivot_longer(-cause, 
                 names_to  = "col",
                 values_to = "val")     %>% 
    unite(cause_col, c("cause", "col")) %>% 
    group_by(cause_col)                 %>% 
    mutate(row = row_number())          %>% 
    pivot_wider(names_from  = "cause_col", 
                values_from = "val")    %>%
    dplyr::select(-row)                 %>% 
    ungroup()
  
  return(data)
}

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

lt_create <- function(data) {
  
  data1 <- data                  %>% 
    as_tibble()                  %>%
    rowwise()                    %>% 
    group_split()                %>%
    map(~ LT_abrige_mx(mx = .x)) 
  
  return(data1)
}
