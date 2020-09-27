CoDa_multi_para_wei_both <- function(dx,  # dx life table deaths stacked horizontial with dimensions time x age for each stackked matrix 
                                     nr_rank, # ????
                                     x, 
                                     w_base, 
                                     age_weight_2) {
  
  # -------------------------------------------------------------------------------------------------- #
  x           <- (1:nrow(dx))
  # weighted rows (rows = years)
  w           <- w_base * (1 - w_base) ^ x
  w_scaled    <- w / sum(w)
  # -------------------------------------------------------------------------------------------------- #
  
  close_dx <- acomp(dx)
  ax       <- geometricmeanCol(close_dx)
  dx_cent  <-close_dx - ax
  clr_cent <- clr(dx_cent)
  
  # SVD: bx and kt ?? (similar to LC?)
  par <- svd.triplet(clr_cent,
                     ncp   = nr_rank,
                     row.w = rev(w_scaled),
                     col.w = age_weight_2)
  U   <- par$U
  V   <- par$V
  S   <- diag(par$vs)
  
  if(nr_rank == 1) {
    bx <- V[ ,1]
    kt <- S[1,1] * U[ ,1]
    return(list(ax       = ax,
                bx       = bx,
                kt       = kt,
                dx       = dx,
                par      = par,
                clr_cent = clr_cent))
  } else if(nr_rank == 2) {
    bx <- V[ ,1]
    kt <- S[1,1] * U[ ,1] 
    
    bx2 <- V[ ,2]
    kt2 <- S[2,2] * U[ ,2] 
    
    return(list(ax  = ax,
                bx  = bx,
                kt  = kt,
                bx2 = bx2,
                kt2 = kt2,
                dx  = dx,
                par = par))
  } else {
    bx  <- V[ ,1]
    kt  <- S[1,1] * U[ ,1] 
    
    bx2 <- V[ ,2]
    kt2 <- S[2,2] * U[ ,2]  
    
    bx3 <- V[ ,3]
    kt3 <- S[3,3] * U[ ,3]
    return(list(ax  = ax,
                bx  = bx,
                kt  = kt,
                bx2 = bx2,
                kt2 = kt2,
                bx3 = bx3,
                kt3 = kt3,
                dx  = dx,
                par = par))  
  }
}





CoDa_multi_para <- function(dx, nr_rank) {
  
  close_dx <- acomp(dx)
  ax       <- geometricmeanCol(close_dx)
  dx_cent  <- close_dx - ax
  clr_cent <- clr(dx_cent)
  
  # SVD: bx and kt
  par <- svd(clr_cent, 
             nu = nr_rank, 
             nv = nr_rank)
  
  U <- par$u
  V <- par$v
  S <- diag(par$d)
  
  if(nr_rank == 1) {
    
    bx <- V[ ,1]
    kt <- S[1,1] * U[ ,1]
    
    return(list(ax = ax,
                bx = bx,
                kt = kt,
                dx = dx,
                par = par))
  } else if(nr_rank == 2) {
    
    bx  <- V[ ,1]
    kt  <- S[1,1] * U[,1] 
    bx2 <- V[ ,2]
    kt2 <- S[2,2] * U[ ,2] 
    
    return(list(ax = ax,
                bx = bx,
                kt = kt,
                bx2 = bx2,
                kt2 = kt2,
                dx  = dx,
                par = par)) 
  } else {
    bx  <- V[ ,1]
    kt  <- S[1,1] * U[ ,1] 
    bx2 <- V[ ,2]
    kt2 <- S[2,2] * U[ ,2]
    bx3 <- V[ ,3]
    kt3 <- S[3,3] * U[ ,3]
    
    return(list(ax       = ax,
                bx       = bx,
                kt       = kt,
                bx2      = bx2,
                kt2      = kt2,
                bx3      = bx3,
                kt3      = kt3,
                dx       = dx,
                par      = par,
                clr_cent = clr_cent))
  }
}



LC_no_jump <- function(mx, t) {
  
  ln_mx     <- log(mx)
  ax        <- colMeans(ln_mx, na.rm = TRUE)
  
  # simply subtracts the mean from columns can keep this wat  
  lnmx_cent <- sweep(ln_mx, 2, ax, FUN = "-")
  
  # SVD: bx and kt
  par <- svd(lnmx_cent)
  U   <- par$u
  V   <- par$v
  S   <- diag(par$d)
  R2  <- par$d[1] ^ 2 / sum(par$d ^ 2)
  bx  <- V[ ,1] / sum(V[ ,1])
  kt  <- S[1,1] * U[ ,1] * sum(V[ ,1])
  
  # Pick an ARIMA model (order=c())
  
  fcts_arima <- compose(partial(forecast, h = t), 
                        partial(Arima, 
                                order         = c(0, 1, 0), 
                                include.drift = TRUE))
  
  fcst   <- fcts_arima(kt)
  kt_fit <- fcst$fitted
  kt_for <- fcst$mean 
  
  variability <- cumsum(par$d ^ 2 /
                          sum(par$d[1:length(par$d)] ^ 2))
  
  # projections
  logmx_proj <- matrix(c(kt_fit, kt_for), (nrow(ln_mx) + t), 1) %*% 
    t(bx)
  
  lnmx_proj2 <- sweep(logmx_proj, 2, ax, FUN = "+")
  proj       <- exp(lnmx_proj2)
  
  output <- list(mx,
                 logmx_proj,
                 bx         = bx,
                 variability,
                 kt         = c(kt, kt_for),
                 forefit_mx = proj,
                 R2         = R2)
  return(output)
}


LT_abrige_mx = function(mx, age = ages_kept, nax = 2.5) {
  
  mx   <- unlist(mx)
  n    <- c(diff(as.numeric(age)), 9999)  
  qx   <- (n * mx) / (1 + (n - nax) * mx) 
  qx   <- c(qx[-(length(qx))], 1)
  qx   <- ifelse(qx > 1, 0.998, qx)
  nage <- length(age)
  npx  <- 1 - qx
  l0   <- 100000
  lx   <- round(cumprod(c(l0, npx)))
  ndx  <- -diff(lx)
  lxpn <- lx[-1]
  nLx  <- n * lxpn + ndx * nax
  nLx  <- c(nLx[-(length(qx))], lx[length(qx)] / mx[length(qx)])
  Tx   <- c(rev(cumsum(rev(nLx[-length(nLx)]))), 0)
  Tx   <- c(Tx[-(length(qx))], nLx[length(qx)])
  lx   <- lx[1:length(age)]
  ex   <- Tx / lx
  
  data <- tibble(mx, nax, qx, lx, 
                 ndx, nLx, Tx, ex)
  
  return(data)
}
