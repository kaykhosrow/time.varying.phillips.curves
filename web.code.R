############################################################
## Time-Varying Phillips Curve – Euro Area
## Multi-bandwidth, multi-kernel, OLS + IV export for website
############################################################

############################################################
## 0. PACKAGES
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(grid)   # for unit()

############################################################
## 1. SETTINGS
############################################################

## Bandwidth grid: h1 in {0, 0.1, ..., 1.0}
H1_GRID <- seq(0, 1, by = 0.1)

## Kernels to use
KERNELS <- c("gaussian", "epanechnikov", "rolling")

## Confidence level for CIs
CI_LEVEL <- 0.90
Z_CRIT   <- qnorm(0.5 + CI_LEVEL / 2)

############################################################
## 2. READ DATA AND BUILD PHILLIPS VARIABLES
############################################################

## Assumptions:
## - input.csv in working directory
## - Column 1: DATE (dd/mm/yyyy)
## - Column 2: ignored
## - Column 3: Unemployment rate
## - Column 4: CPI index

df_raw <- read.csv("input.csv", stringsAsFactors = FALSE)

date_raw <- df_raw[[1]]
U_raw    <- df_raw[[3]]
CPI_raw  <- df_raw[[4]]

date <- dmy(date_raw)
U    <- as.numeric(U_raw)
CPI  <- as.numeric(CPI_raw)

## Helper lag/lead (simple vector versions)
lag_vec  <- function(x, k = 1) c(rep(NA, k), x[1:(length(x) - k)])
lead_vec <- function(x, k = 1) c(x[(k + 1):length(x)], rep(NA, k))

## Inflation and changes
pi        <- 100 * (log(CPI) - lag_vec(log(CPI), 1))  # inflation
dpi       <- pi - lag_vec(pi, 1)                      # Δπ
dpi_lag1  <- lag_vec(dpi, 1)
dpi_lead1 <- lead_vec(dpi, 1)
dpi_l2    <- lag_vec(dpi, 2)
dpi_l3    <- lag_vec(dpi, 3)
dpi_l4    <- lag_vec(dpi, 4)

## Unemployment changes
du    <- U - lag_vec(U, 1)
du_l1 <- lag_vec(du, 1)
du_l2 <- lag_vec(du, 2)
du_l3 <- lag_vec(du, 3)
du_l4 <- lag_vec(du, 4)

## Assemble and drop NA (due to lags/leads)
df <- data.frame(
  date      = date,
  U         = U,
  CPI       = CPI,
  pi        = pi,
  dpi       = dpi,
  dpi_lag1  = dpi_lag1,
  dpi_lead1 = dpi_lead1,
  dpi_l2    = dpi_l2,
  dpi_l3    = dpi_l3,
  dpi_l4    = dpi_l4,
  du        = du,
  du_l1     = du_l1,
  du_l2     = du_l2,
  du_l3     = du_l3,
  du_l4     = du_l4
)

df <- na.omit(df)
T_full <- nrow(df)
cat("Effective sample after lags/leads:", T_full, "observations\n\n")

############################################################
## 3. BASIC MATRIX HELPERS
############################################################

vech_lower <- function(M) {
  M <- as.matrix(M)
  M[lower.tri(M, diag = TRUE)]
}

unvech_sym <- function(v, k) {
  M <- matrix(0, k, k)
  idx <- 1
  for (j in 1:k) {
    for (i in j:k) {
      M[i, j] <- v[idx]
      M[j, i] <- v[idx]
      idx <- idx + 1
    }
  }
  M
}

hdprod_rows <- function(X) {
  X <- as.matrix(X)
  T <- nrow(X)
  k <- ncol(X)
  m <- k * (k + 1) / 2
  out <- matrix(NA_real_, nrow = T, ncol = m)
  for (t in 1:T) {
    xt   <- X[t, , drop = FALSE]
    out[t, ] <- vech_lower(t(xt) %*% xt)
  }
  out
}

## SAFE symmetric matrix square root
mat_sqrt_sym <- function(A) {
  A <- as.matrix(A)
  if (any(!is.finite(A))) {
    return(matrix(0, nrow(A), ncol(A)))
  }
  A <- 0.5 * (A + t(A))
  ev <- eigen(A, symmetric = TRUE)
  vals <- pmax(ev$values, 0)
  if (all(vals == 0)) {
    return(matrix(0, nrow(A), ncol(A)))
  }
  P <- ev$vectors
  P %*% (diag(sqrt(vals), length(vals))) %*% t(P)
}

############################################################
## 4. GENERIC KERNEL WEIGHTS (GAUSSIAN, EPANECHNIKOV, ROLLING)
############################################################

kernel_weights <- function(T, H, kernel = c("gaussian", "epanechnikov", "rolling")) {
  kernel <- match.arg(kernel)
  
  seq_vals <- -T:(T - 1)
  
  if (kernel == "gaussian") {
    base <- exp(-0.5 * (abs(seq_vals) / H)^2)
  } else if (kernel == "epanechnikov") {
    u <- seq_vals / H
    base <- ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
  } else if (kernel == "rolling") {
    base <- ifelse(abs(seq_vals) <= H, 1, 0)
  } else {
    stop("Unknown kernel type: ", kernel)
  }
  
  W_mat <- matrix(NA_real_, nrow = T, ncol = T)
  ini <- T + 1
  fin <- 2 * T
  for (t in 1:T) {
    W_mat[t, ] <- base[ini:fin]
    ini <- ini - 1
    fin <- fin - 1
  }
  W_mat
}

############################################################
## 5. KERNEL LOOP (TIME-VARYING LS) – GENERIC KERNEL
############################################################

kernel_loop_tv <- function(y, X, Z, H, kernel = c("gaussian", "epanechnikov", "rolling")) {
  kernel <- match.arg(kernel)
  
  y <- as.matrix(y)  # T x q
  X <- as.matrix(X)  # T x k
  Z <- as.matrix(Z)  # T x p
  
  T <- nrow(y)
  q <- ncol(y)
  k <- ncol(X)
  
  W_mat <- kernel_weights(T, H, kernel)
  ret   <- matrix(NA_real_, nrow = T, ncol = q * k)
  
  for (t in 1:T) {
    W   <- W_mat[t, ]
    PX  <- Z * W
    num <- crossprod(PX, y)
    den <- crossprod(PX, X)
    
    B <- tryCatch(
      solve(den, num),
      error = function(e) matrix(NA_real_, nrow = k, ncol = q)
    )
    
    ret[t, ] <- as.vector(B)
  }
  
  ret
}

############################################################
## 6. TV-OLS (GENERIC KERNEL, SAFE)
############################################################

tv_OLS_tv <- function(y, X, h1, kernel = c("gaussian", "epanechnikov", "rolling")) {
  kernel <- match.arg(kernel)
  
  y <- as.numeric(y)
  X <- as.matrix(X)
  T <- length(y)
  k <- ncol(X)
  
  H  <- T^h1
  y_mat <- matrix(y, ncol = 1)
  
  bethat <- kernel_loop_tv(
    y      = y_mat,
    X      = X,
    Z      = X,
    H      = H,
    kernel = kernel
  )
  
  y_hat <- rowSums(bethat * X)
  u     <- y - y_hat
  
  xx    <- hdprod_rows(X)
  u2    <- u^2
  XXu   <- xx * u2
  W_mat <- kernel_weights(T, H, kernel)
  
  SE   <- matrix(NA_real_, nrow = T, ncol = k)
  VCV  <- matrix(NA_real_, nrow = T, ncol = k * (k + 1) / 2)
  
  for (t in 1:T) {
    W  <- W_mat[t, ]
    K  <- sum(W)
    K2 <- sum(W^2)
    
    if (!is.finite(K) || K <= 0 || !is.finite(K2) || K2 <= 0) next
    
    XXuW <- as.vector(crossprod(XXu, W))
    xxW  <- as.vector(crossprod(xx,  W))
    
    if (any(!is.finite(XXuW)) || any(!is.finite(xxW))) next
    
    Szu <- unvech_sym(XXuW / K, k)
    Szz <- unvech_sym(xxW  / K, k)
    
    iSzz <- tryCatch(
      solve(Szz),
      error = function(e) NULL
    )
    if (is.null(iSzz)) next
    
    COV_t <- (K2 / (K^2)) * iSzz %*% Szu %*% iSzz
    if (any(!is.finite(COV_t))) next
    
    VCV[t, ] <- vech_lower(COV_t)
    SE[t, ]  <- sqrt(pmax(diag(COV_t), 0))
  }
  
  list(
    coeff   = bethat,
    stderr  = SE,
    vcv     = VCV,
    yhat    = y_hat,
    uhat    = u,
    H       = H,
    h1      = h1,
    kernel  = kernel
  )
}

############################################################
## 7. TV-IV WITH ROBUST SE + TIME-VARYING HAUSMAN (GENERIC, SAFE)
############################################################

build_Pt_init <- function(covar_names, inst_names) {
  q <- length(inst_names)
  k <- length(covar_names)
  Pt <- matrix(0, nrow = q, ncol = k)
  rownames(Pt) <- inst_names
  colnames(Pt) <- covar_names
  for (j in 1:k) {
    pos <- which(inst_names == covar_names[j])
    if (length(pos) > 0) {
      Pt[pos[1], j] <- 1
    }
  }
  Pt
}

tv_IV_tv <- function(y_vec, COVAR, INST,
                     endog_idx,
                     h1, h2 = h1,
                     kernel = c("gaussian", "epanechnikov", "rolling")) {
  kernel <- match.arg(kernel)
  
  y_vec  <- as.numeric(y_vec)
  X_cov  <- as.matrix(COVAR)
  Z_inst <- as.matrix(INST)
  
  T  <- length(y_vec)
  k  <- ncol(X_cov)
  q  <- ncol(Z_inst)
  n_endo <- length(endog_idx)
  if (n_endo == 0) stop("No endogenous variables supplied (endog_idx empty).")
  
  H  <- T^h1
  H2 <- T^h2
  
  covar_names <- colnames(X_cov)
  inst_names  <- colnames(Z_inst)
  if (is.null(covar_names)) {
    covar_names <- paste0("x", seq_len(k))
    colnames(X_cov) <- covar_names
  }
  if (is.null(inst_names)) {
    inst_names <- paste0("z", seq_len(q))
    colnames(Z_inst) <- inst_names
  }
  Pt_init <- build_Pt_init(covar_names, inst_names)
  
  Xe <- X_cov[, endog_idx, drop = FALSE]
  
  psihat_full <- kernel_loop_tv(
    y      = Xe,
    X      = Z_inst,
    Z      = Z_inst,
    H      = H2,
    kernel = kernel
  )
  
  xfithat <- X_cov
  for (j in seq_len(n_endo)) {
    cols_j <- ((j - 1) * q + 1):(j * q)
    psi_j  <- psihat_full[, cols_j, drop = FALSE]
    xfithat[, endog_idx[j]] <- rowSums(psi_j * Z_inst)
  }
  
  y_mat <- matrix(y_vec, ncol = 1)
  beta_iv <- kernel_loop_tv(
    y      = y_mat,
    X      = X_cov,
    Z      = xfithat,
    H      = H,
    kernel = kernel
  )
  
  y_hat_iv <- rowSums(beta_iv * X_cov)
  u        <- y_vec - y_hat_iv
  u2       <- u^2
  
  beta_ols <- kernel_loop_tv(
    y      = y_mat,
    X      = X_cov,
    Z      = X_cov,
    H      = H,
    kernel = kernel
  )
  
  b_e_ols <- beta_ols[, endog_idx, drop = FALSE]
  b_e_iv  <- beta_iv[,  endog_idx, drop = FALSE]
  
  Xe   <- X_cov[, endog_idx, drop = FALSE]
  Xfe  <- xfithat[, endog_idx, drop = FALSE]
  Vvec <- Xe - Xfe
  
  xx   <- hdprod_rows(Xe)
  xpxp <- hdprod_rows(Xfe)
  vv   <- hdprod_rows(Vvec)
  
  zz  <- hdprod_rows(Z_inst)
  zzu <- zz * u2
  
  W_all <- kernel_weights(T, H, kernel)
  
  SE_iv   <- matrix(NA_real_, nrow = T, ncol = k)
  VCV_iv  <- matrix(NA_real_, nrow = T, ncol = k * (k + 1) / 2)
  H_stat  <- rep(NA_real_, T)
  p_val   <- rep(NA_real_, T)
  
  for (t in 1:T) {
    W  <- W_all[t, ]
    K  <- sum(W)
    K2 <- sum(W^2)
    
    if (!is.finite(K) || K <= 0 || !is.finite(K2) || K2 <= 0) next
    
    xxW   <- as.vector(crossprod(xx,   W))
    xpxpW <- as.vector(crossprod(xpxp, W))
    vvW   <- as.vector(crossprod(vv,   W))
    
    if (any(!is.finite(xxW)) || any(!is.finite(xpxpW)) || any(!is.finite(vvW))) {
      next
    }
    
    Sxx   <- unvech_sym(xxW   / K, n_endo)
    Sxpxp <- unvech_sym(xpxpW / K, n_endo)
    Svv   <- unvech_sym(vvW   / K, n_endo)
    
    s2u <- sum(W * u2) / K
    if (!is.finite(s2u) || s2u <= 0) {
      next
    }
    
    hSxx   <- mat_sqrt_sym(Sxx)
    hSxpxp <- mat_sqrt_sym(Sxpxp)
    
    diff_b <- b_e_ols[t, ] - b_e_iv[t, ]
    Vt     <- hSxx %*% hSxpxp %*% matrix(diff_b, ncol = 1)
    
    Svv_sym <- 0.5 * (Svv + t(Svv))
    Svv_sym <- Svv_sym + diag(1e-10, n_endo)
    iSvv <- tryCatch(
      solve(Svv_sym),
      error = function(e) NULL
    )
    if (!is.null(iSvv)) {
      H_t <- ((K^2) / (K2 * s2u)) * as.numeric(t(Vt) %*% iSvv %*% Vt)
      if (is.finite(H_t) && H_t >= 0) {
        H_stat[t] <- H_t
        p_val[t]  <- stats::pchisq(H_t, df = n_endo, lower.tail = FALSE)
      }
    }
    
    m_z <- ncol(zz)
    zzuW <- as.vector(crossprod(zzu, W))
    zzW  <- as.vector(crossprod(zz,  W))
    
    if (any(!is.finite(zzuW)) || any(!is.finite(zzW))) next
    
    Szu <- unvech_sym(zzuW / K, q)
    Szz <- unvech_sym(zzW  / K, q)
    
    Pt_t <- Pt_init
    psi_t_vec <- psihat_full[t, ]
    psi_t_mat <- matrix(psi_t_vec, nrow = q, ncol = n_endo)
    for (j in seq_len(n_endo)) {
      Pt_t[, endog_idx[j]] <- psi_t_mat[, j]
    }
    
    PSzuP <- t(Pt_t) %*% Szu %*% Pt_t
    PSzzP <- t(Pt_t) %*% Szz %*% Pt_t
    
    iPSzzP <- tryCatch(
      solve(PSzzP),
      error = function(e) NULL
    )
    if (!is.null(iPSzzP)) {
      COV_t <- (K2 / (K^2)) * iPSzzP %*% PSzuP %*% iPSzzP
      if (any(is.finite(COV_t))) {
        VCV_iv[t, ] <- vech_lower(COV_t)
        SE_iv[t, ]  <- sqrt(pmax(diag(COV_t), 0))
      }
    }
  }
  
  list(
    coeff_iv   = beta_iv,
    stderr_iv  = SE_iv,
    vcv_iv     = VCV_iv,
    coeff_ols  = beta_ols,
    H_stat     = H_stat,
    p_val      = p_val,
    yhat_iv    = y_hat_iv,
    uhat_iv    = u,
    H          = H,
    h1         = h1,
    kernel     = kernel
  )
}

############################################################
## 8. BUILD CI DATA FRAMES
############################################################

build_ci_df <- function(beta_mat, se_mat, dates, coef_names) {
  beta <- as.matrix(beta_mat)
  se   <- as.matrix(se_mat)
  
  colnames(beta) <- coef_names
  colnames(se)   <- coef_names
  
  df_beta <- as.data.frame(beta)
  df_se   <- as.data.frame(se)
  
  df_beta$date <- dates
  df_se$date   <- dates
  
  df_long <- df_beta %>%
    tidyr::pivot_longer(
      cols      = -date,
      names_to  = "coef",
      values_to = "estimate"
    )
  
  df_se_long <- df_se %>%
    tidyr::pivot_longer(
      cols      = -date,
      names_to  = "coef",
      values_to = "se"
    )
  
  df_full <- df_long %>%
    dplyr::left_join(df_se_long, by = c("date", "coef")) %>%
    dplyr::mutate(
      lower = estimate - Z_CRIT * se,
      upper = estimate + Z_CRIT * se
    )
  
  df_full
}

############################################################
## 9. DEFINE STANDARD AND NK PHILLIPS CURVE REGRESSORS
############################################################

y <- df$dpi

## Standard PC: Δπ_t ~ Δπ_{t-1} + Δu_t
X_standard <- cbind(
  const    = 1,
  dpi_lag1 = df$dpi_lag1,
  du       = df$du
)

Z_standard <- cbind(
  const    = 1,
  dpi_lag1 = df$dpi_lag1,
  du_l1    = df$du_l1,
  du_l2    = df$du_l2,
  du_l3    = df$du_l3,
  du_l4    = df$du_l4,
  dpi_l2   = df$dpi_l2
)

colnames(X_standard) <- c("const", "dpi_lag1", "du")
colnames(Z_standard) <- c("const", "dpi_lag1",
                          "du_l1", "du_l2", "du_l3", "du_l4",
                          "dpi_l2")

## NK PC: Δπ_t ~ Δπ_{t+1} + Δπ_{t-1} + Δu_t
X_nk <- cbind(
  const     = 1,
  dpi_lead1 = df$dpi_lead1,
  dpi_lag1  = df$dpi_lag1,
  du        = df$du
)

Z_nk <- cbind(
  const    = 1,
  dpi_lag1 = df$dpi_lag1,
  dpi_l2   = df$dpi_l2,
  dpi_l3   = df$dpi_l3,
  dpi_l4   = df$dpi_l4,
  du_l1    = df$du_l1,
  du_l2    = df$du_l2,
  du_l3    = df$du_l3,
  du_l4    = df$du_l4
)

colnames(X_nk) <- c("const", "dpi_lead1", "dpi_lag1", "du")
colnames(Z_nk) <- c("const", "dpi_lag1",
                    "dpi_l2", "dpi_l3", "dpi_l4",
                    "du_l1", "du_l2", "du_l3", "du_l4")

coef_names_standard <- c("(Intercept)", "dpi_lag1", "du")
coef_names_nk       <- c("(Intercept)", "dpi_lead1", "dpi_lag1", "du")

############################################################
## 10. LOOP OVER BANDWIDTHS, KERNELS, ESTIMATORS
############################################################

results_list <- list()
res_idx <- 1L

for (kernel_type in KERNELS) {
  cat("Kernel:", kernel_type, "\n")
  
  for (h1 in H1_GRID) {
    cat("  h1 =", h1, "\n")
    
    ## STANDARD PC
    ols_standard <- tv_OLS_tv(
      y      = y,
      X      = X_standard,
      h1     = h1,
      kernel = kernel_type
    )
    
    iv_standard <- tv_IV_tv(
      y_vec     = y,
      COVAR     = X_standard,
      INST      = Z_standard,
      endog_idx = 3L,
      h1        = h1,
      h2        = h1,
      kernel    = kernel_type
    )
    
    df_ols_standard_ci <- build_ci_df(
      beta_mat   = ols_standard$coeff,
      se_mat     = ols_standard$stderr,
      dates      = df$date,
      coef_names = coef_names_standard
    ) %>%
      mutate(
        h1        = h1,
        kernel    = kernel_type,
        pc_type   = "standard",
        estimator = "ols"
      )
    
    df_iv_standard_ci <- build_ci_df(
      beta_mat   = iv_standard$coeff_iv,
      se_mat     = iv_standard$stderr_iv,
      dates      = df$date,
      coef_names = coef_names_standard
    ) %>%
      mutate(
        h1        = h1,
        kernel    = kernel_type,
        pc_type   = "standard",
        estimator = "iv"
      )
    
    ## NK PC
    endog_idx_nk <- c(2L, 4L)  # dpi_lead1 and du endogenous
    
    ols_nk <- tv_OLS_tv(
      y      = y,
      X      = X_nk,
      h1     = h1,
      kernel = kernel_type
    )
    
    iv_nk <- tv_IV_tv(
      y_vec     = y,
      COVAR     = X_nk,
      INST      = Z_nk,
      endog_idx = endog_idx_nk,
      h1        = h1,
      h2        = h1,
      kernel    = kernel_type
    )
    
    df_ols_nk_ci <- build_ci_df(
      beta_mat   = ols_nk$coeff,
      se_mat     = ols_nk$stderr,
      dates      = df$date,
      coef_names = coef_names_nk
    ) %>%
      mutate(
        h1        = h1,
        kernel    = kernel_type,
        pc_type   = "nk",
        estimator = "ols"
      )
    
    df_iv_nk_ci <- build_ci_df(
      beta_mat   = iv_nk$coeff_iv,
      se_mat     = iv_nk$stderr_iv,
      dates      = df$date,
      coef_names = coef_names_nk
    ) %>%
      mutate(
        h1        = h1,
        kernel    = kernel_type,
        pc_type   = "nk",
        estimator = "iv"
      )
    
    results_list[[res_idx]] <- df_ols_standard_ci; res_idx <- res_idx + 1L
    results_list[[res_idx]] <- df_iv_standard_ci;  res_idx <- res_idx + 1L
    results_list[[res_idx]] <- df_ols_nk_ci;       res_idx <- res_idx + 1L
    results_list[[res_idx]] <- df_iv_nk_ci;        res_idx <- res_idx + 1L
  }
}

df_all <- bind_rows(results_list)

df_export <- df_all %>%
  transmute(
    date      = format(date, "%Y-%m-%d"),
    h1        = h1,
    kernel    = kernel,
    pc_type   = pc_type,
    estimator = estimator,
    coef      = coef,
    estimate  = estimate,
    lower     = lower,
    upper     = upper
  )

############################################################
## 11. WRITE CSV TO REPOSITORY (data/)
############################################################

dir.create("data", showWarnings = FALSE)
out_file <- file.path("data", "euro_area_tvpc_full.csv")
write.csv(df_export, out_file, row.names = FALSE)

cat("Exported TV Phillips curve results to", out_file, "\n")
