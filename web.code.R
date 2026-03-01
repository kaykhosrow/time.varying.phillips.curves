############################################################
## Time-Varying Phillips Curves - Multi-country Export
## Exports one CSV for all countries to upload to GitHub
## Output: data/tvpc_all_countries.csv
############################################################

############################################################
## 0. PACKAGES
############################################################

library(dplyr)
library(tidyr)
library(lubridate)

############################################################
## 1. SETTINGS
############################################################

H1_GRID  <- seq(0, 1, by = 0.1)
KERNELS  <- c("gaussian", "epanechnikov", "rolling")
CI_LEVEL <- 0.90
Z_CRIT   <- qnorm(0.5 + CI_LEVEL / 2)

############################################################
## 2. HELPERS
############################################################

lag_vec  <- function(x, k = 1) c(rep(NA, k), x[1:(length(x) - k)])
lead_vec <- function(x, k = 1) c(x[(k + 1):length(x)], rep(NA, k))

clean_num <- function(x) as.numeric(gsub("^#N/A$|^$", NA, trimws(x)))

############################################################
## 3. DATA BUILDER
############################################################

build_tvpc_df_from_csv <- function(file_path) {
  df_raw  <- read.csv(file_path, stringsAsFactors = FALSE)
  df_raw  <- df_raw[trimws(df_raw[[1]]) != "" & !is.na(df_raw[[1]]), ]
  
  date <- lubridate::dmy(df_raw[[1]])
  U    <- clean_num(df_raw[[3]])
  CPI  <- clean_num(df_raw[[4]])
  
  pi        <- 100 * (log(CPI) - lag_vec(log(CPI), 1))
  dpi       <- pi - lag_vec(pi, 1)
  dpi_lag1  <- lag_vec(dpi, 1)
  dpi_lead1 <- lead_vec(dpi, 1)
  dpi_l2    <- lag_vec(dpi, 2)
  dpi_l3    <- lag_vec(dpi, 3)
  dpi_l4    <- lag_vec(dpi, 4)
  du        <- U - lag_vec(U, 1)
  du_l1     <- lag_vec(du, 1)
  du_l2     <- lag_vec(du, 2)
  du_l3     <- lag_vec(du, 3)
  du_l4     <- lag_vec(du, 4)
  
  df <- data.frame(
    date = date, U = U, CPI = CPI, pi = pi, dpi = dpi,
    dpi_lag1 = dpi_lag1, dpi_lead1 = dpi_lead1,
    dpi_l2 = dpi_l2, dpi_l3 = dpi_l3, dpi_l4 = dpi_l4,
    du = du, du_l1 = du_l1, du_l2 = du_l2,
    du_l3 = du_l3, du_l4 = du_l4
  )
  stats::na.omit(df)
}

############################################################
## 4. MATRIX HELPERS
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
    xt       <- X[t, , drop = FALSE]
    out[t, ] <- vech_lower(t(xt) %*% xt)
  }
  out
}

mat_sqrt_sym <- function(A) {
  A <- as.matrix(A)
  if (any(!is.finite(A))) return(matrix(0, nrow(A), ncol(A)))
  A <- 0.5 * (A + t(A))
  ev <- eigen(A, symmetric = TRUE)
  vals <- pmax(ev$values, 0)
  if (all(vals == 0)) return(matrix(0, nrow(A), ncol(A)))
  P <- ev$vectors
  P %*% diag(sqrt(vals), length(vals)) %*% t(P)
}

############################################################
## 5. KERNEL WEIGHTS
############################################################

kernel_weights <- function(T, H, kernel = c("gaussian", "epanechnikov", "rolling")) {
  kernel    <- match.arg(kernel)
  seq_vals  <- -T:(T - 1)
  
  if (kernel == "gaussian") {
    base <- exp(-0.5 * (abs(seq_vals) / H)^2)
  } else if (kernel == "epanechnikov") {
    u    <- seq_vals / H
    base <- ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
  } else {
    base <- ifelse(abs(seq_vals) <= H, 1, 0)
  }
  
  W_mat <- matrix(NA_real_, nrow = T, ncol = T)
  ini   <- T + 1
  fin   <- 2 * T
  for (t in 1:T) {
    W_mat[t, ] <- base[ini:fin]
    ini <- ini - 1
    fin <- fin - 1
  }
  W_mat
}

############################################################
## 6. KERNEL LOOP
############################################################

kernel_loop_tv <- function(y, X, Z, H, kernel = "gaussian") {
  y <- as.matrix(y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  T <- nrow(y); q <- ncol(y); k <- ncol(X)
  W_mat <- kernel_weights(T, H, kernel)
  ret   <- matrix(NA_real_, nrow = T, ncol = q * k)
  for (t in 1:T) {
    W   <- W_mat[t, ]
    PX  <- Z * W
    num <- crossprod(PX, y)
    den <- crossprod(PX, X)
    B   <- tryCatch(solve(den, num),
                    error = function(e) matrix(NA_real_, k, q))
    ret[t, ] <- as.vector(B)
  }
  ret
}

############################################################
## 7. TV-OLS
############################################################

tv_OLS_tv <- function(y, X, h1, kernel = "gaussian") {
  y     <- as.numeric(y)
  X     <- as.matrix(X)
  T     <- length(y)
  k     <- ncol(X)
  H     <- T^h1
  y_mat <- matrix(y, ncol = 1)
  
  bethat <- kernel_loop_tv(y_mat, X, X, H, kernel)
  y_hat  <- rowSums(bethat * X)
  u      <- y - y_hat
  xx     <- hdprod_rows(X)
  u2     <- u^2
  XXu    <- xx * u2
  W_mat  <- kernel_weights(T, H, kernel)
  
  SE  <- matrix(NA_real_, nrow = T, ncol = k)
  VCV <- matrix(NA_real_, nrow = T, ncol = k * (k + 1) / 2)
  
  for (t in 1:T) {
    W  <- W_mat[t, ]
    K  <- sum(W);  K2 <- sum(W^2)
    if (!is.finite(K) || K <= 0) next
    XXuW <- as.vector(crossprod(XXu, W))
    xxW  <- as.vector(crossprod(xx,  W))
    if (any(!is.finite(XXuW)) || any(!is.finite(xxW))) next
    Szu  <- unvech_sym(XXuW / K, k)
    Szz  <- unvech_sym(xxW  / K, k)
    iSzz <- tryCatch(solve(Szz), error = function(e) NULL)
    if (is.null(iSzz)) next
    COV_t    <- (K2 / K^2) * iSzz %*% Szu %*% iSzz
    if (any(!is.finite(COV_t))) next
    VCV[t, ] <- vech_lower(COV_t)
    SE[t, ]  <- sqrt(pmax(diag(COV_t), 0))
  }
  
  list(coeff = bethat, stderr = SE)
}

############################################################
## 8. TV-IV
############################################################

build_Pt_init <- function(covar_names, inst_names) {
  q  <- length(inst_names); k <- length(covar_names)
  Pt <- matrix(0, nrow = q, ncol = k)
  rownames(Pt) <- inst_names; colnames(Pt) <- covar_names
  for (j in 1:k) {
    pos <- which(inst_names == covar_names[j])
    if (length(pos) > 0) Pt[pos[1], j] <- 1
  }
  Pt
}

tv_IV_tv <- function(y_vec, COVAR, INST, endog_idx, h1, h2 = h1, kernel = "gaussian") {
  y_vec  <- as.numeric(y_vec)
  X_cov  <- as.matrix(COVAR)
  Z_inst <- as.matrix(INST)
  T      <- length(y_vec); k <- ncol(X_cov); q <- ncol(Z_inst)
  n_endo <- length(endog_idx)
  H      <- T^h1; H2 <- T^h2
  
  covar_names <- colnames(X_cov)
  inst_names  <- colnames(Z_inst)
  if (is.null(covar_names)) { covar_names <- paste0("x", seq_len(k)); colnames(X_cov) <- covar_names }
  if (is.null(inst_names))  { inst_names  <- paste0("z", seq_len(q)); colnames(Z_inst) <- inst_names }
  Pt_init <- build_Pt_init(covar_names, inst_names)
  
  ## First stage - one regression per endogenous variable
  psihat_list <- vector("list", n_endo)
  xfithat     <- X_cov
  for (j in seq_len(n_endo)) {
    xe_j            <- X_cov[, endog_idx[j], drop = FALSE]
    psi_j           <- kernel_loop_tv(xe_j, Z_inst, Z_inst, H2, kernel)
    psihat_list[[j]] <- psi_j
    xfithat[, endog_idx[j]] <- rowSums(psi_j * Z_inst)
  }
  psihat_full <- do.call(cbind, psihat_list)
  
  ## Second stage
  y_mat    <- matrix(y_vec, ncol = 1)
  beta_iv  <- kernel_loop_tv(y_mat, X_cov, xfithat, H, kernel)
  y_hat_iv <- rowSums(beta_iv * X_cov)
  u        <- y_vec - y_hat_iv
  u2       <- u^2
  
  ## OLS for Hausman
  beta_ols <- kernel_loop_tv(y_mat, X_cov, X_cov, H, kernel)
  b_e_ols  <- beta_ols[, endog_idx, drop = FALSE]
  b_e_iv   <- beta_iv[,  endog_idx, drop = FALSE]
  
  Xe   <- X_cov[, endog_idx, drop = FALSE]
  Xfe  <- xfithat[, endog_idx, drop = FALSE]
  Vvec <- Xe - Xfe
  xx   <- hdprod_rows(Xe)
  xpxp <- hdprod_rows(Xfe)
  vv   <- hdprod_rows(Vvec)
  zz   <- hdprod_rows(Z_inst)
  zzu  <- zz * u2
  W_all <- kernel_weights(T, H, kernel)
  
  SE_iv  <- matrix(NA_real_, nrow = T, ncol = k)
  VCV_iv <- matrix(NA_real_, nrow = T, ncol = k * (k + 1) / 2)
  H_stat <- rep(NA_real_, T)
  p_val  <- rep(NA_real_, T)
  
  for (t in 1:T) {
    W  <- W_all[t, ]; K <- sum(W); K2 <- sum(W^2)
    if (!is.finite(K) || K <= 0) next
    
    xxW   <- as.vector(crossprod(xx,   W))
    xpxpW <- as.vector(crossprod(xpxp, W))
    vvW   <- as.vector(crossprod(vv,   W))
    if (any(!is.finite(xxW)) || any(!is.finite(xpxpW)) || any(!is.finite(vvW))) next
    
    Sxx   <- unvech_sym(xxW   / K, n_endo)
    Sxpxp <- unvech_sym(xpxpW / K, n_endo)
    Svv   <- unvech_sym(vvW   / K, n_endo)
    s2u   <- sum(W * u2) / K
    if (!is.finite(s2u) || s2u <= 0) next
    
    hSxx   <- mat_sqrt_sym(Sxx)
    hSxpxp <- mat_sqrt_sym(Sxpxp)
    diff_b <- b_e_ols[t, ] - b_e_iv[t, ]
    Vt     <- hSxx %*% hSxpxp %*% matrix(diff_b, ncol = 1)
    
    Svv_sym <- 0.5 * (Svv + t(Svv)) + diag(1e-10, n_endo)
    iSvv <- tryCatch(solve(Svv_sym), error = function(e) NULL)
    if (!is.null(iSvv)) {
      H_t <- ((K^2) / (K2 * s2u)) * as.numeric(t(Vt) %*% iSvv %*% Vt)
      if (is.finite(H_t) && H_t >= 0) {
        H_stat[t] <- H_t
        p_val[t]  <- stats::pchisq(H_t, df = n_endo, lower.tail = FALSE)
      }
    }
    
    zzuW <- as.vector(crossprod(zzu, W))
    zzW  <- as.vector(crossprod(zz,  W))
    if (any(!is.finite(zzuW)) || any(!is.finite(zzW))) next
    
    Szu <- unvech_sym(zzuW / K, q)
    Szz <- unvech_sym(zzW  / K, q)
    
    Pt_t      <- Pt_init
    psi_t_mat <- matrix(psihat_full[t, ], nrow = q, ncol = n_endo)
    for (j in seq_len(n_endo)) Pt_t[, endog_idx[j]] <- psi_t_mat[, j]
    
    PSzuP  <- t(Pt_t) %*% Szu %*% Pt_t
    PSzzP  <- t(Pt_t) %*% Szz %*% Pt_t
    iPSzzP <- tryCatch(solve(PSzzP), error = function(e) NULL)
    if (!is.null(iPSzzP)) {
      COV_t <- (K2 / K^2) * iPSzzP %*% PSzuP %*% iPSzzP
      if (any(is.finite(COV_t))) {
        VCV_iv[t, ] <- vech_lower(COV_t)
        SE_iv[t, ]  <- sqrt(pmax(diag(COV_t), 0))
      }
    }
  }
  
  colnames(beta_iv) <- covar_names
  colnames(SE_iv)   <- covar_names
  
  list(beta_iv = beta_iv, SE_iv = SE_iv, H_stat = H_stat, p_val = p_val)
}

############################################################
## 9. CI BUILDER
############################################################

build_ci_df <- function(beta_mat, se_mat, dates, coef_names) {
  beta <- as.matrix(beta_mat); se <- as.matrix(se_mat)
  colnames(beta) <- coef_names; colnames(se) <- coef_names
  df_b <- as.data.frame(beta); df_s <- as.data.frame(se)
  df_b$date <- dates;          df_s$date <- dates
  
  df_long <- df_b %>% tidyr::pivot_longer(-date, names_to = "coef", values_to = "estimate")
  df_se   <- df_s %>% tidyr::pivot_longer(-date, names_to = "coef", values_to = "se")
  
  df_long %>%
    dplyr::left_join(df_se, by = c("date", "coef")) %>%
    dplyr::mutate(lower = estimate - Z_CRIT * se,
                  upper = estimate + Z_CRIT * se)
}

############################################################
## 10. COUNTRY LIST
############################################################

country_specs <- list(
  list(code = "ALL", label = "Euro Area",      region = "EA", file = "input.eur.csv"),
  list(code = "AUS", label = "Austria",        region = "EA", file = "input.aus.csv"),
  list(code = "BEL", label = "Belgium",        region = "EA", file = "input.bel.csv"),
  list(code = "CYP", label = "Cyprus",         region = "EA", file = "input.cyp.csv"),
  list(code = "EST", label = "Estonia",        region = "EA", file = "input.est.csv"),
  list(code = "FIN", label = "Finland",        region = "EA", file = "input.fin.csv"),
  list(code = "BUL", label = "Bulgaria",       region = "EA", file = "input.bul.csv"),
  list(code = "CRO", label = "Croatia",        region = "EA", file = "input.cro.csv"),
  list(code = "FRA", label = "France",         region = "EA", file = "input.fra.csv"),
  list(code = "GER", label = "Germany",        region = "EA", file = "input.ger.csv"),
  list(code = "GRE", label = "Greece",         region = "EA", file = "input.gre.csv"),
  list(code = "IRE", label = "Ireland",        region = "EA", file = "input.ire.csv"),
  list(code = "ITA", label = "Italy",          region = "EA", file = "input.ita.csv"),
  list(code = "LAT", label = "Latvia",         region = "EA", file = "input.lat.csv"),
  list(code = "LIT", label = "Lithuania",      region = "EA", file = "input.lit.csv"),
  list(code = "LUX", label = "Luxembourg",     region = "EA", file = "input.lux.csv"),
  list(code = "MAL", label = "Malta",          region = "EA", file = "input.mal.csv"),
  list(code = "NET", label = "Netherlands",    region = "EA", file = "input.net.csv"),
  list(code = "POR", label = "Portugal",       region = "EA", file = "input.por.csv"),
  list(code = "SLK", label = "Slovakia",       region = "EA", file = "input.slk.csv"),
  list(code = "SLN", label = "Slovenia",       region = "EA", file = "input.sln.csv"),
  list(code = "SPA", label = "Spain",          region = "EA", file = "input.spa.csv"),
  list(code = "UK",  label = "United Kingdom", region = "UK", file = "input.uk.csv"),
  list(code = "US",  label = "United States",  region = "US", file = "input.us.csv")
)

############################################################
## 11. MAIN EXPORT LOOP
############################################################

all_results <- list()

for (spec in country_specs) {
  data_path <- file.path("Data", spec$file)
  cat("============================================================\n")
  cat("Processing:", spec$label, "(", spec$code, ")\n")
  
  if (!file.exists(data_path)) {
    cat("  FILE NOT FOUND, skipping.\n")
    next
  }
  
  df <- build_tvpc_df_from_csv(data_path)
  y  <- df$dpi
  
  ## Regressors
  X_standard <- cbind(const = 1, dpi_lag1 = df$dpi_lag1, du = df$du)
  Z_standard <- cbind(const = 1, dpi_lag1 = df$dpi_lag1,
                      du_l1 = df$du_l1, du_l2 = df$du_l2,
                      du_l3 = df$du_l3, du_l4 = df$du_l4,
                      dpi_l2 = df$dpi_l2)
  
  X_nk <- cbind(const = 1, dpi_lead1 = df$dpi_lead1,
                dpi_lag1 = df$dpi_lag1, du = df$du)
  Z_nk <- cbind(const = 1, dpi_lag1 = df$dpi_lag1,
                dpi_l2 = df$dpi_l2, dpi_l3 = df$dpi_l3, dpi_l4 = df$dpi_l4,
                du_l1 = df$du_l1, du_l2 = df$du_l2,
                du_l3 = df$du_l3, du_l4 = df$du_l4)
  
  coef_names_standard <- c("const", "dpi_lag1", "du")
  coef_names_nk       <- c("const", "dpi_lead1", "dpi_lag1", "du")
  
  country_results <- list()
  idx <- 1L
  
  for (kernel_type in KERNELS) {
    cat("  Kernel:", kernel_type, "\n")
    
    for (h1 in H1_GRID) {
      
      ## --- Standard PC ---
      ols_s <- tv_OLS_tv(y, X_standard, h1, kernel_type)
      iv_s  <- tv_IV_tv(y, X_standard, Z_standard, endog_idx = 3L,
                        h1 = h1, kernel = kernel_type)
      
      df_ols_s <- build_ci_df(ols_s$coeff,  ols_s$stderr,
                              df$date, coef_names_standard) %>%
        dplyr::mutate(estimator = "ols", pc_type = "standard",
                      h1 = h1, kernel = kernel_type,
                      hausman_p = NA_real_)
      
      ## Hausman p-values are per time point - attach as a separate coef row
      df_iv_s <- build_ci_df(iv_s$beta_iv, iv_s$SE_iv,
                             df$date, coef_names_standard) %>%
        dplyr::mutate(estimator = "iv", pc_type = "standard",
                      h1 = h1, kernel = kernel_type,
                      hausman_p = NA_real_)
      
      ## Hausman p-value rows (one row per date, coef = "hausman")
      df_hausman_s <- data.frame(
        date      = df$date,
        coef      = "hausman",
        estimate  = NA_real_,
        se        = NA_real_,
        lower     = NA_real_,
        upper     = NA_real_,
        estimator = "iv",
        pc_type   = "standard",
        h1        = h1,
        kernel    = kernel_type,
        hausman_p = iv_s$p_val
      )
      
      ## --- NK PC ---
      ols_nk <- tv_OLS_tv(y, X_nk, h1, kernel_type)
      iv_nk  <- tv_IV_tv(y, X_nk, Z_nk, endog_idx = c(2L, 4L),
                         h1 = h1, kernel = kernel_type)
      
      df_ols_nk <- build_ci_df(ols_nk$coeff,  ols_nk$stderr,
                               df$date, coef_names_nk) %>%
        dplyr::mutate(estimator = "ols", pc_type = "nk",
                      h1 = h1, kernel = kernel_type,
                      hausman_p = NA_real_)
      
      df_iv_nk <- build_ci_df(iv_nk$beta_iv, iv_nk$SE_iv,
                              df$date, coef_names_nk) %>%
        dplyr::mutate(estimator = "iv", pc_type = "nk",
                      h1 = h1, kernel = kernel_type,
                      hausman_p = NA_real_)
      
      df_hausman_nk <- data.frame(
        date      = df$date,
        coef      = "hausman",
        estimate  = NA_real_,
        se        = NA_real_,
        lower     = NA_real_,
        upper     = NA_real_,
        estimator = "iv",
        pc_type   = "nk",
        h1        = h1,
        kernel    = kernel_type,
        hausman_p = iv_nk$p_val
      )
      
      country_results[[idx]]   <- df_ols_s;      idx <- idx + 1L
      country_results[[idx]]   <- df_iv_s;       idx <- idx + 1L
      country_results[[idx]]   <- df_hausman_s;  idx <- idx + 1L
      country_results[[idx]]   <- df_ols_nk;     idx <- idx + 1L
      country_results[[idx]]   <- df_iv_nk;      idx <- idx + 1L
      country_results[[idx]]   <- df_hausman_nk; idx <- idx + 1L
    }
  }
  
  df_country <- dplyr::bind_rows(country_results) %>%
    dplyr::mutate(country_code  = spec$code,
                  country_label = spec$label,
                  region        = spec$region)
  
  ## Export one CSV per country
  df_export <- df_country %>%
    dplyr::transmute(
      date      = format(date, "%Y-%m-%d"),
      pc_type   = pc_type,
      estimator = estimator,
      kernel    = kernel,
      h1        = h1,
      coef      = coef,
      estimate  = estimate,
      lower     = lower,
      upper     = upper,
      hausman_p = hausman_p
    )
  
  dir.create("data", showWarnings = FALSE)
  out_file <- file.path("data", paste0("tvpc_", tolower(spec$code), ".csv"))
  write.csv(df_export, out_file, row.names = FALSE)
  cat("  Exported:", out_file, "(", nrow(df_export), "rows )\n")
}

cat("\n============================================================\n")
cat("All countries exported to data/tvpc_<code>.csv\n")
cat("============================================================\n")