library(foreach)

# functions ---------------------------------------------------------------
na_interval <- function() {
  list(point = NA_real_, lb = NA_real_, ub = NA_real_)
}

safe_numeric <- function(x) {
  if (length(x) == 0) return(NA_real_)
  x <- as.numeric(x[1])
  if (!is.finite(x)) return(NA_real_)
  x
}

safe_mean_stat <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_var_stat <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  var(x)
}

clean_interval <- function(point, se) {
  point <- safe_numeric(point)
  se <- safe_numeric(se)
  if (is.na(point) || is.na(se) || se < 0) return(na_interval())
  list(point = point, lb = point - 1.96 * se, ub = point + 1.96 * se)
}

# the proposed method (CIMA)
sim_CIMA <- function(data) {
  out_na <- list(CIMA = na_interval())

  tryCatch({
    dt_trial <- data[data$R == 1, , drop = FALSE]
    dt_target <- data[data$R == 0, , drop = FALSE]

    if (nrow(dt_trial) < 2 || nrow(dt_target) < 2) return(out_na)

    overall <- data.frame(
      trial = "S1",
      subgroup = "overall",
      stratum = "TRUE",
      stringsAsFactors = FALSE
    )
    dataTE <- expand.grid(
      trial = "S1",
      subgroup = c("X1", "X2"),
      stratum = c("x == 0", "x == 1"),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    dataTE <- dataTE[order(dataTE$trial, dataTE$subgroup, dataTE$stratum), , drop = FALSE]
    X3 <- data.frame(
      trial = c("S1", "S1"),
      subgroup = c("X3", "X3"),
      stratum = c("x <= 0", "x > 0"),
      stringsAsFactors = FALSE
    )
    X3 <- X3[order(X3$trial, X3$stratum), , drop = FALSE]
    dataTE <- rbind(overall, dataTE, X3)

    AgD <- lapply(seq_len(nrow(dataTE)), function(k) {
      subgroup_name <- dataTE$subgroup[k]
      stratum_rule <- dataTE$stratum[k]

      subgroup_values <- dt_trial[[subgroup_name]]
      subgroup_ind <- if (identical(stratum_rule, "TRUE")) {
        rep(TRUE, nrow(dt_trial))
      } else {
        eval(parse(text = stratum_rule), envir = list(x = subgroup_values))
      }

      if (length(subgroup_ind) != nrow(dt_trial) || all(!subgroup_ind)) {
        return(c(NA_real_, NA_real_, 0))
      }

      subdat <- dt_trial[subgroup_ind, , drop = FALSE]
      y_t <- subdat$Y[subdat$A == 1]
      y_c <- subdat$Y[subdat$A == 0]
      n_t <- length(y_t)
      n_c <- length(y_c)

      if (n_t == 0 || n_c == 0) return(c(NA_real_, NA_real_, n_t + n_c))

      mu_hat <- mean(y_t) - mean(y_c)
      event_t <- sum(y_t)
      event_c <- sum(y_c)
      sigma_hat <- sqrt(event_t * (n_t - event_t) / n_t^3 + event_c * (n_c - event_c) / n_c^3)

      c(mu_hat, sigma_hat, n_t + n_c)
    })

    AgD_mat <- do.call(rbind, AgD)
    dataTE$TE <- AgD_mat[, 1]
    dataTE$seTE <- AgD_mat[, 2]
    dataTE$n <- AgD_mat[, 3]

    n_overall <- dataTE$n[dataTE$trial == "S1" & dataTE$subgroup == "overall"]
    n_overall <- safe_numeric(n_overall)
    dataTE$prop <- if (is.na(n_overall) || n_overall <= 0) NA_real_ else dataTE$n / n_overall

    keep <- is.finite(dataTE$seTE) & dataTE$seTE > 0
    dataTE <- dataTE[keep, , drop = FALSE]
    if (nrow(dataTE) == 0) return(out_na)

    mu_vec <- c(
      mean(dt_trial$X1),
      mean(dt_trial$X2),
      mean(dt_trial$X3),
      mean(dt_trial$X3 * dt_trial$X3)
    )
    if (any(!is.finite(mu_vec))) return(out_na)
    mu <- list("S1" = as.numeric(mu_vec))

    n_target <- nrow(dt_target)
    split_n <- floor(n_target / 2)
    if (split_n < 1 || split_n >= n_target) return(out_na)

    idx <- sample.int(n_target, split_n)
    Q <- dt_target[idx, c("X1", "X2", "X3"), drop = FALSE]
    target <- dt_target[-idx, c("X1", "X2", "X3"), drop = FALSE]

    fit <- try(CIMAgD::CIMA(
      dataTE,
      mu,
      Q,
      target,
      "~ X1 + X2 + X3",
      "~ X1 + X2 + X3 + I(X3^2)"
    ), silent = TRUE)

    if (inherits(fit, "try-error")) return(out_na)

    list(CIMA = clean_interval(fit$ATE$ATE, fit$ATE$se))
  }, error = function(e) {
    out_na
  })
}

# IPD g-formula
sim_IPD <- function(data) {
  out_na <- list(IPD = na_interval())

  tryCatch({
    dt_trial <- data[data$R == 1, , drop = FALSE]
    dt_target <- data[data$R == 0, , drop = FALSE]

    if (nrow(dt_trial) < 2 || nrow(dt_target) < 2) return(out_na)
    if (length(unique(dt_trial$A)) < 2 || length(unique(dt_trial$Y)) < 2) return(out_na)

    fit <- try(glm(Y ~ A * (X1 + X2 + X3), data = dt_trial, family = binomial()), silent = TRUE)
    if (inherits(fit, "try-error")) return(out_na)

    params <- fit$coefficients
    if (any(!is.finite(params))) return(out_na)

    targetX <- dt_target[, c("X1", "X2", "X3"), drop = FALSE]
    g_out <- gformula(params, targetX)
    TE <- g_out$ate

    vcovfit <- try(sandwich::vcovHC(fit), silent = TRUE)
    if (inherits(vcovfit, "try-error")) return(out_na)

    deriv_TE <- try(numDeriv::grad(function(p, X) gformula(p, X)$ate, x = params, X = targetX), silent = TRUE)
    if (inherits(deriv_TE, "try-error")) return(out_na)

    res <- g_out$res
    var1 <- as.numeric(crossprod(res) / length(res)^2)
    var2 <- as.numeric(t(deriv_TE) %*% vcovfit %*% deriv_TE)
    se <- sqrt(max(0, var1 + var2))

    list(IPD = clean_interval(TE, se))
  }, error = function(e) {
    out_na
  })
}

# g-formula function
gformula <- function(params, X) {
  pred_A1 <- cbind(A = 1, X)
  designMat_A1 <- model.matrix(~ A * (X1 + X2 + X3), pred_A1)
  y_hat_A1 <- logitnorm::invlogit(designMat_A1 %*% params)

  pred_A0 <- cbind(A = 0, X)
  designMat_A0 <- model.matrix(~ A * (X1 + X2 + X3), pred_A0)
  y_hat_A0 <- logitnorm::invlogit(designMat_A0 %*% params)

  ate <- mean(y_hat_A1) - mean(y_hat_A0)
  res <- y_hat_A1 - y_hat_A0 - ate

  list(ate = ate, res = res)
}

# simulation data (mixed)
simulateData <- function(N, eta, beta, theta) {
  def <- simstudy::defData(varname = "X1", formula = eta[1], dist = "binary")
  def <- simstudy::defData(def, varname = "X2", formula = eta[2], dist = "binary")
  covariates <- simstudy::genData(N, def)
  covariates <- simstudy::addCorData(covariates, idname = "id", mu = eta[3], sigma = 1, rho = eta[4], corstr = "cs")
  colnames(covariates) <- c("id", "X1", "X2", "X3")
  covariates <- as.data.frame(covariates)

  Pr_R <- plogis(beta[1] + beta[2] * covariates$X1 + beta[3] * covariates$X2 + beta[4] * covariates$X3)
  R <- rbinom(N, 1, Pr_R)

  data_overall <- cbind(data.frame(R = R), covariates)
  data_overall$A <- rbinom(nrow(data_overall), 1, 0.5)

  pY1 <- plogis(theta[1, 1] + theta[1, 2] * data_overall$X1 + theta[1, 3] * data_overall$X2 + theta[1, 4] * data_overall$X3)
  pY0 <- plogis(theta[2, 1] + theta[2, 2] * data_overall$X1 + theta[2, 3] * data_overall$X2 + theta[2, 4] * data_overall$X3)
  data_overall$`Y(1)` <- rbinom(nrow(data_overall), 1, pY1)
  data_overall$`Y(0)` <- rbinom(nrow(data_overall), 1, pY0)
  data_overall$Y <- data_overall$A * data_overall$`Y(1)` + (1 - data_overall$A) * data_overall$`Y(0)`

  data_overall
}

# scenarios ---------------------------------------------------------------
N <- list(
  N1 = 1000,
  N2 = 2000
)

eta <- list(
  eta1 = c(0.3, 0.3, 0, 0.3),
  eta2 = c(0.5, 0.5, 0, 0.5)
)

beta <- list(
  beta1 = c(log(0.8), log(0.5), log(0.5), log(0.5)),
  beta2 = c(log(0.4), log(2), log(2), log(2))
)

theta <- list(
  theta1 = matrix(
    c(log(0.5), log(2), log(0.5), log(1.25),
      log(0.5), log(0.5), log(2), log(0.8)
    ),
    nrow = 2, ncol = 4, byrow = TRUE
  ),
  theta2 = matrix(
    c(log(0.5), log(1.25), log(0.8), log(1.1),
      log(0.5), log(0.8), log(1.25), log(0.9)
    ),
    nrow = 2, ncol = 4, byrow = TRUE
  )
)

scenarios <- expand.grid(
  N = c("N1", "N2"),
  eta = c("eta1", "eta2"),
  beta = c("beta1", "beta2"),
  theta = c("theta1", "theta2"),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# test data generation --------------------------------------------------
# setting <- 1
# 
# dt <- simulateData(
#   N = N[[scenarios$N[setting]]],
#   eta = eta[[scenarios$eta[setting]]],
#   beta = beta[[scenarios$beta[setting]]],
#   theta = theta[[scenarios$theta[setting]]]
# )
# 
# # check the data
# with(dt, table(R))
# aggregate(cbind(X1, X2, X3) ~ R, data = dt, FUN = mean)
# aggregate(Y ~ A + R, data = dt, FUN = mean)
# mean(dt$Y[dt$R == 0 & dt$A == 1]) - mean(dt$Y[dt$R == 0 & dt$A == 0])
# 
# # test models
# system.time(cima <- sim_CIMA(dt))
# cima
# sim_IPD(dt)

# simulation study -------------------------------------------------------
# parallel
doParallel::registerDoParallel(4)

system.time(
  foreach(index = 1:16, .combine = "rbind") %do% {
    foreach(i = 1:1000) %dopar% {
      cat(sprintf("Setting %d, Simulation %d\n", index, i))

      iter_out <- tryCatch({
        dt <- simulateData(
          N = N[[scenarios$N[index]]],
          eta = eta[[scenarios$eta[index]]],
          beta = beta[[scenarios$beta[index]]],
          theta = theta[[scenarios$theta[index]]]
        )

        true_value <- mean(dt$`Y(1)`[dt$R == 0] - dt$`Y(0)`[dt$R == 0])
        if (!is.finite(true_value)) true_value <- NA_real_

        results_CIMA <- sim_CIMA(dt)
        results_IPD <- sim_IPD(dt)

        list(
          true_value = list(true_RD = true_value),
          CIMA = results_CIMA,
          IPD = results_IPD
        )
      }, error = function(e) {
        list(
          true_value = list(true_RD = NA_real_),
          CIMA = list(CIMA = na_interval()),
          IPD = list(IPD = na_interval())
        )
      })

      iter_out
    } -> results

    true_vals <- sapply(results, function(x) safe_numeric(x$true_value$true_RD))
    true_RD_mean <- safe_mean_stat(true_vals)

    cima_point <- sapply(results, function(x) safe_numeric(x$CIMA$CIMA$point))
    ipd_point <- sapply(results, function(x) safe_numeric(x$IPD$IPD$point))
    cima_lb <- sapply(results, function(x) safe_numeric(x$CIMA$CIMA$lb))
    cima_ub <- sapply(results, function(x) safe_numeric(x$CIMA$CIMA$ub))
    ipd_lb <- sapply(results, function(x) safe_numeric(x$IPD$IPD$lb))
    ipd_ub <- sapply(results, function(x) safe_numeric(x$IPD$IPD$ub))

    bias_CIMA <- safe_mean_stat(cima_point - true_RD_mean)
    bias_IPD <- safe_mean_stat(ipd_point - true_RD_mean)

    var_CIMA <- safe_var_stat(cima_point)
    var_IPD <- safe_var_stat(ipd_point)

    coverage_CIMA <- safe_mean_stat(as.numeric(cima_lb <= true_RD_mean & cima_ub >= true_RD_mean))
    coverage_IPD <- safe_mean_stat(as.numeric(ipd_lb <= true_RD_mean & ipd_ub >= true_RD_mean))

    MSE_CIMA <- safe_mean_stat((cima_point - true_RD_mean)^2)
    MSE_IPD <- safe_mean_stat((ipd_point - true_RD_mean)^2)

    MAE_CIMA <- safe_mean_stat(abs(cima_point - true_RD_mean))
    MAE_IPD <- safe_mean_stat(abs(ipd_point - true_RD_mean))

    performance <- c(
      bias_CIMA, bias_IPD,
      var_CIMA, var_IPD,
      coverage_CIMA, coverage_IPD,
      MSE_CIMA, MSE_IPD,
      MAE_CIMA, MAE_IPD
    )
    names(performance) <- c(
      "bias_CIMA", "bias_IPD",
      "var_CIMA", "var_IPD",
      "coverage_CIMA", "coverage_IPD",
      "MSE_CIMA", "MSE_IPD",
      "MAE_CIMA", "MAE_IPD"
    )

    performance
  } -> sim_results
)

sim_results
write.csv(sim_results, "sim_results.csv", row.names = FALSE)
