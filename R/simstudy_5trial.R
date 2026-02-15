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

# standard meta-analysis and meta-regression
sim_meta <- function(data) {
  out_na <- list(meta = na_interval(), metareg = na_interval())

  tryCatch({
    data_R1 <- data[data$R == 1, , drop = FALSE]
    if (nrow(data_R1) == 0) return(out_na)

    split_trials <- split(data_R1, data_R1$S)
    data_R1_sum <- do.call(rbind, lapply(names(split_trials), function(s_name) {
      d <- split_trials[[s_name]]
      data.frame(
        S = s_name,
        Y0_A1 = sum(d$A == 1 & d$Y == 0),
        Y1_A1 = sum(d$A == 1 & d$Y == 1),
        N_A1 = sum(d$A == 1),
        Y0_A0 = sum(d$A == 0 & d$Y == 0),
        Y1_A0 = sum(d$A == 0 & d$Y == 1),
        N_A0 = sum(d$A == 0),
        X1_mean = mean(d$X1),
        X2_mean = mean(d$X2),
        X3_mean = mean(d$X3),
        stringsAsFactors = FALSE
      )
    }))

    data_R1_sum <- data_R1_sum[data_R1_sum$N_A1 > 0 & data_R1_sum$N_A0 > 0, , drop = FALSE]
    if (nrow(data_R1_sum) == 0) return(out_na)

    data_R0 <- data[data$R == 0, , drop = FALSE]
    if (nrow(data_R0) == 0) return(out_na)

    data_R0_mean <- data.frame(
      X1_mean = mean(data_R0$X1),
      X2_mean = mean(data_R0$X2),
      X3_mean = mean(data_R0$X3)
    )

    fit_meta <- try(meta::metabin(
      event.e = Y1_A1,
      n.e = N_A1,
      event.c = Y1_A0,
      n.c = N_A0,
      data = data_R1_sum,
      studlab = S,
      sm = "RD",
      method.tau = "REML"
    ), silent = TRUE)

    if (inherits(fit_meta, "try-error")) return(out_na)

    fit_reg <- try(meta::metareg(fit_meta, ~ X1_mean + X2_mean + X3_mean), silent = TRUE)

    predictions <- tryCatch({
      if (inherits(fit_reg, "try-error")) stop("metareg failed")
      metafor::predict.rma(fit_reg, newmods = as.matrix(data_R0_mean))
    }, error = function(e) {
      list(pred = NA_real_, ci.lb = NA_real_, ci.ub = NA_real_)
    })

    list(
      meta = list(
        point = safe_numeric(fit_meta$TE.random),
        lb = safe_numeric(fit_meta$lower.random),
        ub = safe_numeric(fit_meta$upper.random)
      ),
      metareg = list(
        point = safe_numeric(predictions$pred),
        lb = safe_numeric(predictions$ci.lb),
        ub = safe_numeric(predictions$ci.ub)
      )
    )
  }, error = function(e) {
    out_na
  })
}

# the proposed method (CIMA)
sim_CIMA <- function(data) {
  out_na <- list(CIMA = na_interval())

  tryCatch({
    dt_trial <- data[data$R == 1, , drop = FALSE]
    dt_target <- data[data$R == 0, , drop = FALSE]

    if (nrow(dt_trial) < 10 || nrow(dt_target) < 2) return(out_na)

    trial_labels <- paste0("S", 1:5)

    overall <- data.frame(
      trial = trial_labels,
      subgroup = rep("overall", 5),
      stratum = rep("TRUE", 5),
      stringsAsFactors = FALSE
    )
    dataTE <- expand.grid(
      trial = trial_labels,
      subgroup = c("X1", "X2"),
      stratum = c("x == 0", "x == 1"),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    dataTE <- dataTE[order(dataTE$trial, dataTE$subgroup, dataTE$stratum), , drop = FALSE]
    X3 <- expand.grid(
      trial = trial_labels,
      subgroup = "X3",
      stratum = c("x <= 0", "x > 0"),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    X3 <- X3[order(X3$trial, X3$stratum), , drop = FALSE]
    dataTE <- rbind(overall, dataTE, X3)

    AgD <- lapply(seq_len(nrow(dataTE)), function(k) {
      trial_name <- dataTE$trial[k]
      subgroup_name <- dataTE$subgroup[k]
      stratum_rule <- dataTE$stratum[k]

      trial_ind <- dt_trial$S == trial_name
      if (!any(trial_ind)) return(c(NA_real_, NA_real_, 0))

      d_trial <- dt_trial[trial_ind, , drop = FALSE]
      subgroup_values <- d_trial[[subgroup_name]]
      subgroup_ind <- if (identical(stratum_rule, "TRUE")) {
        rep(TRUE, nrow(d_trial))
      } else {
        eval(parse(text = stratum_rule), envir = list(x = subgroup_values))
      }

      if (length(subgroup_ind) != nrow(d_trial) || all(!subgroup_ind)) {
        return(c(NA_real_, NA_real_, 0))
      }

      subdat <- d_trial[subgroup_ind, , drop = FALSE]
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

    dataTE$prop <- NA_real_
    for (trial_name in trial_labels) {
      i_trial <- dataTE$trial == trial_name
      i_overall <- i_trial & dataTE$subgroup == "overall"
      n_overall <- safe_numeric(dataTE$n[i_overall])
      if (is.finite(n_overall) && n_overall > 0) {
        dataTE$prop[i_trial] <- dataTE$n[i_trial] / n_overall
      }
    }

    keep <- is.finite(dataTE$seTE) & dataTE$seTE > 0
    dataTE <- dataTE[keep, , drop = FALSE]
    if (nrow(dataTE) == 0) return(out_na)

    mu <- setNames(vector("list", length(trial_labels)), trial_labels)
    for (trial_name in trial_labels) {
      d <- dt_trial[dt_trial$S == trial_name, , drop = FALSE]
      if (nrow(d) == 0) return(out_na)
      mu_vec <- c(mean(d$X1), mean(d$X2), mean(d$X3), mean(d$X3 * d$X3))
      if (any(!is.finite(mu_vec))) return(out_na)
      mu[[trial_name]] <- as.numeric(mu_vec)
    }

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
simulateData <- function(N, eta, beta, gamma, theta) {
  def <- simstudy::defData(varname = "X1", formula = eta[1], dist = "binary")
  def <- simstudy::defData(def, varname = "X2", formula = eta[2], dist = "binary")
  covariates <- simstudy::genData(N, def)
  covariates <- simstudy::addCorData(covariates, idname = "id", mu = eta[3], sigma = 1, rho = eta[4], corstr = "cs")
  colnames(covariates) <- c("id", "X1", "X2", "X3")
  covariates <- as.data.frame(covariates)

  Pr_R <- plogis(beta[1] + beta[2] * covariates$X1 + beta[3] * covariates$X2 + beta[4] * covariates$X3)
  R <- rbinom(N, 1, Pr_R)

  data_overall <- cbind(data.frame(R = R), covariates)
  data_overall$S <- NA_character_

  idx_trial <- which(data_overall$R == 1)
  if (length(idx_trial) > 0) {
    x1 <- data_overall$X1[idx_trial]
    x2 <- data_overall$X2[idx_trial]
    x3 <- data_overall$X3[idx_trial]

    lp <- cbind(
      0,
      gamma[1, 1] + gamma[1, 2] * x1 + gamma[1, 3] * x2 + gamma[1, 4] * x3,
      gamma[2, 1] + gamma[2, 2] * x1 + gamma[2, 3] * x2 + gamma[2, 4] * x3,
      gamma[3, 1] + gamma[3, 2] * x1 + gamma[3, 3] * x2 + gamma[3, 4] * x3,
      gamma[4, 1] + gamma[4, 2] * x1 + gamma[4, 3] * x2 + gamma[4, 4] * x3
    )
    probs <- exp(lp)
    probs <- probs / rowSums(probs)

    labels <- c("S1", "S2", "S3", "S4", "S5")
    sampled_S <- apply(probs, 1, function(p) sample(labels, 1, prob = p))
    data_overall$S[idx_trial] <- sampled_S
  }

  data_overall$A <- rbinom(nrow(data_overall), 1, 0.5)

  pY1 <- plogis(theta[1, 1] + theta[1, 2] * data_overall$X1 + theta[1, 3] * data_overall$X2 + theta[1, 4] * data_overall$X3)
  pY0 <- plogis(theta[2, 1] + theta[2, 2] * data_overall$X1 + theta[2, 3] * data_overall$X2 + theta[2, 4] * data_overall$X3)
  data_overall$`Y(1)` <- rbinom(nrow(data_overall), 1, pY1)
  data_overall$`Y(0)` <- rbinom(nrow(data_overall), 1, pY0)
  data_overall$Y <- data_overall$A * data_overall$`Y(1)` + (1 - data_overall$A) * data_overall$`Y(0)`

  data_overall
}

# scenarios --------------------------------------------------------------
N <- list(
  N1 = 5000
)

eta <- list(
  eta1 = c(0.3, 0.3, 0, 0.3),
  eta2 = c(0.5, 0.5, 0, 0.5)
)

beta <- list(
  beta1 = c(log(2), log(0.5), log(0.5), log(0.5)),
  beta2 = c(log(0.8), log(2), log(2), log(2))
)

gamma <- list(
  gamma1 = matrix(
    c(log(2), log(0.5), log(2), log(0.5),
      log(2), log(0.8), log(1.25), log(0.8),
      log(2), log(0.5), log(2), log(0.5),
      log(2), log(0.8), log(1.25), log(0.8)
    ),
    nrow = 4, ncol = 4, byrow = TRUE
  ),
  gamma2 = matrix(
    c(log(2), log(2), log(0.5), log(2),
      log(2), log(1.25), log(0.8), log(1.25),
      log(2), log(2), log(0.5), log(2),
      log(2), log(1.25), log(0.8), log(2)
    ),
    nrow = 4, ncol = 4, byrow = TRUE
  )
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
  N = c("N1"),
  eta = c("eta1", "eta2"),
  beta = c("beta1", "beta2"),
  gamma = c("gamma1", "gamma2"),
  theta = c("theta1", "theta2"),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# test data generation ---------------------------------------
# setting <- 1
# 
# dt <- simulateData(
#   N = N[[scenarios$N[setting]]],
#   eta = eta[[scenarios$eta[setting]]],
#   beta = beta[[scenarios$beta[setting]]],
#   gamma = gamma[[scenarios$gamma[setting]]],
#   theta = theta[[scenarios$theta[setting]]]
# )
# 
# # check the data
# with(dt, table(R))
# with(dt[dt$R == 1, , drop = FALSE], table(S, useNA = "ifany"))
# aggregate(cbind(X1, X2, X3) ~ R + S, data = dt, FUN = mean)
# mean(dt$Y[dt$R == 0 & dt$A == 1]) - mean(dt$Y[dt$R == 0 & dt$A == 0])
# 
# # test models
# sim_meta(dt)
# system.time(cima <- sim_CIMA(dt))
# cima
# sim_IPD(dt)

# simulation study ------------------------------------------------------
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
          gamma = gamma[[scenarios$gamma[index]]],
          theta = theta[[scenarios$theta[index]]]
        )

        true_value <- mean(dt$`Y(1)`[dt$R == 0] - dt$`Y(0)`[dt$R == 0])
        if (!is.finite(true_value)) true_value <- NA_real_

        results_meta <- sim_meta(dt)
        results_CIMA <- sim_CIMA(dt)
        results_IPD <- sim_IPD(dt)

        list(
          true_value = list(true_RD = true_value),
          meta = results_meta,
          CIMA = results_CIMA,
          IPD = results_IPD
        )
      }, error = function(e) {
        list(
          true_value = list(true_RD = NA_real_),
          meta = list(meta = na_interval(), metareg = na_interval()),
          CIMA = list(CIMA = na_interval()),
          IPD = list(IPD = na_interval())
        )
      })

      iter_out
    } -> results

    true_vals <- sapply(results, function(x) safe_numeric(x$true_value$true_RD))
    true_RD_mean <- safe_mean_stat(true_vals)

    meta_point <- sapply(results, function(x) safe_numeric(x$meta$meta$point))
    metareg_point <- sapply(results, function(x) safe_numeric(x$meta$metareg$point))
    cima_point <- sapply(results, function(x) safe_numeric(x$CIMA$CIMA$point))
    ipd_point <- sapply(results, function(x) safe_numeric(x$IPD$IPD$point))

    meta_lb <- sapply(results, function(x) safe_numeric(x$meta$meta$lb))
    meta_ub <- sapply(results, function(x) safe_numeric(x$meta$meta$ub))
    metareg_lb <- sapply(results, function(x) safe_numeric(x$meta$metareg$lb))
    metareg_ub <- sapply(results, function(x) safe_numeric(x$meta$metareg$ub))
    cima_lb <- sapply(results, function(x) safe_numeric(x$CIMA$CIMA$lb))
    cima_ub <- sapply(results, function(x) safe_numeric(x$CIMA$CIMA$ub))
    ipd_lb <- sapply(results, function(x) safe_numeric(x$IPD$IPD$lb))
    ipd_ub <- sapply(results, function(x) safe_numeric(x$IPD$IPD$ub))

    bias_meta <- safe_mean_stat(meta_point - true_RD_mean)
    bias_metareg <- safe_mean_stat(metareg_point - true_RD_mean)
    bias_CIMA <- safe_mean_stat(cima_point - true_RD_mean)
    bias_IPD <- safe_mean_stat(ipd_point - true_RD_mean)

    var_meta <- safe_var_stat(meta_point)
    var_metareg <- safe_var_stat(metareg_point)
    var_CIMA <- safe_var_stat(cima_point)
    var_IPD <- safe_var_stat(ipd_point)

    coverage_meta <- safe_mean_stat(as.numeric(meta_lb <= true_RD_mean & meta_ub >= true_RD_mean))
    coverage_metareg <- safe_mean_stat(as.numeric(metareg_lb <= true_RD_mean & metareg_ub >= true_RD_mean))
    coverage_CIMA <- safe_mean_stat(as.numeric(cima_lb <= true_RD_mean & cima_ub >= true_RD_mean))
    coverage_IPD <- safe_mean_stat(as.numeric(ipd_lb <= true_RD_mean & ipd_ub >= true_RD_mean))

    MSE_meta <- safe_mean_stat((meta_point - true_RD_mean)^2)
    MSE_metareg <- safe_mean_stat((metareg_point - true_RD_mean)^2)
    MSE_CIMA <- safe_mean_stat((cima_point - true_RD_mean)^2)
    MSE_IPD <- safe_mean_stat((ipd_point - true_RD_mean)^2)

    MAE_meta <- safe_mean_stat(abs(meta_point - true_RD_mean))
    MAE_metareg <- safe_mean_stat(abs(metareg_point - true_RD_mean))
    MAE_CIMA <- safe_mean_stat(abs(cima_point - true_RD_mean))
    MAE_IPD <- safe_mean_stat(abs(ipd_point - true_RD_mean))

    performance <- c(
      bias_meta, bias_metareg, bias_CIMA, bias_IPD,
      var_meta, var_metareg, var_CIMA, var_IPD,
      coverage_meta, coverage_metareg, coverage_CIMA, coverage_IPD,
      MSE_meta, MSE_metareg, MSE_CIMA, MSE_IPD,
      MAE_meta, MAE_metareg, MAE_CIMA, MAE_IPD
    )
    names(performance) <- c(
      "bias_meta", "bias_metareg", "bias_CIMA", "bias_IPD",
      "var_meta", "var_metareg", "var_CIMA", "var_IPD",
      "coverage_meta", "coverage_metareg", "coverage_CIMA", "coverage_IPD",
      "MSE_meta", "MSE_metareg", "MSE_CIMA", "MSE_IPD",
      "MAE_meta", "MAE_metareg", "MAE_CIMA", "MAE_IPD"
    )

    performance
  } -> sim_results
)

sim_results
write.csv(sim_results, "sim_results.csv", row.names = FALSE)
