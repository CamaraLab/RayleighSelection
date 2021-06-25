# This file contains helper functions for rayleight_selection. They are not exported.

gpd.approx <- function(samples, observed, n.samp.quar = 250,
                       nextremes = c(seq(50, 250, 25), seq(300, 500, 50), seq(600, 1000, 100)),
                       alpha = 0.15){
  # Computes gpd approximations of p-values and quartiles of these approximations
  #
  # Arguments
  #  samples: Samples as a matrix, the i-th row of the matrix corresponds to i-th observed value.
  #  observed: Observed values as a numeric.
  #  n.samp.quar: Number of samples taken when computing quartiles
  #  nextremes: Vector of integers with the candidate number of extremes for fitting a GDP.
  #  alpha: level of FDR control for choosing the number of extremes.
  #
  # Return
  #  data frame with the natural log of the p-values and the first ant third quartiles of the approximation

  out <- data.frame(log.p = rep(NA, length(observed)),
                    log.q1 = rep(NA, length(observed)),
                    log.q3 = rep(NA, length(observed)))

  # Only consider values in the lowers quantile for fitting a gpd
  num.ext <- sort(nextremes[nextremes*4 < ncol(samples)])
  if(length(num.ext) == 0) return(out)

  for(i in 1:length(observed)){
    # Transform samples and scores
    samp.quant <- quantile(samples[i,], probs = c(0.05, 0.25), names = FALSE)
    samp.loc <- samp.quant[2]
    samp.scale <- samp.quant[2] - samp.quant[1]
    samp <- ( 1 - (samples[i, samples[i,] <= samp.loc] - samp.loc)/samp.scale )^pow
    score <- ( 1 - (observed[i] - samp.loc)/samp.scale)^pow

    # Using ForwardStop to select the number of extreme samples
    gpd.test <- tryCatch(
      eva::gpdSeqTests(samp, nextremes = num.ext, method = "ad"),
      error = function(e) NULL)

    if(is.null(gpd.test)) next

    if(all(gpd.test$ForwardStop > alpha)){
      n.selected <- num.ext[length(num.ext)]
    }else if(any(gpd.test$ForwardStop > alpha)){
      if(gpd.test$ForwardStop[1] <= alpha) next
      n.selected <- num.ext[which(gpd.test$ForwardStop <= alpha)[1] - 1]
    }else next

    # Fitting gpd and computing (log of) p-value
    fit <- eva::gpdFit(samp, nextremes = n.selected)
    loc <- unname(min(fit$exceedances))
    scale <- unname(fit$par.ests[1])
    shape <- unname(fit$par.ests[2])
    if(scale <= 0 | ( (score - loc) / scale )*shape <= -1 ) next
    if(shape == 0){
      log.p <- (score - loc) / scale
    }else{
      log.p <- (-1/shape) * log1p( ( (score - loc) / scale)*shape )
    }
    if(is.infinite(log.p) | is.na(log.p)) next

    # Sampling parameters to compute quartiles of p-value
    coef <- MASS::mvrnorm(n = n.samp.quar, mu = fit$par.ests, Sigma = fit$varcov)
    degenerate <- (coef[,1] <= 0) | ( ( (score - loc)/coef[,1] )*coef[,2] <= -1)
    coef <- coef[!degenerate,]
    coef.1 <- coef[coef[,2] != 0,]
    coef.2 <- coef[coef[,2] == 0,]
    log.p.samples <- (-1/coef.1[,2]) * log1p( ( (score - loc) / coef.1[,1] )*coef.1[,2] )
    log.p.samples <- append(log.p.samples, (score - loc) / coef.2[,1])
    log.p.samples <- append(log.p.samples, rep(-Inf, sum(degenerate)))
    q <- quantile(log.p.samples, probs = c(0.25, 0.75), names = FALSE)
    out$log.p[i] <- log.p
    out$log.q1[i] <- q[1]
    out$log.q3[i] <- q[2]
  }
  return(out)
}

optim.p <- function(observed, func, scorer, dim, use.gpd, min_perms, max_perms, n.cores = 1,
                    nextremes = c(seq(50, 250, 25), seq(300, 500, 50), seq(600, 1000, 100)),
                    alpha = 0.15){
  # Optimizes the calculation of p-values by permutation
  #
  # Arguments
  #  observed: Observed values
  #  func: Functions considered as rows of a matrix
  #  scorer: LaplacianScorer object
  #  dim: Dimension of forms (0 or 1)
  #  use.gpd: Boolean indicating if a generalized pareto is to be used
  #  n.cores: Number of cores to be used in sampling
  #  Other args: Parameters for gpd.approx
  #
  # Return
  #  data frame with p-values and number of samples where convergence was achieved

  n_funcs <- nrow(func)
  idx <- 1:n_funcs # keeps track of the index of non-convergent functions
  p <- rep(NA, n_funcs)
  n.conv <- rep(NA, n_funcs)

  # if the score is positive infinity the p-value is one
  p[is.infinite(observed) & observed > 0] <- 1
  n.conv[is.infinite(observed) & observed > 0] <- 0
  idx<- idx[!(is.infinite(observed) & observed > 0)]

  n_samps <- min_perms
  samp <- scorer$sample_scores(func[idx, ,drop = F], n_samps, dim, n.cores)
  # setting up matrix to keep track of gpd p-values
  if(use.gpd) log.p.gpd <- matrix(rep(NA, 3*length(idx)), nrow = length(idx))
  while(TRUE){
    # if there are enough (at least 10) small the p-value converges
    num.le <- apply(samp <= observed[idx], 1, sum)
    conv <- num.le >= 10
    p[idx[conv]] <- num.le[conv]/n_samps
    n.conv[idx[conv]] <- n_samps
    idx <- idx[!conv]
    samp <- samp[!conv, , drop = F]

    if(use.gpd & any(!conv)){
      # try gpd approximation
      log.p.gpd <- log.p.gpd[!conv, , drop = F]

      gpd.out <- gpd.approx(samp, observed[idx], nextremes = nextremes, alpha = alpha)

      #condition 0: all values are finite
      cond0 <- apply(log.p.gpd, 1, function(x) all(is.finite(x))) &
        is.finite(gpd.out$log.q1) & is.finite(gpd.out$log.q3) &  is.finite(gpd.out$log.p)
      #condition 1: relative variation of log p values is small
      cond1 <- apply(abs(log.p.gpd/ gpd.out$log.p - 1) < 0.1, 1, all)
      #condition 3: quartiles are close to predicted values
      cond2 <- gpd.out$log.q1/gpd.out$log.p < 1.1 & gpd.out$log.q3/gpd.out$log.p > 0.9

      conv <- cond0 & cond1 & cond2
      conv[is.na(conv)] <- FALSE
      p[idx[conv]] <- exp(gpd.out$log.p[conv])
      n.conv[idx[conv]] <- n_samps
      idx <- idx[!conv]
      samp <- samp[!conv, , drop = F]
      log.p.gpd <- cbind(log.p.gpd[!conv, , drop = F], gpd.out$log.p[!conv])
      if(2*n_samps <= max_perms) log.p.gpd <- log.p.gpd[, 2:ncol(log.p.gpd), drop = F]
    }
    if(n_samps < max_perms){
      # get new samples
      num.new.samps = min(n_samps, max_perms - n_samps)
      samp <- cbind(samp, scorer$sample_scores(func[idx, , drop = F], num.new.samps, dim, n.cores))
      n_samps <- n_samps + num.new.samps
    }else{
      break
    }
  }
  p[idx] <- num.le[!conv]/n_samps
  return(data.frame(p = p, n.conv = n.conv))
}
