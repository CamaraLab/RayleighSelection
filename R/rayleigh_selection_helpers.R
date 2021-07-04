# This file contains helper functions for rayleight_selection. They are not exported.

gpd.approx <- function(samples, observed, n.samp.quar = 250, pow = 1,
                       nextremes = c(seq(50, 250, 25), seq(300, 500, 50), seq(600, 1000, 100)),
                       alpha = 0.15){
  # Computes gpd approximations of p-values and quartiles of these approximations
  #
  # Arguments
  #  samples: samples of the null distribution
  #  observed: observed values as a numeric
  #  n.samp.quar: number of samples taken when computing quartiles
  #  nextremes: vector of integers with the candidate number of extremes for fitting a GDP
  #  alpha: level of FDR control for choosing the number of extremes
  #
  # Value
  #  data frame with the natural log of the p-values and the first ant third quartiles of the approximation
  samples <- as.numeric(samples)

  out <- data.frame(log.p = NA, log.q1 = NA, log.q3 = NA)

  # Only consider values in the lowers quantile for fitting a gpd
  num.ext <- sort(nextremes[nextremes*4 < length(samples)])
  if(length(num.ext) == 0) return(out)

  # Transform samples and scores
  samp.quant <- quantile(samples, probs = c(0.05, 0.25), names = FALSE)
  samp.loc <- samp.quant[2]
  samp.scale <- samp.quant[2] - samp.quant[1]
  samp <- ( 1 - (samples[(samples <= samp.loc) & is.finite(samples)] - samp.loc)/samp.scale )^pow
  score <- ( 1 - (observed - samp.loc)/samp.scale)^pow

  # Using ForwardStop to select the number of extreme samples
  gpd.test <- tryCatch(
    eva::gpdSeqTests(samp, nextremes = num.ext, method = "ad"),
    error = function(e) NULL)

  if(is.null(gpd.test)) return(out)

  if(all(gpd.test$ForwardStop > alpha)){
    n.selected <- num.ext[length(num.ext)] #choose largest number is num.ext
  }else if(any(gpd.test$ForwardStop > alpha)){
    if(gpd.test$ForwardStop[1] <= alpha) return(out) #tail is not GPD
    n.selected <- num.ext[which(gpd.test$ForwardStop <= alpha)[1] - 1] #select largest GDP tail
  }else return(out)

  # Fitting gpd and computing (log of) p-value
  fit <- eva::gpdFit(samp, nextremes = n.selected)
  loc <- unname(min(fit$exceedances))
  scale <- unname(fit$par.ests[1])
  shape <- unname(fit$par.ests[2])
  if(scale <= 0 || ( (score - loc) / scale )*shape <= -1 ) return(out)
  if(shape == 0){
    log.p <- (score - loc) / scale
  }else{
    log.p <- (-1/shape) * log1p( ( (score - loc) / scale)*shape )
  }
  if(is.infinite(log.p) || is.na(log.p)) return(out)

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
  out$log.p <- log.p
  out$log.q1 <- q[1]
  out$log.q3 <- q[2]
  return(out)
}


add.samples <- function(scorer, f, n.perm, dim, n.cores = 1, cov = NULL,
                        n.complexes = 1, old.samps = NULL){
  # Adds new samples to old samples
  # Arguments
  #   scorer: Scorer object
  #   f: functions as rows of a matrix
  #   n.perm: number of new permutations to be added
  #   dim: dimension of forms (0 or 1)
  #   n.cores: number of cores used in sampling
  #   cov: covariates as rows of a matrix
  #   n.complexes: number of complexes considered together
  #   old samples: list with old samples
  # Value
  #   List with samples. Each element of the list corresponds to a
  #   simplicial complex, if there are no covariates, the list has a matrix
  #   of size (nrow(f), n.perm + #old samples), where each row corresponds to a
  #   function.
  #      If there are covariates each element of the list is a list with
  #   two entries:
  #     func_scores: matrix with scores of each function in the same form as above
  #     cov_scores: 3-d array with scores of covariates. Position (i, j, k) has the
  #                 score of the j-th covariate corresponding to the k-th sample of
  #                 the i-th function
  if(is.null(cov)){
    new.samps <- scorer$sample_scores(f, n.perm, dim, n.cores)
    if(n.complexes == 1){
      new.samps <- list(new.samps)
    } else {
      new.samps <- asplit(new.samps, 3)
    }
    if(is.null(old.samps)){
      return(new.samps)
    }else{
      #joining old and new samples
      for(i in 1:n.complexes){
        new.samps[[i]] <- cbind(old.samps[[i]], new.samps[[i]])
      }
      return(new.samps)
    }
  }else{
    #there are covariates
    new.samps <- scorer$sample_with_covariate(f, cov, n.perm, dim, n.cores)
    if(n.complexes == 1) new.samps <- list(new.samps)
    if(is.null(old.samps)){
      return(new.samps)
    }else{
      #joining old and new samples
      for(i in 1:n.complexes){
        new.samps[[i]][["func_scores"]] <- cbind(old.samps[[i]][["func_scores"]],
                                                 new.samps[[i]][["func_scores"]])
        new.samps[[i]][["cov_scores"]] <- abind::abind(old.samps[[i]][["cov_scores"]],
                                                       new.samps[[i]][["cov_scores"]],
                                                       along = 3)
      }
    }
    return(new.samps)
  }
}

get.null.samps <- function(f.scores, cov.scores, samps){
  # Get samples from the null distribution and observed value
  # Arguments
  #   f.scores: observed scores of function
  #   cov.scores: observed scores of covariates (or NULL)
  #   samps: samples as returned by add.samples
  # Value
  #   List with two entries:
  #      observed: observed values as a numeric vector (f.scores if cov.scores = NULL,
  #                 residue of linear model relative to observed values else)
  #      sampled: samples from the null distribution as a matrix, each row relative
  #               to a function (scores of shuffled functions if cov.scores = NULL,
  #               residues of linear model else)
  f.scores <- as.numeric(f.scores)
  if(is.null(cov.scores)){
    return(list(observed = f.scores, sampled = samps))
  }else{
    observed <- rep(NA, length(f.scores))
    sampled <- matrix(nrow = length(f.scores), ncol = ncol(samps[["func_scores"]]))
    for(i in seq_along(f.scores)){
      rr <- regresion.residues(f.scores[i], cov.scores, samps[["func_scores"]][i,],
                               samps[["cov_scores"]][i, ,])
      observed[i] <- rr[["observed"]]
      sampled[i,] <- rr[["sampled"]]
    }
    return(list(observed = observed, sampled = sampled))
  }
}


regresion.residues <- function(f.score, cov.scores, f.samps, cov.samps){
  # Accesses the significance of score taking covariates into consideration
  # Arguments
  #   f.score: observed laplacian score of function (double)
  #   cov.scores: observed laplacian score of covariates
  #   f.samps: samples of laplacian scores of function as a vector
  #   cov.samps: samples of laplacian score of covariates as a vector or as a matrix
  #      where each row has the samples of a covariate.
  # Value
  #   list with residue corresponding to obsed value and residues of samples.

  if(is.infinite(f.score)) return(1)

  if(is(cov.samps, "numeric")) cov.samps <- t(as.matrix(cov.samps))
  if(is(cov.scores, "numeric")
     || (is(cov.scores, "array") && length(dim(cov.scores)) == 1)){
    cov.scores <- as.matrix(cov.scores)
  }

  cov.samps <- cov.samps[is.finite(cov.scores), , drop = F]
  cov.scores <- cov.scores[is.finite(cov.scores), ,drop = F]
  if(nrow(cov.samps) == 0) return(list(observed = f.score, sampled = f.samps))

  samples <- data.frame(
    fs = f.samps,
    cs = t(cov.samps)
  )
  observed <- data.frame(
    fs = f.score,
    cs = t(cov.scores)
  )

  #marking rows with infinite values
  finite.rows <- apply(samples, 1, function(x) all(is.finite(x)))

  if(all(!finite.rows)){
    warning("There are no samples after eleminating infite values")
    return(list(observed = f.score, sampled = f.samps))
  }
  if(any(!finite.rows)){
    warning("There are Inf values in the samples")
  }

  #test if any of the covariates equals the function
  func.is.cov <- apply(
    samples[finite.rows, 2:ncol(samples), drop = F], 2,
    function(x) all(abs(x - samples[finite.rows,1]) <= 1e-9*(abs(x) + abs(samples[finite.rows,1]))/2)
  )
  if(any(func.is.cov)) return(list(observed = Inf, sampled = f.samps))

  linear.model <- lm(fs ~ . , data = samples[finite.rows, , drop = F])
  obs.residual <- f.score - predict(linear.model, observed)
  sampled.residuals <- f.samps - predict(linear.model, samples)
  sampled.residuals[is.na(sampled.residuals)] <- -Inf #Na's will be considered as -Inf
  return(list(observed = unname(obs.residual), sampled = unname(sampled.residuals)))
}

combine.p.values <- function(p.vals, samples, method = "KM"){
  # Combines p-values associated to samples using Brown's method
  # Arguments
  #   p.vals: p-values as a numeric
  #   samples: samples as a matrix, each column corresponding to a complex and each row
  #            to a joint sample
  #   method: Method used to combine p-values, can be "KM" for the Kost-McDermott method
  #         or "EBM" for the empirical Brown's method (Gibbs et. al. 16)
  # Value
  #    combined p-value (double)

  if(any(p.vals == 0)) return(0)

  samples <- samples[ apply(is.finite(samples), 1, all),  ]
  if(nrow(samples) == 0) return(1)

  if(method == "EBM"){
    w <- apply(
      samples, 2,
      function(x) -2*log(ecdf(x)(x))
    )
    cov.log.p <- cov(w)
  }
  if(method == "KM"){
    rho <- cor(samples)
    cov.log.p <- 3.263*rho + 0.71*rho^2 + 0.027*rho^3
  }
  expectation <- 2*length(p.vals)
  variance <- 2*expectation + 2*sum(cov.log.p[lower.tri(cov.log.p)])
  f.val <- (expectation^2) / variance
  c.val <- variance / (2*expectation)
  p.combined <- pchisq(-2*sum(log(p.vals)) / c.val, df = 2*f.val, lower.tail = F)
  return(p.combined)
}


compute.p <- function(f, scorer, dim, min.perm, use.gpd = FALSE, cov = NULL,
                      max.perm = min.perm, n.cores = 1, combination.method = "KM",
                      pow = NA, nextremes = NULL, alpha = NA){
  # Compute scores and p-values
  # Arguments
  #   f: functions is a numeric vector, 1d array or rows of a matrix
  #   scorer: LaplacianScorer or ScorerEnsemble object
  #   dim: dimension of forms (0 or 1)
  #   min.perm: minimum number of permutations to be performed by function
  #   use.gpd: boolean determining if GPD approximation will be used
  #   cov: covariates as a numeric vector or as rows of a matrix
  #   max.perm: maximum number of permutation to be performed by function
  #   n.cores: number of cores used in sampling
  #   combination.method: method used to combine p-values when scorer is a
  #                       ScorerEnsemble, must be "KM" or "EBM"
  #   other arguments: parameters for GPD approximation
  #
  # Value
  #   Data frame with scores, p-values and when p-values converged and combined
  #   p-values, in this order

  if(is(f, "numeric") || (is(f, "array") && length(dim(f)) == 1)) f <- t(as.matrix(f))
  if(!is.null(cov) && is(cov, "numeric")) cov <- t(as.matrix(cov))

  out <- list(R = scorer$score(f, dim))
  f.scores <- asplit(out$R,2)

  n.funcs <- nrow(f)
  n.complexes <- length(f.scores)
  idx.nc <- 1:n.funcs # keeps track of the index of non-convergent functions
  p <- matrix(rep(NA, n.funcs*n.complexes), nrow = n.funcs)
  combined.p <- rep(NA, n.funcs)
  n.conv <- rep(NA, n.funcs)
  n.perm <- min.perm
  samps <- NULL

  if(!is.null(cov)){
    cov.scores <- asplit(scorer$score(cov, dim), 2)
  }else{
    cov.scores <- list()
    for(i in 1:n.complexes) cov.scores <- append(cov.scores, list(NULL)) #list of NULL
  }

  if(use.gpd){
    # log.p.gdp has previous approximations
    log.p.gpd <- NULL
    for(i in 1:n.complexes){
      log.p.gpd[[i]] <-  matrix(rep(NA, 3*n.funcs), nrow = n.funcs)
    }
  }

  while(length(idx.nc) > 0){
    #generating new samples
    if(is.null(samps)){
      samps <- add.samples(scorer, f, n.perm, dim, n.cores = n.cores, cov = cov,
                           n.complexes = n.complexes, old.samps = NULL)
    }else{
      n.new <- min(n.perm, max.perm - n.perm)
      n.perm <- n.perm + n.new
      samps <- add.samples(scorer, f, n.new, dim, n.cores = n.cores, cov = cov,
                           n.complexes = n.complexes, old.samps = samps)
    }
    #getting residues if applicable
    null.data <- mapply(get.null.samps, f.scores, cov.scores, samps, SIMPLIFY = F)

    for(i in 1:n.complexes){
      obs <- null.data[[i]][["observed"]]
      null.samps <- null.data[[i]][["sampled"]]

      if(use.gpd) new.log.p <- rep(NA, n.funcs) #vector for p-values obtained by gpd

      #tying to compute p-values
      for(k in 1:n.funcs){
        if(!is.na(p[idx.nc[k],i])) next

        if(is.infinite(obs[k])){
          p[idx.nc[k], i] <- 1
          next
        }
        n.le <- sum(null.samps[k,] <= obs[k])
        if(n.le >= 10){
          p[idx.nc[k], i] <- n.le/n.perm
          next
        }
        if(!use.gpd) next

        log.p.prev <- log.p.gpd[[i]][k,]
        gpd.out <- gpd.approx(null.samps, obs[k], pow = pow, nextremes = nextremes,
                              alpha = alpha)
        log.p <- gpd.out$log.p
        new.log.p[k] <- log.p
        log.q1 <- gpd.out$log.q1
        log.q3 <- gpd.out$log.q3

        if(any(is.na(log.p.prev)) || is.na(log.p)
           || is.na(log.q1) || is.na(log.q3)) next

        if(any(is.infinite(log.p.prev)) || is.infinite(log.p)
           || is.infinite(log.q1) || is.infinite(log.q3)) next

        #condition 1: relative variation of log p values is small
        cond1 <- all(abs(log.p.prev/ log.p - 1) < 0.1)
        #condition 2: quartiles are close to predicted values
        cond2 <- log.q1/log.p < 1.1 && log.q3/log.p > 0.9
        if(cond1 && cond2)  p[idx.nc[k], i] <- exp(log.p)
      }

      if(use.gpd){
        log.p.gpd[[i]] <- unname(cbind(log.p.gpd[[i]], new.log.p))
        # forget oldest approximation if the n.perm in the next iteration will be
        # 10x or more what it was when the orders approximation was done.
        if(n.perm/8 <= max.perm/10){
          log.p.gpd[[i]] <- log.p.gpd[[i]][, 2:ncol(log.p.gpd[[i]]), drop = F]
        }
      }
    }
    # convergent p-values
    conv <- apply(p[idx.nc, , drop = F], 1, function(x) !any(is.na(x)) )
    idx.new.conv <- idx.nc[conv]

    # combining convergent p-values
    if(n.complexes > 1 && sum(conv) > 0){
      for(k in 1:seq_along(idx.new.conv)){
        combined.obs <- matrix(nrow = n.perm, ncol = n.complexes)
        for(i in 1:n.complexes){
          combined.obs[,i] <- null.data[[i]][["sampled"]][which(conv)[k],]
        }
        combined.p[idx.new.conv[k]] <- combine.p.values(p[idx.new.conv[k], ],
                                                        combined.obs,
                                                        method = combination.method)
      }
    }

    # saving where convergence happened
    n.conv[idx.nc[conv]] <- n.perm

    # removing information we no longer need
    idx.nc <- idx.nc[!conv]
    f <- f[!conv, , drop = F]
    n.funcs <- nrow(f)
    for(i in 1:n.complexes){
      f.scores[[i]] <- f.scores[[i]][!conv]
      if(is.null(cov)){
        samps[[i]] <- samps[[i]][!conv, , drop = F]
      }else{
        samps[[i]][["func_scores"]] <- samps[[i]][["func_scores"]][!conv, , drop = F]
        samps[[i]][["cov_scores"]] <- samps[[i]][["cov_scores"]][!conv, , , drop = F]
      }
      if(use.gpd) log.p.gpd[[i]] <- log.p.gpd[[i]][!conv, , drop = F]
    }

    # if this is the last iterations, set non-convergent p-values to be the p-values
    # obtained by permutation
    if(n.perm == max.perm && n.funcs > 0){
      for(k in 1:n.funcs){
        if(n.complexes > 1) combined.obs <- matrix(nrow = n.perm, ncol = n.complexes)
        for(i in 1:n.complexes){
          obs <- null.data[[i]][["observed"]][which(!conv)[k]]
          null.samps <- null.data[[i]][["sampled"]][which(!conv)[k],]
          p[idx.nc[k], i] <- sum(null.samps <= obs) / n.perm
          if(n.complexes > 1) combined.obs[,i] <- null.samps
        }
        if(n.complexes > 1){
          combined.p[idx.nc[k]] <- combine.p.values(p[idx.nc[k],],
                                                    combined.obs,
                                                    method = combination.method)
        }
      }
    idx.nc <- NULL
    }
  }

  out$p <- p
  out$n.conv <- n.conv
  out$combined.p <- combined.p

  return(as.data.frame(out))
}
