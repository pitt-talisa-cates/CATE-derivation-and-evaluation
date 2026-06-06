

# This function takes the cv_k datasets from oest.cv and post-processes for plotting
#
# Arguments:
# 1. cv.obj      : output from oest.cv
# 2. stats       : summary statistics function(s) to use (via match.fun)
# 3. scale       : number of x-axis evaluation points for both outmat and ratelist.
#                  Must be a positive integer. If scale exceeds the number of unique
#                  score values (pooled/final-model) or the minimum fold size
#                  (fold-specific), scale is automatically reduced and a warning issued.
# 4. cores       : number of cores for mclapply (if cores=1, uses lapply)
# 5. direction   : "both", "greater", or "less"
# 6. score       : "cate" or "risk"
# 7. threshold   : method for defining x-axis evaluation points:
#
#   "pooled"        Quantiles of pooled out-of-fold scores. X-axis is raw score.
#                   Conservative under the null; stable cross-fold scale.
#                   Default and equivalent to v3 behavior.
#
#   "fold-specific" Evenly spaced proportion grid on [0,1]. Each fold contributes
#                   ARDs at its own nearest quantile. X-axis is proportion treated.
#                   Consistent with Zhao et al. (2013) ADa(q) estimator.
#                   Most sensitive but potentially anti-conservative under null.
#
#   "final-model"   Quantiles of scores from a final model applied to full data.
#                   Requires final.scores argument. X-axis is raw score.
#                   Closest to Yadlowsky et al. (2025) fixed-rule framework.
#
# 8. final.scores : vector of final-model scores for all patients (required if
#                   threshold="final-model")
# 9. se.method    : inference method for AUTOC/AUQINI SEs: "none" (default),
#                   "bootstrap", or "perturbation" (see rate_se for details)
# 10. B           : number of replicates for se.method (ignored if se.method="none")


# ---- Outcome stat functions (unchanged from v3) ----

stats.bin <- function(data=NULL,
                      ntest=nrow(data),
                      alternative="less",
                      names=FALSE) {
  if (names==FALSE) {
    n.w1 <- length(data[data$W==1,"Y"])
    n.w0 <- length(data[data$W==0,"Y"])
    x.w1 <- sum(data[data$W==1,"Y"])
    x.w0 <- sum(data[data$W==0,"Y"])
    ss   <- mean(data$score)
    if (n.w1>0 & n.w0>0) {
      p1    <- (x.w1/n.w1)
      p0    <- (x.w0/n.w0)
      ad    <- p1-p0
      pp    <- (n.w1*p1 + n.w0*p0)/(n.w1+n.w0)
      denom <- sqrt(pp*(1-pp)*(n.w1^-1 + n.w0^-1))
      x2    <- (ad/denom)^2
      pval  <- pnorm(ad/denom, mean=0, sd=1,
                     lower.tail=ifelse(alternative=="less",TRUE,FALSE))
      c(ad, x2, pval, (pval<0.05)*1, ss)
    } else {
      rep(NA, 5)
    }
  } else {
    c("ad","x2","pval","pval.lt0.05","mean.score")
  }
}


stats.bin.unkp <- function(data=NULL,
                           ntest=nrow(data),
                           alternative="less",
                           names=FALSE) {
  if (names==FALSE) {
    ww  <- data$W
    yy  <- data$Y
    px  <- data$px
    m1x <- data$m1x
    m0x <- data$m0x
    if (length(yy[ww==0])>0 & length(yy[ww==1])>0) {
      p1i <- (ww*yy/px) - (m1x*(ww-px)/px)
      p0i <- ((1-ww)*yy/(1-px)) - m0x*((ww-px)/(1-px))
      ad  <- mean(p1i-p0i)
      c(ad, NA, NA, NA)
    } else {
      rep(NA, 4)
    }
  } else {
    c("ad","x2","pval","pval.lt0.05")
  }
}


stats.norm <- function(data=NULL,
                       ntest=nrow(data),
                       alternative="less",
                       names=FALSE) {
  if (names==FALSE) {
    y1   <- data[data$W==1,"Y"]
    y0   <- data[data$W==0,"Y"]
    n1   <- length(y1)
    n0   <- length(y0)
    x.w1 <- sum(y1)
    x.w0 <- sum(y0)
    if (n1>0 & n0>0) {
      mu1      <- (x.w1/n1)
      mu0      <- (x.w0/n0)
      var1     <- sum((y1-mu1)^2)/(n1-1)
      var0     <- sum((y0-mu0)^2)/(n0-1)
      ad       <- mu1-mu0
      var_satt <- var1/n1 + var0/n0
      df       <- (var_satt^2)/((var1/n1)^2/(n1-1) + (var0/n0)^2/(n0-1))
      t        <- ad/sqrt(var_satt)
      pval     <- pt(t, df, lower.tail=ifelse(alternative=="less",TRUE,FALSE))
      c(ad, t, pval, (pval<0.05)*1)
    } else {
      rep(NA, 4)
    }
  } else {
    c("ad","t","pval","pval.lt0.05")
  }
}


# ---- Main function ----

post.cv.fun <- function(cv.obj,
                        stats        = c("stats.bin"),
                        scale        = 10,
                        cores        = 1,
                        direction    = "both",
                        score        = "cate",
                        threshold    = "pooled",
                        final.scores = NULL,
                        se.method    = "none",
                        B            = 200) {
  # se.method controls how standard errors for AUTOC and AUQINI are computed:
  #   "none"          No SEs computed (default)
  #   "bootstrap"     Half-sample bootstrap (Yadlowsky et al. 2025 Corollary 5):
  #                   draw floor(n/2) obs without replacement within each fold,
  #                   reconstruct averaged TOC curve, compute RATE metrics.
  #                   Valid for all three threshold options.
  #                   Note: for "fold-specific", q values are recomputed within
  #                   each half-sample to preserve proportion scale correctness.
  #   "perturbation"  Perturbation resampling (Zhao et al. 2013 Appendix E):
  #                   multiply each observation's Y and W contributions by iid
  #                   Exp(1) weights, reconstruct averaged TOC curve, compute
  #                   RATE metrics. Valid for all three threshold options.
  #                   Recommended for "fold-specific" as it avoids q recomputation
  #                   and is directly justified by Zhao et al. theory.
  
  # Validate arguments
  if (is.null(scale) || !is.numeric(scale) || length(scale) != 1 ||
      scale < 1 || scale != round(scale))
    stop('scale must be a positive integer (e.g. scale=20).')
  if (!threshold %in% c("pooled","fold-specific","final-model"))
    stop('threshold must be one of "pooled", "fold-specific", or "final-model"')
  if (threshold == "final-model" & is.null(final.scores))
    stop('final.scores must be provided when threshold="final-model"')
  if (!se.method %in% c("none","bootstrap","perturbation"))
    stop('se.method must be one of "none", "bootstrap", or "perturbation"')
  
  # Apply a jitter in case there are tied score values
  cv.obj <- lapply(cv.obj, function(fold) {
    dat <- fold[[1]]
    unique_scores <- sort(unique(dat$score))
    if (length(unique_scores) < nrow(dat)) {
      warning("Non-unique scores found. Jittering to avoid issues with ranking.")
      jitter_bound <- if (length(unique_scores) == 1) 1e-10 else min(diff(unique_scores)) / 2
      dat$score <- dat$score + runif(nrow(dat), -jitter_bound, jitter_bound)
    }
    list(dat)
  })
  
  # ---- Inner functions ----
  
  summarize.outer <- function(cv.ds, drop.high.first=TRUE, stats) {
    dat     <- cv.ds[[1]]
    dat.ord <- dat[order(dat$score, decreasing=drop.high.first),]
    alt     <- ifelse((drop.high.first==TRUE  & score=="cate") |
                        (drop.high.first==FALSE & score=="risk"), "less", "greater")
    new.list <- lapply(1:nrow(dat.ord), summarize.inner, dat.ord, alt, stats)
    new.df   <- as.data.frame(t(matrix(unlist(new.list),, nrow(dat.ord))))
    names(new.df)[1]  <- "q"
    names(new.df)[-1] <- match.fun(stats)(names=TRUE)
    cbind(dat.ord, new.df)
  }
  
  summarize.inner <- function(k, ordered.dat, alternative="less", stats) {
    dat.sub  <- if (k>1) ordered.dat[-c(1:(k-1)),] else ordered.dat
    stat.out <- match.fun(stats)(dat.sub, nrow(dat.sub), alternative)
    c(1-((k-1)/nrow(ordered.dat)), stat.out)
  }
  
  # ---- Build evaluation grid (sq) ----
  # For "pooled" and "final-model": sq is raw score thresholds
  # For "fold-specific": sq is evenly spaced proportions on (0,1]
  # In all cases, scale is capped at the maximum meaningful resolution
  # for the data and a warning issued if it was reduced.
  
  nte       <- nrow(cv.obj[[1]][[1]])
  nte.trunc <- nte - 30
  
  if (threshold == "pooled") {
    all.scores  <- unlist(lapply(1:length(cv.obj), function(x)
      cv.obj[[x]][[1]][31:nte.trunc, "score"]))
    n.unique    <- length(unique(all.scores))
    if (scale > n.unique) {
      warning(paste0("scale=", scale, " exceeds the number of unique pooled scores (",
                     n.unique, "). Reducing scale to ", n.unique, "."))
      scale <- n.unique
    }
    scores.q <- quantile(all.scores, prob=seq(0, 1, 1/scale)[-1])
    
  } else if (threshold == "final-model") {
    n.unique <- length(unique(final.scores))
    if (scale > n.unique) {
      warning(paste0("scale=", scale, " exceeds the number of unique final model scores (",
                     n.unique, "). Reducing scale to ", n.unique, "."))
      scale <- n.unique
    }
    scores.q <- quantile(final.scores, prob=seq(0, 1, 1/scale)[-1])
    
  } else {
    # fold-specific: cap scale at minimum fold size
    min.fold.n <- min(sapply(cv.obj, function(x) nrow(x[[1]])))
    if (scale > min.fold.n) {
      warning(paste0("scale=", scale, " exceeds the minimum fold size (",
                     min.fold.n, "). Reducing scale to ", min.fold.n,
                     " to avoid duplicate proportion points."))
      scale <- min.fold.n
    }
    scores.q <- seq(1/scale, 1, 1/scale)
  }
  
  # ---- Build summary lists ----
  
  if (direction == "both") {
    
    summary.list.less    <- lapply(cv.obj, summarize.outer,
                                   drop.high.first=TRUE,  stats)
    summary.list.greater <- lapply(cv.obj, summarize.outer,
                                   drop.high.first=FALSE, stats)
    
    cols   <- colnames(summary.list.less[[1]])
    remove <- which(cols %in% c("Y","W","score","q"))
    stats  <- cols[-remove]
    
    outlist  <- list()
    ratelist <- list()
    k <- 1
    
    for (dir in c(1,-1)) {
      
      ldat <- if (dir==1) summary.list.less else summary.list.greater
      sq   <- if (dir==1) rev(scores.q) else scores.q
      
      # ---- Build outmat_ ----
      
      outmat_ <- data.frame(sq)
      
      for (st in stats) {
        
        if (threshold == "fold-specific") {
          # Proportion-based lookup: for each target proportion q,
          # find the nearest q value >= q_target within each fold
          # (since q runs from 1 down to 1/n as we drop more patients,
          # we want the smallest q still >= q_target, i.e. the row that
          # has at least q_target proportion of patients remaining)
          cv.mean <- unlist(lapply(sq, function(q_target) {
            out <- unlist(lapply(ldat, function(y) {
              eligible <- y$q[y$q >= q_target]
              if (length(eligible) > 0)
                y[,st][y$q == min(eligible)][1]
              else
                NA
            }))
            mean(out, na.rm=TRUE)
          }))
          
        } else {
          # Raw score lookup: pooled or final-model
          # Identical to v3 logic
          cv.mean <- unlist(lapply(sq, function(x) {
            out <- unlist(lapply(ldat, function(y) {
              if (sum(dir*y$score <= dir*x) > 0)
                y[,st][dir*y$score == max(dir*y$score[dir*y$score <= dir*x])[1]]
              else
                NA
            }))
            mean(out, na.rm=TRUE)
          }))
        }
        
        outmat_ <- cbind(outmat_, cv.mean)
      }
      
      outlist[[k]] <- outmat_
      
      # ---- Compute RATE metrics from outmat_ ----
      
      rate.point    <- rate_from_outmat(outmat_, sq, dir, threshold)
      ratelist[[k]] <- data.frame(
        AUQINI        = rate.point["AUQINI"],
        AUTOC         = rate.point["AUTOC"],
        `Kendall Tau` = rate.point["Kendall Tau"],
        `Adj AUQINI`  = rate.point["Adj AUQINI"],
        check.names   = FALSE
      )
      
      if (se.method != "none") {
        se.rates              <- rate_se(ldat, sq, dir, threshold,
                                         method=se.method, B=B)
        ratelist[[k]]$SE.AUQINI <- se.rates["SE.AUQINI"]
        ratelist[[k]]$SE.AUTOC  <- se.rates["SE.AUTOC"]
      }
      
      names(ratelist)[k] <- if (dir==1) "Drop High First" else "Drop Low First"
      k <- k+1
    }
    
    # Name outmat columns -- first column label reflects threshold type
    xlab.high <- if (threshold=="fold-specific") "q_DropHighFirst"     else "scores.q_DropHighFirst"
    xlab.low  <- if (threshold=="fold-specific") "q_DropLowFirst"      else "scores.q_DropLowFirst"
    
    drop_high_names <- unlist(lapply(stats, function(x) paste0(x,"_DropHighFirst")))
    drop_low_names  <- unlist(lapply(stats, function(x) paste0(x,"_DropLowFirst")))
    
    colnames(outlist[[1]]) <- c(xlab.high, drop_high_names)
    colnames(outlist[[2]]) <- c(xlab.low,  drop_low_names)
    
    outmat <- cbind(outlist[[1]], outlist[[2]])
    
  } else if (direction == "greater") {
    summary.list.less    <- NA
    summary.list.greater <- lapply(cv.obj, summarize.outer,
                                   drop.high.first=FALSE, stats)
    outmat   <- NULL
    ratelist <- NULL
  } else {
    summary.list.less    <- lapply(cv.obj, summarize.outer,
                                   drop.high.first=TRUE, stats)
    summary.list.greater <- NA
    outmat   <- NULL
    ratelist <- NULL
  }
  
  list("outmat"    = outmat,
       "ratelist"  = ratelist,
       "threshold" = threshold)
}


# ---- rate_from_outmat ----
# Compute AUTOC, AUQINI, Kendall Tau, and Adj AUQINI directly from the
# averaged TOC curve in outmat_, using the trapezoidal rule.
#
# For "pooled" and "final-model": sq is in raw score units, so the
# proportion axis u_j is derived from rank position (j/scale).
# For "fold-specific": sq is already on [0,1], so u_j = sq directly.
#
# In both cases AUTOC and AUQINI are computed on the proportion scale,
# consistent with Yadlowsky et al. (2025) definitions.

rate_from_outmat <- function(outmat_, sq, dir, threshold) {
  
  nb <- nrow(outmat_)
  
  u_j <- if (threshold == "fold-specific") as.numeric(sq) else seq(1/nb, 1, 1/nb)
  
  # Append u=1 anchor (TOC = 0 by definition at u=1)
  #gate     <- c(outmat_$cv.mean, outmat_$cv.mean[nb])
  gate     <- c(outmat_$cv.mean, outmat_$cv.mean[1])
  u_full   <- c(u_j, 1)
  ate      <- gate[max(which(is.finite(gate)))]
  toc      <- gate - ate
  
  du       <- u_full[-1] - u_full[-length(u_full)]
  toc_avg  <- (toc[-length(toc)] + toc[-1]) / 2
  
  # AUTOC: area under TOC curve
  autoc  <- sum(du * toc_avg, na.rm=TRUE)
  
  # AUQINI: linear (Qini) weighting -- alpha(u) = u
  u_mid  <- (u_full[-1] + u_full[-length(u_full)]) / 2
  auqini <- sum(u_mid * toc_avg * du, na.rm=TRUE)
  
  # Kendall tau: rank correlation of ARD with evaluation point order
  tau     <- cor(outmat_$cv.mean, 1:nb, method="kendall", use="complete.obs")
  
  adjqini <- tau * auqini
  
  c(AUQINI=auqini, AUTOC=autoc, `Kendall Tau`=tau, `Adj AUQINI`=adjqini)
}


# ---- rate_se ----
# Compute bootstrap or perturbation SEs for AUTOC and AUQINI.
#
# method="bootstrap"
#   Half-sample bootstrap per Yadlowsky et al. (2025) Corollary 5.
#   Draws floor(n/2) observations without replacement within each fold,
#   reconstructs the averaged TOC curve using the same sq grid, computes
#   RATE metrics. For "fold-specific", q values are recomputed within the
#   half-sample so the proportion scale remains correct.
#
# method="perturbation"
#   Perturbation resampling per Zhao et al. (2013) Appendix E.
#   Each observation i receives an iid weight Vi ~ Exp(1). The weighted
#   ARD at each threshold is computed as:
#     ad* = sum(V*Y*W)/sum(V*W) - sum(V*Y*(1-W))/sum(V*(1-W))
#   The averaged TOC curve is reconstructed with these weighted ARDs.
#   The threshold grid sq stays fixed. No rows are removed so q values
#   remain valid for "fold-specific". Directly justified by Zhao et al.
#   theory for the CV-averaged estimator.
#
# In both cases the SE is the SD of AUTOC (or AUQINI) across B replicates.

rate_se <- function(ldat, sq, dir, threshold, method="bootstrap", B=200) {
  
  se.estimates <- matrix(NA, nrow=B, ncol=2,
                         dimnames=list(NULL, c("AUQINI","AUTOC")))
  
  for (b in 1:B) {
    
    if (method == "bootstrap") {
      
      # Draw half-sample without replacement within each fold
      bldat <- lapply(ldat, function(x) {
        nn  <- nrow(x)
        ss  <- sample(1:nn, floor(nn/2), replace=FALSE)
        sub <- x[ss, ]
        # Recompute q within the half-sample so proportion scale is correct
        # q = 1 - (rank - 1) / n_half, where rank is row position in the
        # already-sorted dataset (summarize.outer sorted by score)
        n_half   <- nrow(sub)
        sub$q    <- 1 - (seq_len(n_half) - 1) / n_half
        sub
      })
      
      if (threshold == "fold-specific") {
        cv.mean.b <- unlist(lapply(sq, function(q_target) {
          out <- unlist(lapply(bldat, function(y) {
            eligible <- y$q[y$q >= q_target]
            if (length(eligible) > 0)
              y[,"ad"][y$q == min(eligible)][1]
            else
              NA
          }))
          mean(out, na.rm=TRUE)
        }))
      } else {
        cv.mean.b <- unlist(lapply(sq, function(x) {
          out <- unlist(lapply(bldat, function(y) {
            if (sum(dir*y$score <= dir*x) > 0)
              y[,"ad"][dir*y$score == max(dir*y$score[dir*y$score <= dir*x])[1]]
            else
              NA
          }))
          mean(out, na.rm=TRUE)
        }))
      }
      
    } else {
      # method == "perturbation"
      # Draw iid Exp(1) weights for every observation in every fold.
      # Each fold gets its own weight vector. Recompute weighted ARDs
      # at each threshold point, keeping sq fixed.
      wldat <- lapply(ldat, function(x) {
        x$V <- rexp(nrow(x), rate=1)   # Exp(1): mean=1, var=1
        x
      })
      
      if (threshold == "fold-specific") {
        cv.mean.b <- unlist(lapply(sq, function(q_target) {
          out <- unlist(lapply(wldat, function(y) {
            eligible <- y$q[y$q >= q_target]
            if (length(eligible) > 0) {
              sub <- y[y$q == min(eligible), ]
              # weighted ARD = sum(V*Y*W)/sum(V*W) - sum(V*Y*(1-W))/sum(V*(1-W))
              # But at a single matched row, sub is one row — we need the
              # cumulative subgroup, so subset to q >= q_target
              sub <- y[y$q >= q_target, ]
              n1  <- sum(sub$V * (sub$W == 1))
              n0  <- sum(sub$V * (sub$W == 0))
              if (n1 > 0 & n0 > 0)
                sum(sub$V * sub$Y * (sub$W==1)) / n1 -
                sum(sub$V * sub$Y * (sub$W==0)) / n0
              else NA
            } else NA
          }))
          mean(out, na.rm=TRUE)
        }))
      } else {
        cv.mean.b <- unlist(lapply(sq, function(x) {
          out <- unlist(lapply(wldat, function(y) {
            sub <- y[dir*y$score <= dir*x, ]
            if (nrow(sub) > 0) {
              n1 <- sum(sub$V * (sub$W == 1))
              n0 <- sum(sub$V * (sub$W == 0))
              if (n1 > 0 & n0 > 0)
                sum(sub$V * sub$Y * (sub$W==1)) / n1 -
                sum(sub$V * sub$Y * (sub$W==0)) / n0
              else NA
            } else NA
          }))
          mean(out, na.rm=TRUE)
        }))
      }
    }
    
    outmat_b <- data.frame(sq=sq, cv.mean=cv.mean.b)
    rate.b   <- rate_from_outmat(outmat_b, sq, dir, threshold)
    
    se.estimates[b,"AUQINI"] <- rate.b["AUQINI"]
    se.estimates[b,"AUTOC"]  <- rate.b["AUTOC"]
  }
  
  c(SE.AUQINI = sd(se.estimates[,"AUQINI"], na.rm=TRUE),
    SE.AUTOC  = sd(se.estimates[,"AUTOC"],  na.rm=TRUE))
}