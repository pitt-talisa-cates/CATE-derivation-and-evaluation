
library("ic.infer")
library("scales")
library("sandwich")

## General functions for use

cmat_maker<-function(num.q){
  mat<-matrix(0,num.q-1,num.q)
  for(r in 1:nrow(mat)){
    mat[r,r]<-(-1)
    mat[r,r+1]<-1
  }
  mat
}

# Should modify the test function to output (y,w,score,b,s), where b is rf prediction and s is cf prediction
pval_constrain<-function(gates,type,num.q,vcov.type){
  
  if (type==1){   # at least one strict inequality
    #browser()
    ui1<-cmat_maker(num.q)
    gatesest<-ic.est(coef(gates)[2:(num.q+1)],
                     ui=ui1,
                     Sigma=vcovHC(gates,vcov.type)[2:(num.q+1),2:(num.q+1)])
    ic.test(gatesest)$p.value
  } else if (type==2){   # inequality between ATEs in extreme bins
    ui2<-t(matrix(c(-1,1),2,1))
    
    gatesest2<-ic.est(coef(gates)[c(2,(num.q+1))],
                      ui=ui2,
                      Sigma=vcovHC(gates,vcov.type)[c(2,(num.q+1)),c(2,(num.q+1))])
    ic.test(gatesest2)$p.value
  } 
}


############################################################################################################
# gates.fun fits models to estimate both GATES and BLP
#
# IMPORTANT: gates.fun expects scores to be increasing with E[Y1-Y0|x]; also
#            w=1 codes treatment and w=0 codes control. The user can control
#            whether the modeling estimates / tests correspond to effects of 
#            orientation E[Y0-Y1|x] instead (flips both w coding and scores).
#
# Arguments:
# master = list of length >=1 containing xyz (data.frames with y,w,score,b,etc.)
# pi = probability of w=1, as a scalar
# num.q = number of quantiles for GATEs
# loo = whether doing internal validation (loo=TRUE) or validation in hold-out sample (loo=FALSE)
# effect = Defaults to "Y1-Y0"-type effects. Use "Y0-Y1" if, for example, the target is risk reduction
# hyp.risk.cor.dir = Only applies if the inputted score is a risk score (estimate of E[Y|X]). 
#                    This is the hypothetical direction of correlation between outcome risk and the CATE,
#                    where the CATE scale is specified as in teh effect argument. Takes values "positive" or
#                    "negative". For example, if the user hypothesizes that this
#                    outcome risk is positively correlated with CATE on the "Y0-Y1" scale, set effect="Y0-Y1"
#                    and hyp.risk.cor.dir to "positive". Set to NULL if the score is a CATE estimate.
# blptest = Whether the user wants to test the score as a best linear predictor (BLP) of the CATE
# gates.mod = If "chernozhukov", uses lm to estimate gates, 95% CIs, and test two hypotheses:
#             1) whether GATEs are monotonically increasing; 2) whether first and last GATE are different.
#             If "prop.test", GATEs estimated via simple prop.tests on disjoint strata.
# alpha = desired type 1 error level (used for correcting p-values when loo=TRUE)
# q.breaks = number of quantiles desired (used for GATES analysis)
# vcov.type = type of robust standard errors (internally passed to vcovHC for use in GATEs model) 
############################################################################################################


# Default ("Y1-Y0") assumes: w=1 is treatment, w=0 is control, and score is
# an estimate of E[Y1-Y0|x] (ascending score = ascending Y1-Y0). All theta
# estimates, CIs, BLP coefficients, and monotonicity/extreme p-values are
# then on the Y1-Y0 scale, with quantile 1 = lowest predicted Y1-Y0.
#
# effect="Y0-Y1" reverses w, pi, and score consistently BEFORE any modeling
# (an ARR-style framing when Y is an adverse outcome). Flipping all three
# together -- rather than negating outputs after the fact -- means every
# returned quantity comes out already and correctly oriented to E[Y0-Y1|x],
# with quantile 1 now = lowest Y0-Y1 (i.e., the highest original score).
# No output needs post-hoc sign-flipping under this option.
#
# Note: if you supply q.breaks, supply them already negated to match
# effect="Y0-Y1", since breaks are applied to score after the flip.


gates.fun<-function(master,
                    pi,
                    num.q,
                    loo=FALSE,
                    effect=c("Y1-Y0","Y0-Y1"),
                    hyp.risk.cor.dir = NULL,
                    blptest=TRUE,
                    gates.mod=c("chernozhukov","prop.test"),
                    alpha=0.05,
                    q.breaks=NULL,
                    vcov.type="HC3"){
  
  effect <- match.arg(effect)
  gates.mod <- match.arg(gates.mod)
  
  if (is.null(pi) && (gates.mod=="chernozhukov" || blptest)) {
    stop("pi cannot be NULL when gates.mod='chernozhukov' or blptest=TRUE")
  }
  
  if (effect == "Y0-Y1") {
    warning(
      "gates.fun: effect='Y0-Y1' reverses w (1-w), pi (1-pi), and score (-score) ",
      "before modeling. All returned quantities (theta, CIs, BLP coefficients, ",
      "monotonicity/extreme p-values, and quantile membership) are on the ",
      "E[Y0-Y1|x] scale, with quantile 1 now corresponding to the lowest Y0-Y1 ",
      "(i.e., the highest original score). If you supply q.breaks, supply them ",
      "already negated to match.",
      call. = FALSE
    )
    master <- lapply(master, function(tt){
      tt$w <- 1 - tt$w
      tt$score <- -tt$score
      tt
    })
    pi <- 1 - pi
  }
  
  # Everything above assumes positive correlation between CATE and score
  # If we have a risk score, we may plausibly hypothesize negative correlation.
  # in which case we'd need to flip scoreagain.
  if (!is.null(hyp.risk.cor.dir)) {
    hyp.risk.cor.dir <- match.arg(hyp.risk.cor.dir, c("positive","negative"))
    if (hyp.risk.cor.dir == "negative") {
      master <- lapply(master, function(tt){ 
        tt$score <- -tt$score
        tt })
    }
  }
  
  if(!is.null(q.breaks)){
    if(num.q != (length(q.breaks)-1)) stop("Make sure that length(q.breaks)-1 equals num.q.")
  }
  
  thetalist<-ulist<-llist<-matrix(,length(master),num.q)
  pval_monoton<-pval_extreme<-pval_blptest<-rep(NA_real_,length(master))
  blp_intercept<-blp_slope<-rep(NA_real_,length(master))
  
  bb_all<-matrix(,length(master),num.q+1)
  n0_all<-matrix(,length(master),num.q)
  N0_all<-matrix(,length(master),num.q)
  n1_all<-matrix(,length(master),num.q)
  N1_all<-matrix(,length(master),num.q)
  
  
  for (i in 1:length(master)){

    t<-master[[i]]

    if(length(unique(quantile(t$score,probs=seq(0,1,1/num.q))))<(num.q+1)) {
      #browser()
      t$score<-rnorm(nrow(t),0,0.001) 
    }

    if(!is.null(q.breaks)){
      bb<-q.breaks
    } else {
      bb<-quantile(t$score,probs=seq(0,1,1/num.q))
      bb_all[i,]<-bb
    }

    t$e<-as.numeric(cut(t$score, 
                        breaks=bb,
                        names=FALSE,
                        include.lowest=TRUE))
    
    n0_all[i,]<-tapply(t$y[t$w==0],t$e[t$w==0],sum)
    N0_all[i,]<-tapply(t$y[t$w==0],t$e[t$w==0],length)
    n1_all[i,]<-tapply(t$y[t$w==1],t$e[t$w==1],sum)
    N1_all[i,]<-tapply(t$y[t$w==1],t$e[t$w==1],length)
    
    # If doing internal validation (loo=FALSE), then z quantile should be 2.28
    # Otherwise use the standard 1.96
    if(!loo){
      cov.coef<-qnorm(1-alpha/4)
    } else {
      cov.coef<-qnorm(1-alpha/2)
    }

    if(gates.mod=="chernozhukov" & !is.null(pi)){
      eis<-c()
      for(enum in 1:num.q){
        t$new.e<-ifelse(t$e==enum,1,0)*(t$w-pi)
        ei<-which(colnames(t)=="new.e")
        eis<-c(eis,ei)
        colnames(t)[ei]<-paste0("e",enum)
      }
      iy<-which(colnames(t)=="y")
      tg<-t[,c(iy,eis)]
      
      #Effects coding suggested in Chernozhukov 2018. Note slightly different betas 
      #   compared  with usual dummy coding. These differences completely go away if, 
      #   instead of subtracting 0.5 in the terms above, subtract the estimated propensity
      #   That said, the substantive results of the analysis do not change over the 100 splits.
      gates<-lm(y~.,data=tg)  
      #if(i %in% c(20,40,60,80)) browser()
      thetalist[i,]<-coef(gates)[2:(1+num.q)]
      ulist[i,]<-coef(gates)[2:(1+num.q)]+cov.coef*sqrt(diag(vcovHC(gates,vcov.type))[-1])
      llist[i,]<-coef(gates)[2:(1+num.q)]-cov.coef*sqrt(diag(vcovHC(gates,vcov.type))[-1])
      
      pval_monoton[i]<-pval_constrain(gates,1,num.q,vcov.type)
      pval_extreme[i]<-pval_constrain(gates,2,num.q,vcov.type)
      
    } else {  # If I haven't supplied Pr(w=1), or if i just want to estimate gates via simple prop.test
      #browser()
      ptlist<-lapply(1:num.q,function(x){
        dat<-t[t$e==x,]
        yy<-dat$y
        ww<-dat$w
        pt<-prop.test(c(sum(yy[ww==0]),sum(yy[ww==1])),
                  c(sum(ww==0),sum(ww==1)))
        list(pt$estimate[2]-pt$estimate[1],
             -pt$conf.int[2],
             -pt$conf.int[1])
      })
      
      mat<-matrix(unlist(ptlist),3,num.q)
      
      thetalist[i,]<-mat[1,]
      ulist[i,]<-mat[3,]
      llist[i,]<-mat[2,]
      
      # would need to modify the pval_constrain code to admit the above lm output
      pval_monoton[i]<-NA
      pval_extreme[i]<-NA
      
    }

    # Also do linear test of calibration?
    if(blptest){
      #browser()
      newt<-t #[!(t$e %in% c(2,3)),]
      score <- newt$score
      mean.score <- mean(score)
      # DF <- data.frame(y = newt$y-newt$b,  
      #                  main.score.term = (newt$w - pi) * mean.score, 
      #                  interact.score.term = (newt$w - pi) * (score - mean.score))
      
      # fit.blp <- lm(y ~ b + main.score.term + interact.score.term , 
      #               data = DF)
      DF <- data.frame(y = newt$y,
                       b = newt$b,
                       main.score.term = (newt$w - pi),
                       interact.score.term = (newt$w - pi) * (score-mean.score))
      
      fit.blp <- lm(y ~ b + main.score.term + interact.score.term , 
                    data = DF)
      
      #browser()
      #fit.blp<-glm(y~w+score+w*score,data=t,family="binomial")
      #fit.blp<-lm(y~w+score+w*score,data=t)
      blp.summary <- lmtest::coeftest(fit.blp, 
                                      vcov = sandwich::vcovCL, 
                                      type = "HC3")
      #browser()
      dimnames(blp.summary)[[2]][4] <- gsub("[|]", "", dimnames(blp.summary)[[2]][4])
      
      blp.summary[, 4] <- ifelse(blp.summary[, 3] < 0, 1 - blp.summary[,4]/2, blp.summary[,4]/2)
      #browser()
      pval_blptest[i]<-blp.summary[4,4]
      blp_intercept[i] <- blp.summary[3,1]
      blp_slope[i]     <- blp.summary[4,1]
      
    }
    
    #browser()
  }

  # The original VEIN result is predicated on having a set of p-values which are random conditional on the training set
  # if we generate a single p-value from the LOO estimates, we have a non-random p-value
  
  #browser()
  
  list("Theta Median" = apply(thetalist,2,median),
       "Theta Upper Median" = apply(ulist,2,median),
       "Theta Lower Median" =apply(llist,2,median),
       "Pval Monoton" = median(pval_monoton)*(2-loo*1),
       "Pval Extreme" = median(pval_extreme)*(2-loo*1),
       "Theta P25" = apply(thetalist,2,quantile,probs=0.25),
       "Theta P75" = apply(thetalist,2,quantile,probs=0.75),
       "BLP slope pval" = median(pval_blptest)*(2-loo*1),
       "BLP intercept coef"= median(blp_intercept),
       "BLP slope coef" = median(blp_slope),
       "N quantile" = table(t$e),
       "QUantile Thresholds All" = bb_all,
       "Mean prediction in quantiles" = tapply(t$score,t$e,mean),
       "n0 quantile" = n0_all,
       "N0 quantile" = N0_all,
       "n1 quantile" = n1_all,
       "N1 quantile" = N1_all)
}


