

# Michael's Q test for 2xK designs with proportions

# Michael GA. 2007. A significance test of interaction in 2 x K designs with proportions. Tutorials in # Quantitative Methods for Psychology 3(1): 1-7. doi: 10.20982/tqmp.03.1.p001

# ps is a 2xK matrix of proportions
# ns is a 2xK matrix of sample sizes for each of the proportions
# alpha is the desired Type I error of the test

michaelq<-function(ps,ns,alpha){
  
  z<-qnorm(1-alpha/2)
  d<-ps[1,]-ps[2,]
  vars<-(z^2+4*ns*ps*(1-ps))/(4*(ns+z^2)^2)
  
  D<-vars[1,]+vars[2,]
  d0<-sum(d/D)/sum(1/D)
  
  Qp<-sum(((d-d0)^2)/D)
  nu<-length(d)+1
  
  pval<-pchisq(Qp,nu,lower.tail=FALSE)
  
  list("Difference"=d,
       "Variance"=vars,
       "Q"=Qp,
       "pval"=pval)
  
}