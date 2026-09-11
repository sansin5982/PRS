# Original course simulation; all people, loci and outcomes are fictional.
# No packages beyond base R are required. These functions never download files.
make_course_data <- function() {
  set.seed(7102026)
  m <- 120L
  n_by_role <- c(Discovery=5000L, Development=2000L, Tuning=1000L,
                 Test=2000L, LD_reference=1000L)
  role <- rep(names(n_by_role),n_by_role); n <- length(role)
  # Two haplotypes. Within each six-variant block, copy a shared allele
  # with probability .75; otherwise draw independently at frequency .30.
  # This creates local correlation, not a realistic genome-wide LD model.
  G <- matrix(0L,n,m)
  for (copy in 1:2) for (b in seq(1L,m,by=6L)) {
    shared <- rbinom(n,1,.3)
    for (j in b:(b+5L))
      G[,j] <- G[,j] + ifelse(runif(n)<.75,shared,rbinom(n,1,.3))
  }
  ids <- sprintf("S%05d",seq_len(n)); variants <- sprintf("sim%03d",seq_len(m))
  dimnames(G) <- list(ids,variants)
  truth <- rep(0,m); truth[seq(1,m,by=6)] <- rep(c(.16,-.08,.12,.06),5)
  fixed_weight <- truth + rnorm(m,0,.025)
  # Fixed numerical scale in the data-generating mechanism, not test-fitted.
  genetic <- as.numeric((G %*% truth - .6*sum(truth))/.5)
  age <- runif(n,40,75); SBP <- rnorm(n,130,15)
  probability <- plogis(-2.8+.045*(age-55)+.015*(SBP-120)+.8*genetic)
  CAD5 <- rbinom(n,1,probability)
  CAD5[role=="LD_reference"] <- NA_integer_ # LD outcomes are not supplied.
  d <- data.frame(IID=ids,role=role,age=age,SBP=SBP,CAD5=CAD5,
                  PRS=as.numeric(G %*% fixed_weight))
  weights <- data.frame(ID=variants,CHR=1L,POS=seq(100000L,by=10000L,length.out=m),
                        A1="A",A2="G",weight=fixed_weight)
  list(people=d,G=G,weights=weights,truth=truth,
       metadata=list(build="SIMULATION_ONLY",outcome="Fully observed five-year CAD",
       limitations="No deaths, censoring, population structure, relatedness or batch effects"))
}
auc <- function(y,p) {
  stopifnot(length(y)==length(p),!anyNA(y),all(y %in% c(0,1)),all(is.finite(p)))
  n1 <- sum(y==1); n0 <- sum(y==0)
  if(n1==0 || n0==0) return(NA_real_)
  (sum(rank(p,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0)
}
fit_course_models <- function(demo=make_course_data()) {
  d <- demo$people
  train <- d[d$role=="Development",]; test <- d[d$role=="Test",]
  center <- mean(train$PRS); spread <- sd(train$PRS)
  stopifnot(is.finite(spread),spread>0)
  transform_people <- function(x) {
    x$age10 <- (x$age-60)/10; x$sbp10 <- (x$SBP-120)/10
    x$PRS_z <- (x$PRS-center)/spread; x
  }
  train <- transform_people(train); test <- transform_people(test)
  base <- glm(CAD5~age10+sbp10,data=train,family=binomial())
  extended <- glm(CAD5~age10+sbp10+PRS_z,data=train,family=binomial())
  stopifnot(base$converged,extended$converged,
            all(is.finite(coef(base))),all(is.finite(coef(extended))))
  test$p_base <- predict(base,newdata=test,type="response")
  test$p_prs <- predict(extended,newdata=test,type="response")
  list(train=train,test=test,base=base,extended=extended,
       reference=c(mean=center,sd=spread))
}
summary_gwas <- function(demo) {
  ix <- which(demo$people$role=="Discovery")
  d <- demo$people[ix,]; G <- demo$G[ix,,drop=FALSE]
  estimates <- vapply(seq_len(ncol(G)),function(j) {
    fit <- glm(d$CAD5~G[,j]+d$age+d$SBP,family=binomial())
    if(!fit$converged) stop("Discovery model failed to converge")
    tab <- coef(summary(fit)); c(beta=tab[2,1],se=tab[2,2],p=tab[2,4])
  },numeric(3))
  cases <- sum(d$CAD5==1); controls <- sum(d$CAD5==0)
  stopifnot(cases>0,controls>0,all(is.finite(estimates)),all(estimates["se",]>0))
  data.frame(demo$weights[c("ID","CHR","POS","A1","A2")],
             beta=estimates["beta",],se=estimates["se",],p=estimates["p",],
             n_eff=4/(1/cases+1/controls),N=nrow(d),row.names=NULL)
}
paired_auc_bootstrap <- function(y,p0,p1,B=500L,seed=1606L) {
  stopifnot(length(y)==length(p0),length(y)==length(p1))
  set.seed(seed)
  delta <- replicate(B,{
    ix <- sample.int(length(y),replace=TRUE)
    auc(y[ix],p1[ix])-auc(y[ix],p0[ix])
  })
  valid <- is.finite(delta)
  if(sum(valid)<.95*B) stop("Too many undefined AUC resamples; revisit event counts")
  list(delta=auc(y,p1)-auc(y,p0),CI=quantile(delta[valid],c(.025,.975)),
       failed=sum(!valid),draws=delta)
}
net_benefit <- function(y,p,threshold) {
  stopifnot(threshold>0,threshold<1,length(y)==length(p),
            !anyNA(y),all(y %in% c(0,1)),all(is.finite(p)))
  positive <- p>=threshold
  mean(positive & y==1)-mean(positive & y==0)*threshold/(1-threshold)
}
