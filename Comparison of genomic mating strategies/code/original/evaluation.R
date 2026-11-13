# User-specified common evaluator. No software objective or simulated offspring.
# Z dosage=0/1/2, dominance W=I(Z=1), and NO intercept b.
pair_moments <- function(xf,xm,gamma,delta) {
  stopifnot(is.matrix(xf),identical(dim(xf),dim(xm)),
    ncol(xf)==length(gamma),length(gamma)==length(delta),
    all(xf %in% 0:2),all(xm %in% 0:2))
  u<-xf/2;v<-xm/2;s<-u+v;h<-s-2*u*v
  varZ<-u*(1-u)+v*(1-v);varW<-h*(1-h);covZW<-h*(1-s)
  meanGV<-drop(s%*%gamma+h%*%delta)
  partZ<-drop(varZ%*%(gamma^2));partW<-drop(varW%*%(delta^2))
  partZW<-drop(covZW%*%(2*gamma*delta));variance<-partZ+partW+partZW
  stopifnot(all(is.finite(meanGV)),all(is.finite(variance)),all(variance>=-1e-10))
  data.frame(ExpectedGV=meanGV,PairProgenyVarGV=pmax(variance,0),
    AdditiveScoreVar=partZ,DominanceScoreVar=partW,TwiceADCov=partZW)
}
evaluate_plan <- function(ix,p,sp,cfg,dir=NULL,seed=NULL) {
  t<-nowsec();check_plan(ix,p,sp);z<-decode(ix,p$N);C<-nrow(z)
  f<-pair_moments(p$Z[z$f,,drop=FALSE],p$Z[p$N+z$m,,drop=FALSE],p$a,p$d)
  F<-p$K[cbind(z$f,p$N+z$m)]/2
  # Native OSNB retains ORIGINAL scores. Evaluation uses the requested coding.
  offset<-p$intercept-sum(p$a)
  error<-max(abs(f$ExpectedGV+offset-p$M[cbind(z$f,z$m)]))
  stopifnot(error<1e-7,all(is.finite(F)))
  families<-cbind(data.frame(Family=seq_len(C),Dam=p$ids[z$f],Sire=p$ids[p$N+z$m],
    Weight=rep(1/C,C)),f,ExpectedF_GRM=F)
  list(families=families,
    metrics=c(ExpectedTotalGV=mean(f$ExpectedGV),
      MeanPairProgenyVarGV=mean(f$PairProgenyVarGV),MeanExpectedF_GRM=mean(F)),
    EvaluationSec=nowsec()-t,MaxMeanIdentityError=error,
    Formula='sum(E(Z)*gamma+E(W)*delta); no intercept',
    VarianceDefinition='mean per-pair progeny variance; NO between-family term; cross-locus covariance omitted')
}
