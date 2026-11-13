# Adapted from the original simulation.R: runMacs GENERIC, 10x1 Morgan, AD trait.
# Observed chip excludes QTL; oracle causal effects are supplied equally as relevant.
generate_data <- function(cfg,s,r,dir) {
  path<-file.path(dir,'problem.rds')
  if(file.exists(path)) {
    x<-readRDS(path);stopifnot(identical(x$cfg,cfg),x$s==s,x$r==r);return(x)
  }
  stopifnot(requireNamespace('AlphaSimR',quietly=TRUE))
  spec<-cfg$scales[s,];N<-spec$N;ns<-cfg$nsnp/cfg$nchr;nq<-cfg$nqtl/cfg$nchr
  stopifnot(ns==as.integer(ns),nq==as.integer(nq))
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  stream_seed(seed_for(cfg,s,r))
  ckpath<-file.path(dir,'founder_checkpoint.rds')
  if(file.exists(ckpath)) {
    ck<-readRDS(ckpath);stopifnot(identical(ck$cfg,cfg),ck$s==s,ck$r==r)
    founder<-do.call(methods::new,c(list(Class='MapPop'),ck$slots))
    assign('.Random.seed',ck$rng,envir=.GlobalEnv)
  } else {
    founder<-AlphaSimR::runMacs(nInd=2L*N,nChr=cfg$nchr,segSites=ns+nq,
      species='GENERIC',manualGenLen=1,ploidy=2L,inbred=FALSE,nThreads=1L)
    stopifnot(all(founder@nLoci==ns+nq))
    slots<-setNames(lapply(methods::slotNames(founder),function(k)methods::slot(founder,k)),methods::slotNames(founder))
    atomic(list(cfg=cfg,s=s,r=r,slots=slots,rng=.Random.seed),ckpath)
  }
  SP<-AlphaSimR::SimParam$new(founder);SP$nThreads<-1L
  SP$restrSegSites(overlap=TRUE)
  Zall<-AlphaSimR::pullSegSiteGeno(founder,simParam=SP)
  markers<-unlist(founder@genMap,use.names=FALSE) # names are defined on the map entries
  marker_names<-unlist(lapply(founder@genMap,names),use.names=FALSE)
  stopifnot(ncol(Zall)==cfg$nsnp+cfg$nqtl,length(marker_names)==ncol(Zall))
  chr<-rep(seq_len(cfg$nchr),founder@nLoci)
  maf<-pmin(colMeans(Zall)/2,1-colMeans(Zall)/2)
  # Same best-available MAF policy as original; no redraws.
  exclude<-unlist(lapply(seq_len(cfg$nchr),function(k){ix<-which(chr==k)
    if(sum(maf[ix]>=.05)>=nq)ix[maf[ix]<.05] else setdiff(ix,ix[order(-maf[ix],ix)][seq_len(nq)])}))
  if(length(exclude))SP$restrSegSites(excludeQtl=marker_names[exclude],overlap=TRUE)
  tp<-file.path(dir,'trait_checkpoint.rds')
  if(file.exists(tp)) {
    tc<-readRDS(tp);stopifnot(identical(tc$cfg,cfg),tc$s==s,tc$r==r)
    SP$manAddTrait(do.call(methods::new,c(list(Class='TraitAD'),tc$slots)))
    assign('.Random.seed',tc$rng,envir=.GlobalEnv)
  } else {
    SP$addTraitAD(nQtlPerChr=nq,mean=0,var=1,meanDD=cfg$meanDD,varDD=cfg$varDD,useVarA=TRUE)
    tr<-SP$traits[[1]];slots<-setNames(lapply(methods::slotNames(tr),function(k)methods::slot(tr,k)),methods::slotNames(tr))
    atomic(list(cfg=cfg,s=s,r=r,slots=slots,rng=.Random.seed),tp)
  }
  tr<-SP$traits[[1]]
  q<-as.integer(rep(c(0L,head(cumsum(founder@nLoci),-1L)),times=tr@lociPerChr)+tr@lociLoc)
  chip<-setdiff(seq_len(ncol(Zall)),q)
  stopifnot(length(q)==cfg$nqtl,length(chip)==cfg$nsnp,
    all(tabulate(chr[q],cfg$nchr)==nq),all(tabulate(chr[chip],cfg$nchr)==ns))
  SP$addSnpChipByName(marker_names[chip],name='observed_noncausal_50k')
  ids<-c(sprintf('F%04d',seq_len(N)),sprintf('M%04d',seq_len(N)))
  pop<-AlphaSimR::newPop(founder,id=ids,simParam=SP);pop@sex<-rep(c('F','M'),each=N)
  Z<-Zall[,q,drop=FALSE];G<-Zall[,chip,drop=FALSE];rm(Zall)
  rownames(Z)<-rownames(G)<-ids
  colnames(Z)<-marker_names[q];colnames(G)<-marker_names[chip]
  a<-as.numeric(tr@addEff);d<-as.numeric(tr@domEff);intercept<-as.numeric(tr@intercept)
  gv<-drop((Z-1)%*%a+(Z==1)%*%d)+intercept
  stopifnot(max(abs(gv-as.numeric(pop@gv)))<1e-8)
  freq<-colMeans(Z)/2;alpha<-a+d*(1-2*freq)
  bv<-drop(sweep(Z,2,2*freq,'-')%*%alpha)
  mt<-nowsec();pf<-Z[seq_len(N),,drop=FALSE]/2;pm<-Z[N+seq_len(N),,drop=FALSE]/2
  M<-outer(drop(pf%*%a),drop(pm%*%a),'+')-sum(a)+
    outer(drop(pf%*%d),drop(pm%*%d),'+')-2*tcrossprod(sweep(pf,2,d,'*'),pm)+intercept
  dimnames(M)<-list(ids[seq_len(N)],ids[N+seq_len(N)])
  mtime<-nowsec()-mt
  # Reference frequency is fixed before selection, identical for all methods.
  pobs<-colMeans(G)/2;den<-2*sum(pobs*(1-pobs))
  gc();t<-nowsec();W<-sweep(G,2,2*pobs,'-');K<-tcrossprod(W)/den;rm(W);gc()
  ktime<-nowsec()-t;dimnames(K)<-list(ids,ids)
  gp<-AlphaSimR::genParam(pop,simParam=SP)
  x<-list(cfg=cfg,s=s,r=r,N=N,Z=Z,G=G,M=M,K=K,a=a,d=d,alpha=alpha,qfreq=freq,
    gv=gv,bv=bv,intercept=intercept,ids=ids,qindex=q,chipindex=chip,
    maf_fallback_chr=sum(vapply(seq_len(cfg$nchr),function(k)sum(maf[chr==k]>=.05)<nq,logical(1))),
    MPreparationSec=mtime,KPreparationSec=ktime,realized=c(varA=gp$varA,varD=gp$varD,varG=gp$varG),
    map=data.frame(ID=marker_names,Chr=chr,PositionM=markers,QTL=seq_along(chr)%in%q),
    observed_frequency=pobs)
  stopifnot(identical(dim(M),c(N,N)),max(abs(K-t(K)))<1e-10,
    abs(mean(gv))<1e-8,length(intersect(q,chip))==0L)
  atomic(x,path)
  csv(data.frame(ID=ids,Sex=pop@sex,TotalGV=gv,ReferenceBV=bv),file.path(dir,'candidates.csv'))
  csv(data.frame(Marker=colnames(Z),a=a,d=d,alpha=alpha,ReferenceFrequency=freq),file.path(dir,'qtl_effects.csv'))
  csv(x$map,file.path(dir,'map.csv'));x
}
restore_population <- function(p,dir) {
  requireNamespace('AlphaSimR',quietly=TRUE)
  ck<-readRDS(file.path(dir,'founder_checkpoint.rds'));tc<-readRDS(file.path(dir,'trait_checkpoint.rds'))
  founder<-do.call(methods::new,c(list(Class='MapPop'),ck$slots))
  SP<-AlphaSimR::SimParam$new(founder);SP$nThreads<-1L
  SP$manAddTrait(do.call(methods::new,c(list(Class='TraitAD'),tc$slots)))
  SP$addSnpChipByName(p$map$ID[p$chipindex],name='observed_noncausal_50k')
  pop<-AlphaSimR::newPop(founder,id=p$ids,simParam=SP);pop@sex<-rep(c('F','M'),each=p$N)
  list(pop=pop,SP=SP)
}

# New experiment: generate the maximum population ONCE per replicate using the
# original generator above. Worker flock protects concurrent scale tasks.
prepare_replicate <- function(cfg,r,dir) {
  problem<-file.path(dir,'problem.rds');orderfile<-file.path(dir,'nested_order.rds')
  if(file.exists(problem)&&file.exists(orderfile))return(invisible(NULL))
  base<-generate_data(cfg,nrow(cfg$scales),r,dir)
  if(!file.exists(orderfile)) {
    stream_seed(seed_for(cfg,0L,r,2L))
    atomic(list(F=sample.int(base$N),M=sample.int(base$N),cfg=cfg),orderfile)
  }
  invisible(NULL)
}
nested_problem <- function(base,ordering,N) {
  stopifnot(N<=base$N,identical(base$cfg,ordering$cfg))
  f<-ordering$F[seq_len(N)];m<-ordering$M[seq_len(N)]
  ix<-c(f,base$N+m);p<-base
  p$N<-N;p$ids<-base$ids[ix]
  p$Z<-base$Z[ix,,drop=FALSE];p$G<-base$G[ix,,drop=FALSE]
  p$K<-base$K[ix,ix,drop=FALSE];p$M<-base$M[f,m,drop=FALSE]
  p$gv<-base$gv[ix];p$bv<-base$bv[ix]
  p$ReferencePopulationN<-length(base$ids);p$CandidateIndices<-ix
  # qfreq and observed_frequency deliberately remain on the SAME maximum-
  # population reference for all scales in this independent replicate.
  p
}
build_original_value_matrix <- function(p) {
  pf<-p$Z[seq_len(p$N),,drop=FALSE]/2
  pm<-p$Z[p$N+seq_len(p$N),,drop=FALSE]/2
  M<-outer(drop(pf%*%p$a),drop(pm%*%p$a),'+')-sum(p$a)+
    outer(drop(pf%*%p$d),drop(pm%*%p$d),'+')-
    2*tcrossprod(sweep(pf,2,p$d,'*'),pm)+p$intercept
  dimnames(M)<-list(p$ids[seq_len(p$N)],p$ids[p$N+seq_len(p$N)])
  M
}
# OSNB-specific score and parent-score construction is INCLUDED in its timer.
# Original run_osnb, original nested functions and embedded C++ are untouched.
run_osnb_timed <- function(p,sp,cfg,api,seed) {
  t<-nowsec();p$M<-build_original_value_matrix(p)
  p$gv<-drop((p$Z-1)%*%p$a+(p$Z==1)%*%p$d)+p$intercept
  prep<-nowsec()-t;ans<-run_osnb(p,sp,cfg,api,seed)
  ans$ScorePreparationSec<-prep;ans$RuntimeSec<-nowsec()-t;ans
}
