scales <- data.frame(Scale=paste0('N',seq(100L,600L,100L)),N=seq(100L,600L,100L),
  nf=seq(12L,72L,12L),nm=seq(6L,36L,6L),pool=rep(200L,6L))
settings <- function(smoke=FALSE) {
  profile <- Sys.getenv('GM_PROFILE','inbreeding_restored')
  if(profile!='inbreeding_restored')stop('This package requires GM_PROFILE=inbreeding_restored; unset an old GM_PROFILE')
  reuse_osnb<-Sys.getenv('REUSE_OSNB','0')
  if(!reuse_osnb %in% c('0','1'))stop('REUSE_OSNB must be 0 or 1')
  list(schema='six-scales-inbreeding-restored-v1',seed=20261028L,
    profile=profile,smoke=smoke,nrep=if(smoke)1L else 100L,
    reuse_previous_osnb=(!smoke && reuse_osnb=='1'),
    scales=if(smoke)data.frame(Scale='smoke',N=20L,nf=6L,nm=3L,pool=200L) else scales,
    nsnp=if(smoke)100L else 50000L,nqtl=if(smoke)20L else 200L,nchr=10L,
    meanDD=.5,varDD=.1,osnb_starts=2L,osnb_mode='best',
    gm_pop=if(smoke)12L else 100L,gm_iter=if(smoke)2L else 300L,
    gm_starts=2L,gm_elite=if(smoke)2L else 10L,gm_mutation=.8,
    gm_var=.2,gm_inbreed=.7,
    alpha_solutions=if(smoke)20L else 100L,
    alpha_iter=if(smoke)2000L else 100000L,
    alpha_stagnation=if(smoke)30L else 500L,
    alpha_degree=30,alpha_min_inbreed_pct=75,alpha_tolerance=1e-5,
    evaluation='parental analytic AD; Z=0/1/2; W=I(Z=1); NO intercept',
    variance='mean pair variance; cross-locus covariances omitted',
    grm_reference='maximum candidate population within replicate')
}
nowsec <- function() unname(proc.time()[3])
atomic <- function(x,p) {
  dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
  if(file.exists(p))stop('Checkpoint exists: ',p)
  tmp<-tempfile('pending_',dirname(p));saveRDS(x,tmp,version=3)
  if(!file.rename(tmp,p))stop('Cannot commit ',p)
}
csv <- function(x,p) {dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);write.csv(x,p,row.names=FALSE)}
seed_for <- function(cfg,s,r,stage=0L) cfg$seed+s*10000000L+r*10000L+stage
stream_seed <- function(seed) {RNGkind("L'Ecuyer-CMRG","Inversion","Rejection");set.seed(seed)}
