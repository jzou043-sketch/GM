encode <- function(f,m,N) as.integer((f-1L)*N+m)
decode <- function(ix,N) data.frame(f=(as.integer(ix)-1L)%/%N+1L,m=(as.integer(ix)-1L)%%N+1L)
check_plan <- function(ix,p,sp) {
  stopifnot(length(ix)==sp$nf,!anyNA(ix),all(ix>=1L&ix<=p$N^2))
  z<-decode(ix,p$N);stopifnot(length(unique(z$f))==sp$nf,
    length(unique(z$m))==sp$nm,all(table(z$m)==2L));invisible(TRUE)
}
prepare_osnb <- function(root) {
  e<-new.env(parent=globalenv());sys.source(file.path(root,'original','osnb_adapter.R'),e)
  e$load_original(file.path(root,'original'))
}
osnb_frames <- function(p) list(f=data.frame(ID=rownames(p$M),GeneticValue=p$gv[1:p$N]),
  m=data.frame(ID=colnames(p$M),GeneticValue=p$gv[p$N+1:p$N]))
run_osnb <- function(p,sp,cfg,api,seed) {
  t<-nowsec();e<-api$build();v<-osnb_frames(p)
  vs<-e$build_value_space(p$M,sp$pool,sp$pool)
  ps<-e$build_parent_value_space(v$f,v$m,sp$pool,sp$pool)
  fpool<-union(vs$female_space,ps$female_space);mpool<-union(vs$male_space,ps$male_space)
  ans<-e$gm_multistart_all_runs(p$M,v$f,v$m,sp$nf,sp$nm,sp$pool,sp$pool,
    cfg$osnb_starts,seed,TRUE,cfg$osnb_mode,FALSE)
  ij<-which(ans$assignment==1,arr.ind=TRUE)
  ix<-encode(match(rownames(ans$assignment)[ij[,1]],rownames(p$M)),
    match(colnames(ans$assignment)[ij[,2]],colnames(p$M)),p$N)
  check_plan(ix,p,sp)
  stopifnot(abs(sum(p$M[cbind(decode(ix,p$N)$f,decode(ix,p$N)$m)])-ans$objective)<1e-7)
  list(plan=ix,answer=ans,RuntimeSec=nowsec()-t,PoolF=length(fpool),PoolM=length(mpool),
       CandidatePoolF=fpool,CandidatePoolM=mpool,NativeObjective=ans$objective)
}
repair_factory <- function(p,sp,seed) {
  N<-p$N;set.seed(seed);pf<-sample.int(N);pm<-sample.int(N);calls<-repaired<-0L
  repair<-function(ix) {
    calls<<-calls+1L;z<-decode(ix,N)
    stopifnot(nrow(z)==sp$nf,all(ix>=1&ix<=N^2))
    if(length(unique(z$f))==sp$nf&&length(unique(z$m))==sp$nm&&all(table(z$m)==2L))
      return(encode(z$f[order(z$f)],z$m[order(z$f)],N))
    repaired<<-repaired+1L
    fs<-order(-tabulate(z$f,N),pf)[1:sp$nf];ms<-order(-tabulate(z$m,N),pm)[1:sp$nm]
    loads<-setNames(integer(length(ms)),ms);ans<-integer(length(fs))
    for(k in seq_along(fs)) {
      available<-ms[loads<2];pref<-z$m[z$f==fs[k]&z$m%in%available]
      m<-if(length(pref))pref[1] else available[which.min(pm[available])]
      ans[k]<-m;loads[as.character(m)]<-loads[as.character(m)]+1L
    }
    out<-encode(fs,ans,N);check_plan(out,p,sp);out[order(fs)]
  }
  list(repair=repair,audit=function()list(calls=calls,repaired=repaired,priorityF=pf,priorityM=pm))
}
compact_stats <- function(Markers,K,markereffects,P) {
  # Exact removal of all-zero contribution columns; native package statistics unchanged.
  active<-which(colSums(P)>0)
  GenomicMating::getstatsM1(Markers[active,,drop=FALSE],K[active,active,drop=FALSE],
    markereffects,P[,active,drop=FALSE])
}
adapt_ga <- function(repair,guide=NULL) {
  native<-GenomicMating::getGaSolutions;f<-native;counts<-c(initial=0L,offspring=0L,pairs=0L)
  first<-TRUE
  init<-function(Candidates,ntoselect) {
    if(first&&!is.null(guide)){first<<-FALSE;return(guide)}
    repair(sample(Candidates,ntoselect,replace=TRUE))
  }
  walk<-function(x) {
    if(missing(x))return(quote(expr=))
    if(!is.call(x))return(x)
    if(identical(x,quote(return(sample(Candidates, ntoselect, replace = TRUE))))) {
      counts['initial']<<-counts['initial']+1L;return(quote(return(init(Candidates,ntoselect))))
    }
    if(identical(x,quote(return(sort(cross))))) {
      counts['offspring']<<-counts['offspring']+1L;return(quote(return(repair(sort(cross)))))
    }
    # Avoid quadratic-memory/time repeated rbind while retaining identical pair order.
    if(identical(x[[1]],as.name('for'))&&identical(x[[2]],as.name('i'))&&
       identical(x[[3]],quote(1:N1))) {
      counts['pairs']<<-counts['pairs']+1L
      return(quote(allcombs<-cbind(rep(seq_len(N1),each=N2),rep(N1+seq_len(N2),times=N1))))
    }
    for(i in seq_along(x))if(i>1L)x[[i]]<-walk(x[[i]])
    x
  }
  body(f)<-walk(body(f));stopifnot(identical(unname(counts),c(1L,1L,1L)))
  env<-new.env(parent=environment(native));env$repair<-repair;env$init<-init
  env$getstatsM1<-compact_stats;environment(f)<-env
  list(fun=f,edits=counts)
}
run_gm <- function(p,sp,cfg,seed) {
  t<-nowsec();runs<-list()
  for(start in seq_len(cfg$gm_starts)) {
    rep<-repair_factory(p,sp,seed+start);fun<-adapt_ga(rep$repair)
    set.seed(seed+start)
    z<-fun$fun(Markers=p$Z[1:p$N,,drop=FALSE]-1,Markers2=p$Z[p$N+1:p$N,,drop=FALSE]-1,
      K=p$K,markereffects=p$a,nmates=sp$nf,minparents=sp$nf+sp$nm,
      impvar=cfg$gm_var,impforinbreed=cfg$gm_inbreed,npopGA=cfg$gm_pop,nitGA=cfg$gm_iter,
      nelite=cfg$gm_elite,mutprob=cfg$gm_mutation,mc.cores=1L,miniters=cfg$gm_iter,
      minitbefstop=min(80L,cfg$gm_iter-1L),tolparconv=1e-6,
      noself=TRUE,method=1L,plotiters=FALSE,plotMates=FALSE)
    ix<-encode(match(z$Mates[,1],rownames(p$M)),match(z$Mates[,2],colnames(p$M)),p$N)
    check_plan(ix,p,sp);st<-as.numeric(unlist(z$Stats))
    cost<- -(1-cfg$gm_var-cfg$gm_inbreed)*st[2]-cfg$gm_var*st[3]+cfg$gm_inbreed*st[1]
    runs[[start]]<-list(plan=ix,answer=z,cost=cost,repair=rep$audit(),ASTEdits=fun$edits)
  }
  best<-which.min(vapply(runs,`[[`,numeric(1),'cost'))
  list(plan=runs[[best]]$plan,starts=runs,RuntimeSec=nowsec()-t,PoolF=p$N,PoolM=p$N,
    NativeObjective=-runs[[best]]$cost)
}
alpha_binary <- function(root) {
  override<-Sys.getenv('ALPHAMATE_BIN')
  path<-if(nzchar(override))override else file.path(root,'vendor',if(.Platform$OS.type=='windows')'AlphaMate.exe' else 'AlphaMate_Unix')
  normalizePath(path,mustWork=TRUE)
}
run_alpha <- function(p,sp,cfg,seed,dir,root) {
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  # Work files kept in a fresh attempt directory; no stale output may be parsed.
  dir<-tempfile('attempt_',dir);dir.create(dir)
  dir<-normalizePath(dir,mustWork=TRUE)
  bin<-alpha_binary(root);if(.Platform$OS.type!='windows')Sys.chmod(bin,'0755')
  t<-nowsec()
  write.table(cbind(p$ids,ifelse(seq_along(p$ids)<=p$N,2,1)),file.path(dir,'Gender.txt'),
    quote=FALSE,row.names=FALSE,col.names=FALSE)
  alpha<-p$a+p$d*(1-2*p$qfreq)
  criterion<-drop(sweep(p$Z,2,2*p$qfreq,'-')%*%alpha)
  write.table(cbind(p$ids,criterion),file.path(dir,'Criterion.txt'),quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table(cbind(p$ids,p$K),file.path(dir,'Nrm.txt'),quote=FALSE,row.names=FALSE,col.names=FALSE)
  spec<-c('Seed'=seed,'GenderFile'='Gender.txt','NrmMatrixFile'='Nrm.txt',
    'SelCriterionFile'='Criterion.txt','NumberOfMatings'=sp$nf,'NumberOfMaleParents'=sp$nm,
    'NumberOfFemaleParents'=sp$nf,'EqualizeMaleContributions'='Yes','EqualizeFemaleContributions'='Yes',
    'LimitMaleContributionsMin'=2,'LimitMaleContributionsMax'=2,
    'LimitFemaleContributionsMin'=1,'LimitFemaleContributionsMax'=1,
    'AllowRepeatedMatings'='No','MateAllocation'='Yes',
    'PreselectMales'='No','PreselectFemales'='No',
    'TargetDegree'=cfg$alpha_degree,'TargetMinInbreedingPct'=cfg$alpha_min_inbreed_pct,
    'EvolAlgNumberOfSolutions'=cfg$alpha_solutions,'EvolAlgNumberOfIterations'=cfg$alpha_iter,
    'EvolAlgNumberOfIterationsStop'=cfg$alpha_stagnation,'NumberOfThreads'=1)
  spec<-spec[!is.na(spec)]
  spec<-c(spec,EvolAlgStopTolerance=cfg$alpha_tolerance,
    EvolAlgStopToleranceCoancestry=cfg$alpha_tolerance)
  writeLines(c(paste(names(spec),spec,sep=' , '),'Stop'),file.path(dir,'AlphaMateSpec.txt'))
  prep<-nowsec()-t;t<-nowsec();old<-getwd();setwd(dir);on.exit(setwd(old))
  status<-system2(bin,'AlphaMateSpec.txt',stdout='stdout.txt',stderr='stderr.txt')
  search<-nowsec()-t
  if(status!=0L)stop('AlphaMate exit ',status,'; see ',dir)
  file<-'MatingPlanModeOptTarget1.txt';if(!file.exists(file))stop('AlphaMate did not produce target mating plan: ',dir)
  a<-read.table(file,header=TRUE,stringsAsFactors=FALSE)
  one<-as.character(a$Parent1);two<-as.character(a$Parent2)
  flip<-startsWith(one,'M');f<-ifelse(flip,two,one);m<-ifelse(flip,one,two)
  ix<-encode(match(f,rownames(p$M)),match(m,colnames(p$M)),p$N);check_plan(ix,p,sp)
  log<-read.table('OptimisationLogModeOptTarget1.txt',header=TRUE)
  last<-tail(log,1)
  list(plan=ix,RuntimeSec=nowsec()-t+prep,InputExportSec=prep,ExecutableSec=search,
    NativeObjective=last$Objective,NativeLog=last,Directory=normalizePath(dir),PoolF=p$N,PoolM=p$N)
}
