run_unit_tests <- function(root,out) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE);records<-list()
  test<-function(name,expr) {
    err<-tryCatch({force(expr);''},error=function(e)conditionMessage(e))
    records[[length(records)+1L]]<<-data.frame(Test=name,Status=if(nzchar(err))'FAIL' else 'PASS',Message=err)
    cat(if(nzchar(err))'FAIL' else 'PASS',name,err,'\n');flush.console()
  }
  cfg<-settings(TRUE);sp<-cfg$scales[1,];N<-sp$N;set.seed(45)
  p<-list(N=N,Z=matrix(sample(0:2,2*N*20,TRUE),2*N,20),a=rnorm(20),d=rnorm(20),intercept=.17)
  p$ids<-c(sprintf('F%04d',seq_len(N)),sprintf('M%04d',seq_len(N)));rownames(p$Z)<-p$ids
  p$qfreq<-colMeans(p$Z)/2;p$alpha<-p$a+p$d*(1-2*p$qfreq)
  p$G<-matrix(sample(0:2,2*N*100,TRUE),2*N,100);p$observed_frequency<-colMeans(p$G)/2
  B<-sweep(p$G,2,2*p$observed_frequency,'-')
  p$K<-tcrossprod(B)/(2*sum(p$observed_frequency*(1-p$observed_frequency)))
  p$bv<-drop(sweep(p$Z,2,2*p$qfreq,'-')%*%p$alpha)
  p$M<-build_original_value_matrix(p)
  p$gv<-drop((p$Z-1)%*%p$a+(p$Z==1)%*%p$d)+p$intercept
  repair<-repair_factory(p,sp,11L)
  test('User formula moments agree with enumeration for all 9 parental genotype pairs',{
    for(x in 0:2)for(y in 0:2) {
      u<-x/2;v<-y/2;prob<-c((1-u)*(1-v),u+v-2*u*v,u*v)
      g<-c(0,.73-.41,2*.73);mu<-sum(prob*g);variance<-sum(prob*(g-mu)^2)
      z<-pair_moments(matrix(x,1,1),matrix(y,1,1),.73,-.41)
      stopifnot(abs(z$ExpectedGV-mu)<1e-12,abs(z$PairProgenyVarGV-variance)<1e-12)
    }
  })
  test('Variance is averaged per pair; NO between-family term or intercept',{
    ix<-repair$repair(sample.int(N^2,sp$nf,TRUE));z<-decode(ix,N)
    ev<-evaluate_plan(ix,p,sp,cfg)
    direct<-vapply(seq_len(nrow(z)),function(k) {
      xf<-p$Z[z$f[k],];xm<-p$Z[N+z$m[k],];s<-xf/2+xm/2;h<-s-xf*xm/2
      sum(s*p$a+h*p$d)
    },numeric(1))
    stopifnot(max(abs(ev$families$ExpectedGV-direct))<1e-10,
      abs(ev$metrics['ExpectedTotalGV']-mean(direct))<1e-10,
      abs(ev$metrics['MeanPairProgenyVarGV']-mean(ev$families$PairProgenyVarGV))<1e-12,
      abs(ev$metrics['MeanExpectedF_GRM']-mean(p$K[cbind(z$f,N+z$m)])/2)<1e-12,
      !any(grepl('Pooled|Between|GroupCoancestry',names(ev$metrics))))
  })
  test('1000 repaired GA plans obey proportional counts and exactly 2 dams per sire',{
    for(i in 1:1000)check_plan(repair$repair(sample.int(N^2,sp$nf,TRUE)),p,sp)
    ix<-repair$repair(rep(1L,sp$nf));check_plan(ix,p,sp)
    stopifnot(setequal(ix,repair$repair(ix)),identical(repair$repair(ix),repair$repair(repair$repair(ix))))
  })
  test('Native compact M1 and original full M1 return the same statistics',{
    for(i in 1:20) {
      ix<-repair$repair(sample.int(N^2,sp$nf,TRUE));z<-decode(ix,N)
      P<-matrix(0,sp$nf,2*N);P[cbind(seq_len(sp$nf),z$f)]<-.5;P[cbind(seq_len(sp$nf),N+z$m)]<-.5
      stopifnot(max(abs(compact_stats(p$Z-1,p$K,p$a,P)-
        GenomicMating::getstatsM1(p$Z-1,p$K,p$a,P)))<1e-8)
    }
  })
  test('Native GA executes; only the 3 documented adapter changes are present',{
    ad<-adapt_ga(repair$repair);stopifnot(identical(unname(ad$edits),c(1L,1L,1L)))
    old<-do.call(rbind,lapply(seq_len(N),function(i)cbind(i,N+seq_len(N))))
    new<-cbind(rep(seq_len(N),each=N),rep(N+seq_len(N),times=N))
    stopifnot(identical(unname(old),unname(new)))
    fit<-run_gm(p,sp,cfg,44L);check_plan(fit$plan,p,sp)
  })
  test('Original OSNB file hashes and all 8 nested function bodies unchanged',{
    stopifnot(unname(tools::md5sum(file.path(root,'original/GM_Rcpp_single.R')))=='be63a35bd7f872d59aec21abbc10c6cb',
      unname(tools::md5sum(file.path(root,'original/osnb_adapter.R')))=='f9f5cd63e29ee1153dec443e133f6493')
    api<-prepare_osnb(root);e<-api$build()
    for(n in names(api$definitions))stopifnot(identical(body(e[[n]]),api$definitions[[n]][[3]][[3]]))
    fit<-run_osnb_timed(p,sp,cfg,api,44L);check_plan(fit$plan,p,sp)
    vs<-e$build_value_space(p$M,200L,200L)
    ps<-e$build_parent_value_space(osnb_frames(p)$f,osnb_frames(p)$m,200L,200L)
    stopifnot(setequal(fit$CandidatePoolF,union(vs$female_space,ps$female_space)),
      setequal(fit$CandidatePoolM,union(vs$male_space,ps$male_space)),
      fit$answer$search_mode=='best')
  })
  test('Actual AlphaMate executable runs and outputs a feasible native plan',{
    fit<-run_alpha(p,sp,cfg,19L,file.path(out,'alpha_native'),root);check_plan(fit$plan,p,sp)
    stopifnot(is.finite(fit$NativeObjective))
    txt<-readLines(file.path(fit$Directory,'AlphaMateSpec.txt'))
    stopifnot('PreselectMales , No'%in%txt,'PreselectFemales , No'%in%txt)
    stopifnot('TargetMinInbreedingPct , 75'%in%txt,'TargetDegree , 30'%in%txt,
      cfg$gm_inbreed>0,cfg$gm_var>0,1-cfg$gm_inbreed-cfg$gm_var>0)
  })
  test('Six production scales x 100; 50k SNP; 200 QTL; pool200; proportional pair counts',{
    prod<-settings(FALSE)
    stopifnot(prod$nsnp==50000L,prod$nqtl==200L,prod$nrep==100L,
      identical(prod$scales$N,seq(100L,600L,100L)),all(prod$scales$pool==200L),
      all(prod$scales$nf==.12*prod$scales$N),all(prod$scales$nm==.06*prod$scales$N),
      all(prod$scales$nf==2L*prod$scales$nm))
    stopifnot(prod$profile=='inbreeding_restored',prod$gm_iter==300L,prod$gm_var==.2,
      prod$gm_inbreed==.7,prod$gm_mutation==.8,prod$alpha_degree==30,
      prod$alpha_min_inbreed_pct==75,prod$alpha_stagnation==500L,prod$alpha_tolerance==1e-5)
  })
  result<-do.call(rbind,records);csv(result,file.path(out,'unit_results.csv'))
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  if(any(result$Status!='PASS'))stop('Validation failed; inspect unit_results.csv and AlphaMate stderr')
  invisible(result)
}
