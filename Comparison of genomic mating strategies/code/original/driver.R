aa<-commandArgs(FALSE);root<-dirname(normalizePath(sub('^--file=','',aa[grepl('^--file=',aa)])))
for(f in c('config.R','simulation.R','methods.R','evaluation.R','reuse_data.R'))source(file.path(root,f))
args<-commandArgs(TRUE);if(!length(args))stop('Stage required');stage<-args[1]
out<-Sys.getenv('OUT_ROOT');if(!nzchar(out))stop('OUT_ROOT required')
dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out)
cfgfile<-file.path(out,'configuration.rds')
pkgs<-c('AlphaSimR','GenomicMating','Rcpp','RcppArmadillo')
versions<-vapply(pkgs,function(x)as.character(packageVersion(x)),character(1))
stopifnot(versions[['GenomicMating']]=='2.0')
files<-c(file.path(root,c('config.R','simulation.R','methods.R','evaluation.R','driver.R',
  'validate.R','plot_results.R','reevaluate_saved.R','reuse_data.R','setup.R','worker.sbatch','submit_all.sh')),
  file.path(root,'original',c('GM_Rcpp_single.R','osnb_adapter.R')),
  file.path(root,'vendor',c('MatingPrograms.R','MatingPrograms.cpp')),alpha_binary(root))
hashes<-tools::md5sum(files);names(hashes)<-basename(files)
native<-new.env(parent=globalenv());sys.source(file.path(root,'vendor','MatingPrograms.R'),native)
stopifnot(identical(body(native$getGaSolutions),body(GenomicMating::getGaSolutions)),
  unname(tools::md5sum(file.path(root,'original/GM_Rcpp_single.R')))=='be63a35bd7f872d59aec21abbc10c6cb',
  unname(tools::md5sum(file.path(root,'original/osnb_adapter.R')))=='f9f5cd63e29ee1153dec443e133f6493')
if(stage=='init') {
  if(file.exists(cfgfile))stop('Configuration already exists; choose a new output directory')
  cfg<-settings(Sys.getenv('SMOKE','0')=='1')
  reuse<-initialize_reuse(cfg,out,root)
  plan<-expand.grid(Replicate=seq_len(cfg$nrep),ScaleIndex=seq_len(nrow(cfg$scales)))
  plan$Task<-seq_len(nrow(plan));plan$Scale<-cfg$scales$Scale[plan$ScaleIndex]
  csv(plan,file.path(out,'task_map.csv'));csv(cfg$scales,file.path(out,'scales.csv'))
  writeLines(capture.output(dput(cfg)),file.path(out,'configuration.txt'))
  dir.create(file.path(out,'source'),showWarnings=FALSE)
  file.copy(files,file.path(out,'source'),overwrite=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  source(file.path(root,'validate.R'));run_unit_tests(root,file.path(out,'validation'))
  atomic(list(cfg=cfg,versions=versions,hashes=hashes,R=as.character(getRversion()),tasks=plan,
    reuse=reuse),cfgfile)
  cat('INIT SUCCESS; profile=',cfg$profile,'; tasks=',nrow(plan),'\n');quit(status=0)
}
if(!file.exists(cfgfile))stop('Run init first')
frozen<-readRDS(cfgfile)
stopifnot(identical(frozen$hashes,hashes),identical(frozen$versions,versions),
  frozen$R==as.character(getRversion()))
cfg<-frozen$cfg;tasks<-frozen$tasks
if(stage=='prepare') {
  stopifnot(length(args)==2L);r<-as.integer(args[2]);stopifnot(r%in%seq_len(cfg$nrep))
  dp<-replicate_directory(frozen,out,r)
  if(frozen$reuse$enabled) {
    stopifnot(file.exists(file.path(dp,'problem.rds')),file.exists(file.path(dp,'nested_order.rds')))
  } else prepare_replicate(cfg,r,dp)
  cat('BASE READY',r,'\n');quit(status=0)
}
if(stage=='run') {
  stopifnot(length(args)==2L);id<-as.integer(args[2]);stopifnot(id%in%tasks$Task)
  task<-tasks[tasks$Task==id,];s<-task$ScaleIndex;r<-task$Replicate;sp<-cfg$scales[s,]
  rd<-file.path(out,'runs',sp$Scale,sprintf('rep_%03d',r));dir.create(rd,recursive=TRUE,showWarnings=FALSE)
  if(file.exists(file.path(rd,'COMPLETE.rds'))){cat('Already complete\n');quit(status=0)}
  dp<-replicate_directory(frozen,out,r)
  if(!file.exists(file.path(dp,'nested_order.rds')))stop('Run prepare for this replicate first')
  t<-nowsec();base<-readRDS(file.path(dp,'problem.rds'))
  ordering<-readRDS(file.path(dp,'nested_order.rds'))
  assert_data_compatible(base$cfg,cfg);p<-nested_problem(base,ordering,sp$N)
  rm(base);gc();loadtime<-nowsec()-t
  inputhash<-unname(tools::md5sum(file.path(dp,'problem.rds')))
  if(frozen$reuse$enabled)stopifnot(inputhash==unname(frozen$reuse$InputMD5[r]))
  csv(data.frame(LocalIndex=seq_along(p$ids),ID=p$ids,Sex=rep(c('F','M'),each=p$N),
    MaximumPopulationIndex=p$CandidateIndices),file.path(rd,'candidate_subset.csv'))
  api<-prepare_osnb(root)
  invisible(api$core$gm_assignment_cpp(matrix(0,sp$nf,sp$nm),TRUE))
  invisible(GenomicMating::getstatsM1(matrix(0,2,1),diag(2),1,matrix(.5,1,2)))
  stream_seed(seed_for(cfg,s,r,1L));order<-sample(c('OSNB','GenomicMating','AlphaMate'))
  csv(data.frame(Order=seq_along(order),Method=order),file.path(rd,'execution_order.csv'))
  rows<-list()
  for(method in order) {
    path<-file.path(rd,paste0(method,'.rds'));seed<-seed_for(cfg,s,r,100L)
    if(file.exists(path)) {
      saved<-readRDS(path);stopifnot(identical(saved$frozen,frozen),saved$inputhash==inputhash)
    } else {
      cat('START',sp$Scale,r,method,format(Sys.time()),'\n');flush.console()
      old_osnb<-if(method=='OSNB' && cfg$reuse_previous_osnb)
        load_previous_osnb(frozen,sp,r,p,cfg,inputhash,root) else NULL
      result<-if(!is.null(old_osnb))old_osnb$result else tryCatch(switch(method,OSNB=run_osnb_timed(p,sp,cfg,api,seed),
        GenomicMating=run_gm(p,sp,cfg,seed+100L),
        AlphaMate=run_alpha(p,sp,cfg,seed+200L,file.path(rd,'alphamate'),root)),
        error=function(e){writeLines(conditionMessage(e),file.path(rd,paste0(method,'_FAILED.txt')));stop(e)})
      check_plan(result$plan,p,sp)
      eval<-evaluate_plan(result$plan,p,sp,cfg)
      stopifnot(all(is.finite(eval$metrics)),is.finite(result$RuntimeSec),result$RuntimeSec>=0)
      row<-cbind(data.frame(Scale=sp$Scale,Replicate=r,Method=method,Profile=cfg$profile,
        NFemale=p$N,NMale=p$N,SelectedFemale=sp$nf,SelectedMale=sp$nm,Pairs=sp$nf,
        PoolParameter=if(method=='OSNB')sp$pool else NA,
        ActualSearchPoolF=result$PoolF,ActualSearchPoolM=result$PoolM,
        SearchPipelineSec=result$RuntimeSec,EvaluationSec=eval$EvaluationSec,
        TimingOrigin=if(!is.null(old_osnb))'reused_previous_run' else 'current_run',
        DataOrigin=if(frozen$reuse$enabled)'previous_exact_genotypes' else 'generated_current_run',
        GmVarianceWeight=cfg$gm_var,GmInbreedingWeight=cfg$gm_inbreed,
        AlphaDegree=cfg$alpha_degree,AlphaMinInbreedingPct=cfg$alpha_min_inbreed_pct,
        DataLoadSec=loadtime,NativeObjective=result$NativeObjective,
        ReferencePopulationN=p$ReferencePopulationN,BaseInputMD5=inputhash,
        GVCoding='Z012_Whet_no_intercept',VarianceDefinition='mean_pair_variance_no_between_term'),
        as.data.frame(as.list(eval$metrics)))
      saved<-list(frozen=frozen,inputhash=inputhash,candidate_indices=p$CandidateIndices,
        result=result,evaluation=eval,row=row)
      atomic(saved,path)
    }
    csv(saved$evaluation$families,file.path(rd,paste0(method,'_mating_plan_and_metrics.csv')))
    rows[[method]]<-saved$row
    cat('DONE',sp$Scale,r,method,'\n');flush.console()
  }
  stopifnot(unname(tools::md5sum(file.path(dp,'problem.rds')))==inputhash)
  rows<-do.call(rbind,rows);csv(rows,file.path(rd,'metrics.csv'))
  atomic(list(frozen=frozen,inputhash=inputhash,rows=rows),file.path(rd,'COMPLETE.rds'))
} else if(stage=='collect') {
  paths<-file.path(out,'runs',tasks$Scale,sprintf('rep_%03d',tasks$Replicate),'COMPLETE.rds')
  ok<-file.exists(paths);csv(data.frame(tasks,Complete=ok),file.path(out,'summary','completion.csv'))
  if(any(!ok))stop('Missing task IDs: ',paste(tasks$Task[!ok],collapse=','))
  all<-do.call(rbind,lapply(paths,function(p){x<-readRDS(p);stopifnot(identical(x$frozen,frozen));x$rows}))
  stopifnot(nrow(all)==nrow(tasks)*3,!anyDuplicated(all[c('Scale','Replicate','Method')]))
  csv(all,file.path(out,'summary','all_replicates.csv'))
  metrics<-c('ExpectedTotalGV','MeanPairProgenyVarGV','MeanExpectedF_GRM','SearchPipelineSec')
  summary<-paired<-list()
  for(sc in cfg$scales$Scale)for(m in c('OSNB','GenomicMating','AlphaMate')) {
    d<-all[all$Scale==sc&all$Method==m,];N<-cfg$scales$N[match(sc,cfg$scales$Scale)]
    for(v in metrics) {
      n<-nrow(d);se<-if(n>1)sd(d[[v]])/sqrt(n) else NA_real_
      h<-if(n>1)qt(.975,n-1)*se else NA_real_
      summary[[length(summary)+1L]]<-data.frame(Scale=sc,NPerSex=N,Method=m,Metric=v,N=n,
        Mean=mean(d[[v]]),SD=sd(d[[v]]),SE=se,Lower95CI=mean(d[[v]])-h,
        Upper95CI=mean(d[[v]])+h,Median=median(d[[v]]),Q25=quantile(d[[v]],.25),Q75=quantile(d[[v]],.75))
    }
    if(m!='OSNB') {
      os<-all[all$Scale==sc&all$Method=='OSNB',];os<-os[match(d$Replicate,os$Replicate),]
      stopifnot(identical(os$Replicate,d$Replicate))
      for(v in metrics) {
        delta<-os[[v]]-d[[v]];h<-if(length(delta)>1)qt(.975,length(delta)-1)*sd(delta)/sqrt(length(delta)) else NA_real_
        paired[[length(paired)+1L]]<-data.frame(Scale=sc,NPerSex=N,
          Comparison=paste('OSNB minus',m),Metric=v,N=length(delta),
          Difference=mean(delta),Lower95CI=mean(delta)-h,Upper95CI=mean(delta)+h)
      }
    }
  }
  csv(do.call(rbind,summary),file.path(out,'summary','method_summary.csv'))
  csv(do.call(rbind,paired),file.path(out,'summary','paired_differences.csv'))
  writeLines(if(cfg$reuse_previous_osnb)
    'OSNB timings are historical; compare genetic metrics directly but do not describe timings as a contemporaneous rerun.' else
    'All three methods timed anew in the current run; same-input fixed-configuration comparison.',
    file.path(out,'summary','TIMING_SOURCE.txt'))
  source(file.path(root,'plot_results.R'));make_plots(out)
  if(!file.exists(file.path(out,'ALL_COMPLETE.rds')))atomic(list(Finished=Sys.time(),frozen=frozen),file.path(out,'ALL_COMPLETE.rds'))
  cat('ALL COMPLETE:',nrow(tasks),'scale/replicate tasks;',nrow(all),'method results\n')
} else stop('Unknown stage')
