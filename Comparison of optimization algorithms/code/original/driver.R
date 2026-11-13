aa<-commandArgs(FALSE);folder<-dirname(normalizePath(sub("^--file=","",aa[grepl("^--file=",aa)])))
for(f in c("common.R","simulation.R","osnb_adapter.R","competitors.R"))source(file.path(folder,f))
args<-commandArgs(TRUE);if(!length(args))stop("Mode required")
mode<-args[1];out<-Sys.getenv("OUT_ROOT",file.path(folder,"output"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
RNGkind("Mersenne-Twister","Inversion","Rejection");invisible(Sys.setlocale("LC_COLLATE","C"))
cfg<-get_config(folder,out,create=mode=="init")
if(mode=="init"){cat("INITIALIZED",out,"\n");quit(status=0)}
frozen_path<-file.path(out,"configuration/frozen_parameters.rds")
load_complete<-function(phase,n) {
  paths<-file.path(out,"optimizer_results",phase,sprintf("rep_%04d",seq_len(n)),"COMPLETE.rds")
  csv(data.frame(Replicate=seq_len(n),Complete=file.exists(paths)),file.path(out,"summary",paste0(phase,"completion.csv")))
  if(any(!file.exists(paths)))stop("Missing ",phase," repeats: ",paste(which(!file.exists(paths)),collapse=","))
  lapply(paths,function(p){x<-readRDS(p);stopifnot(identical(x$cfg,cfg));x})
}
run_rep<-function(id,phase,api) {
  rd<-file.path(out,"optimizer_results",phase,sprintf("rep_%04d",id))
  dir.create(rd,recursive=TRUE,showWarnings=FALSE)
  t<-clock_sec();problem<-generate_problem(id,phase,cfg,out);input_elapsed<-clock_sec()-t
  inputfile<-file.path(out,"simulation_data",phase,sprintf("rep_%04d",id),"problem.rds")
  inputhash<-unname(tools::md5sum(inputfile))
  prep_path<-file.path(rd,"shared_preparation.rds")
  if(!file.exists(prep_path))save_atomic(list(SharedInputPreparationSec=input_elapsed,
    includes="generation or resume-loading, oracle matrix and input exports",inputhash=inputhash),prep_path)
  preparation<-readRDS(prep_path);stopifnot(identical(preparation$inputhash,inputhash))
  seed<-cfg$seed+id*10000L+if(phase=="pilot")1000000L else 2000000L
  execute<-function(alg,ci,params,budget=NULL) {
    path<-file.path(rd,sprintf("%s_config%02d.rds",alg,ci))
    if(file.exists(path)) {
      saved<-readRDS(path)
      stopifnot(identical(saved$cfg,cfg),identical(saved$inputhash,inputhash),
        identical(saved$result$parameters,params),identical(saved$reference_budget,budget))
    } else {
      cat("START",phase,id,alg,ci,now(),"\n");flush.console()
      result<-tryCatch(if(alg=="OSNB")run_original_osnb(problem,seed,api) else
        run_competitor(alg,problem,params,seed,budget,api),error=function(e){
          writeLines(c(now(),conditionMessage(e)),paste0(path,"_FAILED.txt"));stop(e)
        })
      saved<-list(cfg=cfg,inputhash=inputhash,reference_budget=budget,result=result,completed=now())
      save_atomic(saved,path)
    }
    pair_export(saved$result,problem,sub("[.]rds$","",path))
    saved$result
  }
  # OSNB MUST run first to measure its untruncated work, not to pick its result
  # opportunistically. The five competing methods then receive that allowance.
  os<-execute("OSNB",1L,list(pool=200L,n_restart=2L,mode="best"))
  budget<-os$summary$BudgetLimit
  if(phase=="pilot")params<-cfg$grid else {
    frozen<-readRDS(frozen_path);stopifnot(identical(frozen$cfg,cfg))
    params<-lapply(frozen$selected,function(x)list(x))
  }
  tasks<-do.call(rbind,lapply(names(params),function(a)data.frame(Algorithm=a,Config=seq_along(params[[a]]))))
  set.seed(seed+500L);tasks<-tasks[sample.int(nrow(tasks)),]
  csv(rbind(data.frame(Algorithm="OSNB",Config=1L),tasks),file.path(rd,"execution_order.csv"))
  answers<-list(cbind(Replicate=id,Config=1L,os$summary));audits<-list()
  for(k in seq_len(nrow(tasks))) {
    a<-tasks$Algorithm[k];ci<-tasks$Config[k]
    res<-execute(a,ci,params[[a]][[ci]],budget)
    for(s in 1:2) {
      same<-identical(res$starts[[s]]$initial_state,os$starts[[s]]$initial_state)
      stopifnot(same,abs(res$starts[[s]]$summary$InitialObjective-os$starts[[s]]$summary$InitialObjective)<1e-8)
      audits[[length(audits)+1L]]<-data.frame(Replicate=id,Algorithm=a,Config=ci,Start=s,
        OriginalGuideMatched=same,InitialObjective=res$starts[[s]]$summary$InitialObjective,
        PoolFemale=res$starts[[s]]$summary$PoolFemale,PoolMale=res$starts[[s]]$summary$PoolMale)
    }
    answers[[length(answers)+1L]]<-cbind(Replicate=id,Config=ci,res$summary)
  }
  rows<-do.call(rbind,answers)
  stopifnot(all(rows$BudgetLimit==budget),all(rows$Requests<=budget),
    all(rows$Requests==rows$LowerSolves+rows$CacheHits))
  rows$SharedInputPreparationSec<-preparation$SharedInputPreparationSec
  rows$SharedInputPlusMethodSec<-rows$SharedInputPreparationSec+rows$RuntimeSec
  csv(rows,file.path(rd,"replicate_summary.csv"))
  csv(do.call(rbind,audits),file.path(rd,"initialization_audit.csv"))
  complete<-file.path(rd,"COMPLETE.rds")
  if(!file.exists(complete))save_atomic(list(cfg=cfg,inputhash=inputhash,summary=rows),complete)
  cat("COMPLETE",phase,id,"matched budget",budget,"\n")
}
if(mode%in%c("pilot","test")) {
  if(length(args)!=2L)stop("Task index required")
  index<-as.integer(args[2]);stopifnot(is.finite(index),index>=1)
  if(mode=="test"&&!file.exists(frozen_path))stop("Pilot selection must finish first")
  ids<-if(mode=="pilot")index else seq.int((index-1L)*cfg$chunk+1L,min(index*cfg$chunk,cfg$nrep))
  stopifnot(all(ids>=1L),all(ids<=if(mode=="pilot")cfg$npilot else cfg$nrep))
  t<-clock_sec();invisible(lapply(PACKAGES,requireNamespace,quietly=TRUE));api<-load_original(folder)
  # Common small solver warmup is outside the method comparison and recorded.
  invisible(api$core$gm_assignment_cpp(matrix(0,12,6),TRUE))
  csv(data.frame(Stage=mode,Task=index,CompilerAndWarmupSec=clock_sec()-t,Recorded=now()),
    file.path(out,"configuration",paste0("task_setup_",mode,"_",index,".csv")))
  for(id in ids)run_rep(id,mode,api)
} else if(mode=="select") {
  rows<-do.call(rbind,lapply(load_complete("pilot",cfg$npilot),`[[`,"summary"))
  stopifnot(nrow(rows)==cfg$npilot*16L)
  best<-aggregate(BestObjective~Replicate,rows,max);names(best)[2]<-"PilotBest"
  rows<-merge(rows,best,by="Replicate");rows$Loss<-pmax(0,rows$PilotBest-rows$BestObjective)
  selected<-list();tables<-list()
  for(a in names(cfg$grid)) {
    d<-subset(rows,Algorithm==a)
    stopifnot(nrow(d)==cfg$npilot*3L,!anyDuplicated(d[c("Replicate","Config")]))
    tb<-aggregate(cbind(Loss,RuntimeSec)~Config,d,mean)
    ix<-min(tb$Config[tb$Loss<=min(tb$Loss)+1e-10]);selected[[a]]<-cfg$grid[[a]][[ix]]
    tb$Algorithm<-a;tb$Chosen<-tb$Config==ix;tables[[a]]<-tb
  }
  csv(do.call(rbind,tables),file.path(out,"summary/pilot_parameter_comparison.csv"))
  audit<-aggregate(cbind(Requests,LowerSolves,CacheHits,RuntimeSec,BudgetLimit)~Algorithm,rows,sum)
  audit$Configurations<-ifelse(audit$Algorithm=="OSNB",1L,3L)
  audit$Tuned<-audit$Algorithm!="OSNB"
  stopifnot(length(unique(audit$BudgetLimit[audit$Tuned]))==1L)
  csv(audit,file.path(out,"summary/pilot_cost_and_budget_audit.csv"))
  value<-list(cfg=cfg,selected=selected)
  if(file.exists(frozen_path))stopifnot(identical(readRDS(frozen_path),value)) else save_atomic(value,frozen_path)
  csv(data.frame(Algorithm=names(selected),Parameters=vapply(selected,function(p)
    paste(names(p),unlist(p),sep="=",collapse=";"),character(1))),file.path(out,"configuration/frozen_parameters.csv"))
  cat("FROZEN: five competitors tuned on pilot only. Original OSNB fixed.\n")
} else if(mode=="collect") {
  rows<-do.call(rbind,lapply(load_complete("test",cfg$nrep),`[[`,"summary"))
  stopifnot(nrow(rows)==cfg$nrep*6L,all(table(rows$Replicate,rows$Algorithm)==1L),
    setequal(rows$Algorithm,ALGORITHMS))
  best<-aggregate(BestObjective~Replicate,rows,max);names(best)[2]<-"BestObserved"
  rows<-merge(rows,best,by="Replicate");rows$ObjectiveLoss<-pmax(0,rows$BestObserved-rows$BestObjective)
  rows$BestObservedHit<-rows$ObjectiveLoss<=1e-8+1e-10*pmax(1,abs(rows$BestObserved))
  csv(rows,file.path(out,"summary/all_replicates_results.csv"))
  summary<-do.call(rbind,lapply(ALGORITHMS,function(a) {
    d<-subset(rows,Algorithm==a)
    data.frame(Algorithm=a,N=nrow(d),MeanBestObjective=mean(d$BestObjective),SD=sd(d$BestObjective),
      MeanRuntimeSec=mean(d$RuntimeSec),MedianRuntimeSec=median(d$RuntimeSec),
      BestObservedHitRatePercent=100*mean(d$BestObservedHit),MeanObjectiveLoss=mean(d$ObjectiveLoss),
      MeanRequests=mean(d$Requests),MeanLowerSolves=mean(d$LowerSolves),MeanCacheHits=mean(d$CacheHits),
      MeanBudgetLimit=mean(d$BudgetLimit),MeanSharedInputPlusMethodSec=mean(d$SharedInputPlusMethodSec))
  }))
  csv(summary,file.path(out,"summary/algorithm_summary.csv"))
  paired<-do.call(rbind,lapply(setdiff(ALGORITHMS,"OSNB"),function(a) {
    os<-subset(rows,Algorithm=="OSNB");d<-subset(rows,Algorithm==a)
    delta<-os$BestObjective-d$BestObjective[match(os$Replicate,d$Replicate)]
    half<-if(length(delta)>1L)qt(.975,length(delta)-1)*sd(delta)/sqrt(length(delta)) else NA_real_
    data.frame(Comparison=paste("OSNB minus",a),N=length(delta),MeanDifference=mean(delta),
      Lower95CI=mean(delta)-half,Upper95CI=mean(delta)+half)
  }))
  csv(paired,file.path(out,"summary/paired_objective_differences.csv"))
  audit<-rows[c("Replicate","Algorithm","BudgetLimit","Requests","LowerSolves","CacheHits","RuntimeSec")]
  csv(audit,file.path(out,"summary/formal_budget_audit.csv"))
  if(!file.exists(file.path(out,"ALL_COMPLETE.rds")))save_atomic(list(cfg=cfg,finished=now()),file.path(out,"ALL_COMPLETE.rds"))
  print(summary);cat("ALL COMPLETE",cfg$nrep,"formal populations; 6 algorithms.\n")
} else stop("Unknown mode")
