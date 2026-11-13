args<-commandArgs(TRUE);if(length(args)!=1L)stop("Usage: Rscript verify_outputs.R OUTPUT_DIRECTORY")
out<-normalizePath(args[1]);cfg<-readRDS(file.path(out,"configuration/config.rds"))
stopifnot(identical(readRDS(file.path(out,"ALL_COMPLETE.rds"))$cfg,cfg))
algorithms<-c("OSNB","PSO","Tabu","SA","GA","DE")
rows<-read.csv(file.path(out,"summary/all_replicates_results.csv"))
stopifnot(nrow(rows)==cfg$nrep*6L,setequal(rows$Algorithm,algorithms),
  all(table(rows$Replicate,rows$Algorithm)==1L),all(is.finite(rows$BestObjective)),
  all(rows$Requests<=rows$BudgetLimit),all(rows$Requests==rows$LowerSolves+rows$CacheHits))
count<-0L
for(phase in c("pilot","test"))for(id in seq_len(if(phase=="pilot")cfg$npilot else cfg$nrep)){
  rd<-file.path(out,"optimizer_results",phase,sprintf("rep_%04d",id))
  rs<-read.csv(file.path(rd,"replicate_summary.csv"));os<-subset(rs,Algorithm=="OSNB")
  stopifnot(nrow(os)==1L,os$CacheHits==0,os$Requests==os$LowerSolves,
    all(rs$BudgetLimit==2*ceiling(os$Requests/2)),all(rs$Requests<=rs$BudgetLimit))
  audit<-read.csv(file.path(rd,"initialization_audit.csv"))
  stopifnot(all(audit$OriginalGuideMatched),all(audit$PoolFemale==100L),all(audit$PoolMale==100L))
  files<-list.files(rd,pattern="_start[12]_pairs[.]csv$",full.names=TRUE)
  stopifnot(length(files)==if(phase=="pilot")32L else 12L)
  problem<-readRDS(file.path(out,"simulation_data",phase,sprintf("rep_%04d",id),"problem.rds"))
  for(p in files){
    d<-read.csv(p)
    start<-as.integer(sub(".*_start([12])_pairs[.]csv$","\\1",p))
    saved<-readRDS(sub("_start[12]_pairs[.]csv$",".rds",p))$result
    if(saved$algorithm=="SA")stopifnot(identical(saved$starts[[start]]$implementation,'stats::optim(method="SANN")'))
    if(saved$algorithm=="Tabu"){
      stopifnot(identical(saved$starts[[start]]$implementation,"tabuSearch::tabuSearch"),
        length(saved$starts[[start]]$binary_adapter$priority)==200L,
        saved$starts[[start]]$binary_repair_requests<=saved$starts[[start]]$summary$Requests)
    }
    values<-problem$M[cbind(match(d$Dam,rownames(problem$M)),match(d$Sire,colnames(problem$M)))]
    stopifnot(nrow(d)==12L,length(unique(d$Dam))==12L,length(unique(d$Sire))==6L,
      all(table(d$Sire)==2L),all(is.finite(values)),max(abs(values-d$ExpectedValue))<1e-8,
      abs(sum(d$ExpectedValue)-saved$starts[[start]]$objective)<1e-7)
    count<-count+1L
  }
}
tuning<-read.csv(file.path(out,"summary/pilot_cost_and_budget_audit.csv"))
stopifnot(all(tuning$Configurations[tuning$Tuned]==3L),
  length(unique(tuning$BudgetLimit[tuning$Tuned]))==1L,
  all(tuning$Requests<=tuning$BudgetLimit))
cat("PASS:",cfg$nrep,"formal and",cfg$npilot,"pilot populations; original-OSNB matched budgets;",
    count,"saved mating plans; no HC.\n")
