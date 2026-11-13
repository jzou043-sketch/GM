args0<-commandArgs(FALSE)
folder<-dirname(normalizePath(sub("^--file=","",args0[grepl("^--file=",args0)])))
source(file.path(folder,"oracle_lib.R"))
stopifnot(requireNamespace("AlphaSimR",quietly=TRUE),requireNamespace("Rcpp",quietly=TRUE))
RNGkind("Mersenne-Twister","Inversion","Rejection")
out<-Sys.getenv("GM_OUT_ROOT",file.path(folder,"results"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
factors<-c(0,1,2,4,8,16)
SCENARIOS<-paste0("D",sprintf("%02d",factors))
source_files<-c("run_large_check.R","oracle_lib.R","GM_Rcpp_single.R",
                "base_population.rds","reference/qtl_effects_by_scenario.csv")
metadata<-list(design="fixed-a-dominance-six-level-pool200-oracle-v1",multipliers=factors,
  method_seed_base=20261026,residual_seed=20261027,H2=0.70,
  n_ind=500,n_female=250,n_male=250,n_snp=10000,n_qtl=2000,
  n_selected_female=40,n_selected_male=20,n_restart=2,pool_female=200,pool_male=200,
  search_mode="best",reference_effect_scenario="S1_weak",R=as.character(getRversion()),
  packages=vapply(c("AlphaSimR","Rcpp"),function(p) as.character(packageVersion(p)),character(1)),
  md5=setNames(as.character(tools::md5sum(file.path(folder,source_files))),source_files))
meta_path<-file.path(out,"configuration","metadata.rds")
if(file.exists(meta_path)) {
  if(!identical(readRDS(meta_path),metadata)) stop("Source/config/software changed. Choose a new output directory.")
} else commit(metadata,meta_path)
snapshot<-file.path(out,"configuration","source_code")
dir.create(snapshot,recursive=TRUE,showWarnings=FALSE)
for(p in list.files(folder,pattern="[.](R|py|sh|sbatch|md)$",full.names=TRUE)) {
  target<-file.path(snapshot,basename(p))
  if(!file.exists(target)) file.copy(p,target)
}
writeLines(capture.output(sessionInfo()),file.path(out,"configuration","sessionInfo.txt"))
input_file<-file.path(out,"simulation_data","simulation_input.rds")
if(file.exists(input_file)) {
  x<-readRDS(input_file)
  stopifnot(identical(x$experiment_metadata,metadata))
} else {
  x<-readRDS(file.path(folder,"base_population.rds"))
  stopifnot(x$config$n==500,x$config$snp==10000,x$config$qtl==2000,
    all(table(x$sex)==250),ncol(x$z)==10000,nrow(x$z)==500)
  reference<-read.csv(file.path(folder,"reference","qtl_effects_by_scenario.csv"))
  reference<-subset(reference,Scenario=="S1_weak" & IsQTL)
  reference<-reference[match(colnames(x$z)[1:2000],reference$SNP),]
  a<-x$a[,"S1_weak"]; d<-x$d[,"S1_weak"]
  stopifnot(identical(unname(a),reference$AdditiveEffect),
            identical(unname(d),reference$DominanceEffect))
  x$a<-matrix(rep(a,length(factors)),nrow=2000,dimnames=list(colnames(x$z)[1:2000],SCENARIOS))
  x$d<-outer(d,factors); dimnames(x$d)<-dimnames(x$a)
  q<-seq_len(2000); zq<-x$z[,q,drop=FALSE]
  x$raw_g<-zq%*%x$a+(zq==1)%*%x$d
  sp<-AlphaSimR::SimParam$new(x$founder); sp$nThreads<-1L
  sp$importTrait(x$map$AlphaSimRSite[q],addEff=x$a,domEff=x$d,
                 intercept=colSums(x$a),name=SCENARIOS)
  pop<-AlphaSimR::newPop(x$founder,simParam=sp)
  stopifnot(max(abs(pop@gv-x$raw_g))<1e-8)
  gp<-AlphaSimR::genParam(pop,simParam=sp)
  vg<-diag(gp$varG); ve<-vg*(1-0.70)/0.70
  set.seed(metadata$residual_seed); noise<-rnorm(500)
  x$phenotype<-50+sweep(x$raw_g,2,colMeans(x$raw_g),"-")+outer(noise,sqrt(ve))
  x$config$nf<-40; x$config$nm<-20; x$config$poolf<-200; x$config$poolm<-200
  x$config$restarts<-2; x$config$seed<-metadata$method_seed_base
  x$config$source_effects<-"Fixed S1_weak a,d; only d multiplied by c(0,1,2,4,8,16)."
  x$config$experiment<-metadata$design
  x$config$md5<-metadata$md5
  x$experiment_metadata<-metadata
  x$stats<-data.frame(Scenario=SCENARIOS,DominanceMultiplier=factors,
    MeanA=colMeans(x$a),SDA=apply(x$a,2,sd),MeanD=colMeans(x$d),SDD=apply(x$d,2,sd),
    VA=diag(gp$varA),VD=diag(gp$varD),VD_over_VA=diag(gp$varD)/diag(gp$varA),
    VG=vg,CovAD=gp$covAD_L,VE=ve,ModelH2=vg/(vg+ve),
    BaseGeneticMean=colMeans(x$raw_g),PhenotypeMean=colMeans(x$phenotype),
    PhenotypeSD=apply(x$phenotype,2,sd),PhenotypeIntercept=50-colMeans(x$raw_g),
    MeanAbsA=colMeans(abs(x$a)),MeanAbsD=colMeans(abs(x$d)),
    FractionQTLAbsDGreaterThanAbsA=colMeans(abs(x$d)>abs(x$a)))
  x$stats$VarianceDescription<-ifelse(x$stats$VD_over_VA>=1,
    "dominance_variance_at_least_additive_stress_test", "additive_variance_larger")
  write_input(x,out)
  commit(x$founder,file.path(out,"simulation_data","founder_haplotypes.rds"))
  print(x$stats[,c("DominanceMultiplier","VA","VD","VD_over_VA")])
}

collect<-function() {
  all_results<-list(); comparisons<-list()
  status<-data.frame(Scenario=SCENARIOS,Complete=FALSE)
  for(i in seq_along(SCENARIOS)) {
    path<-file.path(out,"results",SCENARIOS[i],"COMPLETE.rds")
    if(!file.exists(path)) next
    ck<-readRDS(path)
    expected_cfg<-x$config; expected_cfg$seed<-metadata$method_seed_base-i
    stopifnot(identical(ck$config,expected_cfg),nrow(ck$summary)==4)
    ans<-ck$summary; all_results[[i]]<-ans; status$Complete[i]<-TRUE
    val<-setNames(ans$TotalObjective,ans$Strategy)
    gain<-setNames(ans$ExpectedGain,ans$Strategy)
    comparisons[[i]]<-data.frame(DominanceMultiplier=factors[i],VD_over_VA=x$stats$VD_over_VA[i],
      PS=unname(val["PS"]),GS=unname(val["GS"]),GS_GM=unname(val["GS+GM"]),GM=unname(val["GM"]),
      GM_minus_GS=unname(val["GM"]-val["GS"]),GM_minus_GS_GM=unname(val["GM"]-val["GS+GM"]),
      GM_vs_GS_GM_percent=if(val["GS+GM"]>1e-8) unname(100*(val["GM"]/val["GS+GM"]-1)) else NA_real_,
      Gain_PS=unname(gain["PS"]),Gain_GS=unname(gain["GS"]),
      Gain_GS_GM=unname(gain["GS+GM"]),Gain_GM=unname(gain["GM"]))
  }
  csv(status,file.path(out,"summary","completion_status.csv"))
  if(length(all_results)) {
    csv(do.call(rbind,all_results),file.path(out,"summary","all_strategy_results.csv"))
    tab<-do.call(rbind,comparisons)
    csv(tab,file.path(out,"summary","dominance_comparison.csv"))
    if(all(status$Complete)) {
      mono<-data.frame(Metric=c("GM objective","GM minus GS","GM minus GS+GM","Percent GM vs GS+GM"),
        Nondecreasing=c(all(diff(tab$GM)>=-1e-8),all(diff(tab$GM_minus_GS)>=-1e-8),
          all(diff(tab$GM_minus_GS_GM)>=-1e-8),all(diff(tab$GM_vs_GS_GM_percent)>=-1e-8)))
      csv(mono,file.path(out,"summary","monotonicity_check.csv"))
      print(tab,row.names=FALSE,digits=7); print(mono,row.names=FALSE)
    }
  }
  writeLines(c("One fixed AlphaSimR population; six prespecified dominance levels; oracle analysis.",
    paste("Completed levels:",sum(status$Complete),"of 6"),
    paste("Strategy result rows:",4*sum(status$Complete),"of 24"),
    "No BGLR, no effect estimation, no global-optimum certificate at this size.",
    "Use ExpectedGain for baseline-adjusted genetic improvement. Ratios depend on trait origin.",
    "All levels are retained, including stress-test levels with very large dominance."),
    file.path(out,"summary","README_results.txt"))
  all(status$Complete)
}
if(collect()) {
  commit(metadata,file.path(out,"ALL_COMPLETE.rds")); quit(status=0)
}
gm<-load_oracle_gm(folder); failed<-FALSE
for(i in seq_along(SCENARIOS)) {
  path<-file.path(out,"results",SCENARIOS[i],"COMPLETE.rds")
  if(file.exists(path)) { cat("SKIP",SCENARIOS[i],"\n"); next }
  cat("START",SCENARIOS[i],as.character(Sys.time()),"\n"); flush.console()
  xi<-x
  # The original run_scenario adds i: all levels thus use identical method seeds.
  xi$config$seed<-metadata$method_seed_base-i
  tryCatch({
    ans<-run_scenario(xi,i,gm,out)
    print(ans[,c("Scenario","Strategy","TotalObjective","ExpectedGain","DecisionSeconds")])
  },error=function(e) {
    failed<<-TRUE; dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    writeLines(conditionMessage(e),file.path(dirname(path),"last_error.txt"))
    message("ERROR ",SCENARIOS[i],": ",conditionMessage(e))
  })
  collect(); flush.console()
}
if(failed || !collect()) quit(status=1)
commit(metadata,file.path(out,"ALL_COMPLETE.rds"))
cat("SUCCESS: one 500-parent population, 6 levels, 4 methods, 24 rows.\n")
