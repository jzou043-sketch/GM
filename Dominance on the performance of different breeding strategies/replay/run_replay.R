args <- commandArgs(FALSE)
folder <- dirname(normalizePath(sub("^--file=","",args[grepl("^--file=",args)])))
source(file.path(folder,"oracle_lib.R"))
source(file.path(folder,"replay_common.R"))
check_package(folder)
if(!requireNamespace("Rcpp",quietly=TRUE)) stop("Rcpp is missing. Run preflight.R first.")
RNGkind("Mersenne-Twister","Inversion","Rejection")
Sys.setlocale("LC_COLLATE","C")
SCENARIOS <- mapping$Scenario
out <- Sys.getenv("GM_OUT_ROOT",file.path(folder,"server_results"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
x <- readRDS(file.path(folder,"input/frozen_population.rds"))
check_input(x)
meta <- list(version="frozen-oracle-0248-v1",R=as.character(getRversion()),
  platform=R.version$platform,Rcpp=as.character(packageVersion("Rcpp")),
  manifest=read.csv(file.path(folder,"manifest_md5.csv")),
  method_seed_base=20261026,mapping=mapping)
cfgfile <- file.path(out,"configuration/replay_config.rds")
if(file.exists(cfgfile)) {
  if(!identical(readRDS(cfgfile),meta)) stop("Output configuration differs. Use a NEW GM_OUT_ROOT.")
} else {
  if(length(list.files(out,all.files=TRUE,no..=TRUE)[!list.files(out,all.files=TRUE,no..=TRUE) %in% c(".run.lock")]))
    stop("Nonempty output folder without matching configuration. Choose a new GM_OUT_ROOT.")
  commit(meta,cfgfile)
}
writeLines(capture.output(sessionInfo()),file.path(out,"configuration/sessionInfo.txt"))
x$config$replay_metadata <- meta
write_input(x,out)
csv(mapping,file.path(out,"simulation_data/scenario_mapping.csv"))
snapshot <- file.path(out,"configuration/source_code")
dir.create(snapshot,recursive=TRUE,showWarnings=FALSE)
for(f in list.files(folder,pattern="[.](R|sh|sbatch)$",full.names=TRUE))
  file.copy(f,file.path(snapshot,basename(f)),overwrite=TRUE)
collect <- function() {
  rows <- list(); status <- data.frame(Scenario=SCENARIOS,Complete=FALSE)
  for(i in seq_along(SCENARIOS)) {
    p <- file.path(out,"results",SCENARIOS[i],"COMPLETE.rds")
    if(!file.exists(p)) next
    ck <- readRDS(p); cfg <- x$config; cfg$seed <- meta$method_seed_base-i
    stopifnot(identical(ck$config,cfg),nrow(ck$summary)==4)
    ans <- ck$summary; ans$ScenarioLabel <- mapping$ScenarioLabel[i]
    ans$DominanceMultiplier <- mapping$DominanceMultiplier[i]
    rows[[i]] <- ans; status$Complete[i] <- TRUE
  }
  csv(status,file.path(out,"summary/completion_status.csv"))
  if(length(rows)) csv(do.call(rbind,rows),file.path(out,"summary/four_scenarios_four_strategies.csv"))
  all(status$Complete)
}
if(!collect()) {
  gm <- load_oracle_gm(folder)
  for(i in seq_along(SCENARIOS)) {
    path <- file.path(out,"results",SCENARIOS[i],"COMPLETE.rds")
    if(file.exists(path)) next
    cat("START",SCENARIOS[i],mapping$ScenarioLabel[i],as.character(Sys.time()),"\n"); flush.console()
    xi <- x
    # run_scenario adds its index to each seed: preserve the same original seeds
    # when dropping D01 and D16 from the original six-scenario experiment.
    xi$config$seed <- meta$method_seed_base-i
    ans <- run_scenario(xi,i,gm,out)
    print(ans[,c("Scenario","Strategy","TotalObjective","ExpectedGain","DecisionSeconds")])
    collect(); flush.console()
  }
}
stopifnot(collect())
verify_results(folder,out)
commit(list(metadata=meta,finished=as.character(Sys.time())),file.path(out,"ALL_COMPLETE.rds"))
cat("SUCCESS: frozen population, 4 scenarios, 4 methods, 16 verified results.\n")
