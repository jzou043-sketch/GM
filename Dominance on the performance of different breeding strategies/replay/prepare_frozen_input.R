# Packaging utility only. Server users run submit.sh, not this script.
args <- commandArgs(TRUE)
if(length(args)!=2) stop("Usage: prepare_frozen_input.R source_results package_directory")
src <- normalizePath(args[1]); folder <- normalizePath(args[2])
x <- readRDS(file.path(src,"simulation_data/simulation_input.rds"))
codes <- c("D00","D02","D04","D08")
i <- match(codes,colnames(x$a)); stopifnot(!anyNA(i))
y <- x[c("config","z","id","sex","map","adjacent_ld","generation_seconds")]
for(n in c("a","d","raw_g","phenotype")) y[[n]] <- x[[n]][,i,drop=FALSE]
y$stats <- x$stats[match(codes,x$stats$Scenario),]
y$source_experiment_metadata <- x$experiment_metadata
y$scenario_mapping <- data.frame(Scenario=codes,DominanceMultiplier=c(0,2,4,8),
  ScenarioLabel=c("No dominance","Weak dominance","Moderate dominance","Strong dominance"))
stopifnot(nrow(y$z)==500,ncol(y$z)==10000,all(table(y$sex)==250),
  y$config$nf==40,y$config$nm==20,y$config$poolf==200,y$config$poolm==200,y$config$restarts==2)
for(k in 1:4) {
  stopifnot(identical(y$a[,k],x$a[,1]),
    max(abs(y$d[,k]-x$d[,2]*y$scenario_mapping$DominanceMultiplier[k]))<1e-12)
}
# Plain R matrices and lists are sufficient for replay; no AlphaSimR S4 object is needed.
saveRDS(y,file.path(folder,"input/frozen_population.rds"),version=3)
write.csv(y$scenario_mapping,file.path(folder,"input/scenario_mapping.csv"),row.names=FALSE)
r <- read.csv(file.path(src,"summary/all_strategy_results.csv"))
r <- r[r$Scenario %in% codes,]
stopifnot(nrow(r)==16)
write.csv(r,file.path(folder,"reference/expected_results.csv"),row.names=FALSE)
p <- do.call(rbind,lapply(codes,function(s) read.csv(file.path(src,"results",s,"all_pairs.csv"))))
write.csv(p,file.path(folder,"reference/expected_pairs.csv"),row.names=FALSE)
write.csv(y$stats,file.path(folder,"reference/genetic_parameters.csv"),row.names=FALSE)
cat("Packaged the exact saved input for four scenarios; no SNP or phenotype regenerated.\n")
