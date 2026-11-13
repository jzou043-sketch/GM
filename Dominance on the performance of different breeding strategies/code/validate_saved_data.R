# Independent audit of the archived Figure 2 data; base R only.
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)][1])
root <- dirname(dirname(normalizePath(script, winslash="/", mustWork=TRUE)))
out <- file.path(root, "validation")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
x <- readRDS(file.path(root, "data/simulation_input.rds"))
frozen <- readRDS(file.path(root, "replay/input/frozen_population.rds"))
for (nm in c("z", "id", "sex", "a", "d", "raw_g", "phenotype", "map")) {
  stopifnot(identical(x[[nm]], frozen[[nm]]))
}
stopifnot(identical(dim(x$z), c(500L,10000L)), all(x$z %in% 0:2),
          all(table(x$sex)==250), !anyDuplicated(x$id), sum(x$map$IsQTL)==2000)
exported <- read.csv(gzfile(file.path(root,"data/genotypes_500_10000.csv.gz")),
                     check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(identical(exported$ID,x$id),identical(exported$Sex,x$sex),
  identical(names(exported)[-(1:2)],colnames(x$z)),
  all(as.matrix(exported[,-(1:2),drop=FALSE])==x$z))
ph <- read.csv(file.path(root,"data/phenotypes_true_values.csv"),stringsAsFactors=FALSE)
stopifnot(nrow(ph)==2000L,!anyDuplicated(paste(ph$Scenario,ph$ID)))
for (k in seq_len(4)) {
  pp <- ph[ph$Scenario==colnames(x$a)[k],];jj <- match(x$id,pp$ID)
  stopifnot(!anyNA(jj),identical(pp$Sex[jj],x$sex),
    max(abs(pp$TrueGeneticValue[jj]-x$raw_g[,k]))<1e-8,
    max(abs(pp$Phenotype[jj]-x$phenotype[,k]))<1e-8)
}
q <- which(x$map$IsQTL)
stopifnot(identical(q,1:2000), all(x$map$MAF[q]>=.05),
          all(table(x$map$Chromosome[q])==200L))
raw <- x$z[,q,drop=FALSE] %*% x$a + (x$z[,q,drop=FALSE]==1) %*% x$d
stopifnot(max(abs(raw-x$raw_g))<1e-8)
for (k in 1:4) stopifnot(identical(x$a[,k],x$a[,1]))
stopifnot(all(x$d[,1]==0), max(abs(x$d[,3]-2*x$d[,2]))<1e-12,
          max(abs(x$d[,4]-4*x$d[,2]))<1e-12)
residual <- sweep(x$phenotype,2,50-colMeans(x$raw_g),"-")-x$raw_g
scaled <- sweep(residual,2,sqrt(x$stats$VE),"/")
stopifnot(max(abs(scaled-scaled[,1]))<1e-10,
          max(abs(x$stats$VG/(x$stats$VG+x$stats$VE)-.7))<1e-10)
arch <- file.path(root,"published_results")
sumtab <- read.csv(file.path(arch,"summary/four_scenarios_four_strategies.csv"))
key <- paste(sumtab$Scenario,sumtab$Strategy)
stopifnot(nrow(sumtab)==16L,!anyDuplicated(key))
checks <- list()
for (sc in c("D00","D02","D04","D08")) {
  k <- match(sc,colnames(x$a))
  pairs <- read.csv(file.path(arch,"results",sc,"all_pairs.csv"))
  mat <- as.matrix(read.csv(gzfile(file.path(arch,"results",sc,"true_value_matrix.csv.gz")),
                            row.names=1,check.names=FALSE))
  # Recompute the entire 250 x 250 matrix independently from parental genotypes.
  fi <- which(x$sex=="F"); mi <- which(x$sex=="M")
  uf <- x$z[fi,q,drop=FALSE]/2; vm <- x$z[mi,q,drop=FALSE]/2
  ad <- x$a[,k]+x$d[,k]
  independent <- outer(drop(uf%*%ad),drop(vm%*%ad),"+") -
    2*sweep(uf,2,x$d[,k],"*")%*%t(vm)
  ord <- independent[match(rownames(mat),x$id[fi]),match(colnames(mat),x$id[mi])]
  matrix_error <- max(abs(ord-mat))
  stopifnot(is.finite(matrix_error),matrix_error<1e-8)
  for (method in c("PS","GS","GS+GM","GM")) {
    pp <- pairs[pairs$Strategy==method,]
    f <- match(pp$Dam,x$id); m <- match(pp$Sire,x$id)
    stopifnot(nrow(pp)==40L,!anyNA(f),!anyNA(m),!anyDuplicated(f),
              all(x$sex[f]=="F"),all(x$sex[m]=="M"),
              length(unique(m))==20L,all(table(m)==2L))
    u <- x$z[f,q,drop=FALSE]/2; v <- x$z[m,q,drop=FALSE]/2
    mu <- drop((u+v)%*%x$a[,k]+(u+v-2*u*v)%*%x$d[,k])
    rr <- sumtab[sumtab$Scenario==sc & sumtab$Strategy==method,]
    stopifnot(nrow(rr)==1L,rr$NPairs==40L,rr$NSires==20L)
    err <- max(abs(mu-pp$TrueExpectedProgenyValue),
      abs(sum(mu)-rr$TotalObjective), abs(mean(mu)-rr$MeanExpectedProgenyValue),
      abs(mean(x$raw_g[,k])-rr$BaseGeneticMean),
      abs(mean(mu)-mean(x$raw_g[,k])-rr$ExpectedGain),
      abs(mean(mu)-mean(x$raw_g[,k])+50-rr$ExpectedTraitMean))
    stopifnot(is.finite(err),err<1e-8)
    checks[[length(checks)+1L]] <- data.frame(Scenario=sc,Strategy=method,
      Pairs=nrow(pp),Constraints="PASS",IndependentMeansAndGains="PASS",
      MaxPairAndSummaryError=err,MaxFullMatrixError=matrix_error)
  }
  gs <- pairs[pairs$Strategy=="GS",]; tw <- pairs[pairs$Strategy=="GS+GM",]
  stopifnot(setequal(gs$Dam,tw$Dam),setequal(gs$Sire,tw$Sire))
  top <- function(v,sex,n) x$id[which(x$sex==sex)[
    order(-v[x$sex==sex],x$id[x$sex==sex])][seq_len(n)]]
  ps <- pairs[pairs$Strategy=="PS",]
  stopifnot(setequal(ps$Dam,top(x$phenotype[,k],"F",40)),
            setequal(ps$Sire,top(x$phenotype[,k],"M",20)),
            setequal(gs$Dam,top(x$raw_g[,k],"F",40)),
            setequal(gs$Sire,top(x$raw_g[,k],"M",20)))
}
# The original replay sources and input files must retain their original checksums.
manifest <- read.csv(file.path(root,"replay/manifest_md5.csv"))
got <- unname(tools::md5sum(file.path(root,"replay",manifest$File)))
stopifnot(!anyNA(got),all(tolower(got)==tolower(manifest$MD5)))
report <- do.call(rbind,checks)
write.csv(report,file.path(out,"independent_formula_checks.csv"),row.names=FALSE)
writeLines(c("PASS: frozen input matches the archived server input.",
  "PASS: exported genotype and phenotype CSV files match the actual RDS input.",
  "PASS: population dimensions, QTL counts, dominance scaling and shared residual draw.",
  "PASS: all four complete 250 x 250 cross-value matrices independently reconstructed.",
  "PASS: 16 strategy/scenario results and all 640 selected pair rows independently verified.",
  "PASS: 40 distinct dams, 20 sires, two matings per sire, and PS/GS truncation selections.",
  "PASS: original replay manifest MD5 checksums.",
  "Absolute tolerance for independently computed genetic values and gains: 1e-8.",
  "This audit evaluates saved data; it does not rerun the optimizer.",
  capture.output(sessionInfo())),file.path(out,"SAVED_DATA_VALIDATION.txt"))
print(report,row.names=FALSE)
cat("PASS: all saved data, formulas, mating constraints, and original checksums.\n")
