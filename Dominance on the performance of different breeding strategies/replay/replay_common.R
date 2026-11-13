mapping <- data.frame(Scenario=c("D00","D02","D04","D08"),
  DominanceMultiplier=c(0,2,4,8),
  ScenarioLabel=c("No dominance","Weak dominance","Moderate dominance","Strong dominance"))
check_package <- function(folder) {
  m <- read.csv(file.path(folder,"manifest_md5.csv"),stringsAsFactors=FALSE)
  got <- unname(tools::md5sum(file.path(folder,m$File)))
  if(anyNA(got) || any(tolower(got)!=tolower(m$MD5)))
    stop("Package checksum mismatch. Re-upload the complete unchanged ZIP.")
  invisible(m)
}
check_input <- function(x) {
  stopifnot(identical(colnames(x$a),mapping$Scenario),identical(colnames(x$d),mapping$Scenario),
    identical(dim(x$z),c(500L,10000L)),all(x$z %in% 0:2),
    all(table(x$sex)==250),!anyDuplicated(x$id),sum(x$map$IsQTL)==2000,
    x$config$nf==40,x$config$nm==20,x$config$poolf==200,x$config$poolm==200,
    x$config$restarts==2,x$config$search=="best",all(abs(x$stats$ModelH2-.7)<1e-10))
  g <- x$z[,1:2000,drop=FALSE]%*%x$a+(x$z[,1:2000,drop=FALSE]==1)%*%x$d
  stopifnot(max(abs(g-x$raw_g))<1e-8,all(is.finite(x$phenotype)))
  for(i in 1:4) stopifnot(identical(x$a[,i],x$a[,1]))
  stopifnot(max(abs(x$d[,3]-2*x$d[,2]))<1e-12,
    max(abs(x$d[,4]-4*x$d[,2]))<1e-12,all(x$d[,1]==0))
}
verify_results <- function(folder,out) {
  ref <- read.csv(file.path(folder,"reference/expected_results.csv"))
  ans <- read.csv(file.path(out,"summary/four_scenarios_four_strategies.csv"))
  key <- function(d) paste(d$Scenario,d$Strategy,sep="/")
  stopifnot(nrow(ans)==16,!anyDuplicated(key(ans)),setequal(key(ans),key(ref)))
  j <- match(key(ans),key(ref))
  checks <- data.frame(Scenario=ans$Scenario,Strategy=ans$Strategy,
    ReferenceObjective=ref$TotalObjective[j],RecomputedObjective=ans$TotalObjective,
    AbsoluteObjectiveDifference=abs(ans$TotalObjective-ref$TotalObjective[j]),
    AbsoluteGainDifference=abs(ans$ExpectedGain-ref$ExpectedGain[j]))
  checks$ObjectiveMatches <- checks$AbsoluteObjectiveDifference<=1e-6
  checks$GainMatches <- checks$AbsoluteGainDifference<=1e-8
  refp <- read.csv(file.path(folder,"reference/expected_pairs.csv"))
  pair_reports <- list()
  for(sc in mapping$Scenario) {
    p <- read.csv(file.path(out,"results",sc,"all_pairs.csv"))
    mat <- as.matrix(read.csv(gzfile(file.path(out,"results",sc,"true_value_matrix.csv.gz")),row.names=1,check.names=FALSE))
    for(method in c("PS","GS","GS+GM","GM")) {
      q <- p[p$Strategy==method,]; rp <- refp[refp$Scenario==sc & refp$Strategy==method,]
      stopifnot(nrow(q)==40,!anyDuplicated(q$Dam),length(unique(q$Sire))==20,
        all(table(q$Sire)==2),all(grepl("^F",q$Dam)),all(grepl("^M",q$Sire)))
      vals <- mat[cbind(match(q$Dam,rownames(mat)),match(q$Sire,colnames(mat)))]
      row <- ans[ans$Scenario==sc & ans$Strategy==method,]
      stopifnot(max(abs(vals-q$TrueExpectedProgenyValue))<1e-8,
        abs(sum(vals)-row$TotalObjective)<1e-7,
        abs(row$ExpectedGain-(row$TotalObjective/40-row$BaseGeneticMean))<1e-8)
      pair_reports[[length(pair_reports)+1L]] <- data.frame(Scenario=sc,Strategy=method,
        SameParentSets=setequal(q$Dam,rp$Dam)&&setequal(q$Sire,rp$Sire),
        SamePairs=setequal(paste(q$Dam,q$Sire),paste(rp$Dam,rp$Sire)))
    }
    gs <- subset(p,Strategy=="GS"); two <- subset(p,Strategy=="GS+GM")
    stopifnot(setequal(gs$Dam,two$Dam),setequal(gs$Sire,two$Sire))
  }
  csv(checks,file.path(out,"summary/reproduction_check.csv"))
  csv(do.call(rbind,pair_reports),file.path(out,"summary/pair_identity_check.csv"))
  ok <- all(checks$ObjectiveMatches & checks$GainMatches)
  writeLines(c(if(ok) "PASS: all 16 objectives and genetic gains reproduce the reference." else
    "FAIL: values differ from reference. Recomputed results are retained, not replaced.",
    "Objective absolute tolerance: 1e-6; expected gain absolute tolerance: 1e-8.",
    "All mating constraints and pair sums validated.",
    "Pair identities are reported separately: tied solutions may have identical objectives but different pairs.",
    "Runtime is measured anew and is not expected to match another computer."),
    file.path(out,"summary/REPRODUCTION_STATUS.txt"))
  if(!ok) stop("Reference comparison failed. See summary/reproduction_check.csv.")
  cat("REPRODUCTION PASS: all 16 results agree within declared tolerances.\n")
  invisible(checks)
}
