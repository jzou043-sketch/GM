SCENARIOS <- c("S0_no_dominance", "S1_weak", "S2_moderate", "S3_strong")
METHODS <- c("PS", "GS", "GS+GM", "GM")
popvar <- function(x) mean((x-mean(x))^2)
csv <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  con <- if(grepl("[.]gz$", path)) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(con)); write.csv(x, con, row.names=FALSE, na="NA")
}
commit <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- tempfile("pending_", dirname(path)); saveRDS(x, tmp)
  if(file.exists(path) && .Platform$OS.type=="windows") unlink(path)
  if(!file.rename(tmp,path)) stop("Could not save checkpoint: ",path)
}
numenv <- function(name, default) {
  x <- as.numeric(Sys.getenv(name, as.character(default)))
  if(length(x)!=1 || !is.finite(x)) stop("Invalid parameter: ",name)
  x
}
config <- function(folder) {
  stopifnot(requireNamespace("AlphaSimR",quietly=TRUE),requireNamespace("Rcpp",quietly=TRUE))
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  cfg <- list(version="legacy-oracle-v1", n=numenv("N_IND",500), chr=numenv("N_CHR",10),
    snp=numenv("N_SNP",10000), qtl=numenv("N_QTL",2000), nf=numenv("N_FEMALE_SELECT",40),
    nm=numenv("N_MALE_SELECT",20), poolf=numenv("FEMALE_POOL_SIZE",200),
    poolm=numenv("MALE_POOL_SIZE",200), restarts=numenv("N_RESTART",2),
    seed=numenv("SEED_VALUE",20261026), h2=0.70, phenotypic_mean=50,
    search="best", source_effects="saved exact coefficients; no resampling or scaling")
  counts <- unlist(cfg[c("n","chr","snp","qtl","nf","nm","poolf","poolm","restarts","seed")])
  stopifnot(all(counts>=1),all(counts==floor(counts)),cfg$n%%2==0,
    cfg$snp%%cfg$chr==0,cfg$qtl%%cfg$chr==0,cfg$qtl<=2000,cfg$qtl<=cfg$snp,
    cfg$nf<=cfg$n/2,cfg$nm<=cfg$n/2,cfg$nf%%cfg$nm==0,
    cfg$poolf>=cfg$nf,cfg$poolm>=cfg$nm,cfg$poolf<=cfg$n/2,cfg$poolm<=cfg$n/2)
  files <- c("GM_Rcpp_single.R","oracle_lib.R","run_oracle.R",
    "legacy_reference/qtl_effects_by_scenario.csv","legacy_reference/simulation_design.csv",
    "legacy_reference/scenario_effect_design.csv")
  cfg$md5 <- setNames(as.character(tools::md5sum(file.path(folder,files))),files)
  cfg$R <- getRversion() |> as.character()
  cfg$packages <- vapply(c("AlphaSimR","Rcpp"),function(p) as.character(packageVersion(p)),character(1))
  cfg
}

# Extract only the original nested function definitions, not the effect-fitting
# statements. Closures share explicit backend settings. The GM file is unchanged.
load_oracle_gm <- function(folder) {
  src <- new.env(parent=globalenv()); sys.source(file.path(folder,"GM_Rcpp_single.R"),envir=src)
  env <- new.env(parent=src)
  env$.gm_cpp <- src$.gm_cpp
  env$backend <- "rcpp"; env$search_tolerance <- 1e-12
  env$maxit <- 200; env$threshold <- 1e-6
  names_loaded <- character()
  for(expr in as.list(body(src$run_genomic_mating_pipeline))[-1]) {
    if(is.call(expr) && identical(expr[[1]],as.name("<-")) &&
       is.symbol(expr[[2]]) && is.call(expr[[3]]) &&
       identical(expr[[3]][[1]],as.name("function"))) {
      eval(expr,envir=env); names_loaded <- c(names_loaded,as.character(expr[[2]]))
    }
  }
  required <- c("gm_multistart_all_runs","solve_mating_ilp_general",
                "extract_pair_result","build_value_space","build_parent_value_space")
  if(!all(required %in% names_loaded)) stop("GM source API changed; oracle adapter needs review.")
  env
}

make_input <- function(cfg, folder, out) {
  start <- proc.time()[[3]]
  founder_file <- file.path(out,"simulation_data","founder_haplotypes.rds")
  if(file.exists(founder_file)) {
    founder <- readRDS(founder_file)
  } else {
    set.seed(cfg$seed)
    founder <- AlphaSimR::runMacs(nInd=cfg$n,nChr=cfg$chr,segSites=cfg$snp/cfg$chr,
                                 inbred=FALSE,species="GENERIC",nThreads=1)
    commit(founder,founder_file)
  }
  physical_z <- AlphaSimR::pullSegSiteGeno(founder)
  stopifnot(nrow(physical_z)==cfg$n,ncol(physical_z)==cfg$snp)
  site <- colnames(physical_z)
  chr <- rep(seq_len(cfg$chr),each=cfg$snp/cfg$chr)
  maf <- pmin(colMeans(physical_z)/2,1-colMeans(physical_z)/2)
  set.seed(cfg$seed+1L)
  qloc <- sort(unlist(lapply(seq_len(cfg$chr),function(ch) {
    eligible <- which(chr==ch & maf>=0.05)
    if(length(eligible)<cfg$qtl/cfg$chr) stop("Not enough QTL with MAF >= 0.05. Saved founder retained; no redraw.")
    sample(eligible,cfg$qtl/cfg$chr)
  })))
  original <- read.csv(file.path(folder,"legacy_reference","qtl_effects_by_scenario.csv"),
                       stringsAsFactors=FALSE)
  if(anyDuplicated(original[c("Scenario","SNP")])) stop("Duplicated saved effects.")
  old_names <- sprintf("SNP%05d",seq_len(cfg$qtl))
  a <- d <- matrix(NA_real_,cfg$qtl,4,dimnames=list(old_names,SCENARIOS))
  for(s in seq_along(SCENARIOS)) {
    eff <- original[original$Scenario==SCENARIOS[s],]
    j <- match(old_names,eff$SNP)
    stopifnot(!anyNA(j),all(eff$IsQTL[j]),all(eff$AdditiveEffect[!eff$IsQTL]==0),
              all(eff$DominanceEffect[!eff$IsQTL]==0))
    a[,s] <- eff$AdditiveEffect[j]; d[,s] <- eff$DominanceEffect[j]
  }
  stopifnot(all(is.finite(a)),all(is.finite(d)),all(d[,1]==0))
  # Keep all original SNP names: QTL columns first, but retain true physical map.
  export_order <- c(qloc,setdiff(seq_len(cfg$snp),qloc))
  z <- physical_z[,export_order,drop=FALSE]
  colnames(z) <- sprintf("SNP%05d",seq_len(cfg$snp))
  id <- c(sprintf("F%03d",seq_len(cfg$n/2)),sprintf("M%03d",seq_len(cfg$n/2)))
  sex <- rep(c("F","M"),each=cfg$n/2)
  qz <- z[,seq_len(cfg$qtl),drop=FALSE]
  raw_g <- qz%*%a + (qz==1)%*%d
  sp <- AlphaSimR::SimParam$new(founder); sp$nThreads <- 1L
  # AlphaSimR uses Z-1, whereas saved GM effects use Z. Offset preserves Za+Hd.
  sp$importTrait(site[qloc],addEff=a,domEff=d,intercept=colSums(a),name=SCENARIOS)
  pop <- AlphaSimR::newPop(founder,simParam=sp)
  stopifnot(max(abs(raw_g-pop@gv))<1e-8)
  gp <- AlphaSimR::genParam(pop,simParam=sp)
  vg <- diag(gp$varG); ve <- vg*(1-cfg$h2)/cfg$h2
  set.seed(cfg$seed+2L); noise <- rnorm(cfg$n)
  centered_g <- sweep(raw_g,2,colMeans(raw_g),"-")
  phenotype <- cfg$phenotypic_mean + centered_g + outer(noise,sqrt(ve))
  genetic <- data.frame(Scenario=SCENARIOS,MeanA=colMeans(a),SDA=apply(a,2,sd),
     MeanD=colMeans(d),SDD=apply(d,2,sd),VA=diag(gp$varA),VD=diag(gp$varD),
     VD_over_VA=diag(gp$varD)/diag(gp$varA),VG=vg,CovAD=gp$covAD_L,
     VE=ve,ModelH2=vg/(vg+ve),RawGeneticMean=colMeans(raw_g),
     PhenotypeIntercept=cfg$phenotypic_mean-colMeans(raw_g),PhenotypeMean=colMeans(phenotype),
     PhenotypeSD=apply(phenotype,2,sd),SampleResidualVariance=ve*popvar(noise))
  map <- data.frame(SNP=colnames(z),IsQTL=seq_len(cfg$snp)<=cfg$qtl,
    OriginalEffectSNP=ifelse(seq_len(cfg$snp)<=cfg$qtl,colnames(z),NA_character_),
    AlphaSimRSite=site[export_order],Chromosome=chr[export_order],
    PhysicalColumnIndex=export_order,GeneticPositionMorgan=unlist(founder@genMap)[export_order],
    AlleleFrequency=colMeans(z)/2,MAF=maf[export_order])
  diagnostics <- lapply(seq_len(cfg$chr),function(ch) {
    ix <- which(chr==ch); ix <- ix[maf[ix]>0]
    left <- ix[-length(ix)]; right <- ix[-1]
    r2 <- vapply(seq_along(left),function(k) cor(physical_z[,left[k]],physical_z[,right[k]])^2,numeric(1))
    data.frame(Chromosome=ch,LeftSite=site[left],RightSite=site[right],AdjacentR2=r2,
      DistanceMorgan=diff(founder@genMap[[ch]][match(ix,which(chr==ch))]))
  })
  input <- list(config=cfg,founder=founder,z=z,id=id,sex=sex,a=a,d=d,raw_g=raw_g,
    phenotype=phenotype,map=map,stats=genetic,adjacent_ld=do.call(rbind,diagnostics),
    generation_seconds=proc.time()[[3]]-start)
  input
}

write_input <- function(x,out) {
  dir <- file.path(out,"simulation_data"); cfg <- x$config
  csv(data.frame(ID=x$id,Sex=x$sex,x$z,check.names=FALSE),file.path(dir,"genotypes_500_10000.csv.gz"))
  csv(x$map,file.path(dir,"marker_map.csv")); csv(x$stats,file.path(dir,"genetic_parameters.csv"))
  csv(x$adjacent_ld,file.path(dir,"adjacent_marker_ld.csv.gz"))
  csv(data.frame(SNP=x$map$SNP,ObservedHeterozygosity=colMeans(x$z==1),
                 AlleleFrequency=x$map$AlleleFrequency),file.path(dir,"marker_diagnostics.csv"))
  true_values <- effects <- list()
  for(s in seq_along(SCENARIOS)) {
    csv(data.frame(ID=x$id,Sex=x$sex,Phenotype=x$phenotype[,s],x$z,check.names=FALSE),
      file.path(dir,paste0(SCENARIOS[s],"_gm_dat.csv.gz")))
    true_values[[s]] <- data.frame(Scenario=SCENARIOS[s],ID=x$id,Sex=x$sex,
      TrueGeneticValue=x$raw_g[,s],Phenotype=x$phenotype[,s])
    effects[[s]] <- data.frame(Scenario=SCENARIOS[s],SNP=colnames(x$z),IsQTL=x$map$IsQTL,
      AdditiveEffect=c(x$a[,s],rep(0,cfg$snp-cfg$qtl)),
      DominanceEffect=c(x$d[,s],rep(0,cfg$snp-cfg$qtl)))
  }
  csv(do.call(rbind,true_values),file.path(dir,"phenotypes_true_values.csv"))
  csv(do.call(rbind,effects),file.path(dir,"true_effects_used.csv.gz"))
  commit(x,file.path(dir,"simulation_input.rds"))
}

run_scenario <- function(x,s,gm,out) {
  cfg <- x$config; sc <- SCENARIOS[s]
  start <- proc.time()[[3]]; f <- which(x$sex=="F"); m <- which(x$sex=="M")
  q <- seq_len(cfg$qtl)
  mat <- gm$.gm_cpp$gm_expected_cpp(x$z[f,q,drop=FALSE],x$z[m,q,drop=FALSE],x$a[,s],x$d[,s])
  dimnames(mat) <- list(x$id[f],x$id[m]); mat_seconds <- proc.time()[[3]]-start
  # Independent probability-based cross-value check for the complete matrix.
  pf <- x$z[f,q,drop=FALSE]/2; pm <- x$z[m,q,drop=FALSE]/2
  independent <- outer(drop(pf%*%(x$a[,s]+x$d[,s])),drop(pm%*%(x$a[,s]+x$d[,s])),"+") -
    2*sweep(pf,2,x$d[,s],"*")%*%t(pm)
  stopifnot(max(abs(mat-independent))<1e-8)
  pv <- data.frame(ID=x$id,Sex=x$sex,GeneticValue=x$raw_g[,s])
  top <- function(v) list(f=f[order(-v[f],x$id[f])][seq_len(cfg$nf)],
                         m=m[order(-v[m],x$id[m])][seq_len(cfg$nm)])
  random_plan <- function(sel,seed) {
    set.seed(seed); slots <- rep(sel$m,each=cfg$nf/cfg$nm)
    cbind(f=sel$f[sample.int(length(sel$f))],m=slots[sample.int(length(slots))])
  }
  plans <- list(); elapsed <- numeric(4); names(elapsed)<-METHODS
  tt<-proc.time()[[3]]; ps<-top(x$phenotype[,s]); plans$PS<-random_plan(ps,cfg$seed+100+s)
  elapsed["PS"]<-proc.time()[[3]]-tt
  tt<-proc.time()[[3]]; gs<-top(x$raw_g[,s]); plans$GS<-random_plan(gs,cfg$seed+200+s)
  elapsed["GS"]<-proc.time()[[3]]-tt
  tt<-proc.time()[[3]]
  fit2 <- gm$solve_mating_ilp_general(mat[x$id[gs$f],x$id[gs$m],drop=FALSE],TRUE)
  ij <- which(fit2$assignment==1,arr.ind=TRUE)
  plans[["GS+GM"]] <- cbind(f=gs$f[ij[,1]],m=gs$m[ij[,2]])
  elapsed["GS+GM"]<-elapsed["GS"]+proc.time()[[3]]-tt
  tt<-proc.time()[[3]]
  fit <- gm$gm_multistart_all_runs(value_mat=mat,female_value_df=pv[f,],male_value_df=pv[m,],
    n_female_select=cfg$nf,n_male_select=cfg$nm,female_pool_size=cfg$poolf,male_pool_size=cfg$poolm,
    n_restart=cfg$restarts,seed=cfg$seed+300+s,force_use_all_males=TRUE,
    search_mode="best",verbose=TRUE)
  elapsed["GM"]<-proc.time()[[3]]-tt
  gp <- gm$extract_pair_result(fit,pv)
  plans$GM <- cbind(f=match(gp$Dam,x$id),m=match(gp$Sire,x$id))
  dir <- file.path(out,"results",sc); dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  pairs <- summaries <- list()
  for(method in METHODS) {
    plan<-plans[[method]]; load<-table(plan[,2])
    stopifnot(nrow(plan)==cfg$nf,!anyDuplicated(plan[,1]),length(load)==cfg$nm,
      all(load==cfg$nf/cfg$nm),all(x$sex[plan[,1]]=="F"),all(x$sex[plan[,2]]=="M"))
    values<-mat[cbind(match(plan[,1],f),match(plan[,2],m))]
    pairs[[method]]<-data.frame(Scenario=sc,Strategy=method,Dam=x$id[plan[,1]],Sire=x$id[plan[,2]],
      TrueExpectedProgenyValue=values,DamGeneticValue=x$raw_g[plan[,1],s],
      SireGeneticValue=x$raw_g[plan[,2],s],
      PairAdvantage=values-(x$raw_g[plan[,1],s]+x$raw_g[plan[,2],s])/2)
    summaries[[method]]<-data.frame(Replicate=1,Scenario=sc,Strategy=method,NPairs=cfg$nf,
      NSires=cfg$nm,TotalObjective=sum(values),MeanExpectedProgenyValue=mean(values),
      BaseGeneticMean=mean(x$raw_g[,s]),ExpectedGain=mean(values)-mean(x$raw_g[,s]),
      ExpectedTraitMean=cfg$phenotypic_mean+mean(values)-mean(x$raw_g[,s]),
      MeanPairAdvantage=mean(pairs[[method]]$PairAdvantage),DecisionSeconds=elapsed[method],
      SharedMatrixSeconds=mat_seconds,Evaluation="oracle_true_effects")
    csv(pairs[[method]],file.path(dir,paste0(gsub("+","_plus_",method,fixed=TRUE),"_pairs.csv")))
  }
  summary<-do.call(rbind,summaries); rownames(summary)<-NULL
  stopifnot(setequal(plans$GS[,1],plans[["GS+GM"]][,1]),
            setequal(plans$GS[,2],plans[["GS+GM"]][,2]),
            summary$TotalObjective[3]>=summary$TotalObjective[2]-1e-7)
  if(s==1) stopifnot(abs(summary$TotalObjective[2]-summary$TotalObjective[3])<1e-7)
  ps_total<-summary$TotalObjective[1]
  valid_scale<-is.finite(ps_total)&&ps_total>1e-8
  summary$PS1000Display<-if(valid_scale) 1000*summary$TotalObjective/ps_total else NA_real_
  summary$ImprovementVsPSPercent<-if(valid_scale) 100*(summary$TotalObjective/ps_total-1) else NA_real_
  summary$DisplayScaleStatus<-if(valid_scale) "valid_positive_PS" else "undefined_nonpositive_PS"
  vs<-gm$build_value_space(mat,cfg$poolf,cfg$poolm)
  ps_space<-gm$build_parent_value_space(pv[f,],pv[m,],cfg$poolf,cfg$poolm)
  pool<-list(value=vs,parent=ps_space,
    union_female=union(vs$female_space,ps_space$female_space),
    union_male=union(vs$male_space,ps_space$male_space))
  summary$GMActualFemaleSearchPool<-length(pool$union_female)
  summary$GMActualMaleSearchPool<-length(pool$union_male)
  csv(summary,file.path(dir,"summary.csv")); csv(do.call(rbind,pairs),file.path(dir,"all_pairs.csv"))
  csv(fit$history,file.path(dir,"GM_restart_history.csv"))
  csv(data.frame(Dam=rownames(mat),mat,check.names=FALSE),file.path(dir,"true_value_matrix.csv.gz"))
  commit(list(fit=fit,pools=pool),file.path(dir,"GM_details.rds"))
  commit(list(config=cfg,summary=summary,completed_at=as.character(Sys.time())),file.path(dir,"COMPLETE.rds"))
  summary
}

summarize_results <- function(out,cfg) {
  tables<-list(); status<-data.frame(Scenario=SCENARIOS,Complete=FALSE)
  for(s in seq_along(SCENARIOS)) {
    path<-file.path(out,"results",SCENARIOS[s],"COMPLETE.rds")
    if(file.exists(path)) {
      r<-readRDS(path); stopifnot(identical(cfg,r$config),nrow(r$summary)==4)
      tables[[s]]<-r$summary; status$Complete[s]<-TRUE
    }
  }
  csv(status,file.path(out,"summary","completion_status.csv"))
  if(length(tables)) csv(do.call(rbind,tables),file.path(out,"summary","four_scenarios_four_strategies.csv"))
  writeLines(c("Oracle simulation; one independent population; no BGLR or effect estimation.",
    paste("Complete scenarios:",sum(status$Complete),"of 4"),
    paste("Strategy result rows:",4*sum(status$Complete),"of 16"),
    "Use raw TotalObjective/ExpectedGain as primary results. PS1000Display is an optional within-scenario ratio.",
    "Do not report single-population results as replicate means, confidence intervals or significance tests."),
    file.path(out,"summary","README_results.txt"))
  all(status$Complete)
}
