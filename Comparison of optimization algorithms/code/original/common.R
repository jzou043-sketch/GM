ALGORITHMS <- c("OSNB", "PSO", "Tabu", "SA", "GA", "DE")
PACKAGES <- c("AlphaSimR", "Rcpp", "pso", "GA", "DEoptim", "lpSolve", "stats", "tabuSearch")
now <- function() as.character(Sys.time())
clock_sec <- function() unname(proc.time()[3])
intenv <- function(name, default) {
  x <- as.numeric(Sys.getenv(name, default))
  if(length(x)!=1 || !is.finite(x) || x<1 || x!=floor(x) || x>.Machine$integer.max) stop("Invalid ",name)
  as.integer(x)
}
save_atomic <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  if(file.exists(path)) stop("Refusing to overwrite committed checkpoint: ",path)
  tmp <- tempfile("pending_", dirname(path)); saveRDS(x,tmp,version=3)
  if(!file.rename(tmp,path)) stop("Checkpoint rename failed: ",path)
}
csv <- function(x, path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  con <- if(grepl("[.]gz$",path)) gzfile(path,"wt") else file(path,"wt")
  on.exit(close(con)); write.csv(x,con,row.names=FALSE)
}
grid <- function() list(
  PSO=lapply(c(.4,1/(2*log(2)),.9),function(w) list(pop=40L,w=w)),
  Tabu=lapply(c(5L,10L,20L),function(t) list(listSize=t)),
  SA=lapply(c(.1,1,10),function(t) list(temp=t,tmax=10L)),
  GA=lapply(c(.05,.1,.2),function(m) list(pop=40L,mutation=m)),
  DE=list(list(pop=40L,F=.5,CR=.5),list(pop=40L,F=.8,CR=.5),list(pop=40L,F=.8,CR=.9)))
get_config <- function(folder,out,create=FALSE) {
  files <- c("common.R","osnb_adapter.R","competitors.R","simulation.R","driver.R","GM_Rcpp_single.R")
  hashes <- tools::md5sum(file.path(folder,files)); names(hashes)<-files
  versions<-vapply(PACKAGES,function(p) as.character(packageVersion(p)),character(1))
  if(anyNA(hashes)) stop("Missing source")
  file<-file.path(out,"configuration/config.rds")
  if(file.exists(file)) {
    cfg<-readRDS(file)
    if(!identical(cfg$hashes,hashes)||!identical(cfg$versions,versions)||
       cfg$R!=as.character(getRversion())) stop("Source/software changed. Use a new output directory.")
    return(cfg)
  }
  if(!create) stop("Run setup first")
  cfg<-list(schema="original-osnb-all-package-competitors-v2", nrep=intenv("N_REP",1000),
    npilot=intenv("N_PILOT",30),chunk=intenv("REPS_PER_TASK",10),seed=intenv("SEED_VALUE",20261027),
    meanDD=.5,varDD=.1,n_restart=2L,pool=200L,mode="best",grid=grid(),
    hashes=hashes,versions=versions,R=as.character(getRversion()))
  if(cfg$seed+max(cfg$nrep,cfg$npilot)*10000+2000002>.Machine$integer.max) stop("Seed overflow")
  save_atomic(cfg,file)
  dir.create(file.path(out,"configuration/source"),recursive=TRUE,showWarnings=FALSE)
  file.copy(list.files(folder,full.names=TRUE,pattern="[.](R|sh|sbatch|txt)$"),
            file.path(out,"configuration/source"),overwrite=FALSE)
  csv(data.frame(Algorithm=rep(names(cfg$grid),lengths(cfg$grid)),
    Config=unlist(lapply(cfg$grid,seq_along)),
    Parameters=unlist(lapply(cfg$grid,function(g) vapply(g,function(p)
      paste(names(p),unlist(p),sep="=",collapse=";"),character(1))))),
    file.path(out,"configuration/tuning_grid.csv"))
  writeLines(capture.output(sessionInfo()),file.path(out,"configuration/sessionInfo.txt"))
  cfg
}
canonical <- function(st) {
  valid<-function(x,n) is.numeric(x)&&length(x)==n&&all(is.finite(x))&&
    all(x==floor(x))&&all(x>=1 & x<=100)&&!anyDuplicated(x)
  if(!is.list(st)||!valid(st$f,12L)||!valid(st$m,6L)) stop("Invalid selected parent set")
  list(f=sort(as.integer(st$f)),m=sort(as.integer(st$m)))
}
state_key <- function(st) paste(paste(st$f,collapse=","),paste(st$m,collapse=","),sep="|")
decode_keys <- function(x) {
  if(length(x)!=200L||any(!is.finite(x))) stop("Invalid random-key vector")
  canonical(list(f=order(-x[1:100],1:100)[1:12],m=order(-x[101:200],1:100)[1:6]))
}
pair_export <- function(result,problem,path) {
  for(i in seq_along(result$starts)) {
    s<-result$starts[[i]]; ij<-which(s$assignment==1,arr.ind=TRUE)
    stopifnot(nrow(ij)==12L,all(rowSums(s$assignment)==1L),all(colSums(s$assignment)==2L))
    pairs<-data.frame(Dam=rownames(problem$M)[s$best_state$f[ij[,1]]],
      Sire=colnames(problem$M)[s$best_state$m[ij[,2]]],
      ExpectedValue=problem$M[cbind(s$best_state$f[ij[,1]],s$best_state$m[ij[,2]])])
    stopifnot(abs(sum(pairs$ExpectedValue)-s$objective)<1e-7)
    csv(pairs,paste0(path,"_start",i,"_pairs.csv"))
    csv(s$trace,paste0(path,"_start",i,"_trace.csv.gz"))
  }
  csv(do.call(rbind,lapply(result$starts,`[[`,"summary")),paste0(path,"_starts.csv"))
  implementation<-do.call(rbind,lapply(seq_along(result$starts),function(i){
    s<-result$starts[[i]]
    data.frame(Start=i,Implementation=if(is.null(s$implementation))"original GM Rcpp" else s$implementation,
      BinaryCardinalityRepairRequests=if(is.null(s$binary_repair_requests))0L else s$binary_repair_requests)
  }))
  csv(implementation,paste0(path,"_implementation.csv"))
}
