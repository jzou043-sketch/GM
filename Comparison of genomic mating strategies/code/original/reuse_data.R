assert_data_compatible <- function(old,cfg) {
  keys<-c('seed','nsnp','nqtl','nchr','meanDD','varDD','scales')
  for(k in keys)if(!identical(old[[k]],cfg[[k]]))
    stop('Previous data are incompatible: ',k,'. Choose matching data or unset PREVIOUS_RESULTS.')
  if(isTRUE(old$smoke)!=isTRUE(cfg$smoke))stop('Cannot mix smoke and production data')
  invisible(TRUE)
}
initialize_reuse <- function(cfg,out,root) {
  previous<-Sys.getenv('PREVIOUS_RESULTS')
  if(cfg$smoke && nzchar(previous)) {
    cat('SMOKE ignores PREVIOUS_RESULTS and generates separate tiny validation data.\n')
    previous<-''
  }
  if(!nzchar(previous)) {
    if(cfg$reuse_previous_osnb)stop('REUSE_OSNB=1 requires PREVIOUS_RESULTS')
    return(list(enabled=FALSE,DataRoot=file.path(out,'simulation')))
  }
  previous<-normalizePath(previous,mustWork=TRUE)
  cf<-file.path(previous,'configuration.rds')
  if(!file.exists(cf))stop('PREVIOUS_RESULTS must point to the FULL old server result directory, not the core-results ZIP')
  old<-readRDS(cf);assert_data_compatible(old$cfg,cfg)
  if(old$cfg$nrep<cfg$nrep)stop('Previous study has too few replicates')
  for(f in c('simulation.R','evaluation.R'))
    if(unname(old$hashes[f])!=unname(tools::md5sum(file.path(root,f))))
      stop('Previous ',f,' differs; refusing to assume genotype/evaluation compatibility')
  dataroot<-normalizePath(file.path(previous,'simulation'),mustWork=TRUE)
  paths<-file.path(dataroot,sprintf('rep_%03d',seq_len(cfg$nrep)),'problem.rds')
  orderfiles<-file.path(dirname(paths),'nested_order.rds')
  if(any(!file.exists(paths)) || any(!file.exists(orderfiles)))stop('Missing old simulation/problem.rds or nested_order.rds files')
  hashes<-tools::md5sum(paths)
  csv(data.frame(Replicate=seq_len(cfg$nrep),InputPath=paths,InputMD5=unname(hashes)),
    file.path(out,'reused_inputs.csv'))
  writeLines(c('Original simulation files are read only; never modified by this experiment.',
    paste('Previous result directory:',previous),
    paste('OSNB:',if(cfg$reuse_previous_osnb)'reuse old plan and historical time' else 'rerun unchanged algorithm with current timing')),
    file.path(out,'DATA_REUSE.txt'))
  list(enabled=TRUE,PreviousRoot=previous,PreviousConfigurationMD5=unname(tools::md5sum(cf)),
    DataRoot=dataroot,InputMD5=hashes)
}
replicate_directory <- function(frozen,out,r) {
  file.path(frozen$reuse$DataRoot,sprintf('rep_%03d',r))
}
load_previous_osnb <- function(frozen,sp,r,p,cfg,inputhash,root) {
  reuse<-frozen$reuse
  if(!reuse$enabled)stop('Historical OSNB needs a previous result directory')
  path<-file.path(reuse$PreviousRoot,'runs',sp$Scale,sprintf('rep_%03d',r),'OSNB.rds')
  if(!file.exists(path))stop('Missing original OSNB checkpoint: ',path)
  saved<-readRDS(path);old<-saved$frozen
  assert_data_compatible(old$cfg,cfg)
  for(k in c('osnb_starts','osnb_mode'))stopifnot(identical(old$cfg[[k]],cfg[[k]]))
  required<-c('GM_Rcpp_single.R','osnb_adapter.R','simulation.R','evaluation.R')
  if(any(unname(old$hashes[required])!=unname(frozen$hashes[required])))
    stop('Historical OSNB source mismatch; set REUSE_OSNB=0')
  # Compare OSNB wrapper definitions, not unrelated AlphaMate code in the same
  # methods.R file. Parsing inspects expressions without executing old code.
  oldmethods<-file.path(reuse$PreviousRoot,'source','methods.R')
  if(file.exists(oldmethods)) {
    definition<-function(path,name) {
      hits<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&
        is.symbol(x[[2]])&&as.character(x[[2]])==name,as.list(parse(path)))
      if(length(hits)!=1L)stop('Cannot identify historical function: ',name)
      hits[[1]][[3]]
    }
    for(n in c('encode','decode','check_plan','prepare_osnb','osnb_frames','run_osnb'))
      if(!identical(definition(oldmethods,n),definition(file.path(root,'methods.R'),n)))
        stop('Historical OSNB wrapper differs: ',n)
  } else if(unname(old$hashes['methods.R'])!=unname(frozen$hashes['methods.R']))
    stop('Missing old source snapshot; cannot verify historical OSNB wrapper')
  stopifnot(saved$inputhash==inputhash,identical(saved$candidate_indices,p$CandidateIndices))
  check_plan(saved$result$plan,p,sp)
  list(result=saved$result,SourcePath=path)
}
