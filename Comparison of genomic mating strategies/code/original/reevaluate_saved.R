# Reevaluate only SAVED mating plans; does not optimize or generate offspring.
aa<-commandArgs(FALSE);root<-dirname(normalizePath(sub('^--file=','',aa[grepl('^--file=',aa)])))
for(f in c('config.R','simulation.R','methods.R','evaluation.R','reuse_data.R'))source(file.path(root,f))
args<-commandArgs(TRUE);if(length(args)!=1L)stop('Rscript reevaluate_saved.R /absolute/results')
out<-normalizePath(args[1]);frozen<-readRDS(file.path(out,'configuration.rds'));cfg<-frozen$cfg
checkfiles<-c('config.R','simulation.R','methods.R','evaluation.R')
stopifnot(identical(unname(tools::md5sum(file.path(root,checkfiles))),unname(frozen$hashes[checkfiles])))
rows<-audit<-list()
for(r in seq_len(cfg$nrep)) {
  dp<-replicate_directory(frozen,out,r)
  base<-readRDS(file.path(dp,'problem.rds'));ordering<-readRDS(file.path(dp,'nested_order.rds'))
  md5<-unname(tools::md5sum(file.path(dp,'problem.rds')))
  for(s in seq_len(nrow(cfg$scales))) {
    sp<-cfg$scales[s,];p<-nested_problem(base,ordering,sp$N)
    rd<-file.path(out,'runs',sp$Scale,sprintf('rep_%03d',r))
    for(method in c('OSNB','GenomicMating','AlphaMate')) {
      file<-file.path(rd,paste0(method,'.rds'))
      if(!file.exists(file))stop('Missing saved method: ',file)
      saved<-readRDS(file);stopifnot(identical(saved$frozen,frozen),saved$inputhash==md5,
        identical(saved$candidate_indices,p$CandidateIndices))
      ev<-evaluate_plan(saved$result$plan,p,sp,cfg)
      error<-max(abs(ev$metrics-unlist(saved$row[names(ev$metrics)])))
      stopifnot(error<1e-8)
      rows[[length(rows)+1L]]<-cbind(saved$row[c('Scale','Replicate','Method','NFemale','Pairs','SearchPipelineSec')],
        as.data.frame(as.list(ev$metrics)))
      audit[[length(audit)+1L]]<-data.frame(Scale=sp$Scale,Replicate=r,Method=method,
        Constraints='PASS',MaxMetricError=error,BaseInputMD5=md5)
    }
  }
  rm(base);gc()
}
csv(do.call(rbind,rows),file.path(out,'summary','reevaluated_metrics.csv'))
csv(do.call(rbind,audit),file.path(out,'summary','independent_formula_audit.csv'))
cat('PASS: all saved mating plans independently reevaluated with common parental formulas.\n')
