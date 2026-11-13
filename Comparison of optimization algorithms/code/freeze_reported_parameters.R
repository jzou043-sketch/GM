# Install the settings transcribed from the original server's displayed frozen
# parameter CSV into a NEW initialized run. Does not repeat pilot selection.
aa<-commandArgs(FALSE);script<-sub('^--file=','',aa[grepl('^--file=',aa)][1])
root<-dirname(dirname(normalizePath(script,winslash='/')))
out<-Sys.getenv('OUT_ROOT');if(!nzchar(out))stop('Set OUT_ROOT to a new initialized run')
cfgfile<-file.path(out,'configuration/config.rds')
if(!file.exists(cfgfile))stop('Run code/original/driver.R init first')
cfg<-readRDS(cfgfile)
selected<-list(PSO=list(pop=40L,w=1/(2*log(2))),Tabu=list(listSize=20L),
  SA=list(temp=.1,tmax=10L),GA=list(pop=40L,mutation=.2),DE=list(pop=40L,F=.8,CR=.9))
value<-list(cfg=cfg,selected=selected)
target<-file.path(out,'configuration/frozen_parameters.rds')
if(file.exists(target))stop('Frozen parameter file already exists; refusing to overwrite')
saveRDS(value,target,version=3)
write.csv(data.frame(Algorithm=names(selected),Parameters=vapply(selected,function(p)
  paste(names(p),unlist(p),sep='=',collapse=';'),character(1))),
  file.path(out,'configuration/frozen_parameters.csv'),row.names=FALSE)
writeLines('Fixed parameter reproduction: settings transcribed from server output; pilot not rerun.',
  file.path(out,'configuration/FIXED_PARAMETER_PROVENANCE.txt'))
cat('FROZEN: reported settings installed without new pilot tuning.\n')
