make_plots <- function(out) {
  d<-read.csv(file.path(out,'summary','method_summary.csv'))
  metrics<-c('ExpectedTotalGV','MeanPairProgenyVarGV','MeanExpectedF_GRM','SearchPipelineSec')
  labels<-c('Expected progeny genetic score (no intercept)',
    'Mean pair progeny genetic variance','Mean paired genomic coancestry (Kij/2)',
    'Search pipeline time (seconds)')
  methods<-c('OSNB','GenomicMating','AlphaMate');cols<-c('#0072B2','#D55E00','#009E73')
  dir.create(file.path(out,'figures'),showWarnings=FALSE)
  for(i in seq_along(metrics)) {
    metric<-metrics[i];z<-d[d$Metric==metric,];logarithmic<-metric=='SearchPipelineSec'
    yy<-c(z$Mean,z$Lower95CI,z$Upper95CI);yy<-yy[is.finite(yy)]
    if(logarithmic)yy<-yy[yy>0]
    if(!length(yy))next
    limits<-range(yy);if(diff(limits)==0)limits<-limits+c(-1,1)*max(abs(limits[1])*.1,1e-6)
    if(logarithmic)limits[1]<-max(limits[1],1e-5)
    for(ext in c('png','pdf')) {
      path<-file.path(out,'figures',paste0(metric,'.',ext))
      if(ext=='png')png(path,width=1500,height=1000,res=150) else pdf(path,width=10,height=6.7)
      plot(range(z$NPerSex)+c(-5,5),limits,type='n',log=if(logarithmic)'y' else '',
        xlab='Candidates per sex',ylab=labels[i],main='Six scales: mean and 95% CI')
      for(j in seq_along(methods)) {
        a<-z[z$Method==methods[j],];a<-a[order(a$NPerSex),]
        lines(a$NPerSex,a$Mean,type='b',col=cols[j],pch=15+j,lwd=2)
        ok<-is.finite(a$Lower95CI)&is.finite(a$Upper95CI)
        if(logarithmic)ok<-ok&a$Lower95CI>0
        if(any(ok))arrows(a$NPerSex[ok],a$Lower95CI[ok],a$NPerSex[ok],a$Upper95CI[ok],
          angle=90,code=3,length=.035,col=cols[j])
      }
      legend('topleft',legend=methods,col=cols,lty=1,pch=16:18,bty='n');dev.off()
    }
  }
  invisible(NULL)
}
if(sys.nframe()==0L) {
  args<-commandArgs(TRUE);if(length(args)!=1L)stop('Rscript plot_results.R /absolute/results')
  make_plots(normalizePath(args[1]))
}
