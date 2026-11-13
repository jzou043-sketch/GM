aa<-commandArgs(FALSE);script<-sub('^--file=','',aa[grepl('^--file=',aa)][1])
root<-dirname(dirname(normalizePath(script,winslash='/')))
d<-read.csv(file.path(root,'results/summary/all_replicates.csv'))
ref<-read.csv(file.path(root,'results/summary/method_summary.csv'))
metrics<-c('ExpectedTotalGV','MeanPairProgenyVarGV','MeanExpectedF_GRM','SearchPipelineSec')
stopifnot(nrow(d)==1800L,!anyDuplicated(d[c('Scale','Replicate','Method')]),
  setequal(d$Method,c('OSNB','GenomicMating','AlphaMate')),
  setequal(d$NFemale,seq(100,600,100)),all(table(d$Scale,d$Method)==100L),
  all(d$NFemale==d$NMale),all(d$SelectedFemale==.12*d$NFemale),
  all(d$SelectedMale==.06*d$NMale),all(d$Pairs==d$SelectedFemale),
  all(d$Profile=='inbreeding_restored'),all(d$TimingOrigin=='current_run'),
  all(d$GmVarianceWeight==.2),all(d$GmInbreedingWeight==.7),
  all(d$AlphaDegree==30),all(d$AlphaMinInbreedingPct==75),
  all(d$ReferencePopulationN==1200),all(d$SearchPipelineSec>0),
  all(d$MeanPairProgenyVarGV>=0),
  all(vapply(metrics,function(m)all(is.finite(d[[m]])),logical(1))),
  all(vapply(split(d$BaseInputMD5,d$Replicate),function(x)length(unique(x))==1L,logical(1))))
audit<-list()
for(sc in unique(d$Scale))for(method in unique(d$Method))for(metric in metrics){
  z<-d[d$Scale==sc&d$Method==method,];r<-ref[ref$Scale==sc&ref$Method==method&ref$Metric==metric,]
  stopifnot(nrow(z)==100L,nrow(r)==1L,r$N==100L)
  actual<-c(mean(z[[metric]]),sd(z[[metric]]),sd(z[[metric]])/10,
    mean(z[[metric]])-qt(.975,99)*sd(z[[metric]])/10,
    mean(z[[metric]])+qt(.975,99)*sd(z[[metric]])/10)
  err<-max(abs(actual-as.numeric(r[1,c('Mean','SD','SE','Lower95CI','Upper95CI')])))
  stopifnot(err<1e-8)
  audit[[length(audit)+1L]]<-data.frame(Scale=sc,Method=method,Metric=metric,
    N=100L,Mean=actual[1],SD=actual[2],MaxSummaryDifference=err)
}
out<-file.path(root,'validation');dir.create(out,showWarnings=FALSE)
write.csv(do.call(rbind,audit),file.path(out,'core_result_audit.csv'),row.names=FALSE)
nzero<-sum(d$Method=='AlphaMate'&d$NativeObjective==0)
writeLines(c('PASS: 1800 rows, 6 nested scales x 100 replicates x 3 methods.',
  'PASS: parameter/timing labels, shared input MD5, means, SD and pointwise CIs.',
  paste('AlphaMate zero native objectives:',nzero,'of 600; raw native target status remains unverified.'),
  'This audit does not reconstruct unseen genotype matrices or mating plans.',
  capture.output(sessionInfo())),file.path(out,'VALIDATION.txt'))
cat('PASS: GM strategy core result audit. AlphaMate objective zeros:',nzero,'/600.\n')
