aa<-commandArgs(FALSE);script<-sub('^--file=','',aa[grepl('^--file=',aa)][1])
root<-dirname(dirname(normalizePath(script,winslash='/')))
d<-read.csv(file.path(root,'results/all_replicates_results.csv'))
s<-read.csv(file.path(root,'results/algorithm_summary.csv'))
alg<-c('OSNB','PSO','Tabu','SA','GA','DE')
stopifnot(nrow(d)==6000L,!anyDuplicated(d[c('Replicate','Algorithm')]),
  setequal(d$Algorithm,alg),identical(sort(unique(d$Replicate)),1:1000),
  all(table(d$Replicate,d$Algorithm)==1L),all(is.finite(d$BestObjective)),
  all(is.finite(d$RuntimeSec)),all(d$RuntimeSec>0),all(d$Requests<=d$BudgetLimit),
  all(d$Requests==d$LowerSolves+d$CacheHits))
best<-ave(d$BestObjective,d$Replicate,FUN=max)
loss<-pmax(0,best-d$BestObjective)
stopifnot(max(abs(best-d$BestObserved))<1e-8,max(abs(loss-d$ObjectiveLoss))<1e-8,
  identical(loss<=1e-8+1e-10*pmax(1,abs(best)),d$BestObservedHit),
  all(vapply(split(d$BudgetLimit,d$Replicate),function(x)length(unique(x))==1L,logical(1))))
audit<-do.call(rbind,lapply(alg,function(a){
  z<-d[d$Algorithm==a,];r<-s[s$Algorithm==a,]
  actual<-c(mean(z$BestObjective),sd(z$BestObjective),mean(z$RuntimeSec),
    median(z$RuntimeSec),100*mean(z$BestObservedHit),mean(z$ObjectiveLoss),
    mean(z$Requests),mean(z$LowerSolves),mean(z$CacheHits),mean(z$BudgetLimit),
    mean(z$SharedInputPlusMethodSec))
  fields<-c('MeanBestObjective','SD','MeanRuntimeSec','MedianRuntimeSec',
    'BestObservedHitRatePercent','MeanObjectiveLoss','MeanRequests','MeanLowerSolves',
    'MeanCacheHits','MeanBudgetLimit','MeanSharedInputPlusMethodSec')
  stopifnot(nrow(r)==1L,r$N==1000L)
  err<-max(abs(actual-as.numeric(r[1,fields])))
  stopifnot(err<1e-8)
  data.frame(Algorithm=a,N=nrow(z),MeanTotalExpectedProgenyValue=mean(z$BestObjective),
    SDTotalExpectedProgenyValue=sd(z$BestObjective),MeanRuntimeSec=mean(z$RuntimeSec),
    SDRuntimeSec=sd(z$RuntimeSec),MaxSummaryDifference=err)
}))
out<-file.path(root,'validation');dir.create(out,showWarnings=FALSE)
write.csv(audit,file.path(out,'core_result_audit.csv'),row.names=FALSE)
writeLines(c('PASS: 6000 rows, 1000 independent replicate IDs, 6 algorithms, no duplicate keys.',
  'PASS: summary means/SD, observed-best values/losses/hits and request accounting.',
  'This validates core tables; it does not verify unseen raw genotypes or mating plans.',
  capture.output(sessionInfo())),file.path(out,'VALIDATION.txt'))
print(audit,row.names=FALSE);cat('PASS: optimizer core result audit.\n')
