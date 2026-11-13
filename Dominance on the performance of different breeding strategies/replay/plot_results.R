args <- commandArgs(FALSE)
folder <- dirname(normalizePath(sub("^--file=","",args[grepl("^--file=",args)])))
out <- Sys.getenv("GM_OUT_ROOT",file.path(folder,"server_results"))
root <- file.path(out,"figures"); dir.create(root,recursive=TRUE,showWarnings=FALSE)
raw <- read.csv(file.path(out,"summary/four_scenarios_four_strategies.csv"))
raw <- raw[,!names(raw) %in% c("ScenarioLabel","DominanceMultiplier")]
gen <- read.csv(file.path(out,"simulation_data/genetic_parameters.csv"))
mapping <- data.frame(Scenario=c("D00","D02","D04","D08"),
  ScenarioLabel=c("No dominance","Weak dominance","Moderate dominance","Strong dominance"),
  DominanceMultiplier=c(0,2,4,8))
methods <- c("PS","GS","GS+GM","GM")
tab <- merge(mapping,raw,by="Scenario",sort=FALSE)
tab <- tab[order(match(tab$Scenario,mapping$Scenario),match(tab$Strategy,methods)),]
stopifnot(nrow(tab)==16,all(table(tab$Scenario,tab$Strategy)==1),
  all(tab$NPairs==40),all(tab$NSires==20),
  max(abs(tab$TotalObjective/40-tab$MeanExpectedProgenyValue))<1e-8,
  max(abs(tab$MeanExpectedProgenyValue-tab$BaseGeneticMean-tab$ExpectedGain))<1e-8)
gen <- merge(mapping,gen,by=c("Scenario","DominanceMultiplier"),sort=FALSE)
gen <- gen[match(mapping$Scenario,gen$Scenario),]
write.csv(tab,file.path(root,"four_scenarios_four_strategies.csv"),row.names=FALSE)
write.csv(gen,file.path(root,"four_scenarios_genetic_parameters.csv"),row.names=FALSE)
write.csv(mapping,file.path(root,"scenario_mapping.csv"),row.names=FALSE)
colors <- c("#087BB5","#FF7F0E","#169B26","#D7191C")
labels <- paste0(mapping$ScenarioLabel,"\n",mapping$DominanceMultiplier,
                 "x; VD/VA = ",sprintf("%.3f",gen$VD_over_VA))
draw <- function(metric,ylab) {
  mat <- sapply(mapping$Scenario,function(s) {
    d <- tab[tab$Scenario==s,]; d[[metric]][match(methods,d$Strategy)]
  })
  top <- max(mat)*1.20
  par(mar=c(5.4,5.4,2.5,1.0),oma=c(1.8,0,0,0),las=1,family="sans",cex=.95)
  pos <- barplot(mat,beside=TRUE,col=colors,border=NA,space=c(.18,1.05),
    ylim=c(0,top),axes=FALSE,names.arg=labels,cex.names=.86,ylab=ylab)
  ticks <- axTicks(2); abline(h=ticks,col="#dedede",lty=3,lwd=.7)
  barplot(mat,beside=TRUE,col=colors,border=NA,space=c(.18,1.05),add=TRUE,
          axes=FALSE,axisnames=FALSE)
  axis(2,at=ticks,las=1); box(bty="l")
  text(as.vector(pos),as.vector(mat)+top*.025,sprintf("%.2f",as.vector(mat)),cex=.78)
  legend("top",legend=methods,fill=colors,border=NA,horiz=TRUE,bty="n",inset=c(0,-.13),xpd=NA)
  mtext("Dominance scenario",side=1,line=3.9)
  mtext("One population; true effect evaluation. Scenario names denote relative levels within this experiment.",
    side=1,outer=TRUE,line=.3,cex=.76,col="#555555")
}
specs <- list(c("ExpectedGain","Expected genetic gain relative to base population","genetic_gain"),
              c("TotalObjective","Total expected progeny genetic value","total_expected_progeny_value"))
for(spec in specs) {
  png(file.path(root,paste0(spec[3],".png")),width=10,height=5.6,units="in",res=180)
  draw(spec[1],spec[2]); dev.off()
  tiff(file.path(root,paste0(spec[3],"_1000dpi.tiff")),width=10,height=5.6,units="in",res=1000,compression="lzw")
  draw(spec[1],spec[2]); dev.off()
}
gains <- reshape(tab[,c("ScenarioLabel","Strategy","ExpectedGain")],
  idvar="ScenarioLabel",timevar="Strategy",direction="wide")
write.csv(gains,file.path(root,"expected_genetic_gain_table.csv"),row.names=FALSE)
writeLines(c("PASS: exactly the requested 0, 2, 4 and 8 multipliers; no new simulations.",
  "PASS: all 16 source result rows retained without changing objective values.",
  "PASS: 40 pairs and 20 sires in every result.",
  "PASS: genetic gain = total objective / 40 - scenario-specific base mean.",
  "NOTE: scenario labels are reassigned for presentation, not effect changes.",
  "NOTE: one population, no replicate uncertainty or significance tests."),
  file.path(root,"validation_report.txt"))
print(gains,row.names=FALSE)
