# Portable Figure 2 rendering from the archived numerical results; base R only.
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", args[grepl("^--file=",args)][1])
root <- dirname(dirname(normalizePath(script,winslash="/",mustWork=TRUE)))
extra <- commandArgs(trailingOnly=TRUE)
input <- if(length(extra)>=1) extra[1] else
  file.path(root,"published_results/summary/four_scenarios_four_strategies.csv")
out <- if(length(extra)>=2) extra[2] else file.path(root,"figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(input,stringsAsFactors=FALSE)
scenarios <- c("D00","D02","D04","D08")
methods <- c("PS","GS","GS+GM","GM")
stopifnot(nrow(d)==16L,!anyDuplicated(paste(d$Scenario,d$Strategy)),
          all(d$NPairs==40),
          max(abs(d$ExpectedGain-(d$MeanExpectedProgenyValue-d$BaseGeneticMean)))<1e-8)
y <- sapply(scenarios,function(sc) {
  z <- d[d$Scenario==sc,]; z$ExpectedGain[match(methods,z$Strategy)]
})
stopifnot(all(is.finite(y)))
rownames(y) <- methods
cols <- c("#087BB5","#FF7F0E","#169B26","#D7191C")
draw <- function() {
  par(mar=c(3.4,4.7,2.0,.7),mgp=c(2.7,.7,0),las=1,
      family="sans",cex=.95,bty="l",xaxs="i",yaxs="i")
  pos <- barplot(y,beside=TRUE,space=c(.04,.75),plot=FALSE)
  top <- ceiling(max(y)+.9)
  plot.new();plot.window(xlim=c(0,max(pos)+.9),ylim=c(0,top),xaxs="i",yaxs="i")
  abline(h=seq(0,top,1),col="#CCCCCC",lty=3,lwd=.7)
  edges <- (pos[4,1:3]+pos[1,2:4])/2
  abline(v=edges,col="#55A5FF",lty=2,lwd=.8)
  barplot(y,beside=TRUE,space=c(.04,.75),col=cols,border=NA,add=TRUE,
          axes=FALSE,axisnames=FALSE)
  axis(2,at=seq(0,top,1));axis(1,at=colMeans(pos),tick=FALSE,
    labels=c("No dominance","Weak dominance","Moderate dominance","Strong dominance"))
  mtext("Expected genetic gain",side=2,line=2.8,las=0)
  box(bty="l")
  legend("top",legend=methods,fill=cols,border=NA,horiz=TRUE,bty="o",
         bg="white",box.col="#BBBBBB",title="Method",inset=.0,cex=.95)
}
pdf(file.path(out,"Figure2.pdf"),width=9,height=5.1,useDingbats=FALSE)
draw();dev.off()
png(file.path(out,"Figure2.png"),width=1800,height=1020,res=200)
draw();dev.off()
tiff(file.path(out,"Figure2_1000dpi.tiff"),width=9,height=5.1,units="in",
     res=1000,compression="lzw")
draw();dev.off()
dd <- d[order(match(d$Scenario,scenarios),match(d$Strategy,methods)),
  c("Scenario","ScenarioLabel","DominanceMultiplier","Strategy","NPairs",
    "TotalObjective","MeanExpectedProgenyValue","BaseGeneticMean","ExpectedGain")]
write.csv(dd,file.path(out,"Figure2_plot_data.csv"),row.names=FALSE)
writeLines(c("Figure 2. Performance of four breeding strategies under increasing dominance effects.",
  "PS, phenotype-based selection; GS, selection using true total parental genetic value;",
  "GS+GM, GS followed by optimized mate allocation; GM, joint parental selection and mate allocation.",
  "Expected genetic gain is the mean expected genetic value of progeny across 40 equally weighted",
  "mating pairs minus the mean genetic value of all 500 candidate parents in the same scenario.",
  "All scenarios use one shared base population (250 females and 250 males), the same genotype",
  "matrix and additive effects, and dominance multipliers of 0, 2, 4 and 8.",
  "Each plan uses 40 dams once each and 20 sires twice each. Bars are single-population results;",
  "they are not averages across independent replicates and are not normalized to PS."),
  file.path(out,"Figure2_caption.txt"))
cat("Figure 2 and its numerical plotting data saved in:",normalizePath(out,winslash="/"),"\n")
