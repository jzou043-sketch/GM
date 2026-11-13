library(ggplot2)
library(patchwork)
library(ragg)

aa<-commandArgs(FALSE)
script<-sub('^--file=','',aa[grepl('^--file=',aa)][1])
root<-dirname(dirname(normalizePath(script,winslash='/')))
source_file <- file.path(root,'results/summary/all_replicates.csv')
summary_file <- file.path(root,'results/summary/method_summary.csv')
out <- file.path(root,'figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
input_hash <- unname(tools::md5sum(source_file))
d <- read.csv(source_file,check.names=FALSE)
ref <- read.csv(summary_file,check.names=FALSE)
methods <- c('OSNB','GenomicMating','AlphaMate')
metrics <- c('ExpectedTotalGV','MeanPairProgenyVarGV','MeanExpectedF_GRM','SearchPipelineSec')
stopifnot(nrow(d)==1800L,!anyDuplicated(d[c('Scale','Replicate','Method')]),
          setequal(unique(d$Method),methods),all(table(d$Scale,d$Method)==100L),
          all(d$Profile=='inbreeding_restored'),all(d$TimingOrigin=='current_run'),
          all(d$GmInbreedingWeight==.7),all(d$AlphaDegree==30),
          all(d$AlphaMinInbreedingPct==75),all(d$SearchPipelineSec>0))
stopifnot(all(vapply(split(d$BaseInputMD5,d$Replicate),function(z)length(unique(z))==1L,logical(1))))
stopifnot(all(vapply(metrics,function(m)all(is.finite(d[[m]])),logical(1))))
d$NPerSex <- d$NFemale
d$Method <- factor(d$Method,levels=methods)
d <- d[order(d$Method,d$Replicate,d$NPerSex),]
summaries <- do.call(rbind,lapply(metrics,function(metric) {
  rows <- do.call(rbind,lapply(split(d,list(d$NPerSex,d$Method),drop=TRUE),function(z) {
    data.frame(Scale=z$Scale[1],NPerSex=z$NPerSex[1],Method=as.character(z$Method[1]),
      Metric=metric,N=nrow(z),Mean=mean(z[[metric]]),SD=sd(z[[metric]]),
      Lower95CI=mean(z[[metric]])-qt(.975,nrow(z)-1)*sd(z[[metric]])/sqrt(nrow(z)),
      Upper95CI=mean(z[[metric]])+qt(.975,nrow(z)-1)*sd(z[[metric]])/sqrt(nrow(z)))
  }))
  rows
}))
summaries$Method <- factor(summaries$Method,levels=methods)
for(i in seq_len(nrow(summaries))) {
  z <- summaries[i,]
  r <- ref[ref$Scale==z$Scale & ref$Method==as.character(z$Method) & ref$Metric==z$Metric,]
  stopifnot(nrow(r)==1L,z$N==100L,abs(z$Mean-r$Mean)<1e-10,abs(z$SD-r$SD)<1e-10,
            abs(z$Lower95CI-r$Lower95CI)<1e-8,abs(z$Upper95CI-r$Upper95CI)<1e-8)
}
colors <- c(OSNB='#D55E00',GenomicMating='#0072B2',AlphaMate='#009E73')
linetypes <- c(OSNB='solid',GenomicMating='dashed',AlphaMate='dotdash')
paper_theme <- theme_classic(base_size=10,base_family='Arial') +
  theme(axis.text=element_text(color='grey15',size=9),
    axis.title=element_text(face='plain',size=10),
    axis.title.x=element_text(margin=margin(t=6)),
    axis.title.y=element_text(margin=margin(r=7)),
    axis.line=element_line(color='grey25',linewidth=.35),
    axis.ticks=element_line(color='grey25',linewidth=.3),
    axis.ticks.length=grid::unit(1.1,'mm'),
    panel.grid.major.y=element_line(color='grey93',linewidth=.25),
    plot.tag=element_text(face='bold',size=12,color='black'),
    plot.tag.position=c(.005,1),
    legend.position='top',legend.title=element_blank(),
    legend.text=element_text(size=10),legend.key.width=grid::unit(11,'mm'),
    plot.margin=margin(10,10,6,7))

make_plot <- function(metric,tag,ylab) {
  s <- summaries[summaries$Metric==metric,]
  p <- ggplot() +
    geom_ribbon(data=s,aes(x=NPerSex,ymin=Lower95CI,ymax=Upper95CI,fill=Method,group=Method),
      alpha=.13,colour=NA,show.legend=FALSE) +
    geom_line(data=s,aes(x=NPerSex,y=Mean,color=Method,linetype=Method,group=Method),
      linewidth=.72) +
    geom_point(data=s,aes(x=NPerSex,y=Mean,color=Method),shape=16,
      size=1.85,show.legend=FALSE) +
    scale_color_manual(values=colors,breaks=methods) +
    scale_fill_manual(values=colors,breaks=methods) +
    scale_linetype_manual(values=linetypes,breaks=methods) +
    scale_x_continuous(breaks=seq(100,600,100),expand=expansion(mult=c(.035,.035))) +
    labs(tag=tag,x='Candidates per sex',y=ylab,color=NULL,fill=NULL,linetype=NULL) +
    paper_theme
  if(metric=='SearchPipelineSec') {
    stopifnot(all(s$Lower95CI>0))
    p <- p + scale_y_log10(breaks=c(.1,.5,1,5,10,50,100,500,1000,5000),
      labels=c('0.1','0.5','1','5','10','50','100','500','1000','5000'),
      expand=expansion(mult=c(.06,.08)))
  } else if(metric=='MeanExpectedF_GRM') {
    p <- p + scale_y_continuous(breaks=c(-.020,-.015,-.010,-.005,0),
      labels=function(x)sprintf('%.3f',x),expand=expansion(mult=c(.07,.07)))
  } else {
    p <- p + scale_y_continuous(expand=expansion(mult=c(.05,.07)))
  }
  # Each ribbon and line uses six means based on all 100 replicates per cell.
  built <- ggplot_build(p)$data[[1]]
  stopifnot(nrow(built)==18L,length(unique(built$group))==3L,
            all(table(built$group)==6L))
  p
}
plots <- list(
  A=make_plot('ExpectedTotalGV','A','Mean progeny genetic value'),
  B=make_plot('MeanPairProgenyVarGV','B','Mean pair genetic variance'),
  C=make_plot('MeanExpectedF_GRM','C',expression(bar(F)^GRM)),
  D=make_plot('SearchPipelineSec','D','Runtime (s)')
)
fig <- wrap_plots(G=guide_area(),A=plots$A,B=plots$B,C=plots$C,D=plots$D,
  design='GG\nAB\nCD',heights=c(.13,1,1),guides='collect') &
  theme(legend.position='top',legend.justification='center')
ggsave(file.path(out,'GM6_three_methods_main.png'),fig,width=8,height=5.6,
       units='in',dpi=220,device=ragg::agg_png,bg='white')
ggsave(file.path(out,'GM6_three_methods_main.pdf'),fig,width=8,height=5.6,
       units='in',device=grDevices::cairo_pdf,bg='white')
ggsave(file.path(out,'GM6_three_methods_main_1000dpi.tiff'),fig,width=8,height=5.6,
       units='in',dpi=1000,device=ragg::agg_tiff,compression='lzw',bg='white')
for(letter in names(plots)) {
  p <- plots[[letter]]+labs(tag=NULL)+theme(legend.position='top')
  ggsave(file.path(out,paste0('Panel_',letter,'.pdf')),p,width=6.5,height=3.823,
         units='in',device=grDevices::cairo_pdf,bg='white')
  ggsave(file.path(out,paste0('Panel_',letter,'_1000dpi.tiff')),p,width=6.5,height=3.823,
         units='in',dpi=1000,device=ragg::agg_tiff,compression='lzw',bg='white')
}
write.csv(summaries,file.path(out,'plotted_means_and_SD.csv'),row.names=FALSE)
writeLines(c(
 'Comparison of three genomic mating strategies across six candidate-population sizes.',
 'Lines and points show arithmetic means across 100 independent simulation replicates per scale and method. Shaded bands show pointwise 95% confidence intervals calculated as mean +/- t(0.975,99) times SD/sqrt(100). They are not standard deviations or simultaneous confidence bands.',
 'The x-axis is the number of candidate parents per sex (100, 200, 300, 400, 500, and 600), not generation, optimizer iteration, or replicate index. No artificial zero baseline was added and no smoothing or monotonic fitting was used.',
 '(A) Expected progeny total genetic value (additive plus dominance), averaged across the selected mating pairs; this is not genetic gain relative to a base population or a sum over mating pairs.',
 '(B) Mean per-pair progeny genetic variance, excluding the between-family component and cross-locus covariances.',
 '(C) Mean paired genomic coancestry, calculated as the average Kij/2 relative to the fixed reference population.',
 '(D) Method-specific preparation and optimization runtime in seconds. The axis is logarithmic, with tick labels in original seconds. The thick trajectory is the arithmetic mean runtime, not the geometric mean.',
 'All 1,800 method results are retained. In the current data 476/600 AlphaMate native objective values are zero; target achievement remains to be checked against raw optimization logs. This plotting operation does not change or validate those native target statuses.',
 paste('Source:',source_file),paste('Source MD5:',input_hash)
),file.path(out,'figure_caption_and_provenance.txt'))
stopifnot(identical(unname(tools::md5sum(source_file)),input_hash))
cat('Verified: all 1,800 records used; 100 replicates per point; no observations removed.\n')
cat('Mean, SD and 95% CI values reconciled against the supplied summary.\n')
cat('AlphaMate native objective zeros:',sum(d$Method=='AlphaMate' & d$NativeObjective==0),'of 600.\n')
