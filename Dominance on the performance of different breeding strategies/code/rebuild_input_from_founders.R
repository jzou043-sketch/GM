# Reconstruct the reported simulated inputs from preserved phased founder
# haplotypes, the selected-locus map, and the actual coefficient table.
# This preserves the original MaCS realization rather than drawing new founders.
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", args[grepl("^--file=",args)][1])
root <- dirname(dirname(normalizePath(script,winslash="/",mustWork=TRUE)))
if(!requireNamespace("AlphaSimR",quietly=TRUE)) stop("Install AlphaSimR (original version: 2.1.0).")
RNGkind("Mersenne-Twister","Inversion","Rejection")
founder <- readRDS(file.path(root,"replay/provenance/founder_haplotypes.rds"))
map <- read.csv(file.path(root,"data/marker_map.csv"),stringsAsFactors=FALSE)
eff <- read.csv(gzfile(file.path(root,"data/true_effects_used.csv.gz")),stringsAsFactors=FALSE)
ref <- readRDS(file.path(root,"data/simulation_input.rds"))
physical <- AlphaSimR::pullSegSiteGeno(founder)
stopifnot(identical(dim(physical),c(500L,10000L)),length(founder@genMap)==10L)
z <- physical[,map$PhysicalColumnIndex,drop=FALSE]
colnames(z) <- map$SNP
q <- which(map$IsQTL)
stopifnot(identical(q,1:2000))
# Independently repeat the original QTL selection on the preserved founders.
chr <- rep(1:10,each=1000)
maf <- pmin(colMeans(physical)/2,1-colMeans(physical)/2)
set.seed(20261027)
selected <- sort(unlist(lapply(1:10,function(ch)
  sample(which(chr==ch & maf>=.05),200))))
stopifnot(identical(as.integer(map$PhysicalColumnIndex[q]),as.integer(selected)))
scenarios <- c("D00","D02","D04","D08")
a <- d <- matrix(NA_real_,2000,4,dimnames=list(map$SNP[q],scenarios))
for(k in seq_along(scenarios)) {
  ee <- eff[eff$Scenario==scenarios[k],]
  jj <- match(map$SNP[q],ee$SNP)
  stopifnot(!anyNA(jj),!anyDuplicated(ee$SNP),all(ee$IsQTL[jj]))
  a[,k] <- ee$AdditiveEffect[jj]; d[,k] <- ee$DominanceEffect[jj]
  stopifnot(all(ee$AdditiveEffect[!ee$IsQTL]==0),all(ee$DominanceEffect[!ee$IsQTL]==0))
}
g <- z[,q,drop=FALSE]%*%a+(z[,q,drop=FALSE]==1)%*%d
sp <- AlphaSimR::SimParam$new(founder); sp$nThreads <- 1L
# AlphaSimR uses a centered additive genotype code. Its intercept offsets this
# to match the archived uncentered 0/1/2 additive plus heterozygosity model.
sp$importTrait(map$AlphaSimRSite[q],addEff=a,domEff=d,
               intercept=colSums(a),name=scenarios)
pop <- AlphaSimR::newPop(founder,simParam=sp)
stopifnot(max(abs(pop@gv-g))<1e-8)
gp <- AlphaSimR::genParam(pop,simParam=sp)
vg <- diag(gp$varG); ve <- vg*(1-.70)/.70
set.seed(20261027); noise <- rnorm(500)
phenotype <- 50+sweep(g,2,colMeans(g),"-")+outer(noise,sqrt(ve))
stopifnot(identical(unname(z),unname(ref$z)),
  max(abs(a-ref$a))<1e-12,max(abs(d-ref$d))<1e-12,
  max(abs(g-ref$raw_g))<1e-8,max(abs(vg-ref$stats$VG))<1e-8,
  max(abs(ve-ref$stats$VE))<1e-8,max(abs(phenotype-ref$phenotype))<1e-8)
out <- file.path(root,"recomputed_results/rebuilt_data")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
saveRDS(list(z=z,id=ref$id,sex=ref$sex,map=map,a=a,d=d,raw_g=g,
             phenotype=phenotype,VG=vg,VE=ve),file.path(out,"rebuilt_input.rds"))
writeLines(c("PASS: original phased founders -> original QTL choice and genotype matrix.",
  "PASS: exported effects -> genetic values, AlphaSimR variance and phenotypes.",
  "This is reconstruction of the saved simulation realization, not a new replicate.",
  capture.output(sessionInfo())),file.path(root,"validation/FOUNDER_RECONSTRUCTION.txt"))
cat("PASS: complete simulated input reconstructed from preserved phased founders.\n")
