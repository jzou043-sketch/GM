aa<-commandArgs(FALSE);root<-dirname(normalizePath(sub('^--file=','',aa[grepl('^--file=',aa)])))
lib<-Sys.getenv('R_LIBS_USER');if(nzchar(lib)){dir.create(lib,recursive=TRUE,showWarnings=FALSE);.libPaths(c(lib,.libPaths()))}
deps<-c('AlphaSimR','Rcpp','RcppArmadillo','emoa','scatterplot3d','qtl','SOMbrero','kohonen','plotly','dplyr','magrittr','LowRankQP')
missing<-deps[!vapply(deps,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos='https://cloud.r-project.org',Ncpus=1L)
stopifnot(all(vapply(deps,requireNamespace,logical(1),quietly=TRUE)))
if(!requireNamespace('GenomicMating',quietly=TRUE))install.packages(file.path(root,'vendor','GenomicMating_2.0.tar.gz'),repos=NULL,type='source')
stopifnot(as.character(packageVersion('GenomicMating'))=='2.0')
Rcpp::cppFunction('int gm_compiler_check(){return 1;}');stopifnot(gm_compiler_check()==1)
