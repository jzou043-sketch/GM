needed<-c("AlphaSimR","Rcpp","pso","GA","DEoptim","lpSolve","tabuSearch")
if(getRversion()<"4.1.0") stop("R >= 4.1 required")
lib<-Sys.getenv("R_LIBS_USER")
if(nzchar(lib)) {dir.create(lib,recursive=TRUE,showWarnings=FALSE); .libPaths(c(lib,.libPaths()))}
missing<-needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) {
  if(Sys.getenv("INSTALL_MISSING","1")!="1") stop("Missing packages: ",paste(missing,collapse=","))
  install.packages(missing,repos="https://cloud.r-project.org",lib=if(nzchar(lib))lib else .libPaths()[1],Ncpus=1)
}
stopifnot(all(vapply(needed,requireNamespace,logical(1),quietly=TRUE)))
print(data.frame(Package=needed,Version=vapply(needed,function(p)as.character(packageVersion(p)),character(1))))
Rcpp::cppFunction("int fair_gm_compiler_check() { return 1; }")
stopifnot(fair_gm_compiler_check()==1)
