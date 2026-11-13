# Run the saved-data audit, the unchanged optimizer replay, and Figure 2 plotting.
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", args[grepl("^--file=",args)][1])
root <- dirname(dirname(normalizePath(script,winslash="/",mustWork=TRUE)))
rscript <- file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript")
run <- function(rel,extra=character()) {
  cat("\nRUN:",rel,"\n");flush.console()
  status <- system2(rscript,shQuote(c(file.path(root,rel),extra)))
  if(status!=0L) stop("Failed: ",rel," (exit ",status,")")
}
run("code/validate_saved_data.R")
run("replay/preflight.R")
out <- Sys.getenv("GM_OUT_ROOT",file.path(root,"recomputed_results/four_scenarios"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,winslash="/",mustWork=TRUE)
Sys.setenv(GM_OUT_ROOT=out)
run("replay/run_replay.R")
run("code/plot_figure2.R",c(file.path(out,"summary/four_scenarios_four_strategies.csv"),
                            file.path(out,"figures")))
cat("\nSUCCESS: audited data, independently rerun strategies, and replotted Figure 2.\n",
    "Recomputed output:",out,"\n")
