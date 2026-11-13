# Evaluate the ORIGINAL nested function definitions, not rewritten equivalents.
GM_FUNCTIONS <- c("solve_mating_ilp_general","choose_search_mode","build_value_space",
  "build_parent_value_space","greedy_init_in_space","random_init_in_space",
  "local_search_with_space","gm_multistart_all_runs")
load_original <- function(folder) {
  original<-new.env(parent=globalenv())
  sys.source(file.path(folder,"GM_Rcpp_single.R"),original)
  parts<-as.list(body(original$run_genomic_mating_pipeline))[-1L]
  definitions<-list()
  for(x in parts) if(is.call(x)&&identical(x[[1]],as.name("<-"))&&is.symbol(x[[2]])&&
    as.character(x[[2]]) %in% GM_FUNCTIONS) definitions[[as.character(x[[2]])]]<-x
  stopifnot(setequal(names(definitions),GM_FUNCTIONS))
  build<-function(core=original$.gm_cpp) {
    e<-new.env(parent=original); e$backend<-"rcpp"; e$search_tolerance<-1e-12; e$.gm_cpp<-core
    for(x in definitions) eval(x,e)
    e
  }
  list(core=original$.gm_cpp,build=build,definitions=definitions,source=original)
}
parent_frames <- function(problem) list(
  f=data.frame(ID=rownames(problem$M),GeneticValue=as.numeric(problem$gv[1:100])),
  m=data.frame(ID=colnames(problem$M),GeneticValue=as.numeric(problem$gv[101:200])))
original_call <- function(e,problem,seed) {
  v<-parent_frames(problem)
  e$gm_multistart_all_runs(problem$M,v$f,v$m,12L,6L,200L,200L,2L,seed,TRUE,"best",FALSE)
}
# The same original ranking/greedy helpers provide competitors with the two
# original guides. Their calculation is repeated INSIDE each method timer.
make_guides <- function(api,problem) {
  e<-api$build(); v<-parent_frames(problem)
  vs<-e$build_value_space(problem$M,200L,200L)
  ps<-e$build_parent_value_space(v$f,v$m,200L,200L)
  lapply(list(vs,ps),function(space) {
    z<-e$greedy_init_in_space(problem$M,space$female_space,space$male_space,12L,6L)
    canonical(list(f=match(z$selected_females,rownames(problem$M)),
                   m=match(z$selected_males,colnames(problem$M))))
  })
}
run_original_osnb <- function(problem,seed,api) {
  begin<-clock_sec()
  local_runs<-list(); initial_fits<-list()
  proxy<-new.env(parent=emptyenv())
  proxy$gm_assignment_cpp<-function(value,force=TRUE) {
    z<-api$core$gm_assignment_cpp(value,force)
    initial_fits[[length(initial_fits)+1L]]<<-z
    z
  }
  proxy$gm_local_cpp<-function(value,females,males,female_space,male_space,
                             force=TRUE,mode="best",tolerance=1e-12) {
    t<-clock_sec()
    z<-api$core$gm_local_cpp(value,females,males,female_space,male_space,force,mode,tolerance)
    local_runs[[length(local_runs)+1L]]<<-list(fit=z,elapsed=clock_sec()-t,
      initial=list(f=females,m=males),pool_f=female_space,pool_m=male_space)
    z
  }
  e<-api$build(proxy)
  answer<-original_call(e,problem,seed)
  elapsed<-clock_sec()-begin
  stopifnot(length(local_runs)==2L,length(initial_fits)==2L)
  starts<-lapply(1:2,function(i) {
    s<-local_runs[[i]]; z<-s$fit
    # Preserve original selected-parent order: assignment rows/columns use it.
    state<-list(f=as.integer(z$female_indices),m=as.integer(z$male_indices))
    canonical(state)
    calls<-1+z$evaluations
    trace<-data.frame(AcceptedSwaps=0:(length(z$objective_trace)-1L),
      BestObjective=as.numeric(z$objective_trace))
    list(best_state=state,initial_state=canonical(s$initial),objective=z$objective,
      assignment=z$assignment,trace=trace,
      summary=data.frame(Start=i,Source=c("value-greedy","parent-greedy")[i],
        InitialObjective=initial_fits[[i]]$objective,BestObjective=z$objective,
        Requests=calls,LowerSolves=calls,CacheHits=0,LocalSearchSec=s$elapsed,
        PoolFemale=length(s$pool_f),PoolMale=length(s$pool_m),StopReason="original_local_optimum"),
      original_local=s)
  })
  calls<-sum(vapply(starts,function(x)x$summary$Requests,numeric(1)))
  best<-which.max(vapply(starts,`[[`,numeric(1),"objective"))
  stopifnot(abs(answer$objective-starts[[best]]$objective)<1e-8,all(is.finite(calls)))
  # No evaluation cap is injected into the original algorithm. Its full work
  # defines the competitors' allowance for this replicate, rounded up if odd.
  budget<-2*ceiling(calls/2)
  list(algorithm="OSNB",parameters=list(pool=200L,n_restart=2L,mode="best"),
    starts=starts,original_answer=answer,
    summary=data.frame(Algorithm="OSNB",BestObjective=answer$objective,RuntimeSec=elapsed,
      Requests=calls,LowerSolves=calls,CacheHits=0,BudgetLimit=budget,BestStart=best,
      SetupAndInitializationSec=elapsed-sum(vapply(local_runs,`[[`,numeric(1),"elapsed")),
      PostprocessSec=0,DistinctInitialSets=length(unique(vapply(starts,
        function(s)state_key(s$initial_state),character(1))))))
}
