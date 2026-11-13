run_tests<-function(folder,out) {
  e<-new.env(parent=globalenv())
  for(f in c("common.R","simulation.R","osnb_adapter.R","competitors.R"))sys.source(file.path(folder,f),e)
  api<-e$load_original(folder);records<-list()
  test<-function(name,code) {
    msg<-tryCatch({force(code);""},error=function(x)conditionMessage(x))
    records[[length(records)+1L]]<<-data.frame(Test=name,Status=if(nzchar(msg))"FAIL" else "PASS",Message=msg)
    cat(if(nzchar(msg))"FAIL" else "PASS",name,if(nzchar(msg))msg else "","\n");flush.console()
  }
  problem<-function(M) {dimnames(M)<-list(sprintf("F%03d",1:100),sprintf("M%03d",1:100));
    list(M=M,gv=c(rowMeans(M),colMeans(M)))}
  set.seed(2201);random_problem<-problem(matrix(rnorm(10000),100,100))
  test("six algorithms; OSNB absent from tuning grid",{
    stopifnot(identical(e$ALGORITHMS,c("OSNB","PSO","Tabu","SA","GA","DE")),
      identical(names(e$grid()),c("PSO","Tabu","SA","GA","DE")),all(lengths(e$grid())==3L))
  })
  test("original helper definitions evaluated unchanged",{
    h<-api$build()
    for(n in e$GM_FUNCTIONS)stopifnot(identical(body(h[[n]]),api$definitions[[n]][[3]][[3]]))
  })
  test("original pool cap and guide reproduction",{
    h<-api$build();v<-e$parent_frames(random_problem)
    vs<-h$build_value_space(random_problem$M,200,200)
    ps<-h$build_parent_value_space(v$f,v$m,200,200)
    stopifnot(length(vs$female_space)==100,length(ps$male_space)==100)
    g<-e$make_guides(api,random_problem)
    for(i in 1:2){sp<-list(vs,ps)[[i]];z<-h$greedy_init_in_space(random_problem$M,sp$female_space,sp$male_space,12,6)
      stopifnot(identical(g[[i]],e$canonical(list(f=match(z$selected_females,rownames(random_problem$M)),
        m=match(z$selected_males,colnames(random_problem$M))))))}
  })
  for(kind in c("random","flat","additive"))test(paste("original OSNB exact parity",kind),{
    p<-switch(kind,random=random_problem,flat=problem(matrix(1,100,100)),
      additive=problem(outer(seq_len(100),seq_len(100),"+")))
    raw<-e$original_call(api$build(),p,456L)
    got<-e$run_original_osnb(p,456L,api)
    stopifnot(identical(raw,got$original_answer),got$summary$Requests==got$summary$LowerSolves,
      got$summary$Requests<=got$summary$BudgetLimit,got$summary$CacheHits==0)
    for(s in got$starts){z<-s$original_local$fit
      stopifnot(z$evaluations==1+(z$accepted_swaps+1)*1620,
        s$summary$Requests==z$evaluations+1,all(diff(s$trace$BestObjective)>=0),
        all(rowSums(s$assignment)==1),all(colSums(s$assignment)==2))}
  })
  for(kind in c("positive","negative","mixed","ties"))test(paste("lower solver versus binary LP",kind),{
    set.seed(70);M<-switch(kind,positive=matrix(runif(72),12,6),negative=-matrix(runif(72),12,6),
      mixed=matrix(rnorm(72),12,6),ties=matrix(1,12,6))
    cons<-rbind(kronecker(matrix(1,1,6),diag(12)),kronecker(diag(6),matrix(1,1,12)))
    lp<-lpSolve::lp("max",as.vector(M),cons,rep("=",18),c(rep(1,12),rep(2,6)),all.bin=TRUE)
    z<-api$core$gm_assignment_cpp(M,TRUE)
    stopifnot(lp$status==0,abs(lp$objval-z$objective)<1e-8)
  })
  test("duplicate requests charged; unexpected errors propagate",{
    ev<-e$make_evaluator(random_problem$M,api$core,2)
    st<-list(f=1:12,m=1:6);ev$evaluate(st);ev$evaluate(st)
    stopifnot(ev$n==2L,ev$solves==1L,ev$hits==1L,
      inherits(tryCatch(ev$evaluate(st),error=identity),"gm_budget_exhausted"))
  })
  guides<-e$make_guides(api,random_problem)
  guide_value<-vapply(guides,function(g)api$core$gm_assignment_cpp(random_problem$M[g$f,g$m],TRUE)$objective,numeric(1))
  first_banks<-NULL
  for(a in names(e$grid()))for(ci in 1:3)test(paste("competitor guides, budget, plans",a,ci),{
    res<-e$run_competitor(a,random_problem,e$grid()[[a]][[ci]],990L,1000L,api)
    stopifnot(res$summary$Requests<=1000,res$summary$Requests==res$summary$LowerSolves+res$summary$CacheHits)
    for(s in 1:2){x<-res$starts[[s]];stopifnot(identical(x$initial_state,guides[[s]]),
      abs(x$summary$InitialObjective-guide_value[s])<1e-8,x$summary$Requests<=500,
      all(rowSums(x$assignment)==1),all(colSums(x$assignment)==2),
      all(diff(x$trace$BestObjective)>=0),
      abs(sum(random_problem$M[x$best_state$f,x$best_state$m]*x$assignment)-x$objective)<1e-7)}
    if(a%in%c("PSO","GA","DE")){
      banks<-lapply(res$starts,function(s)s$population_bank$keys)
      if(is.null(first_banks))first_banks<-banks else stopifnot(identical(first_banks,banks))
    }
  })
  for(a in names(e$grid()))test(paste("unexpected lower error propagates",a),{
    fake<-api;fake$core<-new.env();fake$core$gm_assignment_cpp<-function(...)stop("TEST_REAL_FAILURE")
    err<-tryCatch(e$run_competitor(a,random_problem,e$grid()[[a]][[1]],990L,1000L,fake),error=identity)
    stopifnot(inherits(err,"error"),grepl("TEST_REAL_FAILURE",conditionMessage(err)))
  })
  test("HC is rejected",{stopifnot(inherits(tryCatch(e$run_competitor("HC",random_problem,list(),1L,1000L,api),error=identity),"error"))})
  test("Tabu binary decoding preserves cardinality and feasible subsets",{
    set.seed(113L);priority<-c(sample.int(100L),sample.int(100L))
    for(i in 1:100){
      bits<-rbinom(200L,1L,runif(1L));a<-e$binary_decode(bits,priority)
      stopifnot(length(a$f)==12L,length(a$m)==6L,!anyDuplicated(a$f),!anyDuplicated(a$m))
      j<-sample.int(200L,1L);bits[j]<-1L-bits[j];b<-e$binary_decode(bits,priority)
      stopifnot(length(setdiff(a$f,b$f))+length(setdiff(a$m,b$m))<=1L)
      encoded<-integer(200L);encoded[c(a$f,100L+a$m)]<-1L
      stopifnot(identical(a,e$binary_decode(encoded,priority)))
    }
  })
  test("SA and Tabu execute exported package functions",{
    check_call<-function(package,fun,algorithm){
      counter<-new.env();counter$n<-0L
      previous<-getOption("gm_package_call_audit");options(gm_package_call_audit=counter)
      on.exit(options(gm_package_call_audit=previous),add=TRUE)
      suppressMessages(trace(fun,where=asNamespace(package),print=FALSE,
        tracer=quote({audit<-getOption("gm_package_call_audit");audit$n<-audit$n+1L})))
      on.exit(suppressMessages(untrace(fun,where=asNamespace(package))),add=TRUE)
      z<-e$run_competitor(algorithm,random_problem,e$grid()[[algorithm]][[1]],990L,1000L,api)
      stopifnot(counter$n==2L)
      for(s in z$starts)stopifnot(grepl(paste0(package,"::",fun),s$implementation,fixed=TRUE))
    }
    check_call("stats","optim","SA")
    check_call("tabuSearch","tabuSearch","Tabu")
  })
  test("AlphaSimR dimensions, original OSNB parity and checkpoint reuse",{
    scratch<-tempfile("gm_original_validation_");dir.create(scratch)
    cfg<-list(seed=20261027L,meanDD=.5,varDD=.1)
    p<-e$generate_problem(1L,"pilot",cfg,scratch)
    stopifnot(identical(dim(p$M),c(100L,100L)),identical(dim(p$Z),c(200L,1000L)),length(p$qtl_index)==200L)
    raw<-e$original_call(api$build(),p,100L);got<-e$run_original_osnb(p,100L,api)
    stopifnot(identical(raw,got$original_answer))
    file<-file.path(scratch,"simulation_data/pilot/rep_0001/problem.rds");h<-tools::md5sum(file)
    again<-e$generate_problem(1L,"pilot",cfg,scratch)
    stopifnot(identical(h,tools::md5sum(file)),identical(p$M,again$M))
  })
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  table<-do.call(rbind,records);e$csv(table,file.path(out,"test_results.csv"))
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
  if(any(table$Status!="PASS"))stop("Validation failed: ",paste(table$Test[table$Status!="PASS"],collapse="; "))
  cat("ALL",nrow(table),"CHECKS PASSED\n")
  invisible(table)
}
if(sys.nframe()==0L){aa<-commandArgs(FALSE);folder<-dirname(normalizePath(sub("^--file=","",aa[grepl("^--file=",aa)])))
  out<-file.path(Sys.getenv("OUT_ROOT",tempdir()),"validation");run_tests(folder,out)}
