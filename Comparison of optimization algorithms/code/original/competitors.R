make_evaluator <- function(M,core,limit) {
  ev<-new.env(parent=emptyenv()); ev$n<-ev$solves<-ev$hits<-0L
  ev$best<- -Inf; ev$state<-NULL; ev$trace<-list(); ev$first<-NA_real_
  cache<-new.env(hash=TRUE,parent=emptyenv()); start<-clock_sec()
  ev$remaining<-function()limit-ev$n
  ev$evaluate<-function(st) {
    if(ev$n>=limit) stop(structure(list(message="Request budget reached",call=NULL),
      class=c("gm_budget_exhausted","error","condition")))
    st<-canonical(st); key<-state_key(st); ev$n<-ev$n+1L
    if(exists(key,cache,inherits=FALSE)) {value<-get(key,cache,inherits=FALSE);ev$hits<-ev$hits+1L}
    else {
      value<-core$gm_assignment_cpp(M[st$f,st$m,drop=FALSE],TRUE)$objective
      if(length(value)!=1L||!is.finite(value))stop("Invalid lower-solver result")
      ev$solves<-ev$solves+1L; assign(key,value,cache)
    }
    if(ev$n==1L)ev$first<-value
    improved<-is.infinite(ev$best)||value>ev$best+1e-12*max(1,abs(ev$best))
    if(improved){ev$best<-value;ev$state<-st}
    if(improved||ev$n%%100==0||ev$n==limit)ev$trace[[length(ev$trace)+1L]]<-
      data.frame(Request=ev$n,BestObjective=ev$best,ElapsedSec=clock_sec()-start)
    value
  }
  ev
}
population_bank <- function(seed,guide) {
  set.seed(seed)
  keys<-t(matrix(runif(200L*40L),nrow=200L,ncol=40L))
  g<-keys[1,]*.4; selected<-c(guide$f,100L+guide$m);g[selected]<-.6+keys[1,selected]*.4
  keys[1,]<-g
  list(keys=keys,guide=g,states=lapply(1:40,function(i)decode_keys(keys[i,])))
}
binary_decode <- function(bits,priority) {
  bits<-as.numeric(bits)
  if(length(bits)!=200L||anyNA(bits)||any(!bits%in%c(0,1))||length(priority)!=200L)
    stop("Invalid tabu binary encoding")
  # Fixed, seeded tie priorities are independent of genetic/mating values.
  # Ones are preferred; excess ones are trimmed and missing slots filled.
  canonical(list(f=order(-bits[1:100],priority[1:100])[1:12],
                 m=order(-bits[101:200],priority[101:200])[1:6]))
}
run_competitor_kernel <- function(name,params,guide,seed,ev) {
  if(!name %in% c("PSO","Tabu","SA","GA","DE"))stop("Unsupported competitor: ",name)
  ev$implementation<-switch(name,PSO="pso::psoptim",GA="GA::ga",DE="DEoptim::DEoptim",
    SA='stats::optim(method="SANN")',Tabu="tabuSearch::tabuSearch")
  if(name %in% c("PSO","GA","DE")) {
    stopifnot(params$pop==40L)
    bank<-population_bank(seed,guide);ev$bank<-bank;seen<-0L;set.seed(seed)
    objective<-function(keys) {
      st<-decode_keys(keys)
      if(seen<40L && (!identical(st,bank$states[[seen+1L]])||
         max(abs(keys-bank$keys[seen+1L,]))>1e-14))stop("Package initialization differs from saved bank: ",name)
      score<-ev$evaluate(st);seen<<-seen+1L;score
    }
    cap<-as.integer(max(1000,ev$remaining()+2));lower<-rep(0,200);upper<-rep(1,200)
    if(name=="PSO") {
      pso::psoptim(par=bank$guide,fn=function(x)-objective(x),lower=lower,upper=upper,
        control=list(s=40,w=params$w,type="SPSO2007",hybrid=FALSE,maxit=cap,maxf=Inf,trace=0))
    } else if(name=="GA") {
      GA::ga(type="real-valued",fitness=objective,lower=lower,upper=upper,popSize=40,
        suggestions=bank$keys,pmutation=params$mutation,pcrossover=.8,optim=FALSE,
        maxiter=cap,run=cap,monitor=FALSE,parallel=FALSE)
    } else {
      DEoptim::DEoptim(fn=function(x)-objective(x),lower=lower,upper=upper,
        control=DEoptim::DEoptim.control(NP=40,initialpop=bank$keys,F=params$F,CR=params$CR,
          strategy=2,itermax=cap,steptol=cap,trace=FALSE,parallelType="none"))
    }
    return("package_stop")
  }
  set.seed(seed)
  if(name=="Tabu") {
    if(!requireNamespace("tabuSearch",quietly=TRUE))stop("Install package tabuSearch")
    priority<-c(sample.int(100L),sample.int(100L))
    bits<-integer(200L);bits[c(guide$f,100L+guide$m)]<-1L
    stopifnot(identical(binary_decode(bits,priority),guide))
    # A sufficiently long preliminary phase lets the external request limit
    # stop the package before its additional intensification/diversification.
    iterations<-as.integer(ceiling((ev$remaining()-1L)/200L)+2L)
    ev$binary_adapter<-list(priority=priority,initial=bits,neigh=200L,
      iters=iterations,listSize=params$listSize,nRestarts=0L,repeatAll=1L,
      rule="ones-first-with-fixed-random-tie-priority")
    ev$binary_repair_requests<-0L
    objective<-function(x) {
      raw<-as.numeric(x)
      st<-binary_decode(raw,priority)
      value<-ev$evaluate(st)
      if(sum(raw[1:100])!=12L||sum(raw[101:200])!=6L)
        ev$binary_repair_requests<-ev$binary_repair_requests+1L
      value
    }
    tabuSearch::tabuSearch(size=200L,iters=iterations,objFunc=objective,config=bits,
      neigh=200L,listSize=params$listSize,nRestarts=0L,repeatAll=1L,verbose=FALSE)
    return("package_stop")
  }
  unpack<-function(x)canonical(list(f=x[1:12],m=x[13:18]))
  proposal<-function(x) {
    st<-unpack(x);f<-setdiff(1:100,st$f);m<-setdiff(1:100,st$m)
    nf<-12L*length(f);k<-sample.int(nf+6L*length(m),1L)
    if(k<=nf)st$f[(k-1L)%/%length(f)+1L]<-f[(k-1L)%%length(f)+1L] else {
      k<-k-nf;st$m[(k-1L)%/%length(m)+1L]<-m[(k-1L)%%length(m)+1L]
    }
    st<-canonical(st);c(st$f,st$m)
  }
  ev$temperature<-list(implementation="stats::optim SANN",temp=params$temp,tmax=params$tmax,
    schedule="package logarithmic cooling",proposal="uniform feasible single-parent swap")
  # gr is the documented proposal callback for SANN, not an acceptance rule.
  # Metropolis acceptance, cooling and the search loop remain in stats::optim.
  stats::optim(par=c(guide$f,guide$m),fn=function(x)-ev$evaluate(unpack(x)),gr=proposal,
    method="SANN",control=list(maxit=ev$remaining(),temp=params$temp,tmax=params$tmax,trace=FALSE))
  "package_stop"
}
run_competitor <- function(name,problem,params,seed,budget,api) {
  if(!name%in%c("PSO","Tabu","SA","GA","DE"))stop("Unsupported competitor")
  if(budget<82 || budget%%2)stop("Invalid comparison budget")
  begin<-clock_sec();guides<-make_guides(api,problem);prep<-clock_sec()-begin
  starts<-vector("list",2);post_time<-0
  for(i in 1:2) {
    # Reserve one lower solve per start for the final saved plan. The original
    # OSNB's extra initial-report solves also count toward its total allowance.
    ev<-make_evaluator(problem$M,api$core,budget/2-1L);t<-clock_sec()
    reason<-tryCatch(run_competitor_kernel(name,params,guides[[i]],seed+i,ev),
      gm_budget_exhausted=function(e)"budget")
    search_elapsed<-clock_sec()-t
    if(is.null(ev$state)||!is.finite(ev$best))stop("No valid result")
    t<-clock_sec();st<-ev$state
    fit<-api$core$gm_assignment_cpp(problem$M[st$f,st$m,drop=FALSE],TRUE)
    post<-clock_sec()-t;post_time<-post_time+post
    stopifnot(abs(fit$objective-ev$best)<1e-8,ev$n==ev$solves+ev$hits)
    starts[[i]]<-list(initial_state=guides[[i]],best_state=st,objective=ev$best,
      assignment=fit$assignment,trace=do.call(rbind,ev$trace),population_bank=ev$bank,
      temperature=ev$temperature,implementation=ev$implementation,
      binary_adapter=ev$binary_adapter,binary_repair_requests=ev$binary_repair_requests,
      summary=data.frame(Start=i,Source=c("original-value-guide","original-parent-guide")[i],
        InitialObjective=ev$first,BestObjective=ev$best,Requests=ev$n+1L,
        LowerSolves=ev$solves+1L,CacheHits=ev$hits,LocalSearchSec=search_elapsed,
        PoolFemale=100L,PoolMale=100L,StopReason=reason))
  }
  elapsed<-clock_sec()-begin;ss<-do.call(rbind,lapply(starts,`[[`,"summary"));best<-which.max(ss$BestObjective)
  stopifnot(sum(ss$Requests)<=budget)
  list(algorithm=name,parameters=params,starts=starts,
    summary=data.frame(Algorithm=name,BestObjective=ss$BestObjective[best],RuntimeSec=elapsed,
      Requests=sum(ss$Requests),LowerSolves=sum(ss$LowerSolves),CacheHits=sum(ss$CacheHits),
      BudgetLimit=budget,BestStart=best,SetupAndInitializationSec=prep,PostprocessSec=post_time,
      DistinctInitialSets=length(unique(vapply(guides,state_key,character(1))))))
}
