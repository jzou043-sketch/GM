required_pkgs <- c("lpSolve", "DEoptim", "pso", "GA", "parallel")
to_install <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(lpSolve)
  library(DEoptim)
  library(pso)
  library(GA)
  library(parallel)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

get_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) default else as.integer(value)
}

get_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) default else value
}

# =========================================================
# Current simulation settings requested by the user
# =========================================================
n_rep <- get_int("N_REP", 1000)
n_female_total <- get_int("N_FEMALE", 100)
n_male_total <- get_int("N_MALE", 100)
n_snp <- get_int("N_SNP", 1000)
n_major <- get_int("N_QTL", 200)
n_female_select <- get_int("N_FEMALE_SELECT", 12)
n_male_select <- get_int("N_MALE_SELECT", 6)
n_start <- get_int("N_START", 2)

female_pool_size <- get_int("FEMALE_POOL_SIZE", min(50, n_female_total))
male_pool_size <- get_int("MALE_POOL_SIZE", min(50, n_male_total))
max_neighbors_per_iter <- get_int("MAX_NEIGHBORS_PER_ITER", 300)

seed_base <- get_int("SEED_VALUE", 2026)
n_workers <- get_int("N_WORKERS", max(1, min(10, parallel::detectCores() - 1)))
chunk_size <- get_int("CHUNK_SIZE", 20)

out_dir <- get_chr("OUT_DIR", "optimizer_comparison_current_1000rep_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw_output_path <- file.path(out_dir, "benchmark_1000rep_optimizer_current_replicate_best.csv")
start_output_path <- file.path(out_dir, "benchmark_1000rep_optimizer_current_start_level.csv")
summary_output_path <- file.path(out_dir, "benchmark_1000rep_optimizer_current_summary.csv")
design_output_path <- file.path(out_dir, "benchmark_1000rep_optimizer_current_design.csv")

# =========================================================
# Hyperparameters retained from the previous benchmark file
# =========================================================
tabu_max_iter <- 80
tabu_tenure <- 7
sa_max_iter <- 140
sa_temp0 <- 1
sa_cooling <- 0.97
de_itermax <- 15
de_NP <- 70
pso_maxit <- 45
pso_s <- 50
ga_maxiter <- 35
ga_popSize <- 50
ga_run <- 12

algorithm_order <- c(
  "Particle swarm optimization",
  "Tabu search",
  "Proposed method",
  "Hill climbing",
  "Simulated annealing",
  "Genetic algorithm",
  "Differential evolution"
)

write.csv(
  data.frame(
    Parameter = c(
      "N_REP", "N_FEMALE", "N_MALE", "N_SNP", "N_QTL",
      "N_FEMALE_SELECT", "N_MALE_SELECT", "N_START",
      "FEMALE_POOL_SIZE", "MALE_POOL_SIZE", "MAX_NEIGHBORS_PER_ITER",
      "Tabu max_iter", "Tabu tenure",
      "SA max_iter", "SA temp0", "SA cooling",
      "DE itermax", "DE NP",
      "PSO maxit", "PSO swarm size",
      "GA maxiter", "GA popSize", "GA run",
      "Effect scenario", "Value objective"
    ),
    Value = c(
      n_rep, n_female_total, n_male_total, n_snp, n_major,
      n_female_select, n_male_select, n_start,
      female_pool_size, male_pool_size, max_neighbors_per_iter,
      tabu_max_iter, tabu_tenure,
      sa_max_iter, sa_temp0, sa_cooling,
      de_itermax, de_NP,
      pso_maxit, pso_s,
      ga_maxiter, ga_popSize, ga_run,
      "Strong dominance: additive U(0.50,0.80), dominance U(0.90,1.30)",
      "Expected progeny value only"
    )
  ),
  design_output_path,
  row.names = FALSE
)

cat("Output directory:", out_dir, "\n")
cat("Replicates:", n_rep, "\n")
cat("Parents:", n_female_total, "female +", n_male_total, "male\n")
cat("SNP:", n_snp, "| Effective SNP:", n_major, "\n")
cat("Select:", n_female_select, "female +", n_male_select, "male\n")
cat("Starts per algorithm:", n_start, "\n")
cat("Workers:", n_workers, "\n")

# =========================================================
# Simulation and mating-value construction
# =========================================================
simulate_value_problem <- function(rep_id) {
  set.seed(seed_base + rep_id)

  qtl_idx <- sort(sample(seq_len(n_snp), n_major))
  qtl_p <- runif(n_major, 0.4, 0.6)
  alpha_true <- runif(n_major, 0.50, 0.80)
  delta_true <- runif(n_major, 0.90, 1.30)

  G_female <- matrix(0L, nrow = n_female_total, ncol = n_major)
  G_male <- matrix(0L, nrow = n_male_total, ncol = n_major)
  for (k in seq_len(n_major)) {
    G_female[, k] <- rbinom(n_female_total, size = 2, prob = qtl_p[k])
    G_male[, k] <- rbinom(n_male_total, size = 2, prob = qtl_p[k])
  }

  female_ids <- paste0("F", seq_len(n_female_total))
  male_ids <- paste0("M", seq_len(n_male_total))

  rownames(G_female) <- female_ids
  rownames(G_male) <- male_ids

  p_female <- G_female / 2
  p_male <- G_male / 2

  # Expected progeny value:
  # E[Z] = pf + pm
  # E[W] = pf + pm - 2 pf pm
  # E[G] = alpha E[Z] + delta E[W]
  h <- alpha_true + delta_true
  female_main <- as.numeric(p_female %*% h)
  male_main <- as.numeric(p_male %*% h)
  cross_term <- (p_female * matrix(delta_true, nrow = n_female_total, ncol = n_major, byrow = TRUE)) %*% t(p_male)
  value_mat <- outer(female_main, male_main, "+") - 2 * cross_term

  rownames(value_mat) <- female_ids
  colnames(value_mat) <- male_ids

  parent_value_df <- data.frame(
    ID = c(female_ids, male_ids),
    Sex = c(rep("F", n_female_total), rep("M", n_male_total)),
    GeneticValue = c(
      as.numeric(G_female %*% alpha_true + (G_female == 1) %*% delta_true),
      as.numeric(G_male %*% alpha_true + (G_male == 1) %*% delta_true)
    ),
    stringsAsFactors = FALSE
  )

  if (rep_id == 1) {
    full_geno <- matrix(0L, nrow = n_female_total + n_male_total, ncol = n_snp)
    for (k in seq_len(n_snp)) {
      pk <- runif(1, 0.4, 0.6)
      full_geno[, k] <- rbinom(n_female_total + n_male_total, size = 2, prob = pk)
    }
    full_geno[seq_len(n_female_total), qtl_idx] <- G_female
    full_geno[n_female_total + seq_len(n_male_total), qtl_idx] <- G_male
    rownames(full_geno) <- c(female_ids, male_ids)
    colnames(full_geno) <- paste0("SNP", seq_len(n_snp))

    write.csv(full_geno, file.path(out_dir, "example_rep001_parent_genotypes.csv"))
    write.csv(
      data.frame(
        SNP = paste0("SNP", qtl_idx),
        QTLIndex = qtl_idx,
        AlleleFrequency = qtl_p,
        AdditiveEffect = alpha_true,
        DominanceEffect = delta_true,
        stringsAsFactors = FALSE
      ),
      file.path(out_dir, "example_rep001_marker_effects.csv"),
      row.names = FALSE
    )
    write.csv(value_mat, file.path(out_dir, "example_rep001_value_matrix.csv"))
  }

  list(
    value_mat = value_mat,
    parent_value_df = parent_value_df,
    female_ids = female_ids,
    male_ids = male_ids
  )
}

# =========================================================
# Common subset evaluation by lower-level ILP
# =========================================================
solve_mating_ilp_general <- function(value_mat_sub, force_use_all_males = TRUE) {
  nf <- nrow(value_mat_sub)
  nm <- ncol(value_mat_sub)
  obj <- as.vector(value_mat_sub)

  mother_block <- matrix(0, nrow = nf, ncol = nf * nm)
  for (i in seq_len(nf)) {
    mother_block[i, i + (0:(nm - 1)) * nf] <- 1
  }

  sire_block <- matrix(0, nrow = nm, ncol = nf * nm)
  for (j in seq_len(nm)) {
    sire_block[j, ((j - 1) * nf + 1):(j * nf)] <- 1
  }

  lower <- if (force_use_all_males && nf >= nm) floor(nf / nm) else 0
  upper <- ceiling(nf / nm)

  fit <- lpSolve::lp(
    direction = "max",
    objective.in = obj,
    const.mat = rbind(mother_block, sire_block, sire_block),
    const.dir = c(rep("=", nf), rep(">=", nm), rep("<=", nm)),
    const.rhs = c(rep(1, nf), rep(lower, nm), rep(upper, nm)),
    all.bin = TRUE
  )

  if (fit$status != 0) stop("ILP failed for the given subset.")
  list(objective = fit$objval)
}

make_evaluator <- function(value_mat) {
  cache_env <- new.env(hash = TRUE, parent = emptyenv())
  force(value_mat)

  function(female_sel, male_sel) {
    female_sel <- sort(female_sel)
    male_sel <- sort(male_sel)
    key <- paste0(paste(female_sel, collapse = ","), "|", paste(male_sel, collapse = ","))
    if (!exists(key, envir = cache_env, inherits = FALSE)) {
      submat <- value_mat[female_sel, male_sel, drop = FALSE]
      fit <- solve_mating_ilp_general(submat, force_use_all_males = TRUE)
      assign(key, fit$objective, envir = cache_env)
    }
    get(key, envir = cache_env, inherits = FALSE)
  }
}

as_state <- function(female, male) {
  list(female = sort(unique(female)), male = sort(unique(male)))
}

state_signature <- function(state) {
  paste0(paste(sort(state$female), collapse = ","), "|", paste(sort(state$male), collapse = ","))
}

sample_if_needed <- function(x, max_n) {
  if (is.infinite(max_n) || length(x) <= max_n) return(x)
  x[sample.int(length(x), max_n)]
}

all_neighbors <- function(state, female_space, male_space, max_neighbors = max_neighbors_per_iter) {
  female_out <- state$female
  female_in <- setdiff(female_space, state$female)
  male_out <- state$male
  male_in <- setdiff(male_space, state$male)

  neighbors <- vector("list", length(female_out) * length(female_in) + length(male_out) * length(male_in))
  idx <- 1L
  for (f_out in female_out) {
    for (f_in in female_in) {
      neighbors[[idx]] <- as_state(c(setdiff(state$female, f_out), f_in), state$male)
      idx <- idx + 1L
    }
  }
  for (m_out in male_out) {
    for (m_in in male_in) {
      neighbors[[idx]] <- as_state(state$female, c(setdiff(state$male, m_out), m_in))
      idx <- idx + 1L
    }
  }

  sample_if_needed(neighbors, max_neighbors)
}

build_value_space <- function(value_mat, female_ids, male_ids, female_pool_size, male_pool_size) {
  female_score <- rowMeans(value_mat)
  male_score <- colMeans(value_mat)
  list(
    female_space = names(sort(female_score, decreasing = TRUE))[seq_len(min(female_pool_size, length(female_ids)))],
    male_space = names(sort(male_score, decreasing = TRUE))[seq_len(min(male_pool_size, length(male_ids)))]
  )
}

build_parent_value_space <- function(parent_value_df, female_pool_size, male_pool_size) {
  female_value_df <- subset(parent_value_df, Sex == "F")
  male_value_df <- subset(parent_value_df, Sex == "M")
  female_ord <- female_value_df[order(female_value_df$GeneticValue, decreasing = TRUE), ]
  male_ord <- male_value_df[order(male_value_df$GeneticValue, decreasing = TRUE), ]
  list(
    female_space = female_ord$ID[seq_len(min(female_pool_size, nrow(female_ord)))],
    male_space = male_ord$ID[seq_len(min(male_pool_size, nrow(male_ord)))]
  )
}

choose_search_mode <- function(
    n_female_select, n_male_select, female_space_size, male_space_size,
    neighborhood_threshold = 1500, ilp_var_threshold = 200, workload_threshold = 200000
) {
  neighborhood_count <-
    n_female_select * max(0, female_space_size - n_female_select) +
    n_male_select * max(0, male_space_size - n_male_select)
  ilp_vars <- n_female_select * n_male_select
  workload <- neighborhood_count * ilp_vars
  if (neighborhood_count <= neighborhood_threshold &&
      ilp_vars <= ilp_var_threshold &&
      workload <= workload_threshold) "best" else "first"
}

greedy_init_in_space <- function(value_mat, female_space, male_space) {
  submat <- value_mat[female_space, male_space, drop = FALSE]
  male_score <- colMeans(submat)
  selected_males <- names(sort(male_score, decreasing = TRUE))[seq_len(n_male_select)]
  female_score <- rowMeans(submat[, selected_males, drop = FALSE])
  selected_females <- names(sort(female_score, decreasing = TRUE))[seq_len(n_female_select)]
  as_state(selected_females, selected_males)
}

random_init_in_space <- function(female_space, male_space, seed) {
  set.seed(seed)
  as_state(
    sample(female_space, n_female_select, replace = FALSE),
    sample(male_space, n_male_select, replace = FALSE)
  )
}

hill_climb_from_state <- function(initial_state, female_space, male_space, evaluate_subset, mode = "best") {
  state <- as_state(initial_state$female, initial_state$male)
  current_obj <- evaluate_subset(state$female, state$male)
  repeat {
    neighbors <- all_neighbors(state, female_space, male_space)
    if (length(neighbors) == 0) break

    if (mode == "best") {
      best_state <- state
      best_obj <- current_obj
      for (cand in neighbors) {
        obj <- evaluate_subset(cand$female, cand$male)
        if (obj > best_obj) {
          best_obj <- obj
          best_state <- cand
        }
      }
      if (best_obj <= current_obj + 1e-12) break
      state <- best_state
      current_obj <- best_obj
    } else {
      improved <- FALSE
      for (cand in neighbors) {
        obj <- evaluate_subset(cand$female, cand$male)
        if (obj > current_obj + 1e-12) {
          state <- cand
          current_obj <- obj
          improved <- TRUE
          break
        }
      }
      if (!improved) break
    }
  }
  list(objective = current_obj, state = state)
}

tabu_search_from_state <- function(initial_state, female_ids, male_ids, evaluate_subset, max_iter = tabu_max_iter, tenure = tabu_tenure) {
  state <- as_state(initial_state$female, initial_state$male)
  current_obj <- evaluate_subset(state$female, state$male)
  best_state <- state
  best_obj <- current_obj
  tabu_queue <- character(0)

  for (iter in seq_len(max_iter)) {
    neighbors <- all_neighbors(state, female_ids, male_ids)
    signatures <- vapply(neighbors, state_signature, character(1))
    objs <- vapply(neighbors, function(st) evaluate_subset(st$female, st$male), numeric(1))
    admissible <- !(signatures %in% tabu_queue)
    aspiration <- objs > best_obj + 1e-12
    allowed <- admissible | aspiration
    if (!any(allowed)) next
    idx <- which(allowed)[which.max(objs[allowed])]
    state <- neighbors[[idx]]
    current_obj <- objs[idx]
    tabu_queue <- c(tabu_queue, signatures[idx])
    if (length(tabu_queue) > tenure) tabu_queue <- tabu_queue[-1]
    if (current_obj > best_obj + 1e-12) {
      best_obj <- current_obj
      best_state <- state
    }
  }
  list(objective = best_obj, state = best_state)
}

simulated_annealing_from_state <- function(initial_state, female_ids, male_ids, evaluate_subset, max_iter = sa_max_iter, temp0 = sa_temp0, cooling = sa_cooling) {
  state <- as_state(initial_state$female, initial_state$male)
  current_obj <- evaluate_subset(state$female, state$male)
  best_state <- state
  best_obj <- current_obj
  temp <- temp0
  for (iter in seq_len(max_iter)) {
    neighbors <- all_neighbors(state, female_ids, male_ids)
    cand <- neighbors[[sample.int(length(neighbors), 1)]]
    cand_obj <- evaluate_subset(cand$female, cand$male)
    delta <- cand_obj - current_obj
    if (delta >= 0 || runif(1) < exp(delta / max(temp, 1e-8))) {
      state <- cand
      current_obj <- cand_obj
      if (current_obj > best_obj) {
        best_obj <- current_obj
        best_state <- state
      }
    }
    temp <- temp * cooling
  }
  list(objective = best_obj, state = best_state)
}

decode_real_keys <- function(x, female_ids, male_ids) {
  female_score <- x[seq_len(length(female_ids))]
  male_score <- x[length(female_ids) + seq_len(length(male_ids))]
  female_sel <- female_ids[order(female_score, decreasing = TRUE)][seq_len(n_female_select)]
  male_sel <- male_ids[order(male_score, decreasing = TRUE)][seq_len(n_male_select)]
  as_state(female_sel, male_sel)
}

run_deoptim_once <- function(seed, female_ids, male_ids, evaluate_subset, itermax = de_itermax, NP = de_NP) {
  set.seed(seed)
  fn <- function(x) {
    st <- decode_real_keys(x, female_ids, male_ids)
    -evaluate_subset(st$female, st$male)
  }
  fit <- DEoptim::DEoptim(
    fn = fn,
    lower = rep(0, length(female_ids) + length(male_ids)),
    upper = rep(1, length(female_ids) + length(male_ids)),
    control = DEoptim::DEoptim.control(itermax = itermax, NP = NP, trace = FALSE)
  )
  decoded <- decode_real_keys(fit$optim$bestmem, female_ids, male_ids)
  list(objective = evaluate_subset(decoded$female, decoded$male), state = decoded)
}

run_pso_once <- function(seed, female_ids, male_ids, evaluate_subset, maxit = pso_maxit, s = pso_s) {
  set.seed(seed)
  fn <- function(x) {
    st <- decode_real_keys(x, female_ids, male_ids)
    -evaluate_subset(st$female, st$male)
  }
  fit <- pso::psoptim(
    par = rep(0.5, length(female_ids) + length(male_ids)),
    fn = fn,
    lower = rep(0, length(female_ids) + length(male_ids)),
    upper = rep(1, length(female_ids) + length(male_ids)),
    control = list(maxit = maxit, s = s, trace = FALSE)
  )
  decoded <- decode_real_keys(fit$par, female_ids, male_ids)
  list(objective = evaluate_subset(decoded$female, decoded$male), state = decoded)
}

run_ga_once <- function(seed, female_ids, male_ids, evaluate_subset, maxiter = ga_maxiter, popSize = ga_popSize) {
  set.seed(seed)
  fitness_fn <- function(x) {
    st <- decode_real_keys(x, female_ids, male_ids)
    evaluate_subset(st$female, st$male)
  }
  fit <- GA::ga(
    type = "real-valued",
    fitness = fitness_fn,
    lower = rep(0, length(female_ids) + length(male_ids)),
    upper = rep(1, length(female_ids) + length(male_ids)),
    popSize = popSize,
    maxiter = maxiter,
    run = ga_run,
    optim = FALSE,
    monitor = FALSE,
    seed = seed
  )
  decoded <- decode_real_keys(fit@solution[1, ], female_ids, male_ids)
  list(objective = evaluate_subset(decoded$female, decoded$male), state = decoded)
}

benchmark_one_replicate <- function(rep_id) {
  problem <- simulate_value_problem(rep_id)
  value_mat <- problem$value_mat
  parent_value_df <- problem$parent_value_df
  female_ids <- problem$female_ids
  male_ids <- problem$male_ids
  evaluate_subset <- make_evaluator(value_mat)

  value_space <- build_value_space(value_mat, female_ids, male_ids, female_pool_size, male_pool_size)
  parent_space <- build_parent_value_space(parent_value_df, female_pool_size, male_pool_size)
  female_space_union <- union(value_space$female_space, parent_space$female_space)
  male_space_union <- union(value_space$male_space, parent_space$male_space)
  proposed_mode <- choose_search_mode(n_female_select, n_male_select, length(female_space_union), length(male_space_union))

  proposed_starts <- list(
    greedy_init_in_space(value_mat, value_space$female_space, value_space$male_space),
    greedy_init_in_space(value_mat, parent_space$female_space, parent_space$male_space)
  )[seq_len(n_start)]

  fullspace_starts <- list(
    greedy_init_in_space(value_mat, female_ids, male_ids),
    random_init_in_space(female_ids, male_ids, 7200 + rep_id * 10 + 1)
  )[seq_len(n_start)]

  benchmark_from_fixed_starts <- function(name, runner, starts) {
    start_records <- vector("list", length(starts))
    for (s_idx in seq_along(starts)) {
      elapsed <- system.time({
        res <- runner(starts[[s_idx]])
      })
      start_records[[s_idx]] <- data.frame(
        replicate = rep_id,
        Algorithm = name,
        Start = s_idx,
        Objective = res$objective,
        ElapsedSec = unname(elapsed["elapsed"]),
        stringsAsFactors = FALSE
      )
    }
    start_df <- do.call(rbind, start_records)
    best_idx <- which.max(start_df$Objective)
    best_df <- data.frame(
      replicate = rep_id,
      Algorithm = name,
      BestObjective = start_df$Objective[best_idx],
      MeanObjective = mean(start_df$Objective),
      TotalElapsedSec = sum(start_df$ElapsedSec),
      BestStart = start_df$Start[best_idx],
      stringsAsFactors = FALSE
    )
    list(best = best_df, starts = start_df)
  }

  benchmark_from_seed_runs <- function(name, runner, seeds) {
    start_records <- vector("list", length(seeds))
    for (s_idx in seq_along(seeds)) {
      elapsed <- system.time({
        res <- runner(seeds[s_idx])
      })
      start_records[[s_idx]] <- data.frame(
        replicate = rep_id,
        Algorithm = name,
        Start = s_idx,
        Objective = res$objective,
        ElapsedSec = unname(elapsed["elapsed"]),
        stringsAsFactors = FALSE
      )
    }
    start_df <- do.call(rbind, start_records)
    best_idx <- which.max(start_df$Objective)
    best_df <- data.frame(
      replicate = rep_id,
      Algorithm = name,
      BestObjective = start_df$Objective[best_idx],
      MeanObjective = mean(start_df$Objective),
      TotalElapsedSec = sum(start_df$ElapsedSec),
      BestStart = start_df$Start[best_idx],
      stringsAsFactors = FALSE
    )
    list(best = best_df, starts = start_df)
  }

  results <- list(
    benchmark_from_fixed_starts(
      "Proposed method",
      function(st) hill_climb_from_state(st, female_space_union, male_space_union, evaluate_subset, proposed_mode),
      proposed_starts
    ),
    benchmark_from_fixed_starts(
      "Hill climbing",
      function(st) hill_climb_from_state(st, female_ids, male_ids, evaluate_subset, "best"),
      fullspace_starts
    ),
    benchmark_from_fixed_starts(
      "Tabu search",
      function(st) tabu_search_from_state(st, female_ids, male_ids, evaluate_subset),
      fullspace_starts
    ),
    benchmark_from_fixed_starts(
      "Simulated annealing",
      function(st) simulated_annealing_from_state(st, female_ids, male_ids, evaluate_subset),
      fullspace_starts
    ),
    benchmark_from_seed_runs(
      "Genetic algorithm",
      function(sd) run_ga_once(sd, female_ids, male_ids, evaluate_subset),
      8101 + rep_id * 10 + seq_len(n_start) - 1
    ),
    benchmark_from_seed_runs(
      "Differential evolution",
      function(sd) run_deoptim_once(sd, female_ids, male_ids, evaluate_subset),
      8201 + rep_id * 10 + seq_len(n_start) - 1
    ),
    benchmark_from_seed_runs(
      "Particle swarm optimization",
      function(sd) run_pso_once(sd, female_ids, male_ids, evaluate_subset),
      8301 + rep_id * 10 + seq_len(n_start) - 1
    )
  )

  list(
    best = do.call(rbind, lapply(results, `[[`, "best")),
    starts = do.call(rbind, lapply(results, `[[`, "starts"))
  )
}

build_summary <- function(raw_results) {
  replicate_best <- aggregate(BestObjective ~ replicate, data = raw_results, FUN = max)
  names(replicate_best)[2] <- "ReplicateBest"
  raw_aug <- merge(raw_results, replicate_best, by = "replicate", all.x = TRUE)
  raw_aug$TopHit <- raw_aug$BestObjective >= raw_aug$ReplicateBest - 1e-10
  raw_aug$Regret <- raw_aug$ReplicateBest - raw_aug$BestObjective

  summary_results <- do.call(
    rbind,
    lapply(split(raw_aug, raw_aug$Algorithm), function(df) {
      data.frame(
        Algorithm = df$Algorithm[1],
        MeanBestObjective = mean(df$BestObjective),
        SDBestObjective = sd(df$BestObjective),
        MedianBestObjective = median(df$BestObjective),
        IQRBestObjective = IQR(df$BestObjective),
        MeanMeanObjective = mean(df$MeanObjective),
        MeanRuntimeSec = mean(df$TotalElapsedSec),
        MedianRuntimeSec = median(df$TotalElapsedSec),
        TopHitRatePct = 100 * mean(df$TopHit),
        MeanRegret = mean(df$Regret),
        stringsAsFactors = FALSE
      )
    })
  )

  summary_results <- summary_results[match(algorithm_order, summary_results$Algorithm), ]
  row.names(summary_results) <- NULL
  list(raw_aug = raw_aug, summary_results = summary_results)
}

initialize_worker <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  suppressPackageStartupMessages({
    library(lpSolve)
    library(DEoptim)
    library(pso)
    library(GA)
  })
  TRUE
}

cl <- makeCluster(n_workers)
on.exit(stopCluster(cl), add = TRUE)

export_names <- setdiff(ls(envir = .GlobalEnv), "cl")
clusterExport(cl, export_names, envir = .GlobalEnv)
clusterEvalQ(cl, initialize_worker())

if (file.exists(raw_output_path)) {
  existing_raw <- read.csv(raw_output_path, stringsAsFactors = FALSE)
  completed_rep_ids <- sort(unique(existing_raw$replicate))
  message("Resuming from existing replicate-best file. Completed replicates: ", length(completed_rep_ids))
} else {
  existing_raw <- NULL
  completed_rep_ids <- integer(0)
}

if (file.exists(start_output_path)) {
  existing_starts <- read.csv(start_output_path, stringsAsFactors = FALSE)
} else {
  existing_starts <- NULL
}

remaining_rep_ids <- setdiff(seq_len(n_rep), completed_rep_ids)
raw_results <- existing_raw
start_results <- existing_starts

if (length(remaining_rep_ids) > 0) {
  chunks <- split(remaining_rep_ids, ceiling(seq_along(remaining_rep_ids) / chunk_size))
  for (chunk_idx in seq_along(chunks)) {
    rep_ids_chunk <- chunks[[chunk_idx]]
    message(
      "Running chunk ", chunk_idx, "/", length(chunks),
      " | replicates ", min(rep_ids_chunk), "-", max(rep_ids_chunk),
      " | n = ", length(rep_ids_chunk)
    )

    raw_list_chunk <- parLapplyLB(cl, rep_ids_chunk, benchmark_one_replicate)
    chunk_best <- do.call(rbind, lapply(raw_list_chunk, `[[`, "best"))
    chunk_starts <- do.call(rbind, lapply(raw_list_chunk, `[[`, "starts"))

    raw_results <- if (is.null(raw_results)) chunk_best else rbind(raw_results, chunk_best)
    start_results <- if (is.null(start_results)) chunk_starts else rbind(start_results, chunk_starts)
    row.names(raw_results) <- NULL
    row.names(start_results) <- NULL

    built <- build_summary(raw_results)
    write.csv(built$raw_aug, file = raw_output_path, row.names = FALSE)
    write.csv(start_results, file = start_output_path, row.names = FALSE)
    write.csv(built$summary_results, file = summary_output_path, row.names = FALSE)
  }
}

built <- build_summary(read.csv(raw_output_path, stringsAsFactors = FALSE))
raw_results <- built$raw_aug
summary_results <- built$summary_results

write.csv(raw_results, file = raw_output_path, row.names = FALSE)
write.csv(summary_results, file = summary_output_path, row.names = FALSE)

writeLines(
  c(
    "Current optimizer benchmark based on the previous optimizer-only benchmark.",
    "",
    "Requested simulation scale:",
    paste0("- Female parents: ", n_female_total),
    paste0("- Male parents: ", n_male_total),
    paste0("- SNPs: ", n_snp),
    paste0("- Effective SNPs: ", n_major),
    paste0("- Selected females: ", n_female_select),
    paste0("- Selected males: ", n_male_select),
    paste0("- Replicates: ", n_rep),
    paste0("- Starts per algorithm: ", n_start),
    "",
    "Core hyperparameters retained from the previous benchmark:",
    paste0("- Tabu search: max_iter=", tabu_max_iter, ", tenure=", tabu_tenure),
    paste0("- Simulated annealing: max_iter=", sa_max_iter, ", temp0=", sa_temp0, ", cooling=", sa_cooling),
    paste0("- Differential evolution: itermax=", de_itermax, ", NP=", de_NP),
    paste0("- Particle swarm optimization: maxit=", pso_maxit, ", s=", pso_s),
    paste0("- Genetic algorithm: maxiter=", ga_maxiter, ", popSize=", ga_popSize, ", run=", ga_run),
    "",
    "Note: MAX_NEIGHBORS_PER_ITER limits the number of one-swap neighbors evaluated per step so the 100x100 benchmark can finish.",
    "Set MAX_NEIGHBORS_PER_ITER to a very large number if an exact full-neighborhood run is required.",
    "",
    "Output files:",
    "- benchmark_1000rep_optimizer_current_replicate_best.csv",
    "- benchmark_1000rep_optimizer_current_start_level.csv",
    "- benchmark_1000rep_optimizer_current_summary.csv",
    "- benchmark_1000rep_optimizer_current_design.csv"
  ),
  con = file.path(out_dir, "README.txt")
)

tar_file <- paste0(out_dir, ".tar.gz")
if (file.exists(tar_file)) unlink(tar_file)
utils::tar(tarfile = tar_file, files = out_dir, compression = "gzip")

cat("\n=== Current 1000-replicate optimizer-only summary ===\n")
print(summary_results)
cat("\nOutput directory:", out_dir, "\n")
cat("Archive:", tar_file, "\n")
