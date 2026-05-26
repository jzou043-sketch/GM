if (!requireNamespace("lpSolve", quietly = TRUE)) {
  install.packages("lpSolve", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(lpSolve)
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

get_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) default else as.numeric(value)
}

get_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) default else value
}

# ---------------------------------------------------------
# 1. Founder-population generator
# ---------------------------------------------------------
generate_large_gm_population <- function(
    scenario = c("S1_weak", "S2_medium", "S3_strong"),
    seed = 2026,
    n_female = 100,
    n_male = 100,
    m = 1000,
    n_major = 200,
    H2 = 0.50,
    mu = 50
) {
  scenario <- match.arg(scenario)
  set.seed(seed)

  n_ind <- n_female + n_male
  if (n_major >= m) stop("n_major must be smaller than m.")

  ID <- paste0("ID", sprintf("%04d", seq_len(n_ind)))
  Sex <- c(rep("F", n_female), rep("M", n_male))

  major_idx <- sort(sample(seq_len(m), n_major, replace = FALSE))
  bg_idx <- setdiff(seq_len(m), major_idx)

  G <- matrix(0L, nrow = n_ind, ncol = m)
  colnames(G) <- paste0("SNP", sprintf("%04d", seq_len(m)))
  for (k in seq_len(m)) {
    pk <- runif(1, 0.4, 0.6)
    G[, k] <- rbinom(n_ind, size = 2, prob = pk)
  }

  W_true <- ifelse(G == 1, 1, 0)
  alpha_true <- rep(0, m)
  delta_true <- rep(0, m)

  if (scenario == "S1_weak") {
    alpha_true[major_idx] <- runif(n_major, 0.50, 0.80)
    delta_true[major_idx] <- runif(n_major, 0.05, 0.15)
  } else if (scenario == "S2_medium") {
    alpha_true[major_idx] <- runif(n_major, 0.55, 0.85)
    delta_true[major_idx] <- runif(n_major, 0.35, 0.55)
  } else {
    alpha_true[major_idx] <- runif(n_major, 0.50, 0.80)
    delta_true[major_idx] <- runif(n_major, 0.90, 1.30)
  }

  g_true <- as.numeric(G %*% alpha_true + W_true %*% delta_true)
  var_g <- var(g_true)
  sigma_e2 <- var_g * (1 - H2) / H2
  sigma_e <- sqrt(sigma_e2)
  Phenotype <- as.numeric(mu + g_true + rnorm(n_ind, 0, sigma_e))

  dat <- data.frame(
    ID = ID,
    Sex = Sex,
    G,
    Phenotype = Phenotype,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  truth_summary <- data.frame(
    Scenario = scenario,
    N = n_ind,
    Female = n_female,
    Male = n_male,
    SNP = m,
    EffectiveSNP = n_major,
    BackgroundSNP = length(bg_idx),
    MeanAlpha = mean(alpha_true[major_idx]),
    SDAlpha = sd(alpha_true[major_idx]),
    MeanDelta = mean(delta_true[major_idx]),
    SDDelta = sd(delta_true[major_idx]),
    MeanDeltaOverAlpha = mean(delta_true[major_idx] / alpha_true[major_idx]),
    TrueGeneticMean = mean(g_true),
    TrueGeneticSD = sd(g_true),
    PhenotypeMean = mean(Phenotype),
    PhenotypeSD = sd(Phenotype),
    TrueVarG = var_g,
    ResidualVar = sigma_e2,
    TargetH2 = H2,
    stringsAsFactors = FALSE
  )

  list(
    dat = dat,
    G = G,
    W_true = W_true,
    alpha_true = alpha_true,
    delta_true = delta_true,
    g_true = g_true,
    sigma_e = sigma_e,
    sigma_e2 = sigma_e2,
    mu = mu,
    major_idx = major_idx,
    bg_idx = bg_idx,
    truth_summary = truth_summary
  )
}

# ---------------------------------------------------------
# 2. Utility functions
# ---------------------------------------------------------
compute_true_scores <- function(dat, alpha_true, delta_true) {
  snp_cols <- grep("^SNP", names(dat))
  Z <- as.matrix(dat[, snp_cols, drop = FALSE])
  W <- ifelse(Z == 1, 1, 0)
  as.numeric(Z %*% alpha_true + W %*% delta_true)
}

build_value_matrices_from_truth <- function(dat, alpha_true, delta_true) {
  snp_cols <- grep("^SNP", names(dat))
  Z_all <- as.matrix(dat[, snp_cols, drop = FALSE])

  female_idx <- which(dat$Sex == "F")
  male_idx <- which(dat$Sex == "M")
  female_id <- dat$ID[female_idx]
  male_id <- dat$ID[male_idx]
  G_female <- Z_all[female_idx, , drop = FALSE]
  G_male <- Z_all[male_idx, , drop = FALSE]

  pf <- G_female / 2
  pm <- G_male / 2
  h <- alpha_true + delta_true

  female_main <- as.numeric(pf %*% h)
  male_main <- as.numeric(pm %*% h)
  cross_term <- (pf * matrix(delta_true, nrow = nrow(pf), ncol = length(delta_true), byrow = TRUE)) %*% t(pm)

  mean_mat <- outer(female_main, male_main, "+") - 2 * cross_term
  dimnames(mean_mat) <- list(female_id, male_id)

  genetic_value <- compute_true_scores(dat, alpha_true, delta_true)
  parent_value_df <- data.frame(
    ID = dat$ID,
    Sex = dat$Sex,
    Phenotype = dat$Phenotype,
    GeneticValue = genetic_value,
    stringsAsFactors = FALSE
  )

  list(
    mean_mat = mean_mat,
    value_mat = mean_mat,
    parent_value_df = parent_value_df
  )
}

solve_mating_ilp_general <- function(value_mat_sub) {
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

  lower <- floor(nf / nm)
  upper <- ceiling(nf / nm)

  fit <- lpSolve::lp(
    direction = "max",
    objective.in = obj,
    const.mat = rbind(mother_block, sire_block, sire_block),
    const.dir = c(rep("=", nf), rep(">=", nm), rep("<=", nm)),
    const.rhs = c(rep(1, nf), rep(lower, nm), rep(upper, nm)),
    all.bin = TRUE
  )

  if (fit$status != 0) stop("Lower-level ILP failed.")

  assignment <- matrix(
    fit$solution,
    nrow = nf,
    ncol = nm,
    dimnames = list(rownames(value_mat_sub), colnames(value_mat_sub))
  )

  list(
    objective = fit$objval,
    assignment = assignment,
    sire_load = colSums(assignment)
  )
}

as_state <- function(female, male) {
  list(female = sort(unique(female)), male = sort(unique(male)))
}

build_value_space <- function(value_mat, female_pool_size, male_pool_size) {
  female_pool_size <- min(female_pool_size, nrow(value_mat))
  male_pool_size <- min(male_pool_size, ncol(value_mat))

  female_score <- rowMeans(value_mat)
  male_score <- colMeans(value_mat)
  list(
    female_space = names(sort(female_score, decreasing = TRUE))[seq_len(female_pool_size)],
    male_space = names(sort(male_score, decreasing = TRUE))[seq_len(male_pool_size)]
  )
}

build_parent_value_space <- function(female_value_df, male_value_df, female_pool_size, male_pool_size) {
  female_pool_size <- min(female_pool_size, nrow(female_value_df))
  male_pool_size <- min(male_pool_size, nrow(male_value_df))

  female_ord <- female_value_df[order(female_value_df$GeneticValue, decreasing = TRUE), ]
  male_ord <- male_value_df[order(male_value_df$GeneticValue, decreasing = TRUE), ]
  list(
    female_space = female_ord$ID[seq_len(female_pool_size)],
    male_space = male_ord$ID[seq_len(male_pool_size)]
  )
}

greedy_init_in_space <- function(value_mat, female_space, male_space, n_female_select, n_male_select) {
  submat <- value_mat[female_space, male_space, drop = FALSE]
  male_score <- colMeans(submat)
  selected_males <- names(sort(male_score, decreasing = TRUE))[seq_len(n_male_select)]
  female_score <- rowMeans(submat[, selected_males, drop = FALSE])
  selected_females <- names(sort(female_score, decreasing = TRUE))[seq_len(n_female_select)]
  as_state(selected_females, selected_males)
}

random_init_in_space <- function(female_space, male_space, n_female_select, n_male_select, seed) {
  set.seed(seed)
  as_state(
    sample(female_space, n_female_select, replace = FALSE),
    sample(male_space, n_male_select, replace = FALSE)
  )
}

choose_search_mode <- function(n_female_select, n_male_select, female_space_size, male_space_size) {
  neighborhood_count <-
    n_female_select * max(0, female_space_size - n_female_select) +
    n_male_select * max(0, male_space_size - n_male_select)
  ilp_vars <- n_female_select * n_male_select
  workload <- neighborhood_count * ilp_vars
  if (neighborhood_count <= 1500 && ilp_vars <= 200 && workload <= 200000) "best" else "first"
}

run_proposed_optimizer <- function(
    value_mat,
    female_value_df,
    male_value_df,
    n_female_select = 12,
    n_male_select = 6,
    female_pool_size,
    male_pool_size,
    n_restart = 2,
    seed = 2026
) {
  cache_env <- new.env(hash = TRUE, parent = emptyenv())

  evaluate_state <- function(state) {
    key <- paste0(paste(state$female, collapse = ","), "|", paste(state$male, collapse = ","))
    if (!exists(key, envir = cache_env, inherits = FALSE)) {
      submat <- value_mat[state$female, state$male, drop = FALSE]
      fit <- solve_mating_ilp_general(submat)
      assign(key, fit, envir = cache_env)
    }
    get(key, envir = cache_env, inherits = FALSE)
  }

  hill_climb_from_state <- function(initial_state, female_space, male_space, mode) {
    state <- as_state(initial_state$female, initial_state$male)
    current_fit <- evaluate_state(state)

    repeat {
      fixed_females <- state$female
      fixed_males <- state$male
      unselected_females <- setdiff(female_space, fixed_females)
      unselected_males <- setdiff(male_space, fixed_males)

      if (mode == "best") {
        best_state <- state
        best_fit <- current_fit

        for (f_out in fixed_females) {
          for (f_in in unselected_females) {
            cand_state <- as_state(c(setdiff(fixed_females, f_out), f_in), fixed_males)
            cand_fit <- evaluate_state(cand_state)
            if (cand_fit$objective > best_fit$objective + 1e-12) {
              best_state <- cand_state
              best_fit <- cand_fit
            }
          }
        }

        for (m_out in fixed_males) {
          for (m_in in unselected_males) {
            cand_state <- as_state(fixed_females, c(setdiff(fixed_males, m_out), m_in))
            cand_fit <- evaluate_state(cand_state)
            if (cand_fit$objective > best_fit$objective + 1e-12) {
              best_state <- cand_state
              best_fit <- cand_fit
            }
          }
        }

        if (best_fit$objective <= current_fit$objective + 1e-12) break
        state <- best_state
        current_fit <- best_fit
      } else {
        improved <- FALSE

        for (f_out in fixed_females) {
          if (improved) break
          for (f_in in unselected_females) {
            cand_state <- as_state(c(setdiff(fixed_females, f_out), f_in), fixed_males)
            cand_fit <- evaluate_state(cand_state)
            if (cand_fit$objective > current_fit$objective + 1e-12) {
              state <- cand_state
              current_fit <- cand_fit
              improved <- TRUE
              break
            }
          }
        }

        if (!improved) {
          for (m_out in fixed_males) {
            if (improved) break
            for (m_in in unselected_males) {
              cand_state <- as_state(fixed_females, c(setdiff(fixed_males, m_out), m_in))
              cand_fit <- evaluate_state(cand_state)
              if (cand_fit$objective > current_fit$objective + 1e-12) {
                state <- cand_state
                current_fit <- cand_fit
                improved <- TRUE
                break
              }
            }
          }
        }

        if (!improved) break
      }
    }
    list(state = state, fit = current_fit)
  }

  value_space <- build_value_space(value_mat, female_pool_size, male_pool_size)
  parent_space <- build_parent_value_space(female_value_df, male_value_df, female_pool_size, male_pool_size)
  female_space_union <- union(value_space$female_space, parent_space$female_space)
  male_space_union <- union(value_space$male_space, parent_space$male_space)
  search_mode <- choose_search_mode(n_female_select, n_male_select, length(female_space_union), length(male_space_union))

  set.seed(seed)
  init_list <- vector("list", n_restart)
  init_list[[1]] <- greedy_init_in_space(value_mat, value_space$female_space, value_space$male_space, n_female_select, n_male_select)
  if (n_restart >= 2) {
    init_list[[2]] <- greedy_init_in_space(value_mat, parent_space$female_space, parent_space$male_space, n_female_select, n_male_select)
  }
  if (n_restart >= 3) {
    for (r in 3:n_restart) {
      if (r %% 2 == 1) {
        init_list[[r]] <- random_init_in_space(value_space$female_space, value_space$male_space, n_female_select, n_male_select, seed + r)
      } else {
        init_list[[r]] <- random_init_in_space(parent_space$female_space, parent_space$male_space, n_female_select, n_male_select, seed + r)
      }
    }
  }

  best_ans <- NULL
  best_obj <- -Inf
  for (r in seq_len(n_restart)) {
    ans_r <- hill_climb_from_state(init_list[[r]], female_space_union, male_space_union, search_mode)
    if (ans_r$fit$objective > best_obj + 1e-12) {
      best_ans <- ans_r
      best_obj <- ans_r$fit$objective
    }
  }

  selected_assignment <- which(best_ans$fit$assignment == 1, arr.ind = TRUE)
  value_sub <- value_mat[best_ans$state$female, best_ans$state$male, drop = FALSE]
  pair_df <- data.frame(
    Dam = rownames(best_ans$fit$assignment)[selected_assignment[, 1]],
    Sire = colnames(best_ans$fit$assignment)[selected_assignment[, 2]],
    ExpectedProgenyMean = value_sub[selected_assignment],
    stringsAsFactors = FALSE
  )

  list(
    selected_females = best_ans$state$female,
    selected_males = best_ans$state$male,
    assignment = best_ans$fit$assignment,
    pair_df = pair_df,
    objective = best_ans$fit$objective,
    mean_pair_progeny = mean(pair_df$ExpectedProgenyMean)
  )
}

make_balanced_pairs <- function(selected_females, selected_males, female_score, male_score) {
  # PS and GS represent selection-only baselines. After parent selection,
  # mating is randomized while keeping sire usage approximately balanced.
  female_ord <- sample(selected_females, size = length(selected_females), replace = FALSE)
  male_ord <- sample(selected_males, size = length(selected_males), replace = FALSE)
  sire_use <- rep(male_ord, length.out = length(female_ord))
  sire_use <- sample(sire_use, size = length(female_ord), replace = FALSE)
  data.frame(
    Dam = female_ord,
    Sire = sire_use,
    stringsAsFactors = FALSE
  )
}

simulate_offspring_generation <- function(
    current_dat,
    crosses,
    alpha_true,
    delta_true,
    sigma_e,
    mu,
    generation_index,
    next_id_start,
    n_offspring_per_cross = 20
) {
  snp_cols <- grep("^SNP", names(current_dat))
  m <- length(snp_cols)
  offspring_list <- vector("list", nrow(crosses))
  next_id <- next_id_start

  for (c in seq_len(nrow(crosses))) {
    dam_id <- crosses$Dam[c]
    sire_id <- crosses$Sire[c]
    dam_geno <- as.numeric(current_dat[current_dat$ID == dam_id, snp_cols, drop = TRUE])
    sire_geno <- as.numeric(current_dat[current_dat$ID == sire_id, snp_cols, drop = TRUE])
    p_dam <- dam_geno / 2
    p_sire <- sire_geno / 2

    n_off <- n_offspring_per_cross
    maternal <- matrix(rbinom(n_off * m, size = 1, prob = rep(p_dam, each = n_off)), nrow = n_off, ncol = m)
    paternal <- matrix(rbinom(n_off * m, size = 1, prob = rep(p_sire, each = n_off)), nrow = n_off, ncol = m)
    G_off <- maternal + paternal
    W_off <- ifelse(G_off == 1, 1, 0)
    g_true_off <- as.numeric(G_off %*% alpha_true + W_off %*% delta_true)
    phenotype_off <- as.numeric(mu + g_true_off + rnorm(n_off, 0, sigma_e))
    off_ids <- paste0("G", sprintf("%02d", generation_index), "_ID", sprintf("%05d", next_id:(next_id + n_off - 1)))

    n_f <- ceiling(n_off / 2)
    n_m <- n_off - n_f
    off_sex <- c(rep("F", n_f), rep("M", n_m))

    off_df <- data.frame(
      ID = off_ids,
      Sex = off_sex,
      G_off,
      Phenotype = phenotype_off,
      TrueGV = g_true_off,
      Generation = generation_index,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    names(off_df)[3:(2 + m)] <- names(current_dat)[snp_cols]

    offspring_list[[c]] <- off_df
    next_id <- next_id + n_off
  }

  list(
    dat = do.call(rbind, offspring_list),
    next_id = next_id
  )
}

run_multigeneration_strategy <- function(
    founder_dat,
    alpha_true,
    delta_true,
    sigma_e,
    mu,
    strategy = c("PS mating", "GS mating", "GS+GM", "GM"),
    n_generations = 10,
    final_females = 12,
    final_males = 6,
    gs_gm_preselect_female = 36,
    gs_gm_preselect_male = 18,
    gm_female_pool_size = Inf,
    gm_male_pool_size = Inf,
    seed = 2026,
    n_restart = 2,
    n_offspring_per_cross = 20
) {
  strategy <- match.arg(strategy)
  current_dat <- founder_dat
  current_dat$TrueGV <- compute_true_scores(current_dat, alpha_true, delta_true)
  current_dat$Generation <- 0

  history <- data.frame(
    Strategy = strategy,
    Generation = 0,
    MeanTrueGV = mean(current_dat$TrueGV),
    MaxTrueGV = max(current_dat$TrueGV),
    MeanPhenotype = mean(current_dat$Phenotype),
    TotalProgenyMeanObjective = NA_real_,
    MeanSelectedPairProgenyMean = NA_real_,
    stringsAsFactors = FALSE
  )

  pair_history <- data.frame()
  next_id <- 1L

  for (gen in seq_len(n_generations)) {
    matrices <- build_value_matrices_from_truth(current_dat, alpha_true, delta_true)

    parent_value_df <- matrices$parent_value_df
    female_value_df <- subset(parent_value_df, Sex == "F")
    male_value_df <- subset(parent_value_df, Sex == "M")
    female_pheno <- setNames(current_dat$Phenotype[current_dat$Sex == "F"], current_dat$ID[current_dat$Sex == "F"])
    male_pheno <- setNames(current_dat$Phenotype[current_dat$Sex == "M"], current_dat$ID[current_dat$Sex == "M"])
    female_gv <- setNames(female_value_df$GeneticValue, female_value_df$ID)
    male_gv <- setNames(male_value_df$GeneticValue, male_value_df$ID)

    if (strategy == "PS mating") {
      selected_females <- names(sort(female_pheno, decreasing = TRUE))[seq_len(final_females)]
      selected_males <- names(sort(male_pheno, decreasing = TRUE))[seq_len(final_males)]
      crosses <- make_balanced_pairs(selected_females, selected_males, female_pheno, male_pheno)
      pair_mean <- matrices$value_mat[cbind(crosses$Dam, crosses$Sire)]
      objective_value <- sum(pair_mean)
      mean_pair_value <- mean(pair_mean)
    } else if (strategy == "GS mating") {
      selected_females <- names(sort(female_gv, decreasing = TRUE))[seq_len(final_females)]
      selected_males <- names(sort(male_gv, decreasing = TRUE))[seq_len(final_males)]
      crosses <- make_balanced_pairs(selected_females, selected_males, female_gv, male_gv)
      pair_mean <- matrices$value_mat[cbind(crosses$Dam, crosses$Sire)]
      objective_value <- sum(pair_mean)
      mean_pair_value <- mean(pair_mean)
    } else if (strategy == "GS+GM") {
      pre_f <- names(sort(female_gv, decreasing = TRUE))[seq_len(min(gs_gm_preselect_female, length(female_gv)))]
      pre_m <- names(sort(male_gv, decreasing = TRUE))[seq_len(min(gs_gm_preselect_male, length(male_gv)))]
      value_sub <- matrices$value_mat[pre_f, pre_m, drop = FALSE]
      female_val_sub <- female_value_df[female_value_df$ID %in% pre_f, , drop = FALSE]
      male_val_sub <- male_value_df[male_value_df$ID %in% pre_m, , drop = FALSE]
      opt_ans <- run_proposed_optimizer(
        value_mat = value_sub,
        female_value_df = female_val_sub,
        male_value_df = male_val_sub,
        n_female_select = final_females,
        n_male_select = final_males,
        female_pool_size = length(pre_f),
        male_pool_size = length(pre_m),
        n_restart = n_restart,
        seed = seed + gen
      )
      crosses <- opt_ans$pair_df[, c("Dam", "Sire"), drop = FALSE]
      pair_mean <- opt_ans$pair_df$ExpectedProgenyMean
      objective_value <- opt_ans$objective
      mean_pair_value <- opt_ans$mean_pair_progeny
    } else {
      gm_f_pool <- min(if (is.infinite(gm_female_pool_size)) nrow(matrices$value_mat) else gm_female_pool_size, nrow(matrices$value_mat))
      gm_m_pool <- min(if (is.infinite(gm_male_pool_size)) ncol(matrices$value_mat) else gm_male_pool_size, ncol(matrices$value_mat))
      opt_ans <- run_proposed_optimizer(
        value_mat = matrices$value_mat,
        female_value_df = female_value_df,
        male_value_df = male_value_df,
        n_female_select = final_females,
        n_male_select = final_males,
        female_pool_size = gm_f_pool,
        male_pool_size = gm_m_pool,
        n_restart = n_restart,
        seed = seed + gen
      )
      crosses <- opt_ans$pair_df[, c("Dam", "Sire"), drop = FALSE]
      pair_mean <- opt_ans$pair_df$ExpectedProgenyMean
      objective_value <- opt_ans$objective
      mean_pair_value <- opt_ans$mean_pair_progeny
    }

    pair_history <- rbind(
      pair_history,
      data.frame(
        Strategy = strategy,
        Generation = gen,
        Dam = crosses$Dam,
        Sire = crosses$Sire,
        ExpectedProgenyMean = as.numeric(pair_mean),
        stringsAsFactors = FALSE
      )
    )

    offspring <- simulate_offspring_generation(
      current_dat = current_dat,
      crosses = crosses,
      alpha_true = alpha_true,
      delta_true = delta_true,
      sigma_e = sigma_e,
      mu = mu,
      generation_index = gen,
      next_id_start = next_id,
      n_offspring_per_cross = n_offspring_per_cross
    )
    next_id <- offspring$next_id
    current_dat <- offspring$dat

    history <- rbind(
      history,
      data.frame(
        Strategy = strategy,
        Generation = gen,
        MeanTrueGV = mean(current_dat$TrueGV),
        MaxTrueGV = max(current_dat$TrueGV),
        MeanPhenotype = mean(current_dat$Phenotype),
        TotalProgenyMeanObjective = objective_value,
        MeanSelectedPairProgenyMean = mean_pair_value,
        stringsAsFactors = FALSE
      )
    )
  }

  list(
    history = history,
    pair_history = pair_history,
    final_population = current_dat
  )
}

# ---------------------------------------------------------
# 3. Run replicate multi-generation comparisons
# ---------------------------------------------------------
suppressPackageStartupMessages({
  library(parallel)
})

out_dir <- get_chr("OUT_DIR", "multigen_expected_mean_large_results")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

history_out <- file.path(out_dir, "multigen_expected_mean_large_history.csv")
pair_out <- file.path(out_dir, "multigen_expected_mean_large_pair_history.csv")
summary_out <- file.path(out_dir, "multigen_expected_mean_large_final_summary.csv")
gen_summary_out <- file.path(out_dir, "multigen_expected_mean_large_generation_summary.csv")
truth_out <- file.path(out_dir, "multigen_expected_mean_large_truth_summary.csv")
design_out <- file.path(out_dir, "multigen_expected_mean_large_design.csv")
plot_out <- file.path(out_dir, "multigen_expected_mean_large_gain_plot.png")

n_replicates <- get_int("N_REP", 100L)
n_generations <- get_int("N_GENERATIONS", 10L)
n_founder_female <- get_int("N_FEMALE", 100L)
n_founder_male <- get_int("N_MALE", 100L)
n_snp <- get_int("N_SNP", 1000L)
n_major <- get_int("N_QTL", 200L)
final_females <- get_int("N_FEMALE_SELECT", 12L)
final_males <- get_int("N_MALE_SELECT", 6L)
gs_gm_preselect_female <- get_int("GS_GM_PRESELECT_FEMALE", 36L)
gs_gm_preselect_male <- get_int("GS_GM_PRESELECT_MALE", 18L)
n_restart <- get_int("N_RESTART", 2L)
n_offspring_per_cross <- get_int("N_OFFSPRING_PER_CROSS", 20L)
H2 <- get_num("H2", 0.50)
mu <- get_num("MU", 50)
scenario <- get_chr("SCENARIO", "S3_strong")
n_workers <- get_int("N_WORKERS", max(1L, min(6L, detectCores(logical = FALSE))))

write.csv(
  data.frame(
    Parameter = c(
      "Scenario", "N_REP", "N_GENERATIONS", "N_FEMALE", "N_MALE",
      "N_SNP", "N_QTL", "N_FEMALE_SELECT", "N_MALE_SELECT",
      "GS_GM_PRESELECT_FEMALE", "GS_GM_PRESELECT_MALE",
      "N_RESTART", "N_OFFSPRING_PER_CROSS", "H2", "MU",
      "Objective"
    ),
    Value = c(
      scenario, n_replicates, n_generations, n_founder_female, n_founder_male,
      n_snp, n_major, final_females, final_males,
      gs_gm_preselect_female, gs_gm_preselect_male,
      n_restart, n_offspring_per_cross, H2, mu,
      "Expected progeny mean only"
    ),
    stringsAsFactors = FALSE
  ),
  design_out,
  row.names = FALSE
)

run_one_replicate <- function(rep_id) {
  base_seed <- 100000 + rep_id * 100
  founder <- generate_large_gm_population(
    scenario = scenario,
    seed = base_seed,
    n_female = n_founder_female,
    n_male = n_founder_male,
    m = n_snp,
    n_major = n_major,
    H2 = H2,
    mu = mu
  )

  founder_dat <- founder$dat

  common_args <- list(
    founder_dat = founder_dat,
    alpha_true = founder$alpha_true,
    delta_true = founder$delta_true,
    sigma_e = founder$sigma_e,
    mu = founder$mu,
    n_generations = n_generations,
    final_females = final_females,
    final_males = final_males,
    gs_gm_preselect_female = gs_gm_preselect_female,
    gs_gm_preselect_male = gs_gm_preselect_male,
    n_offspring_per_cross = n_offspring_per_cross
  )

  res_ps <- do.call(run_multigeneration_strategy, c(common_args, list(strategy = "PS mating", seed = base_seed + 1, n_restart = 1)))
  res_gs <- do.call(run_multigeneration_strategy, c(common_args, list(strategy = "GS mating", seed = base_seed + 1001, n_restart = 1)))
  res_gsgm <- do.call(run_multigeneration_strategy, c(common_args, list(strategy = "GS+GM", seed = base_seed + 2001, n_restart = n_restart)))
  res_gm <- do.call(run_multigeneration_strategy, c(common_args, list(strategy = "GM", seed = base_seed + 3001, n_restart = n_restart)))

  history_rep <- rbind(
    res_ps$history,
    res_gs$history,
    res_gsgm$history,
    res_gm$history
  )

  pair_rep <- rbind(
    res_ps$pair_history,
    res_gs$pair_history,
    res_gsgm$pair_history,
    res_gm$pair_history
  )

  history_rep$Replicate <- rep_id
  pair_rep$Replicate <- rep_id

  list(
    history = history_rep,
    pair_history = pair_rep,
    truth = founder$truth_summary
  )
}

rep_ids <- seq_len(n_replicates)

worker_objects <- c(
  "generate_large_gm_population",
  "compute_true_scores",
  "build_value_matrices_from_truth",
  "solve_mating_ilp_general",
  "as_state",
  "build_value_space",
  "build_parent_value_space",
  "greedy_init_in_space",
  "random_init_in_space",
  "choose_search_mode",
  "run_proposed_optimizer",
  "make_balanced_pairs",
  "simulate_offspring_generation",
  "run_multigeneration_strategy",
  "run_one_replicate",
  "n_replicates", "n_generations", "n_founder_female", "n_founder_male",
  "n_snp", "n_major", "final_females", "final_males",
  "gs_gm_preselect_female", "gs_gm_preselect_male", "n_restart",
  "n_offspring_per_cross", "H2", "mu", "scenario"
)

if (n_workers > 1L) {
  cl <- makeCluster(n_workers)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(lpSolve))
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      BLIS_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1",
      NUMEXPR_NUM_THREADS = "1"
    )
    NULL
  })
  clusterExport(cl, varlist = worker_objects, envir = environment())
  result_list <- parLapplyLB(cl, rep_ids, run_one_replicate)
} else {
  result_list <- lapply(rep_ids, run_one_replicate)
}

history_all <- do.call(rbind, lapply(result_list, `[[`, "history"))
pair_all <- do.call(rbind, lapply(result_list, `[[`, "pair_history"))
truth_all <- do.call(rbind, lapply(result_list, `[[`, "truth"))
truth_all$Replicate <- rep_ids

baseline <- aggregate(
  MeanTrueGV ~ Strategy + Replicate,
  data = history_all[history_all$Generation == 0, ],
  FUN = mean
)
names(baseline)[names(baseline) == "MeanTrueGV"] <- "BaselineGV"
history_all <- merge(history_all, baseline, by = c("Strategy", "Replicate"), all.x = TRUE)
history_all$GeneticGain <- history_all$MeanTrueGV - history_all$BaselineGV
history_all <- history_all[order(history_all$Replicate, history_all$Strategy, history_all$Generation), ]
row.names(history_all) <- NULL

summary_generation <- do.call(
  rbind,
  lapply(split(history_all, list(history_all$Strategy, history_all$Generation), drop = TRUE), function(df) {
    data.frame(
      Strategy = df$Strategy[1],
      Generation = df$Generation[1],
      MeanGeneticGain = mean(df$GeneticGain),
      SDGeneticGain = sd(df$GeneticGain),
      MeanTrueGV = mean(df$MeanTrueGV),
      SDTrueGV = sd(df$MeanTrueGV),
      MeanPhenotype = mean(df$MeanPhenotype),
      SDPhenotype = sd(df$MeanPhenotype),
      MeanSelectedPairProgenyMean = mean(df$MeanSelectedPairProgenyMean, na.rm = TRUE),
      SDSelectedPairProgenyMean = sd(df$MeanSelectedPairProgenyMean, na.rm = TRUE),
      MeanTotalProgenyMeanObjective = mean(df$TotalProgenyMeanObjective, na.rm = TRUE),
      SDTotalProgenyMeanObjective = sd(df$TotalProgenyMeanObjective, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
summary_generation <- summary_generation[order(summary_generation$Strategy, summary_generation$Generation), ]
row.names(summary_generation) <- NULL

final_df <- history_all[history_all$Generation == max(history_all$Generation), ]
summary_final <- do.call(
  rbind,
  lapply(split(final_df, final_df$Strategy), function(df) {
    data.frame(
      Strategy = df$Strategy[1],
      FinalGeneration = max(df$Generation),
      MeanFinalTrueGV = mean(df$MeanTrueGV),
      SDFinalTrueGV = sd(df$MeanTrueGV),
      MeanFinalMaxTrueGV = mean(df$MaxTrueGV),
      MeanFinalPhenotype = mean(df$MeanPhenotype),
      MeanFinalSelectedPairProgenyMean = mean(df$MeanSelectedPairProgenyMean),
      MeanFinalTotalProgenyMeanObjective = mean(df$TotalProgenyMeanObjective),
      MeanGeneticGain = mean(df$GeneticGain),
      SDGeneticGain = sd(df$GeneticGain),
      stringsAsFactors = FALSE
    )
  })
)
summary_final <- summary_final[order(-summary_final$MeanFinalTrueGV), ]
row.names(summary_final) <- NULL

write.csv(history_all, file = history_out, row.names = FALSE)
write.csv(pair_all, file = pair_out, row.names = FALSE)
write.csv(truth_all, file = truth_out, row.names = FALSE)
write.csv(summary_generation, file = gen_summary_out, row.names = FALSE)
write.csv(summary_final, file = summary_out, row.names = FALSE)

strategy_levels <- c("PS mating", "GS mating", "GS+GM", "GM")
strategy_cols <- c(
  "PS mating" = "gray45",
  "GS mating" = "dodgerblue3",
  "GS+GM" = "forestgreen",
  "GM" = "red3"
)
strategy_lty <- c(
  "PS mating" = 2,
  "GS mating" = 3,
  "GS+GM" = 4,
  "GM" = 1
)

png(plot_out, width = 1200, height = 800, res = 140)
par(mar = c(5, 5, 2, 2))
x_rng <- range(history_all$Generation, na.rm = TRUE)
y_rng <- range(history_all$GeneticGain, na.rm = TRUE)
plot(
  NA, NA,
  xlim = x_rng,
  ylim = y_rng,
  xlab = "Generation",
  ylab = "Genetic gain",
  cex.lab = 1.3,
  cex.axis = 1.1
)

for (st in strategy_levels) {
  sub_m <- summary_generation[summary_generation$Strategy == st, ]
  lines(
    sub_m$Generation,
    sub_m$MeanGeneticGain,
    col = strategy_cols[[st]],
    lwd = 3,
    lty = strategy_lty[[st]]
  )
  points(
    sub_m$Generation,
    sub_m$MeanGeneticGain,
    col = strategy_cols[[st]],
    pch = 16,
    cex = 0.8
  )
}

legend(
  "topleft",
  legend = strategy_levels,
  col = unname(strategy_cols[strategy_levels]),
  lty = unname(strategy_lty[strategy_levels]),
  lwd = 3,
  bty = "n",
  cex = 1.1
)
dev.off()

tar_file <- paste0(out_dir, ".tar.gz")
if (file.exists(tar_file)) unlink(tar_file)
utils::tar(tarfile = tar_file, files = out_dir, compression = "gzip")

cat("\n=== Multi-generation expected-progeny-mean final summary ===\n")
print(summary_final)
cat("\nFiles written to:\n")
cat(history_out, "\n")
cat(pair_out, "\n")
cat(truth_out, "\n")
cat(gen_summary_out, "\n")
cat(summary_out, "\n")
cat(plot_out, "\n")
cat(tar_file, "\n")
