# Genomic mating pipeline under an additive-dominance model.
#
# Required input:
#   dat: data.frame containing ID, Sex, Phenotype, and SNP columns named SNP...
#        Sex must be coded as "F" and "M"; SNP genotypes must be coded 0/1/2.
#
# Marker-effect estimation:
#   effect_estimator = "bglr" uses BGLR::BGLR with two marker-effect components:
#       additive  = SNP genotype matrix Z
#       dominance = heterozygote indicator matrix W
#   effect_estimator = "em_reml" uses the internal EM-REML estimator.
#
# Primary optimization target:
#   Expected progeny value, defined as the expected progeny mean for each
#   female-by-male mating combination.
#
# Default optimization settings are chosen for routine empirical analyses:
#   effect_estimator = "bglr", n_restart = 2,
#   female_pool_size = 200, male_pool_size = 200.
run_genomic_mating_pipeline <- function(
    dat,
    n_female_select,
    n_male_select,
    female_pool_size = 200,
    male_pool_size = 200,
    n_restart = 2,
    seed = 2026,
    force_use_all_males = TRUE,
    search_mode = c("best", "first", "auto"),
    threshold = 1e-6,
    maxit = 200,
    verbose = TRUE,
    effect_estimator = c("bglr", "em_reml"),
    bglr_model = "BRR",
    bglr_nIter = 12000,
    bglr_burnIn = 4000,
    bglr_thin = 5,
    bglr_verbose = FALSE,
    bglr_saveAt = NULL
) {
  # Match and validate the local search mode
  search_mode <- match.arg(search_mode)
  if (length(effect_estimator) == 1 && effect_estimator == "em-reml") {
    effect_estimator <- "em_reml"
  }
  effect_estimator <- match.arg(effect_estimator)
  
  # Basic input checks
  if (!all(c("ID", "Sex", "Phenotype") %in% names(dat))) {
    stop("dat must contain ID, Sex, and Phenotype columns.")
  }
  
  dat$ID <- as.character(dat$ID)
  dat$Sex <- as.character(dat$Sex)
  
  if (anyDuplicated(dat$ID)) {
    stop("ID must be unique.")
  }
  
  if (!all(c("F", "M") %in% unique(dat$Sex))) {
    stop("Sex column must contain both 'F' and 'M'.")
  }
  
  if (n_restart < 1) {
    stop("n_restart must be at least 1.")
  }
  
  if (sum(dat$Sex == "F") < n_female_select) {
    stop("Not enough female candidates.")
  }
  
  if (sum(dat$Sex == "M") < n_male_select) {
    stop("Not enough male candidates.")
  }
  
  if (!is.numeric(dat$Phenotype)) {
    stop("Phenotype must be numeric.")
  }

  # Identify SNP columns by name
  snp_cols <- grep("^SNP", names(dat))
  if (length(snp_cols) == 0) {
    stop("No SNP columns were found. SNP column names must start with 'SNP'.")
  }
  
  snp_check <- as.matrix(dat[, snp_cols, drop = FALSE])
  suppressWarnings(storage.mode(snp_check) <- "numeric")
  if (anyNA(snp_check)) {
    stop("SNP genotype columns must be numeric or coercible to numeric and contain no missing values.")
  }
  if (!all(snp_check %in% c(0, 1, 2))) {
    stop("SNP genotype columns must be coded as 0, 1, or 2.")
  }
  
  # The search space must be at least as large as the final selected set
  if (female_pool_size < n_female_select) {
    stop("female_pool_size must be greater than or equal to n_female_select.")
  }
  
  if (male_pool_size < n_male_select) {
    stop("male_pool_size must be greater than or equal to n_male_select.")
  }
  
  # ---------------------------------------------------------
  # Solve the lower-level mating assignment by integer programming
  # ---------------------------------------------------------
  solve_mating_ilp_general <- function(value_mat_sub, force_use_all_males = TRUE) {
    nf <- nrow(value_mat_sub)
    nm <- ncol(value_mat_sub)
    
    # Objective coefficients
    obj <- as.vector(value_mat_sub)
    
    # Each female must be assigned exactly once
    mother_block <- matrix(0, nrow = nf, ncol = nf * nm)
    for (i in seq_len(nf)) {
      mother_block[i, i + (0:(nm - 1)) * nf] <- 1
    }
    
    # Track how many females are assigned to each male
    sire_block <- matrix(0, nrow = nm, ncol = nf * nm)
    for (j in seq_len(nm)) {
      sire_block[j, ((j - 1) * nf + 1):(j * nf)] <- 1
    }
    
    # Allow balanced sire usage even when nf is not divisible by nm
    lower <- if (force_use_all_males && nf >= nm) floor(nf / nm) else 0
    upper <- ceiling(nf / nm)
    
    const.mat <- rbind(mother_block, sire_block, sire_block)
    const.dir <- c(rep("=", nf), rep(">=", nm), rep("<=", nm))
    const.rhs <- c(rep(1, nf), rep(lower, nm), rep(upper, nm))
    
    fit <- lpSolve::lp(
      direction = "max",
      objective.in = obj,
      const.mat = const.mat,
      const.dir = const.dir,
      const.rhs = const.rhs,
      all.bin = TRUE
    )
    
    if (fit$status != 0) {
      stop("ILP failed. Check female/male counts or constraints.")
    }
    
    assignment <- matrix(
      fit$solution,
      nrow = nf,
      ncol = nm,
      dimnames = list(rownames(value_mat_sub), colnames(value_mat_sub))
    )
    
    list(
      objective = fit$objval,
      assignment = assignment,
      sire_load = colSums(assignment),
      status = fit$status
    )
  }
  
  # ---------------------------------------------------------
  # Automatically choose best-improvement or first-improvement
  # ---------------------------------------------------------
  choose_search_mode <- function(
    n_female_select,
    n_male_select,
    female_space_size,
    male_space_size,
    neighborhood_threshold = 1500,
    ilp_var_threshold = 200,
    workload_threshold = 200000
  ) {
    # Total number of 1-swap neighbors
    neighborhood_count <-
      n_female_select * max(0, female_space_size - n_female_select) +
      n_male_select * max(0, male_space_size - n_male_select)
    
    # ILP size for one evaluation
    ilp_vars <- n_female_select * n_male_select
    
    # Rough computational workload
    workload <- neighborhood_count * ilp_vars
    
    mode <- if (neighborhood_count <= neighborhood_threshold &&
                ilp_vars <= ilp_var_threshold &&
                workload <= workload_threshold) {
      "best"
    } else {
      "first"
    }
    
    list(
      mode = mode,
      neighborhood_count = neighborhood_count,
      ilp_vars = ilp_vars,
      workload = workload
    )
  }
  
  # ---------------------------------------------------------
  # Build value-based candidate space from row/column means of value_mat
  # ---------------------------------------------------------
  build_value_space <- function(value_mat, female_pool_size, male_pool_size) {
    female_pool_size <- min(female_pool_size, nrow(value_mat))
    male_pool_size <- min(male_pool_size, ncol(value_mat))
    
    female_score <- rowMeans(value_mat)
    male_score <- colMeans(value_mat)
    
    female_space <- names(sort(female_score, decreasing = TRUE))[seq_len(female_pool_size)]
    male_space <- names(sort(male_score, decreasing = TRUE))[seq_len(male_pool_size)]
    
    list(
      female_space = female_space,
      male_space = male_space
    )
  }
  
  # ---------------------------------------------------------
  # Build parent-based candidate space from parent genetic values
  # ---------------------------------------------------------
  build_parent_value_space <- function(female_value_df, male_value_df, female_pool_size, male_pool_size) {
    female_pool_size <- min(female_pool_size, nrow(female_value_df))
    male_pool_size <- min(male_pool_size, nrow(male_value_df))
    
    female_ord <- female_value_df[order(female_value_df$GeneticValue, decreasing = TRUE), ]
    male_ord <- male_value_df[order(male_value_df$GeneticValue, decreasing = TRUE), ]
    
    female_space <- female_ord$ID[seq_len(female_pool_size)]
    male_space <- male_ord$ID[seq_len(male_pool_size)]
    
    list(
      female_space = female_space,
      male_space = male_space
    )
  }
  
  # ---------------------------------------------------------
  # Greedy initialization within a candidate space
  # ---------------------------------------------------------
  greedy_init_in_space <- function(value_mat, female_space, male_space, n_female_select, n_male_select) {
    submat <- value_mat[female_space, male_space, drop = FALSE]
    
    male_score <- colMeans(submat)
    selected_males <- names(sort(male_score, decreasing = TRUE))[seq_len(n_male_select)]
    
    female_score <- rowMeans(submat[, selected_males, drop = FALSE])
    selected_females <- names(sort(female_score, decreasing = TRUE))[seq_len(n_female_select)]
    
    list(
      selected_females = selected_females,
      selected_males = selected_males
    )
  }
  
  # ---------------------------------------------------------
  # Random initialization within a candidate space
  # ---------------------------------------------------------
  random_init_in_space <- function(female_space, male_space, n_female_select, n_male_select) {
    list(
      selected_females = sample(female_space, n_female_select, replace = FALSE),
      selected_males = sample(male_space, n_male_select, replace = FALSE)
    )
  }
  
  # ---------------------------------------------------------
  # Local search on parent sets using 1-swap neighborhood
  # ---------------------------------------------------------
  local_search_with_space <- function(
    value_mat,
    current_females,
    current_males,
    female_space,
    male_space,
    force_use_all_males = TRUE,
    mode = c("auto", "best", "first"),
    verbose = FALSE
  ) {
    mode <- match.arg(mode)
    
    auto_info <- choose_search_mode(
      n_female_select = length(current_females),
      n_male_select = length(current_males),
      female_space_size = length(female_space),
      male_space_size = length(male_space)
    )
    
    search_mode_local <- if (mode == "auto") auto_info$mode else mode
    
    current_sub <- value_mat[current_females, current_males, drop = FALSE]
    current_fit <- solve_mating_ilp_general(current_sub, force_use_all_males)
    
    if (verbose) {
      cat(
        "Local search mode =", search_mode_local,
        "| neighborhood =", auto_info$neighborhood_count,
        "| ilp_vars =", auto_info$ilp_vars,
        "| workload =", auto_info$workload, "\n"
      )
    }
    
    repeat {
      fixed_females <- current_females
      fixed_males <- current_males
      
      unselected_females <- setdiff(female_space, fixed_females)
      unselected_males <- setdiff(male_space, fixed_males)
      
      if (search_mode_local == "best") {
        # Best-improvement: evaluate all neighbors and take the best one
        best_obj <- current_fit$objective
        best_females <- fixed_females
        best_males <- fixed_males
        best_fit <- current_fit
        
        for (f_out in fixed_females) {
          for (f_in in unselected_females) {
            candidate_females <- c(setdiff(fixed_females, f_out), f_in)
            candidate_sub <- value_mat[candidate_females, fixed_males, drop = FALSE]
            candidate_fit <- solve_mating_ilp_general(candidate_sub, force_use_all_males)
            
            if (candidate_fit$objective > best_obj) {
              best_obj <- candidate_fit$objective
              best_females <- candidate_females
              best_males <- fixed_males
              best_fit <- candidate_fit
            }
          }
        }
        
        for (m_out in fixed_males) {
          for (m_in in unselected_males) {
            candidate_males <- c(setdiff(fixed_males, m_out), m_in)
            candidate_sub <- value_mat[fixed_females, candidate_males, drop = FALSE]
            candidate_fit <- solve_mating_ilp_general(candidate_sub, force_use_all_males)
            
            if (candidate_fit$objective > best_obj) {
              best_obj <- candidate_fit$objective
              best_females <- fixed_females
              best_males <- candidate_males
              best_fit <- candidate_fit
            }
          }
        }
        
        if (best_obj <= current_fit$objective) {
          break
        }
        
        current_females <- best_females
        current_males <- best_males
        current_fit <- best_fit
        
      } else {
        # First-improvement: accept the first improving neighbor found
        improved <- FALSE
        
        for (f_out in fixed_females) {
          if (improved) break
          for (f_in in unselected_females) {
            candidate_females <- c(setdiff(fixed_females, f_out), f_in)
            candidate_sub <- value_mat[candidate_females, fixed_males, drop = FALSE]
            candidate_fit <- solve_mating_ilp_general(candidate_sub, force_use_all_males)
            
            if (candidate_fit$objective > current_fit$objective) {
              current_females <- candidate_females
              current_males <- fixed_males
              current_fit <- candidate_fit
              improved <- TRUE
              break
            }
          }
        }
        
        if (!improved) {
          for (m_out in fixed_males) {
            if (improved) break
            for (m_in in unselected_males) {
              candidate_males <- c(setdiff(fixed_males, m_out), m_in)
              candidate_sub <- value_mat[fixed_females, candidate_males, drop = FALSE]
              candidate_fit <- solve_mating_ilp_general(candidate_sub, force_use_all_males)
              
              if (candidate_fit$objective > current_fit$objective) {
                current_females <- fixed_females
                current_males <- candidate_males
                current_fit <- candidate_fit
                improved <- TRUE
                break
              }
            }
          }
        }
        
        if (!improved) {
          break
        }
      }
    }
    
    list(
      selected_females = current_females,
      selected_males = current_males,
      objective = current_fit$objective,
      assignment = current_fit$assignment,
      sire_load = current_fit$sire_load,
      sub_value_mat = value_mat[current_females, current_males, drop = FALSE],
      search_mode = search_mode_local,
      auto_info = auto_info
    )
  }
  
  # ---------------------------------------------------------
  # Multistart GM optimization
  # ---------------------------------------------------------
  gm_multistart_all_runs <- function(
    value_mat,
    female_value_df,
    male_value_df,
    n_female_select,
    n_male_select,
    female_pool_size,
    male_pool_size,
    n_restart,
    seed,
    force_use_all_males,
    search_mode,
    verbose
  ) {
    set.seed(seed)
    
    value_space <- build_value_space(value_mat, female_pool_size, male_pool_size)
    parent_space <- build_parent_value_space(female_value_df, male_value_df, female_pool_size, male_pool_size)
    
    female_space_all <- union(value_space$female_space, parent_space$female_space)
    male_space_all <- union(value_space$male_space, parent_space$male_space)
    
    auto_info_global <- choose_search_mode(
      n_female_select = n_female_select,
      n_male_select = n_male_select,
      female_space_size = length(female_space_all),
      male_space_size = length(male_space_all)
    )
    
    actual_mode <- if (search_mode == "auto") auto_info_global$mode else search_mode
    
    if (verbose) {
      cat("\nGlobal search mode =", actual_mode, "\n")
      cat("Neighborhood count =", auto_info_global$neighborhood_count, "\n")
      cat("ILP vars =", auto_info_global$ilp_vars, "\n")
      cat("Workload =", auto_info_global$workload, "\n\n")
    }
    
    history <- data.frame(
      restart = integer(0),
      source = character(0),
      init_obj = numeric(0),
      final_obj = numeric(0),
      search_mode = character(0),
      stringsAsFactors = FALSE
    )
    
    init_list <- vector("list", n_restart)
    
    # Start 1: greedy from value-space
    init_list[[1]] <- greedy_init_in_space(
      value_mat,
      value_space$female_space,
      value_space$male_space,
      n_female_select,
      n_male_select
    )
    
    # Start 2: greedy from parent-space
    init_list[[2]] <- greedy_init_in_space(
      value_mat,
      parent_space$female_space,
      parent_space$male_space,
      n_female_select,
      n_male_select
    )
    
    # Starts 3...n_restart: random starts alternating between spaces
    if (n_restart >= 3) {
      for (r in 3:n_restart) {
        if (r %% 2 == 1) {
          init_list[[r]] <- random_init_in_space(
            value_space$female_space,
            value_space$male_space,
            n_female_select,
            n_male_select
          )
        } else {
          init_list[[r]] <- random_init_in_space(
            parent_space$female_space,
            parent_space$male_space,
            n_female_select,
            n_male_select
          )
        }
      }
    }
    
    best_ans <- NULL
    best_final_obj <- -Inf
    best_run_source <- NULL
    
    for (r in seq_len(n_restart)) {
      init_f <- init_list[[r]]$selected_females
      init_m <- init_list[[r]]$selected_males
      
      init_sub <- value_mat[init_f, init_m, drop = FALSE]
      init_fit <- solve_mating_ilp_general(init_sub, force_use_all_males)
      
      src <- if (r == 1) {
        "value-greedy"
      } else if (r == 2) {
        "parent-greedy"
      } else if (r %% 2 == 1) {
        "value-random"
      } else {
        "parent-random"
      }
      
      ans_r <- local_search_with_space(
        value_mat = value_mat,
        current_females = init_f,
        current_males = init_m,
        female_space = female_space_all,
        male_space = male_space_all,
        force_use_all_males = force_use_all_males,
        mode = actual_mode,
        verbose = FALSE
      )
      
      history <- rbind(history, data.frame(
        restart = r,
        source = src,
        init_obj = init_fit$objective,
        final_obj = ans_r$objective,
        search_mode = ans_r$search_mode,
        stringsAsFactors = FALSE
      ))
      
      if (verbose) {
        cat(
          "Restart", r,
          "source =", src,
          "init objective =", round(init_fit$objective, 4),
          "final objective =", round(ans_r$objective, 4),
          "mode =", ans_r$search_mode, "\n"
        )
      }
      
      if (ans_r$objective > best_final_obj) {
        best_final_obj <- ans_r$objective
        best_ans <- ans_r
        best_run_source <- src
      }
    }
    
    best_ans$history <- history
    best_ans$best_run_source <- best_run_source
    best_ans$best_final_objective <- best_final_obj
    best_ans$global_search_mode <- actual_mode
    best_ans$global_auto_info <- auto_info_global
    
    best_ans
  }
  
  # ---------------------------------------------------------
  # Convert the binary assignment matrix into a readable pair table
  # ---------------------------------------------------------
  extract_pair_result <- function(ans, parent_value_df) {
    idx <- which(ans$assignment == 1, arr.ind = TRUE)
    
    pair_result <- data.frame(
      Dam = rownames(ans$assignment)[idx[, 1]],
      Sire = colnames(ans$assignment)[idx[, 2]],
      ExpectedProgenyValue = ans$sub_value_mat[idx],
      row.names = NULL
    )
    
    pair_result$DamGeneticValue <- parent_value_df$GeneticValue[
      match(pair_result$Dam, parent_value_df$ID)
    ]
    
    pair_result$SireGeneticValue <- parent_value_df$GeneticValue[
      match(pair_result$Sire, parent_value_df$ID)
    ]
    
    pair_result[order(pair_result$ExpectedProgenyValue, decreasing = TRUE), ]
  }
  
  # ---------------------------------------------------------
  # Build phenotype and genotype matrices
  # ---------------------------------------------------------
  y <- dat$Phenotype
  
  X <- matrix(1, nrow = length(y), ncol = 1)
  colnames(X) <- "Intercept"
  
  Z <- as.matrix(dat[, snp_cols])
  mode(Z) <- "numeric"
  W <- ifelse(Z == 1, 1, 0)
  
  if (anyNA(y) || anyNA(Z)) {
    stop("Phenotype and SNP genotype columns must not contain missing values.")
  }
  
  n <- nrow(Z)
  qA <- ncol(Z)
  qD <- ncol(W)
  
  # ---------------------------------------------------------
  # 1-2. Estimate marker additive and dominance effects
  # ---------------------------------------------------------
  if (effect_estimator == "em_reml") {
    K_gamma <- Z %*% t(Z)
    K_delta <- W %*% t(W)
    
    XX <- crossprod(X)
    XZ <- crossprod(X, Z)
    XW <- crossprod(X, W)
    ZX <- t(XZ)
    WX <- t(XW)
    ZZ <- crossprod(Z)
    ZW <- crossprod(Z, W)
    WZ <- t(ZW)
    WW <- crossprod(W)
    
    LHS <- rbind(
      cbind(XX, XZ, XW),
      cbind(ZX, ZZ, ZW),
      cbind(WX, WZ, WW)
    )
    
    Xy <- crossprod(X, y)
    Zy <- crossprod(Z, y)
    Wy <- crossprod(W, y)
    RHS <- rbind(Xy, Zy, Wy)
    
    yy <- crossprod(y)
    N <- length(y)
    rankX <- qr(X)$rank
    nh <- ncol(X)
    
    I_A <- diag(qA)
    I_D <- diag(qD)
    
    kA0 <- 2
    kD0 <- 2
    
    results <- matrix(NA, nrow = 0, ncol = 6)
    colnames(results) <- c("iter", "sigmaE", "sigmaA", "sigmaD", "kA", "kD")
    
    # Estimate variance components by EM-REML
    for (iter in seq_len(maxit)) {
      LHS1 <- LHS
      kA <- kA0
      kD <- kD0
      
      idxA <- (nh + 1):(nh + qA)
      idxD <- (nh + qA + 1):(nh + qA + qD)
      
      LHS1[idxA, idxA] <- LHS1[idxA, idxA] + I_A * kA
      LHS1[idxD, idxD] <- LHS1[idxD, idxD] + I_D * kD
      
      sol <- solve(LHS1, RHS)
      C <- solve(LHS1)
      
      ahat <- sol[idxA, 1, drop = FALSE]
      dhat <- sol[idxD, 1, drop = FALSE]
      
      Caa <- C[idxA, idxA, drop = FALSE]
      Cdd <- C[idxD, idxD, drop = FALSE]
      
      sigmaE <- as.numeric(yy - crossprod(sol, RHS)) / (N - rankX)
      
      aIa <- as.numeric(t(ahat) %*% I_A %*% ahat)
      trCA <- sum(diag(Caa %*% I_A))
      sigmaA <- as.numeric((aIa + trCA * sigmaE) / qA)
      
      dId <- as.numeric(t(dhat) %*% I_D %*% dhat)
      trCD <- sum(diag(Cdd %*% I_D))
      sigmaD <- as.numeric((dId + trCD * sigmaE) / qD)
      
      sigmaA <- max(sigmaA, 1e-8)
      sigmaD <- max(sigmaD, 1e-8)
      
      kA0 <- as.numeric(sigmaE / sigmaA)
      kD0 <- as.numeric(sigmaE / sigmaD)
      
      results <- rbind(results, c(iter, sigmaE, sigmaA, sigmaD, kA0, kD0))
      
      if (max(abs(c(kA - kA0, kD - kD0))) < threshold) {
        break
      }
    }
    
    results <- as.data.frame(results)
    
    V <- sigmaA * K_gamma + sigmaD * K_delta + sigmaE * diag(n)
    V <- V + diag(1e-8, n)
    
    Vinv <- solve(V)
    beta_hat <- solve(t(X) %*% Vinv %*% X, t(X) %*% Vinv %*% y)
    r <- y - X %*% beta_hat
    
    xi_gamma_hat <- sigmaA * K_gamma %*% Vinv %*% r
    xi_delta_hat <- sigmaD * K_delta %*% Vinv %*% r
    
    gamma_hat <- t(Z) %*% MASS::ginv(Z %*% t(Z)) %*% xi_gamma_hat
    delta_hat <- t(W) %*% MASS::ginv(W %*% t(W)) %*% xi_delta_hat
    
    gamma_hat <- as.numeric(gamma_hat)
    delta_hat <- as.numeric(delta_hat)
  } else {
    if (!requireNamespace("BGLR", quietly = TRUE)) {
      stop("effect_estimator = 'bglr' requires the BGLR package. Please run install.packages('BGLR') first.")
    }
    
    if (verbose) {
      cat("\nEstimating marker effects by BGLR...\n")
      cat("BGLR model =", bglr_model, "\n")
      cat("nIter =", bglr_nIter, "burnIn =", bglr_burnIn, "thin =", bglr_thin, "\n")
    }
    
    user_saveAt <- !is.null(bglr_saveAt)
    if (!user_saveAt) {
      bglr_saveAt <- paste0(tempfile(pattern = "GM_BGLR_", tmpdir = tempdir()), "_")
      on.exit(unlink(paste0(bglr_saveAt, "*")), add = TRUE)
    }
    
    set.seed(seed)
    has_dominance_signal <- any(W != 0) && any(apply(W, 2, stats::var) > 0)
    ETA <- list(additive = list(X = Z, model = bglr_model))
    if (has_dominance_signal) {
      ETA$dominance <- list(X = W, model = bglr_model)
    } else if (verbose) {
      cat("Dominance design matrix has no variation; BGLR dominance effects are set to zero.\n")
    }
    
    bglr_fit <- BGLR::BGLR(
      y = y,
      ETA = ETA,
      nIter = bglr_nIter,
      burnIn = bglr_burnIn,
      thin = bglr_thin,
      saveAt = bglr_saveAt,
      verbose = bglr_verbose
    )
    
    gamma_hat <- as.numeric(bglr_fit$ETA[[1]]$b)
    delta_hat <- if (has_dominance_signal) {
      as.numeric(bglr_fit$ETA[[2]]$b)
    } else {
      rep(0, qD)
    }
    
    if (length(gamma_hat) != qA || length(delta_hat) != qD) {
      stop("BGLR did not return marker-effect vectors with the expected length.")
    }
    
    additive_value <- as.numeric(Z %*% gamma_hat)
    dominance_value <- as.numeric(W %*% delta_hat)
    
    sigmaE <- as.numeric(bglr_fit$varE)
    sigmaA <- if (length(unique(additive_value)) > 1) var(additive_value) else 1e-8
    sigmaD <- if (length(unique(dominance_value)) > 1) var(dominance_value) else 1e-8
    sigmaA <- max(sigmaA, 1e-8)
    sigmaD <- max(sigmaD, 1e-8)
    kA0 <- as.numeric(sigmaE / sigmaA)
    kD0 <- as.numeric(sigmaE / sigmaD)
    
    results <- data.frame(
      iter = bglr_nIter,
      sigmaE = sigmaE,
      sigmaA = sigmaA,
      sigmaD = sigmaD,
      kA = kA0,
      kD = kD0
    )
  }
  
  if (verbose) {
    cat("\nMarker-effect estimator =", effect_estimator, "\n")
    cat("Final variance components:\n")
    cat("sigmaE =", sigmaE, "\n")
    cat("sigmaA =", sigmaA, "\n")
    cat("sigmaD =", sigmaD, "\n")
    cat("kA     =", kA0, "\n")
    cat("kD     =", kD0, "\n")
    cat("iterations =", nrow(results), "\n")
  }
  
  names(gamma_hat) <- colnames(Z)
  names(delta_hat) <- colnames(Z)
  
  marker_effects <- data.frame(
    SNP = colnames(Z),
    gamma_hat = gamma_hat,
    delta_hat = delta_hat
  )
  
  # ---------------------------------------------------------
  # 3. Compute parent genetic values
  # ---------------------------------------------------------
  parent_genetic_value <- as.numeric(Z %*% gamma_hat + W %*% delta_hat)
  
  parent_value_df <- data.frame(
    ID = dat$ID,
    Sex = dat$Sex,
    Phenotype = dat$Phenotype,
    GeneticValue = parent_genetic_value
  )
  
  female_value_df <- subset(parent_value_df, Sex == "F")
  male_value_df <- subset(parent_value_df, Sex == "M")
  
  # ---------------------------------------------------------
  # 4. Compute the full expected progeny value matrix
  # ---------------------------------------------------------
  female_idx <- which(dat$Sex == "F")
  male_idx <- which(dat$Sex == "M")
  
  G_female <- Z[female_idx, , drop = FALSE]
  G_male <- Z[male_idx, , drop = FALSE]
  
  female_id <- dat$ID[female_idx]
  male_id <- dat$ID[male_idx]
  
  # Expected additive and dominance coefficients in the progeny
  get_pair_coefficients <- function(sire, dam) {
    if (sire == 0 && dam == 0) return(c(EZ = 0.0, EW = 0.0))
    if (sire == 0 && dam == 1) return(c(EZ = 0.5, EW = 0.5))
    if (sire == 0 && dam == 2) return(c(EZ = 1.0, EW = 1.0))
    if (sire == 1 && dam == 0) return(c(EZ = 0.5, EW = 0.5))
    if (sire == 1 && dam == 1) return(c(EZ = 1.0, EW = 0.5))
    if (sire == 1 && dam == 2) return(c(EZ = 1.5, EW = 0.5))
    if (sire == 2 && dam == 0) return(c(EZ = 1.0, EW = 1.0))
    if (sire == 2 && dam == 1) return(c(EZ = 1.5, EW = 0.5))
    if (sire == 2 && dam == 2) return(c(EZ = 2.0, EW = 0.0))
    stop("Genotypes must be 0, 1, or 2.")
  }
  
  # Expected progeny value for one sire-dam pair
  predict_pair_expected_value <- function(sire_geno, dam_geno, gamma_hat, delta_hat) {
    sire_geno <- as.numeric(sire_geno)
    dam_geno <- as.numeric(dam_geno)
    
    m_local <- length(sire_geno)
    marker_value <- numeric(m_local)
    
    for (k in seq_len(m_local)) {
      tmp <- get_pair_coefficients(sire_geno[k], dam_geno[k])
      marker_value[k] <- tmp["EZ"] * gamma_hat[k] + tmp["EW"] * delta_hat[k]
    }
    
    sum(marker_value)
  }
  
  mean_mat <- matrix(
    NA_real_,
    nrow = nrow(G_female),
    ncol = nrow(G_male),
    dimnames = list(female_id, male_id)
  )

  for (i in seq_len(nrow(G_female))) {
    for (j in seq_len(nrow(G_male))) {
      mean_mat[i, j] <- predict_pair_expected_value(
        sire_geno = G_male[j, ],
        dam_geno = G_female[i, ],
        gamma_hat = gamma_hat,
        delta_hat = delta_hat
      )
    }
  }

  # In this version, genomic mating is optimized only by expected progeny value.
  value_mat <- mean_mat
  
  # ---------------------------------------------------------
  # 5. Run genomic mating optimization
  # ---------------------------------------------------------
  ans_gm <- gm_multistart_all_runs(
    value_mat = value_mat,
    female_value_df = female_value_df,
    male_value_df = male_value_df,
    n_female_select = n_female_select,
    n_male_select = n_male_select,
    female_pool_size = female_pool_size,
    male_pool_size = male_pool_size,
    n_restart = n_restart,
    seed = seed,
    force_use_all_males = force_use_all_males,
    search_mode = search_mode,
    verbose = verbose
  )
  
  pair_result_gm <- extract_pair_result(
    ans_gm,
    parent_value_df
  )
  
  if (verbose) {
    cat("\nGM results:\n")
    cat("Best run source =", ans_gm$best_run_source, "\n")
    cat("Final total expected progeny value =", ans_gm$best_final_objective, "\n")
    cat("GM global search mode =", ans_gm$global_search_mode, "\n")
    cat("GM sire load =\n")
    print(ans_gm$sire_load)
    
    cat("\nGM best pair result:\n")
    print(pair_result_gm)
  }
  
  # ---------------------------------------------------------
  # 6. Return all major outputs
  # ---------------------------------------------------------
  list(
    effect_estimator = effect_estimator,
    bglr_settings = if (effect_estimator == "bglr") {
      list(
        model = bglr_model,
        nIter = bglr_nIter,
        burnIn = bglr_burnIn,
        thin = bglr_thin,
        verbose = bglr_verbose,
        saveAt = bglr_saveAt
      )
    } else {
      NULL
    },
    variance_components = list(
      results = results,
      sigmaE = sigmaE,
      sigmaA = sigmaA,
      sigmaD = sigmaD,
      kA = kA0,
      kD = kD0
    ),
    marker_effects = marker_effects,
    parent_genetic_values = parent_value_df,
    expected_progeny_value_matrix = mean_mat,
    gm_results = list(
      fit_gm = ans_gm,
      pair_result_gm = pair_result_gm
    )
  )
}

