# Optimizer-only AlphaSimR benchmark. No phenotypes or estimated effects.
# The 1,000 observed SNPs include all 200 QTL. Z is allele dosage (0/1/2),
# NOT the centered additive design: gv = (Z - 1) a + I(Z == 1) d + intercept.
# Output is independent of optimizer settings, mating ratios, and call order.
# Founders use the MaCS GENERIC coalescent/demographic model, not independent
# haplotypes. Each simulated chromosome is 1 Morgan; retained marker endpoints
# need not span the whole chromosome. Preserve their irregular MaCS positions.
# Prefer QTL MAF >= 0.05; otherwise use the chromosome's 20 highest-MAF sites
# in the same saved founder, breaking MAF ties by locus order. Never redraw.

.sim_schema <- "AlphaSimR_optimizer_inputs_v3_macs_best_available"

.sim_scalar <- function(x, name, integer = FALSE, lower = -Inf,
                        upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x))) {
    stop(name, " must be a finite ", if (integer) "integer" else "number",
         " in [", lower, ", ", upper, "].", call. = FALSE)
  }
  if (integer) as.integer(x) else as.double(x)
}

validate_simulation_config <- function(cfg = list()) {
  if (!is.list(cfg) || is.data.frame(cfg) ||
      (length(cfg) && (is.null(names(cfg)) || anyNA(names(cfg)) ||
                      any(!nzchar(names(cfg))) || anyDuplicated(names(cfg))))) {
    stop("cfg must be a uniquely named list.", call. = FALSE)
  }
  get <- function(name, default) if (is.null(cfg[[name]])) default else cfg[[name]]
  fixed <- list(
    n_females = c("n_females", "n_female", "nFemale", "nFemales", "nF"),
    n_males = c("n_males", "n_male", "nMale", "nMales", "nM"),
    n_snps = c("n_snps", "n_snp", "nSnp", "nSNP", "nSNPs"),
    n_qtl = c("n_qtl", "n_qtls", "nQtl", "nQTL", "nQTLs"),
    n_chr = c("n_chr", "nChr"),
    n_ind = c("n_ind", "nInd"),
    snps_per_chr = c("snps_per_chr", "nSnpPerChr", "segSites"),
    qtl_per_chr = c("qtl_per_chr", "nQtlPerChr"),
    ploidy = "ploidy", varA = c("varA", "var_A"),
    trait_mean = c("trait_mean", "mean"),
    qtl_maf = c("qtl_maf", "min_qtl_maf", "qtl_min_maf"),
    gen_len = c("gen_len", "genLen"))
  values <- list(n_females = 100L, n_males = 100L, n_snps = 1000L,
                 n_qtl = 200L, n_chr = 10L, n_ind = 200L,
                 snps_per_chr = 100L, qtl_per_chr = 20L, ploidy = 2L,
                 varA = 1, trait_mean = 0, qtl_maf = 0.05, gen_len = 1)
  for (key in names(fixed)) {
    for (alias in fixed[[key]]) {
      if (!is.null(cfg[[alias]])) {
        value <- .sim_scalar(cfg[[alias]], paste0("cfg$", alias))
        if (value != values[[key]]) {
          stop("The benchmark fixes ", alias, " at ", values[[key]], ".",
               call. = FALSE)
        }
      }
    }
  }
  if (!is.null(cfg$useVarA) && !identical(cfg$useVarA, TRUE)) {
    stop("The benchmark requires cfg$useVarA = TRUE.", call. = FALSE)
  }
  fixed_text <- list(founder_method = "runMacs", species = "GENERIC",
                     qtl_maf_policy = "best_available_without_redraw")
  for (key in names(fixed_text)) {
    if (!is.null(cfg[[key]]) && !identical(cfg[[key]], fixed_text[[key]])) {
      stop("The benchmark requires cfg$", key, " = '", fixed_text[[key]], "'.",
           call. = FALSE)
    }
  }
  # Algorithm-specific cfg fields deliberately do not affect generation/hash.
  c(list(schema = .sim_schema,
         seed = .sim_scalar(get("seed", 20261027L), "cfg$seed", TRUE,
                            0, .Machine$integer.max),
         meanDD = .sim_scalar(get("meanDD", 0.5), "cfg$meanDD"),
         varDD = .sim_scalar(get("varDD", 0.1), "cfg$varDD", lower = 0)),
    values, fixed_text, list(useVarA = TRUE,
                 rng = "L'Ecuyer-CMRG/Inversion/Rejection", n_threads = 1L))
}

.sim_hash <- function(object) {
  path <- tempfile("alphasimr_hash_", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, ascii = FALSE, version = 2, compress = FALSE)
  unname(tools::md5sum(path))
}

.sim_slots <- function(object) {
  slots <- methods::slotNames(object)
  stats::setNames(lapply(slots, function(s) methods::slot(object, s)), slots)
}

.sim_restore <- function(class, slots) {
  do.call(methods::new, c(list(Class = class), slots))
}

.sim_publish <- function(path, writer) {
  # Never truncate a checkpoint or replace conflicting data. Failed partials
  # remain beside the destination for diagnosis; only successful ones move.
  partial <- tempfile(paste0(basename(path), ".partial_"), dirname(path))
  writer(partial)
  if (file.exists(path)) {
    hashes <- unname(tools::md5sum(c(path, partial)))
    if (anyNA(hashes) || hashes[1L] != hashes[2L]) {
      stop("Existing artifact differs; preserved both files: ", path,
           " and ", partial, call. = FALSE)
    }
    unlink(partial)
  } else if (!file.rename(partial, path)) {
    stop("Cannot publish artifact; partial preserved: ", partial, call. = FALSE)
  }
  invisible(path)
}

.sim_rds <- function(object, path) {
  .sim_publish(path, function(p) saveRDS(object, p, version = 2, compress = "gzip"))
}

.sim_csv <- function(object, path) {
  .sim_publish(path, function(p) {
    con <- gzfile(p, open = "wt", encoding = "UTF-8", compression = 6)
    on.exit(close(con), add = TRUE)
    utils::write.table(object, con, sep = ",", dec = ".", quote = TRUE,
                       row.names = FALSE, col.names = TRUE, na = "NA",
                       qmethod = "double")
  })
}

.sim_check_header <- function(object, meta, label) {
  if (!is.list(object) || !identical(object$meta, meta)) {
    stop(label, " does not match this generation config, phase, replicate, ",
         "or AlphaSimR version. Existing data were not changed.", call. = FALSE)
  }
  invisible(TRUE)
}

.sim_cross_matrix <- function(Z, a, d, intercept) {
  pf <- Z[seq_len(100L), , drop = FALSE] / 2
  pm <- Z[100L + seq_len(100L), , drop = FALSE] / 2
  additive <- outer(drop(pf %*% a), drop(pm %*% a), "+") - sum(a)
  dominance <- outer(drop(pf %*% d), drop(pm %*% d), "+") -
    2 * tcrossprod(sweep(pf, 2L, d, "*"), pm)
  M <- additive + dominance + intercept
  dimnames(M) <- list(rownames(pf), rownames(pm))
  M
}

.sim_verify <- function(problem, tolerance = 1e-10) {
  Z <- problem$geno
  ids <- c(sprintf("F%03d", seq_len(100L)), sprintf("M%03d", seq_len(100L)))
  sex <- rep(c("F", "M"), each = 100L)
  if (!is.matrix(Z) || !identical(dim(Z), c(200L, 1000L)) ||
      !all(Z %in% 0:2) || !identical(rownames(Z), ids) ||
      !identical(problem$id, ids) || !identical(problem$sex, sex) ||
      !identical(problem$Z, Z) ||
      !identical(dim(problem$M), c(100L, 100L)) ||
      !identical(dimnames(problem$M), list(ids[1:100], ids[101:200])) ||
      length(problem$gv) != 200L || !identical(names(problem$gv), ids) ||
      length(problem$a) != 1000L || length(problem$d) != 1000L ||
      !identical(names(problem$a), colnames(Z)) ||
      !identical(names(problem$d), colnames(Z)) ||
      !all(is.finite(c(problem$M, problem$gv, problem$a, problem$d,
                       problem$trait_intercept)))) {
    stop("Invalid benchmark dimensions, identities, coding, or finite values.",
         call. = FALSE)
  }
  q <- problem$qtl_index
  if (length(q) != 200L || anyNA(q) || anyDuplicated(q) ||
      any(!q %in% seq_len(1000L)) ||
      nrow(problem$map) != 1000L ||
      !identical(problem$map$marker, colnames(Z)) ||
      !identical(which(problem$map$is_qtl), q) ||
      !all(tabulate(problem$map$chr, nbins = 10L) == 100L) ||
      !all(tabulate(problem$map$chr[q], nbins = 10L) == 20L) ||
      any(problem$a[-q] != 0) || any(problem$d[-q] != 0)) {
    stop("Invalid SNP/QTL map or effects alignment.", call. = FALSE)
  }
  allele_frequency <- colMeans(Z) / 2
  maf <- pmin(allele_frequency, 1 - allele_frequency)
  for (k in seq_len(10L)) {
    sites <- which(problem$map$chr == k)
    selected <- q[problem$map$chr[q] == k]
    eligible <- sites[maf[sites] >= 0.05]
    if (length(eligible) >= 20L) {
      valid <- all(selected %in% eligible)
    } else {
      best <- sites[order(-maf[sites], sites)][seq_len(20L)]
      valid <- identical(sort(selected), sort(best))
    }
    if (!valid) {
      stop("QTL selection violates the prespecified best-available-MAF policy ",
           "on chromosome ", k, ". Founder checkpoint retained.", call. = FALSE)
    }
  }
  # Independent genotype-class lookup, not the matrix expression used for M.
  gv_lookup <- rep(problem$trait_intercept, 200L)
  for (j in q) {
    gv_lookup <- gv_lookup + c(-problem$a[j], problem$d[j], problem$a[j])[Z[, j] + 1L]
  }
  gv_matrix <- drop((Z - 1) %*% problem$a + (Z == 1) %*% problem$d) +
    problem$trait_intercept
  gv_error <- max(abs(gv_lookup - problem$gv), abs(gv_matrix - problem$gv))
  # Independently enumerate offspring genotype probabilities for every pair.
  M_reference <- matrix(0, 100L, 100L)
  for (i in seq_len(100L)) {
    pf <- Z[i, q] / 2
    for (j in seq_len(100L)) {
      pm <- Z[100L + j, q] / 2
      prob0 <- (1 - pf) * (1 - pm)
      prob2 <- pf * pm
      prob1 <- pf * (1 - pm) + (1 - pf) * pm
      M_reference[i, j] <- sum(-problem$a[q] * prob0 + problem$d[q] * prob1 +
                                  problem$a[q] * prob2) + problem$trait_intercept
    }
  }
  M_error <- max(abs(problem$M - M_reference))
  if (gv_error > tolerance || M_error > tolerance) {
    stop("True-effect verification failed: max gv error = ", gv_error,
         "; max M error = ", M_error, ". Checkpoints retained.", call. = FALSE)
  }
  list(gv_max_abs_error = gv_error, M_max_abs_error = M_error,
       verification_tolerance = tolerance)
}

.sim_exports <- function(problem, directory, meta) {
  components <- list(M = problem$M, gv = problem$gv, geno = problem$geno,
                     effects = problem$effects, map = problem$map,
                     stats = problem$stats, config = problem$config,
                     id_sex = data.frame(id = problem$id, sex = problem$sex))
  for (name in names(components)) {
    .sim_rds(components[[name]], file.path(directory, paste0(name, ".rds")))
  }
  .sim_csv(data.frame(id = problem$id, sex = problem$sex, problem$geno,
                      check.names = FALSE), file.path(directory, "geno.csv.gz"))
  .sim_csv(data.frame(female_id = rownames(problem$M), problem$M,
                      check.names = FALSE), file.path(directory, "M.csv.gz"))
  .sim_csv(data.frame(id = problem$id, sex = problem$sex, gv = unname(problem$gv)),
           file.path(directory, "gv.csv.gz"))
  .sim_csv(problem$effects, file.path(directory, "effects.csv.gz"))
  .sim_csv(problem$map, file.path(directory, "map.csv.gz"))
  .sim_csv(problem$chromosome_stats, file.path(directory, "chromosome_stats.csv.gz"))
  .sim_csv(data.frame(statistic = names(problem$stats),
                      value = vapply(problem$stats, as.character, character(1))),
           file.path(directory, "stats.csv.gz"))
  .sim_csv(data.frame(setting = names(problem$config),
                      value = vapply(problem$config, as.character, character(1))),
           file.path(directory, "config.csv.gz"))
  filenames <- c("founder_checkpoint.rds", "trait_checkpoint.rds", "problem.rds",
                  paste0(names(components), ".rds"),
                  paste0(c("geno", "M", "gv", "effects", "map",
                           "chromosome_stats", "stats", "config"), ".csv.gz"))
  manifest <- list(meta = meta, status = "complete",
                    files = stats::setNames(unname(tools::md5sum(
                      file.path(directory, filenames))), filenames))
  if (anyNA(manifest$files)) stop("Cannot checksum every artifact.", call. = FALSE)
  .sim_rds(manifest, file.path(directory, "manifest.rds"))
  invisible(manifest)
}

# Returns plain R objects; no R6 pointer, phenotype, optimizer, or ratio state.
# Resume by calling with the same out, phase, rep_id and generation config.
# A different generation config must use a different output directory.
generate_problem <- function(rep_id, phase = c("pilot", "test"), cfg = list(),
                             out = ".") {
  rep_id <- .sim_scalar(rep_id, "rep_id", TRUE, 1, 1000000)
  if (missing(phase)) phase <- "pilot"
  if (!is.character(phase) || length(phase) != 1L || is.na(phase) ||
      !phase %in% c("pilot", "test")) {
    stop("phase must be exactly 'pilot' or 'test'.", call. = FALSE)
  }
  config <- validate_simulation_config(cfg)
  if (!is.character(out) || length(out) != 1L || is.na(out) || !nzchar(out)) {
    stop("out must be one nonempty directory path.", call. = FALSE)
  }
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    stop("AlphaSimR is required to generate or validate these inputs.", call. = FALSE)
  }
  meta <- list(schema = .sim_schema, config_hash = .sim_hash(config),
                phase = phase, rep_id = rep_id,
                AlphaSimR_version = as.character(utils::packageVersion("AlphaSimR")))
  directory <- file.path(out, "simulation_data", phase, sprintf("rep_%04d", rep_id))
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE)) {
    stop("Cannot create replicate directory: ", directory, call. = FALSE)
  }
  paths <- stats::setNames(file.path(directory, paste0(
    c("founder_checkpoint", "trait_checkpoint", "problem", "manifest"), ".rds")),
    c("founder", "trait", "problem", "manifest"))

  # Restore the caller's RNG, including the absence of .Random.seed.
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)

  tryCatch({
    if (file.exists(paths[["manifest"]])) {
      manifest <- readRDS(paths[["manifest"]])
      .sim_check_header(manifest, meta, "Completed manifest")
      observed <- unname(tools::md5sum(file.path(directory, names(manifest$files))))
      if (!identical(manifest$status, "complete") || anyNA(observed) ||
          !identical(observed, unname(manifest$files))) {
        stop("Completed artifacts are missing or changed; no files were replaced.",
             call. = FALSE)
      }
      problem <- readRDS(paths[["problem"]])
      .sim_check_header(problem, meta, "Saved problem")
      .sim_verify(problem)
      return(problem)
    }
    if (file.exists(paths[["problem"]])) {
      problem <- readRDS(paths[["problem"]])
      .sim_check_header(problem, meta, "Saved problem")
      .sim_verify(problem)
      .sim_exports(problem, directory, meta)
      return(problem)
    }
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    set.seed(config$seed)
    stream <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    # Phases use separate 2^127 streams; replicates use 2^76 substreams.
    # This is independent of generation order and never uses a retry seed.
    if (phase == "test") stream <- parallel::nextRNGStream(stream)
    if (rep_id > 1L) for (i in seq_len(rep_id - 1L))
      stream <- parallel::nextRNGSubStream(stream)
    assign(".Random.seed", stream, envir = .GlobalEnv)

    if (file.exists(paths[["founder"]])) {
      checkpoint <- readRDS(paths[["founder"]])
      .sim_check_header(checkpoint, meta, "Founder checkpoint")
      if (!identical(checkpoint$initial_rng, stream) ||
          !identical(.sim_hash(checkpoint$founder), checkpoint$founder_hash)) {
        stop("Founder checkpoint integrity check failed; it will not be redrawn.",
             call. = FALSE)
      }
    } else {
      if (file.exists(paths[["trait"]])) {
        stop("Trait checkpoint exists without its founder; refusing to redraw.",
             call. = FALSE)
      }
      founder <- tryCatch(
        AlphaSimR::runMacs(nInd = 200L, nChr = 10L, segSites = 100L,
                            species = "GENERIC", manualGenLen = 1,
                            ploidy = 2L, inbred = FALSE, nThreads = 1L),
        error = function(e) stop("MaCS founder generation failed: ",
                                 conditionMessage(e),
                                 ". No alternative founder model or seed was used.",
                                 call. = FALSE))
      # Save exact phased packed haplotypes BEFORE inspecting MAF or drawing QTL.
      # Plain slot lists are portable and do not serialize SimParam/R6 objects.
      founder_slots <- .sim_slots(founder)
      checkpoint <- list(meta = meta, config = config, founder_class = "MapPop",
                           founder = founder_slots,
                           founder_hash = .sim_hash(founder_slots),
                           initial_rng = stream,
                           rng_after_founder = get(".Random.seed", envir = .GlobalEnv),
                           R_version = as.character(getRversion()))
      .sim_rds(checkpoint, paths[["founder"]])
      # All downstream work uses the saved representation, even on first run.
      checkpoint <- readRDS(paths[["founder"]])
    }
    founder <- .sim_restore(checkpoint$founder_class, checkpoint$founder)
    assign(".Random.seed", checkpoint$rng_after_founder, envir = .GlobalEnv)
    SP <- AlphaSimR::SimParam$new(founder)
    SP$nThreads <- 1L
    SP$restrSegSites(overlap = TRUE)
    Z <- AlphaSimR::pullSegSiteGeno(founder, simParam = SP)
    ids <- c(sprintf("F%03d", seq_len(100L)), sprintf("M%03d", seq_len(100L)))
    sex <- rep(c("F", "M"), each = 100L)
    markers <- unlist(lapply(founder@genMap, names), use.names = FALSE)
    dimnames(Z) <- list(ids, markers)
    freq <- colMeans(Z) / 2
    maf <- pmin(freq, 1 - freq)
    chr <- rep(seq_len(10L), each = 100L)
    eligible <- tabulate(chr[maf >= 0.05], nbins = 10L)
    feasible <- eligible >= 20L
    # Restrict only QTL eligibility, not the 1,000 observed SNPs or founders.
    excluded <- unlist(lapply(seq_len(10L), function(k) {
      sites <- which(chr == k)
      if (feasible[k]) sites[maf[sites] < 0.05]
      else setdiff(sites, sites[order(-maf[sites], sites)][seq_len(20L)])
    }), use.names = FALSE)
    if (length(excluded)) SP$restrSegSites(excludeQtl = markers[excluded], overlap = TRUE)
    SP$addSnpChipByName(markers, name = "all_1000_SNPs_including_QTL")

    if (file.exists(paths[["trait"]])) {
      trait_checkpoint <- readRDS(paths[["trait"]])
      .sim_check_header(trait_checkpoint, meta, "Trait checkpoint")
      if (!identical(trait_checkpoint$founder_hash, checkpoint$founder_hash) ||
          !identical(.sim_hash(trait_checkpoint$trait), trait_checkpoint$trait_hash)) {
        stop("Trait checkpoint integrity check failed; effects will not be redrawn.",
             call. = FALSE)
      }
      SP$manAddTrait(.sim_restore("TraitAD", trait_checkpoint$trait))
      assign(".Random.seed", trait_checkpoint$rng_after_trait, envir = .GlobalEnv)
    } else {
      SP$addTraitAD(nQtlPerChr = 20L, mean = 0, var = 1,
                     meanDD = config$meanDD, varDD = config$varDD, useVarA = TRUE)
      trait_slots <- .sim_slots(SP$traits[[1L]])
      trait_checkpoint <- list(meta = meta, founder_hash = checkpoint$founder_hash,
                                trait = trait_slots, trait_hash = .sim_hash(trait_slots),
                                rng_after_trait = get(".Random.seed", envir = .GlobalEnv))
      .sim_rds(trait_checkpoint, paths[["trait"]])
    }
    trait <- SP$traits[[1L]]
    qtl_index <- as.integer(rep(c(0L, cumsum(founder@nLoci)[-10L]),
                                times = trait@lociPerChr) + trait@lociLoc)
    a <- d <- stats::setNames(numeric(1000L), markers)
    a[qtl_index] <- trait@addEff
    d[qtl_index] <- trait@domEff
    pop <- AlphaSimR::newPop(founder, id = ids, simParam = SP)
    pop@sex <- sex
    # newPop's default residual variance is NA. No setPheno call or phenotype
    # enters the problem; pop@gv is retained verbatim as the reference values.
    gv <- stats::setNames(as.numeric(pop@gv[, 1L]), ids)
    snp_Z <- AlphaSimR::pullSnpGeno(pop, simParam = SP)
    qtl_Z <- AlphaSimR::pullQtlGeno(pop, simParam = SP)
    if (!identical(unname(Z), unname(snp_Z)) ||
        !identical(unname(Z[, qtl_index, drop = FALSE]), unname(qtl_Z))) {
      stop("AlphaSimR SNP/QTL dosage alignment check failed.", call. = FALSE)
    }
    M <- .sim_cross_matrix(Z, a, d, trait@intercept)
    map <- data.frame(marker = markers, snp_index = seq_len(1000L), chr = chr,
                       locus = rep(seq_len(100L), 10L),
                       position_M = unlist(founder@genMap, use.names = FALSE),
                       allele_frequency = freq, maf = maf,
                       is_qtl = seq_len(1000L) %in% qtl_index, row.names = NULL)
    effects <- data.frame(marker = markers, snp_index = seq_len(1000L),
                           is_qtl = map$is_qtl, a = unname(a), d = unname(d))
    chromosome_stats <- data.frame(chr = seq_len(10L), n_snp = 100L, n_qtl = 20L,
                                    simulated_length_M = 1,
                                    retained_marker_span_M = vapply(founder@genMap,
                                      function(x) diff(range(x)), numeric(1)),
                                    maf_eligible = eligible, maf_constraint_feasible = feasible,
                                    maf_fallback_used = !feasible,
                                    maf_status = ifelse(feasible, "preferred_threshold",
                                                        "best_available_fallback"),
                                    qtl_maf_below_005 = tabulate(
                                      chr[qtl_index[maf[qtl_index] < 0.05]], nbins = 10L),
                                    qtl_maf_min = vapply(seq_len(10L), function(k)
                                      min(maf[qtl_index[chr[qtl_index] == k]]), numeric(1)))
    gp <- AlphaSimR::genParam(pop, simParam = SP)
    problem <- list(meta = meta, config = config, config_hash = meta$config_hash,
                      rep_id = rep_id, phase = phase, seed = config$seed,
                      rng_seed = stream, M = M, gv = gv, id = ids, sex = sex,
                      geno = Z, Z = Z, a = a, d = d, qtl_index = qtl_index,
                      a_qtl = a[qtl_index], d_qtl = d[qtl_index],
                      trait_intercept = as.numeric(trait@intercept),
                      intercept = as.numeric(trait@intercept), M_center_offset = 0,
                      map = map, effects = effects, chromosome_stats = chromosome_stats,
                      founder_hash = checkpoint$founder_hash,
                      stats = list(n_females = 100L, n_males = 100L, n_snps = 1000L,
                                   n_qtl = 200L, n_chr = 10L, meanDD = config$meanDD,
                                   founder_method = "runMacs", species = "GENERIC",
                                   simulated_chr_length_M = 1,
                                   varDD = config$varDD, varA = as.numeric(gp$varA),
                                   varD = as.numeric(gp$varD), varG = as.numeric(gp$varG),
                                   gv_mean = mean(gv), gv_sd = stats::sd(gv),
                                   M_mean = mean(M), M_sd = stats::sd(as.vector(M)),
                                   qtl_maf_min = min(maf[qtl_index]),
                                   qtl_maf_below_005 = sum(maf[qtl_index] < 0.05),
                                   qtl_maf_policy = config$qtl_maf_policy,
                                   maf_fallback_used = any(!feasible),
                                   maf_fallback_chr_count = sum(!feasible),
                                   maf_status = if (all(feasible)) "preferred_threshold"
                                     else "best_available_fallback",
                                   maf_constraint_feasible = all(feasible),
                                   maf_constraint_met = all(maf[qtl_index] >= 0.05),
                                   trait_intercept = as.numeric(trait@intercept)))
    verification <- .sim_verify(problem)
    if (abs(problem$stats$varA - 1) > 1e-10 || abs(mean(gv)) > 1e-10 ||
        any(maf[qtl_index][feasible[chr[qtl_index]]] < 0.05)) {
      stop("Trait scaling, centering, or feasible QTL MAF preference failed.", call. = FALSE)
    }
    problem$stats <- c(problem$stats, verification)
    .sim_rds(problem, paths[["problem"]])
    .sim_exports(problem, directory, meta)
    problem
  }, error = function(e) {
    # Keep every successful checkpoint/export and a separate failure record.
    failure <- list(meta = meta, message = conditionMessage(e),
                      time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                      rng = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                        get(".Random.seed", envir = .GlobalEnv) else NULL)
    try(saveRDS(failure, tempfile("error_", directory, fileext = ".rds"),
                version = 2), silent = TRUE)
    stop(e)
  })
}
