suppressPackageStartupMessages({
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Package 'lpSolve' is required because GM uses ILP.")
  }
})

arg_or <- function(args, i, default) {
  if (length(args) >= i && nzchar(args[[i]])) args[[i]] else default
}

split_arg <- function(x) {
  x <- gsub(";", ",", x)
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

safe_num <- function(x) {
  out <- gsub("-", "m", format(x, trim = TRUE, scientific = FALSE))
  gsub("\\.", "p", out)
}

trait_label_from_prefix <- function(prefix) {
  p <- basename(prefix)
  if (grepl("area", p, ignore.case = TRUE)) return("Area")
  if (grepl("content", p, ignore.case = TRUE)) return("Content")
  if (grepl("den", p, ignore.case = TRUE)) return("Density")
  p
}

read_gm_matrix <- function(file) {
  x <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  rn <- as.character(x[[1]])
  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- rn
  mat
}

find_or_extract_file <- function(prefix, suffix, bundle_file = "") {
  prefix_dir <- dirname(prefix)
  if (identical(prefix_dir, ".")) prefix_dir <- getwd()
  prefix_base <- basename(prefix)
  fname <- paste0(prefix_base, suffix)

  candidates <- c(
    file.path(prefix_dir, fname),
    file.path(prefix_dir, "parameter_sensitivity_work", fname),
    file.path(getwd(), fname),
    file.path(getwd(), "parameter_sensitivity_work", fname)
  )

  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) return(normalizePath(hit[[1]], winslash = "/", mustWork = TRUE))

  bundle_candidates <- c(
    bundle_file,
    file.path(prefix_dir, "pgen1005250_results_bundle.tar.gz"),
    file.path(prefix_dir, "DOWNLOAD_ME_pgen1005250_results_bundle.tar.gz"),
    file.path(getwd(), "pgen1005250_results_bundle.tar.gz"),
    file.path(getwd(), "DOWNLOAD_ME_pgen1005250_results_bundle.tar.gz")
  )
  bundle_candidates <- bundle_candidates[nzchar(bundle_candidates)]
  bundle_hit <- bundle_candidates[file.exists(bundle_candidates)]
  if (length(bundle_hit) == 0) {
    stop("Cannot find input file: ", fname, "\nAlso cannot find pgen1005250 results bundle.")
  }

  extract_dir <- file.path(prefix_dir, "pgen1005250_parameter_input_extract")
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  utils::untar(bundle_hit[[1]], files = fname, exdir = extract_dir)

  extracted <- file.path(extract_dir, fname)
  if (!file.exists(extracted)) {
    stop("Tried to extract ", fname, " from ", bundle_hit[[1]], " but did not find it.")
  }
  normalizePath(extracted, winslash = "/", mustWork = TRUE)
}

assert_same_dimnames <- function(a, b, label_a, label_b) {
  if (!identical(rownames(a), rownames(b)) || !identical(colnames(a), colnames(b))) {
    stop("Matrix dimnames do not match: ", label_a, " vs ", label_b)
  }
}

plot_heatmap_base <- function(summary_df, value_col, title, out_file, annotate = TRUE) {
  c_vals <- sort(unique(summary_df$c_weight))
  l_vals <- sort(unique(summary_df$lambda_weight))
  z <- matrix(
    NA_real_,
    nrow = length(c_vals),
    ncol = length(l_vals),
    dimnames = list(as.character(c_vals), as.character(l_vals))
  )

  for (i in seq_len(nrow(summary_df))) {
    z[as.character(summary_df$c_weight[i]), as.character(summary_df$lambda_weight[i])] <- summary_df[[value_col]][i]
  }

  zlim <- range(z, na.rm = TRUE)
  if (diff(zlim) < .Machine$double.eps) {
    zlim <- zlim + c(-1, 1) * max(1, abs(zlim[1])) * 0.01
  }
  pal <- grDevices::hcl.colors(120, "YlOrRd", rev = FALSE)

  grDevices::tiff(out_file, width = 7.2, height = 6.2, units = "in", res = 600, compression = "lzw")
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)

  layout(matrix(c(1, 2), nrow = 1), widths = c(5.1, 0.55))
  par(mar = c(5, 5, 4, 1))
  image(
    x = c_vals,
    y = l_vals,
    z = z,
    col = pal,
    zlim = zlim,
    xlab = "c weight",
    ylab = "lambda weight",
    main = title,
    axes = FALSE,
    useRaster = TRUE
  )
  axis(1, at = c_vals, labels = format(c_vals, trim = TRUE))
  axis(2, at = l_vals, labels = format(l_vals, trim = TRUE), las = 1)
  box(lwd = 1.2)

  grid <- expand.grid(c_weight = c_vals, lambda_weight = l_vals)
  vals <- mapply(function(cw, lw) z[as.character(cw), as.character(lw)], grid$c_weight, grid$lambda_weight)
  if (annotate && length(vals) <= 121) {
    text(grid$c_weight, grid$lambda_weight, labels = sprintf("%.1f", vals), cex = 0.45)
  }

  par(mar = c(5, 0.5, 4, 3.5))
  yseq <- seq(zlim[1], zlim[2], length.out = length(pal))
  image(
    x = 1,
    y = yseq,
    z = matrix(yseq, nrow = 1),
    col = pal,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  axis(4, las = 1, cex.axis = 0.8)
  mtext(value_col, side = 4, line = 2.5, cex = 0.8)
}

plot_combined_heatmap <- function(all_summary, value_col, out_file) {
  trait_levels <- unique(all_summary$Trait)
  c_vals <- sort(unique(all_summary$c_weight))
  l_vals <- sort(unique(all_summary$lambda_weight))
  pal <- grDevices::hcl.colors(120, "YlOrRd", rev = FALSE)
  zlim <- range(all_summary[[value_col]], na.rm = TRUE)
  if (diff(zlim) < .Machine$double.eps) {
    zlim <- zlim + c(-1, 1) * max(1, abs(zlim[1])) * 0.01
  }

  grDevices::tiff(out_file, width = 12, height = 4.6, units = "in", res = 600, compression = "lzw")
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)

  layout(matrix(c(seq_along(trait_levels), length(trait_levels) + 1), nrow = 1), widths = c(rep(1, length(trait_levels)), 0.12))
  for (tr in trait_levels) {
    d <- all_summary[all_summary$Trait == tr, , drop = FALSE]
    z <- matrix(
      NA_real_,
      nrow = length(c_vals),
      ncol = length(l_vals),
      dimnames = list(as.character(c_vals), as.character(l_vals))
    )
    for (i in seq_len(nrow(d))) {
      z[as.character(d$c_weight[i]), as.character(d$lambda_weight[i])] <- d[[value_col]][i]
    }
    par(mar = c(4.5, 4.5, 3, 0.8))
    image(
      x = c_vals,
      y = l_vals,
      z = z,
      col = pal,
      zlim = zlim,
      xlab = "c weight",
      ylab = "lambda weight",
      main = tr,
      axes = FALSE,
      useRaster = TRUE
    )
    axis(1, at = c_vals, labels = format(c_vals, trim = TRUE), cex.axis = 0.75)
    axis(2, at = l_vals, labels = format(l_vals, trim = TRUE), las = 1, cex.axis = 0.75)
    box(lwd = 1.1)
  }

  par(mar = c(4.5, 0.4, 3, 3.6))
  yseq <- seq(zlim[1], zlim[2], length.out = length(pal))
  image(x = 1, y = yseq, z = matrix(yseq, nrow = 1), col = pal, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = 0.75)
  mtext(value_col, side = 4, line = 2.4, cex = 0.8)
}

run_one_trait_sensitivity <- function(
    prefix,
    out_root,
    n_female_select,
    n_male_select,
    female_pool_size,
    male_pool_size,
    n_restart,
    seed,
    search_mode,
    grid_step,
    baseline_c,
    baseline_lambda,
    bundle_file = "",
    force_use_all_males = TRUE) {

  trait_label <- trait_label_from_prefix(prefix)
  trait_safe <- gsub("[^A-Za-z0-9]+", "_", trait_label)
  out_dir <- file.path(out_root, paste0(basename(prefix), "_parameter_sensitivity_ms_ilp"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pair_dir <- file.path(out_dir, "pair_results_by_parameter")
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)

  mean_file <- find_or_extract_file(prefix, "_progeny_mean_matrix.csv", bundle_file)
  sd_file <- find_or_extract_file(prefix, "_progeny_sd_matrix.csv", bundle_file)
  k_file <- find_or_extract_file(prefix, "_inbreeding_penalty_matrix.csv", bundle_file)
  parent_file <- find_or_extract_file(prefix, "_gm_parent_genetic_values.csv", bundle_file)

  mean_mat <- read_gm_matrix(mean_file)
  sd_mat <- read_gm_matrix(sd_file)
  k_mat <- read_gm_matrix(k_file)

  assert_same_dimnames(mean_mat, sd_mat, paste0(trait_label, " mean"), paste0(trait_label, " sd"))
  assert_same_dimnames(mean_mat, k_mat, paste0(trait_label, " mean"), paste0(trait_label, " K"))

  parent_df <- read.csv(parent_file, check.names = FALSE, stringsAsFactors = FALSE)
  parent_df$ID <- as.character(parent_df$ID)
  parent_df$Sex <- as.character(parent_df$Sex)

  female_value_df <- parent_df[parent_df$Sex == "F" & parent_df$ID %in% rownames(mean_mat), , drop = FALSE]
  male_value_df <- parent_df[parent_df$Sex == "M" & parent_df$ID %in% colnames(mean_mat), , drop = FALSE]

  if (nrow(female_value_df) < n_female_select) {
    stop("Not enough females for ", trait_label, ": available=", nrow(female_value_df), ", requested=", n_female_select)
  }
  if (nrow(male_value_df) < n_male_select) {
    stop("Not enough males for ", trait_label, ": available=", nrow(male_value_df), ", requested=", n_male_select)
  }

  female_pool_size <- min(female_pool_size, nrow(female_value_df))
  male_pool_size <- min(male_pool_size, nrow(male_value_df))

  c_grid <- seq(0, 2, by = grid_step)
  lambda_grid <- seq(0, 2, by = grid_step)
  summary_list <- vector("list", length(c_grid) * length(lambda_grid))

  if (all(abs(k_mat) < 1e-12, na.rm = TRUE)) {
    cat("Note:", trait_label, "K matrix is all zero; lambda will not change the objective.\n")
  }

  counter <- 0L
  for (cw in c_grid) {
    for (lw in lambda_grid) {
      counter <- counter + 1L
      cat(sprintf("[%s %03d/%03d] c=%.3f lambda=%.3f\n", trait_label, counter, length(summary_list), cw, lw))

      value_mat <- mean_mat + cw * sd_mat - lw * k_mat
      ans <- gm_multistart_all_runs(
        value_mat = value_mat,
        female_value_df = female_value_df,
        male_value_df = male_value_df,
        n_female_select = n_female_select,
        n_male_select = n_male_select,
        female_pool_size = female_pool_size,
        male_pool_size = male_pool_size,
        n_restart = n_restart,
        seed = seed + counter,
        force_use_all_males = force_use_all_males,
        search_mode = search_mode,
        verbose = FALSE
      )

      pair_result <- extract_pair_result(
        ans = ans,
        parent_value_df = parent_df,
        mean_mat = mean_mat,
        sd_mat = sd_mat,
        K_mat = k_mat
      )

      pair_result$Trait <- trait_label
      pair_result$c_weight <- cw
      pair_result$lambda_weight <- lw
      pair_result$ValueUnderCurrentParameter <- pair_result$ProgenyMean + cw * pair_result$ProgenySD - lw * pair_result$InbreedingPenalty
      pair_result$ValueUnderCommonBaseline <- pair_result$ProgenyMean + baseline_c * pair_result$ProgenySD - baseline_lambda * pair_result$InbreedingPenalty

      pair_file <- file.path(pair_dir, paste0("pair_result_", trait_safe, "_c", safe_num(cw), "_lambda", safe_num(lw), ".csv"))
      write.csv(pair_result, pair_file, row.names = FALSE)

      summary_list[[counter]] <- data.frame(
        Trait = trait_label,
        Prefix = basename(prefix),
        c_weight = cw,
        lambda_weight = lw,
        OptimizedObjectiveCurrent = ans$best_final_objective,
        TotalCurrentObjective = sum(pair_result$ValueUnderCurrentParameter),
        TotalCommonBaselineObjective = sum(pair_result$ValueUnderCommonBaseline),
        TotalProgenyMean = sum(pair_result$ProgenyMean),
        MeanProgenyMean = mean(pair_result$ProgenyMean),
        TotalProgenySD = sum(pair_result$ProgenySD),
        MeanProgenySD = mean(pair_result$ProgenySD),
        TotalInbreedingPenalty = sum(pair_result$InbreedingPenalty),
        MeanInbreedingPenalty = mean(pair_result$InbreedingPenalty),
        UniqueDams = length(unique(pair_result$Dam)),
        UniqueSires = length(unique(pair_result$Sire)),
        BestRunSource = ans$best_run_source,
        SearchMode = ans$global_search_mode,
        PairFile = pair_file,
        stringsAsFactors = FALSE
      )
    }
  }

  summary_df <- do.call(rbind, summary_list)
  summary_file <- file.path(out_dir, paste0(basename(prefix), "_parameter_sensitivity_multistart_ilp_summary.csv"))
  write.csv(summary_df, summary_file, row.names = FALSE)

  plot_heatmap_base(
    summary_df,
    "TotalCurrentObjective",
    paste0(trait_label, ": raw objective under each c/lambda"),
    file.path(out_dir, paste0(basename(prefix), "_heatmap_raw_current_objective.tiff"))
  )
  plot_heatmap_base(
    summary_df,
    "TotalCommonBaselineObjective",
    paste0(trait_label, ": common baseline objective"),
    file.path(out_dir, paste0(basename(prefix), "_heatmap_common_baseline_objective.tiff"))
  )
  plot_heatmap_base(
    summary_df,
    "TotalProgenyMean",
    paste0(trait_label, ": total progeny mean"),
    file.path(out_dir, paste0(basename(prefix), "_heatmap_total_progeny_mean.tiff"))
  )
  plot_heatmap_base(
    summary_df,
    "TotalProgenySD",
    paste0(trait_label, ": total progeny SD"),
    file.path(out_dir, paste0(basename(prefix), "_heatmap_total_progeny_sd.tiff"))
  )

  summary_df
}

run_chicken_three_traits_parameter_sensitivity_ms_ilp <- function(
    prefixes = c(
      "pgen1005250_X6_true_trab_area",
      "pgen1005250_X6_true_trab_content",
      "pgen1005250_X6_true_trab_den"
    ),
    out_root = "pgen1005250_parameter_sensitivity_ms_ilp",
    n_female_select = 40,
    n_male_select = 20,
    female_pool_size = 200,
    male_pool_size = 200,
    n_restart = 3,
    seed = 2026,
    search_mode = "best",
    grid_step = 0.2,
    baseline_c = 0.1,
    baseline_lambda = 0,
    gm_code_file = "run_chicken_bglr_gm_compare_all.R",
    bundle_file = "") {

  if (!file.exists(gm_code_file)) {
    stop("Cannot find GM code file: ", gm_code_file)
  }
  source(gm_code_file, encoding = "UTF-8")

  required_functions <- c("gm_multistart_all_runs", "extract_pair_result")
  missing_functions <- required_functions[!vapply(required_functions, exists, logical(1), mode = "function")]
  if (length(missing_functions) > 0) {
    stop("GM code file does not define: ", paste(missing_functions, collapse = ", "))
  }

  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  cat("Traits:", paste(prefixes, collapse = ", "), "\n")
  cat("Grid step:", grid_step, "\n")
  cat("Search mode:", search_mode, "\n")
  cat("n_restart:", n_restart, "\n")
  cat("Selected females:", n_female_select, "\n")
  cat("Selected males:", n_male_select, "\n")
  cat("Female pool:", female_pool_size, "\n")
  cat("Male pool:", male_pool_size, "\n")

  all_summary <- do.call(rbind, lapply(seq_along(prefixes), function(i) {
    run_one_trait_sensitivity(
      prefix = prefixes[[i]],
      out_root = out_root,
      n_female_select = n_female_select,
      n_male_select = n_male_select,
      female_pool_size = female_pool_size,
      male_pool_size = male_pool_size,
      n_restart = n_restart,
      seed = seed + i * 10000L,
      search_mode = search_mode,
      grid_step = grid_step,
      baseline_c = baseline_c,
      baseline_lambda = baseline_lambda,
      bundle_file = bundle_file
    )
  }))

  all_summary_file <- file.path(out_root, "pgen1005250_three_traits_parameter_sensitivity_multistart_ilp_summary.csv")
  write.csv(all_summary, all_summary_file, row.names = FALSE)

  plot_combined_heatmap(
    all_summary,
    "TotalCommonBaselineObjective",
    file.path(out_root, "pgen1005250_three_traits_heatmap_common_baseline_objective.tiff")
  )
  plot_combined_heatmap(
    all_summary,
    "TotalProgenyMean",
    file.path(out_root, "pgen1005250_three_traits_heatmap_total_progeny_mean.tiff")
  )
  plot_combined_heatmap(
    all_summary,
    "TotalProgenySD",
    file.path(out_root, "pgen1005250_three_traits_heatmap_total_progeny_sd.tiff")
  )

  cat("Finished chicken three-trait parameter sensitivity with multistart local search + ILP.\n")
  cat("Combined summary:", normalizePath(all_summary_file, winslash = "/", mustWork = TRUE), "\n")
  cat("Output directory:", normalizePath(out_root, winslash = "/", mustWork = TRUE), "\n")

  invisible(all_summary)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)

  prefixes <- split_arg(arg_or(
    args,
    1,
    "pgen1005250_X6_true_trab_area,pgen1005250_X6_true_trab_content,pgen1005250_X6_true_trab_den"
  ))
  out_root <- arg_or(args, 2, "pgen1005250_parameter_sensitivity_ms_ilp")
  n_female_select <- as.integer(arg_or(args, 3, "40"))
  n_male_select <- as.integer(arg_or(args, 4, "20"))
  female_pool_size <- as.integer(arg_or(args, 5, "200"))
  male_pool_size <- as.integer(arg_or(args, 6, "200"))
  n_restart <- as.integer(arg_or(args, 7, "3"))
  seed <- as.integer(arg_or(args, 8, "2026"))
  search_mode <- arg_or(args, 9, "best")
  grid_step <- as.numeric(arg_or(args, 10, "0.2"))
  baseline_c <- as.numeric(arg_or(args, 11, "0.1"))
  baseline_lambda <- as.numeric(arg_or(args, 12, "0"))
  gm_code_file <- arg_or(args, 13, "run_chicken_bglr_gm_compare_all.R")
  bundle_file <- arg_or(args, 14, "")

  run_chicken_three_traits_parameter_sensitivity_ms_ilp(
    prefixes = prefixes,
    out_root = out_root,
    n_female_select = n_female_select,
    n_male_select = n_male_select,
    female_pool_size = female_pool_size,
    male_pool_size = male_pool_size,
    n_restart = n_restart,
    seed = seed,
    search_mode = search_mode,
    grid_step = grid_step,
    baseline_c = baseline_c,
    baseline_lambda = baseline_lambda,
    gm_code_file = gm_code_file,
    bundle_file = bundle_file
  )
}
