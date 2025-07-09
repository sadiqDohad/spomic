#' @export
get_silverman <- function(pp_subset) {
  # Function to use Silverman's rule of thumb for kernel density estimation
  n <- pp_subset$n
  coords <- spatstat.geom::coords(pp_subset)
  sd_x <- sd(coords$x)
  sd_y <- sd(coords$y)
  sd_total <- sqrt(sd_x^2 + sd_y^2)
  sigma_silverman <- sd_total * n^(-1/5)

  return(sigma_silverman)
}

#' @export
calculate_csr_envelopes <- function(spomic, density_smoothing = "silverman", verbose = FALSE) {
  unique_marks <- unique(spatstat.geom::marks(spomic@pp))
  combinations <- expand.grid(unique_marks, unique_marks)
  pairwise_combinations <- apply(combinations, 1, function(x) paste(x, collapse = "_"))

  envelope_results <- vector("list", length(pairwise_combinations))
  names(envelope_results) <- pairwise_combinations

  for (k in seq_along(pairwise_combinations)) {
    i <- as.character(combinations[k, 1])
    j <- as.character(combinations[k, 2])
    pair_label <- paste0(i, "_", j)

    if (verbose) cat("Calculating CSR envelope for pair:", i, "-", j, "\n")

    result <- tryCatch({
      ij_subset <- spatstat.geom::subset.ppp(spomic@pp, marks %in% c(i, j))

      if (spomic@details$hyperparameters$colocalization_type == "Kcross") {
        spatstat.explore::envelope(
          Y = spatstat.geom::rescale(ij_subset),
          fun = spatstat.explore::Kcross,
          i = i,
          j = j,
          correction = "Ripley",
          nsim = spomic@details$hyperparameters$csr_nsim,
          global = FALSE,
          fix.n = TRUE,

          verbose = verbose
        )
      } else if (spomic@details$hyperparameters$colocalization_type == "Kcross.inhom") {
        stop("Kcross.inhom not implemented yet")
      } else if (spomic@details$hyperparameters$colocalization_type == "CLQ") {
        stop("CLQ not implemented yet")
      } else {
        stop("Set a valid spatial statistic hyperparameter!: c('Kcross', 'Kcross.inhom', 'CLQ')")
      }
    }, error = function(e) {
      warning(sprintf("CSR envelope calculation failed for pair '%s': %s", pair_label, e$message))
      NULL
    })

    envelope_results[[pair_label]] <- result
  }

  # Remove NULL results
  envelope_results <- Filter(Negate(is.null), envelope_results)

  spomic@results$csr_envelope <- envelope_results
  return(spomic)
}

safe_approx <- function(x, y, xout) {
  if (sum(!is.na(x) & !is.na(y)) >= 2) {
    return(approx(x, y, xout = xout)$y)
  } else {
    return(rep(NA, length(xout)))
  }
}

#' @export
find_nonrandom_pairs <- function(spomic) {
  envelope_list <- spomic@results$csr_envelope
  nonrandom_pairs <- list()
  pairs <- names(envelope_list)
  for(pair in names(envelope_list)) {
    env <- envelope_list[[pair]]
    # obs <- approx(env$r, env$obs, xout = spomic@details$hyperparameters$r)$y
    obs <- safe_approx(env$r, env$obs, spomic@details$hyperparameters$r)

    theo <- approx(env$r, env$theo, xout = spomic@details$hyperparameters$r)$y
    lo <- safe_approx(env$r, env$lo, spomic@details$hyperparameters$r)
    hi <- safe_approx(env$r, env$hi, spomic@details$hyperparameters$r)
    # mmean <- approx(env$r, env$mmean, xout = spomic@details$hyperparameters$r)$y

    nonrandom_pairs[[pair]] <- data.frame(
      i_j = pair,
      r = spomic@details$hyperparameters$r,
      observed = obs,
      lower = lo,
      upper = hi,
      nonrandom = !dplyr::between(obs, lo, hi) & dplyr::between(theo, lo, hi)
      # nonrandom = !dplyr::between(obs, lo, hi) & dplyr::between(mmean, lo, hi)

    )
  }
  nonrandom_pairs <- dplyr::bind_rows(nonrandom_pairs)
  cell_pairs <- nonrandom_pairs |> dplyr::filter(nonrandom) |> dplyr::pull(i_j)

  spomic@results$nonrandom_pairs <- cell_pairs

  return(spomic)
}

run_lohboot <- function(spomic, i, j) {

  if(spomic@details$hyperparameters$colocalization_type == "Kcross") {
    return(get_kcross(spomic, i, j))
  } else if(spomic@details$hyperparameters$colocalization_type == "Kcross.inhom") {
    return(get_kcross.inhom(spomic, i, j))

  } else if(spomic@details$hyperparameters$colocalization_type == "Lcross") {
    return(get_lcross(spomic, i, j))
  } else if(spomic@details$hyperparameters$colocalization_type == "Lcross.inhom") {
    return(get_lcross.inhom(spomic, i, j))
  } else if(spomic@details$hyperparameters$colocalization_type == "CLQ") {
    return(get_clq(spomic))
  } else {
    warning("Enter a supported colocalization statistic: {Kcross, Kcross.inhom, Lcross, Lcross.inhom, CLQ}.")
  }

}

#' @export
get_kcross <- function(spomic, i, j) {
  if(any(!c(i, j) %in% spatstat.geom::marks(spomic@pp)) | (sum(spatstat.geom::marks(spomic@pp) == i) < 2) | (sum(spatstat.geom::marks(spomic@pp) == j) < 2)) {
    spomic@results$colocalization_bootstrap[[paste0(i, "_", j)]] <- data.frame(sample = spomic@details$sample,
                                                                               i_j = paste0(i, "_", j),
                                                                               colocalization_type = spomic@details$hyperparameters$colocalization_type,
                                                                               r = spomic@details$hyperparameters$r,
                                                                               colocalization_stat = NA,
                                                                               colocalization_var = NA)
    return(spomic)
  }
  ij_subset <- spatstat.geom::subset.ppp(spomic@pp, marks %in% c(i, j))
  i_subset <- spatstat.geom::subset.ppp(spomic@pp, marks == i)
  j_subset <-  spatstat.geom::subset.ppp(spomic@pp, marks == j)

  if (i != j) {
    suppressWarnings({
      invisible(capture.output({ # The lohboot function has some annoying printouts
        loh_bootstrap <- spatstat.explore::lohboot(spatstat.geom::rescale(ij_subset),
                                                   spatstat.explore::Kcross,
                                                   from = i,
                                                   to = j,
                                                   correction = "Ripley",
                                                   global = ifelse(spomic@details$hyperparameters$fixed_distance, FALSE, TRUE),
                                                   nsim = 100)}))
      })
  } else {
    suppressWarnings({
      invisible(capture.output({ # The lohboot function has some annoying printouts
        loh_bootstrap <- spatstat.explore::lohboot(spatstat.geom::rescale(ij_subset),
                                                   spatstat.explore::Kest,
                                                   correction = "Ripley",
                                                   global = FALSE,
                                                   nsim = 100)}))
      })
  }
  spomic@results$colocalization_bootstrap[[paste0(i, "_", j)]] <- loh_bootstrap
  return(spomic)
}

#' @export
get_kcross.inhom <- function(spomic, i, j) {
  if(any(!c(i, j) %in% spatstat.geom::marks(spomic@pp)) | (sum(spatstat.geom::marks(spomic@pp) == i) < 2) | (sum(spatstat.geom::marks(spomic@pp) == j) < 2)) {
    spomic@results$colocalization_bootstrap[[paste0(i, "_", j)]] <- data.frame(sample = spomic@details$sample,
                                                                               i_j = paste0(i, "_", j),
                                                                               colocalization_type = spomic@details$hyperparameters$colocalization_type,
                                                                               r = spomic@details$hyperparameters$r,
                                                                               colocalization_stat = NA,
                                                                               colocalization_var = NA)
    return(spomic)
  }
  ij_subset <- spatstat.geom::subset.ppp(spomic@pp, marks %in% c(i, j))
  i_subset <- spatstat.geom::subset.ppp(spomic@pp, marks == i)
  j_subset <-  spatstat.geom::subset.ppp(spomic@pp, marks == j)

  if (i != j) {
    suppressWarnings({
      invisible(capture.output({ # The lohboot function has some annoying printouts
        loh_bootstrap <- spatstat.explore::lohboot(spatstat.geom::rescale(ij_subset),
                                                   spatstat.explore::Kcross.inhom,
                                                   from = i,
                                                   to = j,
                                                   correction = "Ripley",
                                                   global = ifelse(spomic@details$hyperparameters$fixed_distance, FALSE, TRUE),
                                                   nsim = 100)}))
    })
  } else {
    suppressWarnings({
      invisible(capture.output({ # The lohboot function has some annoying printouts
        loh_bootstrap <- spatstat.explore::lohboot(spatstat.geom::rescale(ij_subset),
                                                   spatstat.explore::Kinhom,
                                                   correction = "Ripley",
                                                   global = ifelse(spomic@details$hyperparameters$fixed_distance, FALSE, TRUE),
                                                   nsim = 100)}))
    })
  }
  spomic@results$colocalization_bootstrap[[paste0(i, "_", j)]] <- loh_bootstrap
  return(spomic)
}

#' @export
get_lcross <- function(spomic, i, j) {
  if(any(!c(i, j) %in% spatstat.geom::marks(spomic@pp)) | (sum(spatstat.geom::marks(spomic@pp) == i) < 2) | (sum(spatstat.geom::marks(spomic@pp) == j) < 2)) {
    spomic@results$colocalization_bootstrap[[paste0(i, "_", j)]] <- data.frame(sample = spomic@details$sample,
                                                                               i_j = paste0(i, "_", j),
                                                                               colocalization_type = spomic@details$hyperparameters$colocalization_type,
                                                                               r = spomic@details$hyperparameters$r,
                                                                               colocalization_stat = NA,
                                                                               colocalization_var = NA)
    return(spomic)
  }
  ij_subset <- spatstat.geom::subset.ppp(spomic@pp, marks %in% c(i, j))
  i_subset <- spatstat.geom::subset.ppp(spomic@pp, marks == i)
  j_subset <-  spatstat.geom::subset.ppp(spomic@pp, marks == j)

  if (i != j) {
    suppressWarnings({
      invisible(capture.output({ # The lohboot function has some annoying printouts
        loh_bootstrap <- spatstat.explore::lohboot(spatstat.geom::rescale(ij_subset),
                                                   spatstat.explore::Lcross,
                                                   from = i,
                                                   to = j,
                                                   correction = "Ripley",
                                                   global = ifelse(spomic@details$hyperparameters$fixed_distance, FALSE, TRUE),
                                                   nsim = 100)}))
    })
  } else {
    suppressWarnings({
      invisible(capture.output({ # The lohboot function has some annoying printouts
        loh_bootstrap <- spatstat.explore::lohboot(spatstat.geom::rescale(ij_subset),
                                                   spatstat.explore::Lest,
                                                   correction = "Ripley",
                                                   global = ifelse(spomic@details$hyperparameters$fixed_distance, FALSE, TRUE),
                                                   nsim = 100)}))
    })
  }
  spomic@results$colocalization_bootstrap[[paste0(i, "_", j)]] <- loh_bootstrap
  return(spomic)
}

get_lcross.inhom <- function(spomic, i, j) {
  if(any(!c(i, j) %in% spatstat.geom::marks(spomic@pp)) | (sum(spatstat.geom::marks(spomic@pp) == i) < 2) | (sum(spatstat.geom::marks(spomic@pp) == j) < 2)) {
    spomic@results$colocalization_bootstrap[[paste0(i, "_", j)]] <- data.frame(sample = spomic@details$sample,
                                                                               i_j = paste0(i, "_", j),
                                                                               colocalization_type = spomic@details$hyperparameters$colocalization_type,
                                                                               r = spomic@details$hyperparameters$r,
                                                                               colocalization_stat = NA,
                                                                               colocalization_var = NA)
    return(spomic)
  }
  ij_subset <- spatstat.geom::subset.ppp(spomic@pp, marks %in% c(i, j))
  i_subset <- spatstat.geom::subset.ppp(spomic@pp, marks == i)
  j_subset <-  spatstat.geom::subset.ppp(spomic@pp, marks == j)

  if (i != j) {
    suppressWarnings({
      invisible(capture.output({ # The lohboot function has some annoying printouts
        loh_bootstrap <- spatstat.explore::lohboot(spatstat.geom::rescale(ij_subset),
                                                   spatstat.explore::Lcross.inhom,
                                                   from = i,
                                                   to = j,
                                                   correction = "Ripley",
                                                   global = ifelse(spomic@details$hyperparameters$fixed_distance, FALSE, TRUE),
                                                   nsim = 100)}))
    })
  } else {
    suppressWarnings({
      invisible(capture.output({ # The lohboot function has some annoying printouts
        loh_bootstrap <- spatstat.explore::lohboot(spatstat.geom::rescale(ij_subset),
                                                   spatstat.explore::Linhom,
                                                   correction = "Ripley",
                                                   global = ifelse(spomic@details$hyperparameters$fixed_distance, FALSE, TRUE),
                                                   nsim = 100)}))
    })
  }
  spomic@results$colocalization_bootstrap[[paste0(i, "_", j)]] <- loh_bootstrap
  return(spomic)
}

get_clq <- function(spomic, n_sim = 100) {
  sf <- ppp_to_sf(spomic@pp)
  marks_vec <- as.character(sf$marks)
  nb <- sfdep::st_dist_band(sf::st_geometry(sf), upper=spomic@details$hyperparameters$r)


  weights <- sfdep::st_kernel_weights(
    nb, sf,
    kernel = "gaussian",
    adaptive = TRUE
  )

  valid <- lengths(nb) > 1
  sf <- sf[valid, ]
  nb <- nb[valid]
  weights <- weights[valid]
  marks_vec <- marks_vec[valid]

  # Compute colocations
  lclq_results <- sfdep::local_colocation(marks_vec, marks_vec, nb, wt = weights, nsim = 5)
  n_col <- ncol(lclq_results)
  lclq_results <- lclq_results[, 1:(n_col/2)]
  lclq_results <- dplyr::mutate(lclq_results,
                                cell_label = marks_vec,
                                cell_ID = seq_along(marks_vec)
  )
  lclq_results[is.na(lclq_results)] <- 0

  cell_types <- unique(marks_vec)
  long_lclq <- tidyr::pivot_longer(lclq_results, cols = all_of(cell_types), names_to = "variable", values_to = "value") |>
    dplyr::mutate(pair_label = paste(cell_label, variable, sep = "_"))

  pair_labels <- expand.grid(i = cell_types, j = cell_types) |>
    dplyr::mutate(pair_label = paste(i, j, sep = "_")) |>
    dplyr::pull(pair_label)

  # Bootstrap with replicate
  boot_results <- replicate(n_sim, {
    long_lclq |>
      dplyr::group_by(pair_label) |>
      dplyr::summarise(clq = mean(sample(value, replace = TRUE), na.rm = TRUE), .groups = "drop")
  }, simplify = FALSE)

  summary_df <- dplyr::bind_rows(boot_results, .id = "sim") |>
    dplyr::group_by(pair_label) |>
    dplyr::summarise(
      colocalization = mean(clq, na.rm = TRUE),
      colocalization_lo = quantile(clq, probs = 0.025, na.rm = TRUE),
      colocalization_hi = quantile(clq, probs = 0.975, na.rm = TRUE),
      colocalization_var = ((colocalization_hi - colocalization_lo) / (2 * 1.96))^2,
      .groups = "drop"
    )

  spomic@results$colocalization_bootstrap <- summary_df
  return(spomic)
}


get_spatial_stats <- function(spomic) {
  cell_types <- unique(spatstat.geom::marks(spomic@pp))
  cell_pairs <- expand.grid(i = cell_types, j = cell_types, stringsAsFactors = FALSE)

  process_one_pair <- function(pair_row) {
    i <- pair_row$i
    j <- pair_row$j

    temp_spomic <- run_lohboot(spomic, i, j)
    key <- paste0(i, "_", j)
    val <- temp_spomic@results$colocalization_bootstrap[[key]]

    list(key = key, value = val)
  }

  # processed_list <- pbapply::pblapply(seq_len(nrow(cell_pairs)), function(idx) {
  #   process_one_pair(cell_pairs[idx, ])
  # })

  if(spomic@details$hyperparameters$colocalization_type != "CLQ") {
    processed_list <- pbapply::pblapply(seq_len(nrow(cell_pairs)), function(idx) {
      process_one_pair(cell_pairs[idx, ])
    })
    results <- setNames(
      lapply(processed_list, `[[`, "value"),
      sapply(processed_list, `[[`, "key")
    )
    spomic@results$colocalization_bootstrap <- results
  } else {
    spomic <- run_lohboot(spomic, i = "", j = "")
  }
  return(spomic)
}


#' @export
get_all_kcross <- function(spomic, lohboot=TRUE) {
  # Get unique cell types
  cell_types <- unique(spatstat.geom::marks(spomic@pp))

  # Generate all (i, j) pairs including (i, i)
  cell_pairs <- expand.grid(i = cell_types, j = cell_types, stringsAsFactors = FALSE)

  # Store results
  results <- list()

  # Iterate over all (i, j) pairs
  for (k in seq_len(nrow(cell_pairs))) {
    i <- cell_pairs$i[k]
    j <- cell_pairs$j[k]

    # cat("Processing:", i, j, "\n")  # Print progress

    # Run get_kcross for this pair
      temp_spomic <- get_kcross(spomic, i, j)
    # if(lohboot) {
    #   temp_spomic <- get_kcross(spomic, i, j)
    # } else {
    #     crossK <- spatstat.explore::Kcross(X=spatstat.geom::rescale(ij_subset),
    #                              i=i,
    #                              j=j,
    #                              correction="Ripley")
    #
    #   }
    results[[paste0(i, "_", j)]] <- temp_spomic@results$colocalization_bootstrap[[paste0(i, "_", j)]]
  }

  spomic@results$colocalization_bootstrap <- results
  return(spomic)
}

#' @export
get_spatial_summary <- function(spomic) {
  colocalization_type <- spomic@details$hyperparameters$colocalization_type
  r_val <- spomic@details$hyperparameters$r
  sample_name <- spomic@details$sample

  if(colocalization_type != "CLQ") {
    keys <- names(spomic@results$colocalization_bootstrap)
    results <- lapply(keys, function(i_j) {
      colocalization_i_j <- spomic@results$colocalization_bootstrap[[i_j]]

      if (length(colocalization_i_j$iso) > 1) {
        colocalization <- approx(colocalization_i_j$r, colocalization_i_j$iso, xout = r_val)$y
        colocalization_theo <- approx(colocalization_i_j$r, colocalization_i_j$theo, xout = r_val)$y
        colocalization_lo <- tryCatch(
          approx(colocalization_i_j$r, colocalization_i_j$lo, xout = r_val)$y,
          error = function(e) NA
        )
        colocalization_hi <- tryCatch(
          approx(colocalization_i_j$r, colocalization_i_j$hi, xout = r_val)$y,
          error = function(e) NA
        )
        # colocalization_lo <- approx(colocalization_i_j$r, colocalization_i_j$lo, xout = r_val)$y
        # colocalization_hi <- approx(colocalization_i_j$r, colocalization_i_j$hi, xout = r_val)$y
        colocalization_var <- ((colocalization_hi - colocalization_lo) / (2 * 1.96))^2

      } else {
        colocalization <- colocalization_lo <- colocalization_hi <- colocalization_var <- NA
      }

      data.frame(
        sample = sample_name,
        i_j = i_j,
        colocalization_type = colocalization_type,
        colocalization_stat = colocalization,
        colocalization_theo = colocalization_theo,
        colocalization_lo = colocalization_lo,
        colocalization_hi = colocalization_hi,
        colocalization_var = colocalization_var)
    })
    return(dplyr::bind_rows(results))
  } else {
    return(spomic@results$colocalization_bootstrap)
  }
#   r_val <- spomic@details$hyperparameters$r
#   sample_name <- spomic@details$sample
#
# keys <- names(spomic@results$colocalization_bootstrap)
# results <- lapply(keys, function(i_j) {
#   colocalization_i_j <- spomic@results$colocalization_bootstrap[[i_j]]
#
#   if (length(colocalization_i_j$iso) > 1) {
#     colocalization <- approx(colocalization_i_j$r, colocalization_i_j$iso, xout = r_val)$y
#     colocalization_theo <- approx(colocalization_i_j$r, colocalization_i_j$theo, xout = r_val)$y
#     colocalization_lo <- approx(colocalization_i_j$r, colocalization_i_j$lo, xout = r_val)$y
#     colocalization_hi <- approx(colocalization_i_j$r, colocalization_i_j$hi, xout = r_val)$y
#     colocalization_var <- ((colocalization_hi - colocalization_lo) / (2 * 1.96))^2
#
#   } else {
#     colocalization <- colocalization_lo <- colocalization_hi <- colocalization_var <- NA
#   }
#
#   data.frame(
#     sample = sample_name,
#     i_j = i_j,
#     colocalization_type = colocalization_type,
#     colocalization_stat = colocalization,
#     colocalization_theo = colocalization_theo,
#     colocalization_lo = colocalization_lo,
#     colocalization_hi = colocalization_hi,
#     colocalization_var = colocalization_var)
# })

return(dplyr::bind_rows(results))
}



#' @export
get_kcross_summary <- function(spomic) {
  # results <- list()
  # for(i_j in names(spomic@results$colocalization_bootstrap)) {
  #   # print(i_j)
  #
  #   kcross_i_j <- spomic@results$colocalization_bootstrap[[i_j]]
  #   kcross <- approx(kcross_i_j$r, kcross_i_j$iso, xout = spomic@details$hyperparameters$r)$y
  #   kcross_lo <- approx(kcross_i_j$r, kcross_i_j$lo, xout = spomic@details$hyperparameters$r)$y
  #   kcross_hi <- approx(kcross_i_j$r, kcross_i_j$hi, xout = spomic@details$hyperparameters$r)$y
  #   kcross_var <- ((kcross_hi - kcross_lo) / (2 * 1.96))^2
  #
  #   results[[i_j]] <- data.frame(
  #     sample = spomic@details$sample,
  #     i_j = i_j,
  #     kcross = kcross,
  #     kcross_var = kcross_var)
  # }
  results <- list()
  for(i_j in names(spomic@results$colocalization_bootstrap)) {
    print(i_j)

    kcross_i_j <- spomic@results$colocalization_bootstrap[[i_j]]
    if(length(kcross_i_j$iso)>1){
      kcross <- approx(kcross_i_j$r, kcross_i_j$iso, xout = spomic@details$hyperparameters$r)$y
      kcross_lo <- approx(kcross_i_j$r, kcross_i_j$lo, xout = spomic@details$hyperparameters$r)$y
      kcross_hi <- approx(kcross_i_j$r, kcross_i_j$hi, xout = spomic@details$hyperparameters$r)$y
      kcross_var <- ((kcross_hi - kcross_lo) / (2 * 1.96))^2


      kcross_theo <- approx(kcross_i_j$r, kcross_i_j$theo, xout = spomic@details$hyperparameters$r)$y
      kcross_diff <- kcross - kcross_theo
      kcross_diff_lo <- kcross_lo - kcross_theo
      kcross_diff_hi <- kcross_hi - kcross_theo
      kcross_diff_var <- ((kcross_diff_hi - kcross_diff_lo) / (2 * 1.96))^2

      lcross <- sqrt(kcross / pi)
      lcross_theo <- spomic@details$hyperparameters$r
      lcross_lo <- sqrt(kcross_lo / pi)
      lcross_hi <- sqrt(kcross_hi / pi)
      lcross_var <- ((lcross_hi - lcross_lo) / (2 * 1.96))^2




    } else {
      kcross <- NA
      kcross_lo <- NA
      kcross_hi <- NA
      kcross_var <- NA

      kcross_theo <- NA
      kcross_diff <- NA
      kcross_diff_lo <- NA
      kcross_diff_hi <- NA
      kcross_diff_var <- NA

      lcross <- NA
      lcross_theo <- NA
      lcross_lo <- NA
      lcross_hi <- NA
      lcross_var <- NA
    }


    results[[i_j]] <- data.frame(
      sample = spomic@details$sample,
      i_j = i_j,
      kcross = kcross,
      kcross_var = kcross_var,
      kcross_diff = kcross_diff,
      kcross_diff_var = kcross_diff_var,
      lcross = lcross,
      lcross_lo = lcross_lo,
      lcross_hi = lcross_hi,
      lcross_var = lcross_var)}
  return(dplyr::bind_rows(results))
}
