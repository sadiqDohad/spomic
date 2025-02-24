#' @export
calculate_csr_envelopes <- function(spomic, verbose = TRUE) {
  envelope_results <- list()
  unique_marks <- unique(spatstat.geom::marks(spomic@pp))
  combinations <- expand.grid(unique_marks, unique_marks)
  pairwise_combinations <- apply(combinations, 1, function(x) paste(x, collapse = "_"))

  envelope_results <- lapply(1:length(pairwise_combinations), function(k) {
    i <- as.character(combinations[k, 1])
    j <- as.character(combinations[k, 2])

    if(verbose) cat("Calculating CSR envelope for pair:", i, "-", j, "\n")

    envelope_result <- spatstat.explore::envelope(
      Y = spatstat.geom::rescale(spomic@pp),
      fun = spatstat.explore::Kcross.inhom,
      i = i,
      j = j,
      fix.marks = TRUE,
      correction = "Ripley",
      nsim = spomic@details$hyperparameters$csr_nsim,
      global = FALSE,
      verbose = verbose
    )

    envelope_results[[paste0(i, "_", j)]] <- envelope_result

  })
  names(envelope_results) <- pairwise_combinations

  spomic@results$csr_envelope <- envelope_results
  return(spomic)
}

#' @export
find_nonrandom_pairs <- function(spomic) {
  envelope_list <- spomic@results$csr_envelope
  nonrandom_pairs <- list()
  pairs <- names(envelope_list)
  for(pair in names(envelope_list)) {
    env <- envelope_list[[pair]]
    obs <- approx(env$r, env$obs, xout = spomic@details$hyperparameters$r)$y
    lo <- approx(env$r, env$lo, xout = spomic@details$hyperparameters$r)$y
    hi <- approx(env$r, env$hi, xout = spomic@details$hyperparameters$r)$y

    nonrandom_pairs[[pair]] <- data.frame(
      i_j = pair,
      r = spomic@details$hyperparameters$r,
      observed = obs,
      lower = lo,
      upper = hi,
      nonrandom = !dplyr::between(obs, lo, hi)
    )
  }
  nonrandom_pairs <- dplyr::bind_rows(nonrandom_pairs)
  cell_pairs <- nonrandom_pairs |> dplyr::filter(nonrandom) |> dplyr::pull(i_j)

  spomic@results$nonrandom_pairs <- cell_pairs

  return(spomic)
}


#' @export
get_kcross <- function(spomic, i, j) {
  if(any(!c(i, j) %in% spatstat.geom::marks(spomic@pp)) | (sum(spatstat.geom::marks(spomic@pp) == i) < 2)) {
    return(data.frame(sample = spomic@details$sample,
                      i_j = paste0(i, "_", j),
                      colocalization_type = spomic@details$hyperparameters$colocalization_type,
                      r = spomic@details$hyperparameters$r,
                      colocalization_stat = NA,
                      colocalization_var = NA))
  }

  i_subset <- spatstat.geom::subset.ppp(spomic@pp, marks == i)
  j_subset <-  spatstat.geom::subset.ppp(spomic@pp, marks == j)

  suppressWarnings({
    invisible(capture.output({ # The lohboot function has some annoying printouts
      loh_bootstrap <- spatstat.explore::lohboot(spatstat.geom::rescale(spomic@pp),
                                                 spatstat.explore::Kcross.inhom,
                                                 from = i,
                                                 to = j,
                                                 lambdaFrom = spatstat.explore::density.ppp(i_subset),
                                                 lambdaTo = spatstat.explore::density.ppp(j_subset),
                                                 correction = "Ripley",
                                                 global = FALSE,
                                                 nsim = 100)
    })

    )
  })
  spomic@results$colocalization_bootstrap[[paste0(i, "_", j)]] <- loh_bootstrap
  return(spomic)
}

#' @export
get_kcross_summary <- function(spomic) {
  results <- list()
  for(i_j in names(spomic@results$colocalization_bootstrap)) {
    # print(i_j)

    kcross_i_j <- spomic@results$colocalization_bootstrap[[i_j]]
    kcross <- approx(kcross_i_j$r, kcross_i_j$iso, xout = spomic@details$hyperparameters$r)$y
    kcross_lo <- approx(kcross_i_j$r, kcross_i_j$lo, xout = spomic@details$hyperparameters$r)$y
    kcross_hi <- approx(kcross_i_j$r, kcross_i_j$hi, xout = spomic@details$hyperparameters$r)$y
    kcross_var <- ((kcross_hi - kcross_lo) / (2 * 1.96))^2

    results[[i_j]] <- data.frame(
      sample = spomic@details$sample,
      i_j = i_j,
      kcross = kcross,
      kcross_var = kcross_var)
  }
  return(dplyr::bind_rows(results))
}
