
# Function to create outcome
createY <- function(Us, As, option = c('linear', 'nonlinear')){
  # Us is a n x nreps matrix of simulated unmeasured confounder
  # As is a n x nreps matrix of simulated exposure
  # option is a string indicating the form of the outcome model
  # returns a vector of outcome
  option <- match.arg(option)
  n <- nrow(Us)
  nreps <- ncol(Us)
  Ys <- matrix(NA, n, nreps)
  # linear outcome model
  if (option == 'linear'){
    for (i in 1:nreps){
      Ys[,i] <- rnorm(n, -2 + (-2)*Us[,i] + As[,i] - 1*As[,i]*Us[,i], 1)
    }
  }
  # nonlinear outcome model
  if (option == 'nonlinear'){
    for (i in 1:nreps){
      eta <- -2 - 2*Us[,i] + As[,i] - 0.4*As[,i]^2 - 0.5*Us[,i]*As[,i]
      Ys[,i] <- rnorm(n, eta, 1)
    }
  }
  return(Ys)
}

# Function to plot the variables titled "names" in the dataframe "df"
plotfunc <- function(df, names, labels=names,
                    xlimits = c(-125, -65), ylimits = c(25, 50)){
  # df is a sf dataframe including columns names
  # names is a vector of strings with the names of the columns to be plotted
  # labels is a vector of strings to title the ggplots
  # returns a list of ggplots

  K <- length(names)
  gs <- list()
  for (k in 1:K){
    # extract the column with names[k] from df
    var <- df[[names[k]]]
    qs <- round(quantile(var, probs = c(0.1, 0.3, 0.5, 0.7, 0.9)),2)
    # Turn qs into a vector of strings
    qschar <- as.character(qs)
    gs[[k]] <- ggplot(df) +
      xlim(xlimits[1],xlimits[2]) + 
      ylim(ylimits[1], ylimits[2]) +
      geom_sf(aes_string(fill = names[k]), color=NA, size = 0.005) +
      scale_fill_gradient2(low = "#1e90ff", 
                           mid = "white", 
                           high = "#8b0000", 
                           midpoint = 0,
                           breaks = qs,
                           labels = qschar,
                           limits = c(min(var), max(var)),
                           na.value = "white") +
      theme_minimal() +
      theme(plot.title = element_text(size = 24 * 2,hjust = 0.5),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            line = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal", 
            legend.text.align = 0.75,
            legend.key.width = unit(80, "points"),
            panel.grid.major = element_line(colour = "transparent"),
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 25)
            ) + 
      ggtitle(labels[k])
  }
  return(gs)
}

# Function that simulates data, estimates truncated exposure effect using different methods, and saves results to csvs
simfunc <- function(nsims,
                   lat,
                   lon,
                   confounding_mechanism,
                   option = c('linear', 'nonlinear'),
                   methods = c(
                     # 'baseline',
                      'oracle',
                     # 'spatialcoord',
                     # 'trueIV',
                     # 'trueIV-spatialcoord',
                     'IV-TPS',
                     'IV-GraphLaplacian'#,
                     # 'IV-TPS-spatialcoord',
                     # 'IV-GraphLaplacian-spatialcoord'
                   ),
                   B_tps_full,
                   B_gl_full,
                   statemat,
                   W = NULL,
                   cutoff = 0,
                   select_basis = TRUE,
                   n_cores = 1L,
                   results_dir = "results_Apr9/",
                   iv_control = list())
{
  # nsims is the number of simulations
  # lat is a vector of latitudes
  # lon is a vector of longitudes
  # option is the form of the outcome model
  # methods are the methods used to estimate truncated exposure effect
  # statemat is the matrix of state-level indicators
  # cutoff is c
  # select_basis: if TRUE, choose n_ac with the sequential stability rule
  # select_basis_dimension() in ../R/basis_selection.R; otherwise fit a
  # single fixed n_ac (iv_control$fixed_n_ac_fraction * ncol(B)).
  # results_dir: directory to write output CSVs
  #
  # No outer cross-fitting: each simulated dataset gets one full-sample basis
  # selection + effect estimate, matching ../data_application/1_run_analysis.R.

  iv_defaults <- list(
    density_trim = 0.01,
    r2_reference = 0.7,      # reference B^c explains at least this share of Var(A)
    fixed_n_ac_fraction = 0.1,
    n_grid = 25L,
    sl_library = c("SL.gam", "SL.glm", "SL.mean", "SL.glm.interaction"),
    bandwidth = NULL,
    bw_seq = NULL,
    constrain = TRUE,
    family_alpha = 0.05,
    selection_seed = 12001L,
    save_diagnostics = TRUE
  )
  if (!is.list(iv_control)) {
    stop("iv_control must be a list.")
  }
  iv_control <- utils::modifyList(iv_defaults, iv_control)
  iv_methods <- c(
    "IV-TPS",
    "IV-GraphLaplacian",
    "IV-TPS-spatialcoord",
    "IV-GraphLaplacian-spatialcoord"
  )
  required_estimator_functions <- "estimate_truncated_effect_full"
  if (any(methods %in% iv_methods)) {
    required_estimator_functions <- c(
      required_estimator_functions,
      "project_Ac_full",
      "select_basis_dimension"
    )
  }
  missing_estimator_functions <- required_estimator_functions[
    !vapply(required_estimator_functions, exists, logical(1), mode = "function")
  ]
  if (length(missing_estimator_functions) > 0L) {
      stop(
        paste0(
          "Source the R modules before simfunc(): ",
          paste(missing_estimator_functions, collapse = ", ")
        )
      )
  }

  confounding_mechanism <- as.integer(confounding_mechanism)
  option <- match.arg(option)
  
  ################# GENERATE DATA #################
  
  # Compute distance matrix
  distmat <- geosphere::distm(cbind(lon, lat), 
                              fun = geosphere::distHaversine)
  distmat <- distmat/1000000 # scale so range (0,2)
  n <- length(lat)
  dat <- generate_mechanism_data(
    confounding_mechanism, nsims, B_tps_full = B_tps_full, B_gl_full = B_gl_full,
    W = W, lat = lat, lon = lon, distmat = distmat
  )
  Ac <- dat$Ac
  Auc <- dat$Auc
  A <- Ac + Auc # all have dimension n x nsims
  if (confounding_mechanism %in% c(3, 7)) {
    # Two confounders: U1 enters as U does elsewhere, plus an A x U2 term.
    U1 <- dat$U1
    U2 <- dat$U2
    Y <- matrix(NA, n, nsims)
    for (i in seq_len(nsims)) {
      eta <- if (option == 'linear') {
        -2 - 2*U1[,i] + A[,i] - 1*A[,i]*U1[,i] - 1.5*A[,i]*U2[,i]
      } else {
        -2 - 2*U1[,i] + A[,i] - 0.4*A[,i]^2 - 0.5*U1[,i]*A[,i] - 0.5*U2[,i]*A[,i]
      }
      Y[, i] <- rnorm(n, eta, 1)
    }
  } else {
    U <- dat$U
    Y <- createY(Us=U, As=A, option = option)
  }
  
  
  ################# FIT MODELS #################

  # Normalize coordinates to [0,1] once; used by all spatialcoord methods so
  # that SuperLearner splines operate on the same scale as the DGP basis functions.
  lat_n_m <- (lat  - min(lat))  / (max(lat)  - min(lat))
  lon_n_m <- (lon - min(lon)) / (max(lon) - min(lon))

  t_simfunc_start <- proc.time()["elapsed"]

  for (method in methods){
    # Create filename for csvs containing estimates
    filename <- paste0(results_dir, 'conf', confounding_mechanism, '_', option, '_', method, '.csv')
  
    t_method_start <- proc.time()["elapsed"]

    # Run sims in parallel (n_cores = 1 → sequential, identical behaviour).
    results_list <- parallel::mclapply(seq_len(nsims), function(sim) {
      t_sim_start <- proc.time()["elapsed"]
      message(paste(method, sim))

      # ---- build xmat ----
      n_ac_sel <- NA_integer_
      n_uc_sel <- NA_integer_
      selection_diagnostics <- NULL
      cor_Ac   <- NA_real_
      y <- Y[, sim]
      a <- A[, sim]

      # The IV methods execute one full-sample basis selection (searching a
      # grid of candidate confounded dimensions n_ac) followed by one
      # full-sample doubly robust effect estimate at the selected dimension.
      # Non-IV methods go straight to the full-sample estimator below.
      if (method %in% iv_methods) {
        B_method <- if (method %in% c("IV-TPS", "IV-TPS-spatialcoord")) {
          B_tps_full
        } else {
          B_gl_full
        }
        # For the reversed-scale mechanisms (4, 8) the instrumental variation
        # is assumed to sit at LARGE spatial scales, so the basis is simply
        # reversed (smallest scale in column 1, largest in column m).
        # project_Ac_full()/select_basis_dimension_full() treat the first
        # n_ac columns as confounded, so after reversal B^c is the n_ac
        # smallest-scale columns and B^uc the m - n_ac largest-scale columns;
        # walking n_ac down moves progressively smaller-scale columns into
        # B^uc. No knowledge of the DGP's confounded block is used here.
        B_for_selection <- if (confounding_mechanism %in% reversed_mechanisms) {
          B_method[, rev(seq_len(ncol(B_method))), drop = FALSE]
        } else {
          B_method
        }
        x_extra <- if (method %in% c(
          "IV-TPS-spatialcoord",
          "IV-GraphLaplacian-spatialcoord"
        )) {
          out <- cbind(lat_n_m, lon_n_m)
          colnames(out) <- c("Latitude", "Longitude")
          out
        } else {
          matrix(1, nrow = n, ncol = 1, dimnames = list(NULL, "Intercept"))
        }
        m <- ncol(B_for_selection)

        effect_args <- list(
          nuisance_args = list(
            sl_library = iv_control$sl_library,
            density_trim = iv_control$density_trim
          ),
          n_grid = iv_control$n_grid,
          bandwidth = iv_control$bandwidth,
          bw_seq = iv_control$bw_seq,
          density_trim = iv_control$density_trim,
          constrain = iv_control$constrain
        )

        iv_fit <- tryCatch(
          if (select_basis) {
            select_basis_dimension(
              y = y,
              a = a,
              x = x_extra,
              B = B_for_selection,
              cutoff = cutoff,
              effect_estimator = estimate_truncated_effect_full,
              effect_args = effect_args,
              r2_reference = iv_control$r2_reference,
              family_alpha = iv_control$family_alpha,
              seed = iv_control$selection_seed
            )
          } else {
            n_ac_fixed <- max(1L, as.integer(floor(iv_control$fixed_n_ac_fraction * m)))
            projection <- project_Ac_full(B_for_selection, a, n_ac_fixed)
            set.seed(iv_control$selection_seed)
            fixed_fit <- do.call(
              estimate_truncated_effect_full,
              c(
                list(y = y, a = a, w = cbind(Ac = projection$Ac, x_extra), cutoff = cutoff),
                effect_args
              )
            )
            list(
              selected_n_ac = n_ac_fixed,
              selected_n_uc = m - n_ac_fixed,
              selected_fit = fixed_fit,
              selected_projection = projection,
              diagnostics = NULL
            )
          },
          error = function(e) {
            message("IV estimator error: ", e$message)
            NULL
          }
        )

        if (is.null(iv_fit)) {
          return(list(
            muest = NA_real_,
            ci = c(NA_real_, NA_real_),
            n_ac = NA_integer_,
            n_uc = NA_integer_,
            selection_diagnostics = NULL,
            cor_Ac = NA_real_,
            time_s = proc.time()["elapsed"] - t_sim_start
          ))
        }

        n_ac_sel <- as.integer(iv_fit$selected_n_ac)
        n_uc_sel <- as.integer(iv_fit$selected_n_uc)
        selection_diagnostics <- iv_fit$diagnostics
        cor_Ac <- stats::cor(iv_fit$selected_projection$Ac, Ac[, sim])
        return(list(
          muest = iv_fit$selected_fit$psi,
          ci = iv_fit$selected_fit$confidence_interval,
          n_ac = n_ac_sel,
          n_uc = n_uc_sel,
          selection_diagnostics = selection_diagnostics,
          cor_Ac = cor_Ac,
          time_s = proc.time()["elapsed"] - t_sim_start
        ))
      }

      if (method == 'baseline'){
        xmat <- matrix(rep(1, n), ncol = 1)
        colnames(xmat) <- 'Intercept'
      }
      if (method == 'oracle'){
        if (!confounding_mechanism %in% c(3, 7)){
          xmat <- matrix(U[, sim], ncol = 1)
          colnames(xmat) <- 'U'
        } else {
          xmat <- cbind(U1[, sim], U2[, sim])
          colnames(xmat) <- c('U1', 'U2')
        }
      }
      if (method == 'spatialcoord'){
        xmat <- cbind(lat_n_m, lon_n_m)
        colnames(xmat) <- c('Latitude', 'Longitude')
      }
      if (method == 'trueIV-spatialcoord'){
        xmat <- cbind(matrix(Ac[, sim], ncol = 1), lat_n_m, lon_n_m)
        colnames(xmat) <- c('Ac_TPS', 'Latitude', 'Longitude')
      }
      if (method == 'trueIV'){
        xmat <- matrix(Ac[, sim], ncol = 1)
        colnames(xmat) <- 'Ac_TPS'
      }

      non_iv_fit <- tryCatch(
        estimate_truncated_effect_full(
          y = y,
          a = a,
          w = xmat,
          cutoff = cutoff,
          nuisance_args = list(
            sl_library = iv_control$sl_library,
            density_trim = iv_control$density_trim
          ),
          n_grid = iv_control$n_grid,
          bandwidth = iv_control$bandwidth,
          bw_seq = iv_control$bw_seq,
          density_trim = iv_control$density_trim,
          constrain = iv_control$constrain
        ),
        error = function(e) {
          message("non-IV estimator error: ", e$message)
          NULL
        }
      )

      muest <- if (is.null(non_iv_fit)) NA_real_ else non_iv_fit$psi
      ci <- if (is.null(non_iv_fit)) {
        c(NA_real_, NA_real_)
      } else {
        non_iv_fit$confidence_interval
      }

      list(muest  = muest,
           ci     = ci,
           n_ac   = n_ac_sel,
           n_uc   = n_uc_sel,
           selection_diagnostics = selection_diagnostics,
           cor_Ac = cor_Ac,
           time_s = proc.time()["elapsed"] - t_sim_start)

    }, mc.cores = n_cores)

    # ---- aggregate results ----
    safe <- function(r, field, default)
      if (is.list(r)) r[[field]] else default

    muests    <- sapply(results_list, safe, "muest",  NA_real_)
    cis       <- do.call(rbind, lapply(results_list, safe, "ci",
                                       c(NA_real_, NA_real_)))
    n_ac_sels <- as.integer(sapply(results_list, safe, "n_ac", NA_integer_))
    n_uc_sels <- as.integer(sapply(results_list, safe, "n_uc", NA_integer_))
    selection_diagnostics <- lapply(
      results_list,
      safe,
      field = "selection_diagnostics",
      default = NULL
    )
    cor_Ac_sels <- sapply(results_list, safe, "cor_Ac", NA_real_)
    sim_times <- sapply(results_list, safe, "time_s", NA_real_)

    t_method_elapsed <- proc.time()["elapsed"] - t_method_start
    t_total_elapsed  <- proc.time()["elapsed"] - t_simfunc_start
    message(sprintf("[%s] method=%s  sims=%d  method_time=%.1fs  total_time=%.1fs  mean_sim_time=%.2fs",
                    format(Sys.time(), "%H:%M:%S"), method, nsims,
                    t_method_elapsed, t_total_elapsed, mean(sim_times, na.rm = TRUE)))

    # Create dataframe whose first column is a.vals and the rest of cols are muests
    df <- muests
    
    # write results to file
    # Check if file exists
    if (file.exists(filename)){
      # write new sims to file as new columns
      olddf <- read.csv(filename)
      newdf <- cbind(olddf, muests)
      write.csv(newdf, filename, row.names = FALSE)
    }
    # if file for estimates does not exist create it and write results
    else{
      write.csv(df, filename, row.names = FALSE)
    }

    # Save confidence intervals as two separate CSVs (lower / upper bounds),
    # matching the estimates format: nsims rows, one column per batch.
    if (any(!is.na(cis))) {
      filename_ci_lower <- paste0(results_dir, 'conf', confounding_mechanism,
                                  '_', option, '_', method, '_ci_lower.csv')
      filename_ci_upper <- paste0(results_dir, 'conf', confounding_mechanism,
                                  '_', option, '_', method, '_ci_upper.csv')
      df_lower <- data.frame(lower = cis[, 1])
      df_upper <- data.frame(upper = cis[, 2])
      if (file.exists(filename_ci_lower)) {
        write.csv(cbind(read.csv(filename_ci_lower), df_lower),
                  filename_ci_lower, row.names = FALSE)
        write.csv(cbind(read.csv(filename_ci_upper), df_upper),
                  filename_ci_upper, row.names = FALSE)
      } else {
        write.csv(df_lower, filename_ci_lower, row.names = FALSE)
        write.csv(df_upper, filename_ci_upper, row.names = FALSE)
      }
    }

    # Save per-sim wall times (seconds)
  #   filename_time <- paste0(results_dir, 'conf', confounding_mechanism,
  #                           '_', option, '_', method, '_time.csv')
  #   df_time <- data.frame(time_s = sim_times)
  #   if (file.exists(filename_time)) {
  #     write.csv(cbind(read.csv(filename_time), df_time), filename_time, row.names = FALSE)
  #   } else {
  #     write.csv(df_time, filename_time, row.names = FALSE)
  #   }
  #
  #   # Save n_ac/n_uc selections for methods that use basis selection
    if (select_basis && any(!is.na(n_ac_sels))) {
      filename_n_ac <- paste0(results_dir, 'conf', confounding_mechanism,
                              '_', option, '_', method, '_n_ac.csv')
      df_n_ac <- data.frame(n_ac = n_ac_sels)
      if (file.exists(filename_n_ac)) {
        write.csv(cbind(read.csv(filename_n_ac), df_n_ac), filename_n_ac, row.names = FALSE)
      } else {
        write.csv(df_n_ac, filename_n_ac, row.names = FALSE)
      }

      filename_n_uc <- paste0(results_dir, 'conf', confounding_mechanism,
                              '_', option, '_', method, '_n_uc.csv')
      df_n_uc <- data.frame(n_uc = n_uc_sels)
      if (file.exists(filename_n_uc)) {
        olddf_n_uc <- read.csv(filename_n_uc)
        write.csv(cbind(olddf_n_uc, df_n_uc), filename_n_uc, row.names = FALSE)
      } else {
        write.csv(df_n_uc, filename_n_uc, row.names = FALSE)
      }

      if (isTRUE(iv_control$save_diagnostics) &&
          any(lengths(selection_diagnostics) > 0L)) {
        diagnostic_tag <- paste0(
          format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid()
        )
        saveRDS(
          selection_diagnostics,
          file = paste0(
            results_dir, 'conf', confounding_mechanism, '_', option, '_', method,
            '_selection_diagnostics_', diagnostic_tag, '.rds'
          )
        )
      }
    }

    # Save per-sim correlation between estimated Ac and true Ac (IV-TPS and IV-GL only)
    if (any(!is.na(cor_Ac_sels))) {
      filename_cor_Ac <- paste0(results_dir, 'conf', confounding_mechanism,
                                '_', option, '_', method, '_cor_Ac.csv')
      df_cor_Ac <- data.frame(cor_Ac = cor_Ac_sels)
      if (file.exists(filename_cor_Ac)) {
        write.csv(cbind(read.csv(filename_cor_Ac), df_cor_Ac),
                  filename_cor_Ac, row.names = FALSE)
      } else {
        write.csv(df_cor_Ac, filename_cor_Ac, row.names = FALSE)
      }
    }
  }

  invisible(filename)
}

# Function that computes the covariance matrix of the GP
compute_Sigma_GP <- function(distmat, 
                        kappa=2, 
                        rangeu, 
                        rangec,
                        rho = 0.95,
                        sigu = 1, 
                        sigc = 1, 
                        sigz = 1){
  # distmat is the distance matrix
  # kappa is the smoothness parameter
  # rangeu is the range of the GP for the unconfounded part of exposure
  # rangec is the range of the GP for the confounded part of exposure
  # rho is the correlation between the exposure and unmeasured confounder
  # sigu, sigc, sigz are the standard deviations of the Auc, Ac, and U
  # returns the covariance matrix of the GP
  
  n <- nrow(distmat)
  phiu <- rangeu/(2*sqrt(kappa)) # to match Paciorek implementation
  phic <- rangec/(2*sqrt(kappa)) 
  Sigmau <- geoR::matern(u=distmat, phi=phiu, kappa=kappa)
  Sigmac <- geoR::matern(u=distmat, phi=phic, kappa=kappa)
  Sigma <- matrix(0, nrow = 3*n, ncol = 3*n)
  Sigma[1:n, 1:n] <- sigu^2*Sigmau
  Sigma[(n+1):(2*n), (n+1):(2*n)] <- sigc^2*Sigmac
  Sigma[(2*n+1):(3*n), (2*n+1):(3*n)] <- sigz^2*Sigmac
  # # Auc is uncorrelated + indep of Ac and U
  # Sigma[1:n, (n+1):(3*n)] <- 0
  # Sigma[(n+1):(3*n), 1:n] <- 0
  # Ac and U are highly dependent
  Sigma[(n+1):(2*n), (2*n+1):(3*n)] <- rho*sigc*sigz*Sigmac
  Sigma[(2*n+1):(3*n), (n+1):(2*n)] <- rho*sigc*sigz*Sigmac
  return(Sigma)
}

# Function that computes the covariance matrix for two confounders
compute_Sigma_GP_2U <- function(distmat,
                                kappa = 2,
                                rangeu, rangec, rangez1, rangez2,
                                rho1 = 0.9, rho2 = 0.7,
                                sigu = 1, sigc = 1, sigz1 = 1, sigz2 = 1) {
  stopifnot(abs(rho1) <= 1, abs(rho2) <= 1)
  n <- nrow(distmat)
  phi <- function(r) r / (2 * sqrt(kappa))
  
  Ku  <- geoR::matern(u = distmat, phi = phi(rangeu),  kappa = kappa)
  Kc  <- geoR::matern(u = distmat, phi = phi(rangec),  kappa = kappa)
  Kz1 <- geoR::matern(u = distmat, phi = phi(rangez1), kappa = kappa)
  Kz2 <- geoR::matern(u = distmat, phi = phi(rangez2), kappa = kappa)
  
  Sigma <- matrix(0, nrow = 4*n, ncol = 4*n)
  iAuc <- 1:n; iAc <- (n+1):(2*n); iU1 <- (2*n+1):(3*n); iU2 <- (3*n+1):(4*n)
  
  # Auc (independent)
  Sigma[iAuc, iAuc] <- sigu^2 * Ku
  
  # Coregionalized part on Kc (Ac, U1, U2 share it)
  Sigma[iAc, iAc]   <- sigc^2 * Kc
  Sigma[iU1, iU1]   <- (rho1^2) * sigz1^2 * Kc
  Sigma[iU2, iU2]   <- (rho2^2) * sigz2^2 * Kc
  
  Sigma[iAc, iU1] <- rho1 * sigc * sigz1 * Kc
  Sigma[iU1, iAc] <- t(Sigma[iAc, iU1])
  
  Sigma[iAc, iU2] <- rho2 * sigc * sigz2 * Kc
  Sigma[iU2, iAc] <- t(Sigma[iAc, iU2])
  
  # U1–U2 correlation induced by sharing Kc
  Sigma[iU1, iU2] <- (rho1 * rho2 * sigz1 * sigz2) * Kc
  Sigma[iU2, iU1] <- t(Sigma[iU1, iU2])
  
  # Idiosyncratic scales for U1, U2 (their own kernels)
  Sigma[iU1, iU1] <- Sigma[iU1, iU1] + (1 - rho1^2) * sigz1^2 * Kz1
  Sigma[iU2, iU2] <- Sigma[iU2, iU2] + (1 - rho2^2) * sigz2^2 * Kz2
  
  return(Sigma)
}


# Baseline GP data-generating process (supplement Section 4.2, mechanism 1).
# Also used for the reversed-scale mechanism (5) via different theta values.
# Generates 4 independent mean-zero GP fields, then assembles:
#   Auc = 0.6*B1_uc + 0.4*B2_uc
#   Ac  = 1.0*B1_c  + 0.8*B2_c
#   U   = rho1*B1_c + rho2*B2_c + Zu
compute_data_GP <- function(nsims, distmat,
                            theta_B1uc = 0.025, theta_B2uc = 0.01,
                            theta_B1c  = 0.20, theta_B2c  = 0.10,
                            theta_u    = 0.30,
                            rho1 = 0.8, rho2 = 0.6,
                            kappa = 2) {
  n   <- nrow(distmat)
  phi <- function(r) r / (2 * sqrt(kappa))
  # Sample nsims independent draws from N(0, K); always return n x nsims matrix.
  gp  <- function(K) {
    s <- MASS::mvrnorm(nsims, rep(0, n), K)
    if (nsims == 1L) matrix(s, n, 1L) else t(s)
  }
  B1uc <- gp(geoR::matern(u = distmat, phi = phi(theta_B1uc), kappa = kappa))
  B2uc <- gp(geoR::matern(u = distmat, phi = phi(theta_B2uc), kappa = kappa))
  B1c  <- gp(geoR::matern(u = distmat, phi = phi(theta_B1c),  kappa = kappa))
  B2c  <- gp(geoR::matern(u = distmat, phi = phi(theta_B2c),  kappa = kappa))
  Zu   <- gp(geoR::matern(u = distmat, phi = phi(theta_u),    kappa = kappa))
  list(
    Auc = 1.1 * B1uc + 0.7 * B2uc,
    Ac  = 1.0 * B1c  + 0.8 * B2c,
    U   = rho1 * B1c + rho2 * B2c + Zu
  )
}

# Two-confounder GP data-generating process (supplement Section 4.2, mechanism 4).
# Same 4 latent fields as mechanism 1; two confounders each linked to B1_c and B2_c.
#   U1 = rho1*B1_c + rho2*B2_c + Zu1  (Zu1 ~ GP(R(theta_u1)))
#   U2 = rho3*B1_c + rho4*B2_c + Zu2  (Zu2 ~ GP(R(theta_u2)))
compute_data_GP_2U <- function(nsims, distmat,
                               theta_B1uc = 0.025, theta_B2uc = 0.01,
                               theta_B1c  = 0.20, theta_B2c  = 0.10,
                               theta_u1   = 0.30, theta_u2   = 0.20,
                               rho1 = 0.8, rho2 = 0.6,
                               rho3 = 0.5, rho4 = 0.4,
                               kappa = 2) {
  n   <- nrow(distmat)
  phi <- function(r) r / (2 * sqrt(kappa))
  gp  <- function(K) {
    s <- MASS::mvrnorm(nsims, rep(0, n), K)
    if (nsims == 1L) matrix(s, n, 1L) else t(s)
  }
  B1uc <- gp(geoR::matern(u = distmat, phi = phi(theta_B1uc), kappa = kappa))
  B2uc <- gp(geoR::matern(u = distmat, phi = phi(theta_B2uc), kappa = kappa))
  B1c  <- gp(geoR::matern(u = distmat, phi = phi(theta_B1c),  kappa = kappa))
  B2c  <- gp(geoR::matern(u = distmat, phi = phi(theta_B2c),  kappa = kappa))
  Zu1  <- gp(geoR::matern(u = distmat, phi = phi(theta_u1),   kappa = kappa))
  Zu2  <- gp(geoR::matern(u = distmat, phi = phi(theta_u2),   kappa = kappa))
  list(
    Auc = 1.1 * B1uc + 0.7 * B2uc,
    Ac  = 1.0 * B1c  + 0.8 * B2c,
    U1  = rho1 * B1c + rho2 * B2c + Zu1,
    U2  = rho3 * B1c + rho4 * B2c + Zu2
  )
}

# =============================================================================
# Data-generating mechanisms (October 2026)
#
# Every mechanism fixes the confounded share s = Var(A^c) / Var(A) in each
# replicate by rescaling A^uc (A^c and U are untouched). The ground truth is
# therefore the same for every mechanism in R^2 units: the confounded
# variation explains s of Var(A). The selection rule's reference
# (r2_reference = 0.7 in ../R/basis_selection.R) is valid whenever s < 0.7.
# =============================================================================

confounded_share <- 0.5

# Rescale each column of Auc so that Var(Ac) / Var(Ac + Auc) = s.
rescale_to_share <- function(Ac, Auc, s) {
  for (i in seq_len(ncol(Auc))) {
    Auc[, i] <- Auc[, i] * sqrt(var(Ac[, i]) * (1 - s) / (s * var(Auc[, i])))
  }
  Auc
}

# Random combinations of basis columns `idx`, one per replicate. Coefficient
# standard deviations decay as (1, 2, 3, ...)^(-decay) from the first
# (largest-scale) column of idx; decay = 0 gives equal weights.
draw_basis_field <- function(B, idx, nsims, decay = 0) {
  sds <- seq_along(idx)^(-decay)
  coef <- matrix(rnorm(length(idx) * nsims, sd = sds), length(idx), nsims)
  B[, idx, drop = FALSE] %*% coef
}

# Basis mechanisms (1: TPS, 2: GL).
# A^c and U live in the n_c largest-scale columns; A^uc lives in the remaining
# columns with coefficient sd decaying as j^(-decay). decay = 0.5 makes the
# exposure about as smooth as PM2.5: B^c reaches 70% of Var(A) at ~25 columns
# (vs ~370 with equal weights), comparable to the data application.
# Coefficients are scaled by sc = sqrt(n / n_c) so A^c and U have average
# per-unit variance ~1 (B is orthonormal, so unscaled fields would have
# variance ~n_c / n and almost no confounding).
#   Ac  = B[, 1:n_c] (sc a1)
#   U   = B[, 1:n_c] (sc (rho1 a1 + rho2 a2))
#   Auc = B[, (n_c+1):m] (decaying coefficients), rescaled to share s
compute_data_basis <- function(B, nsims, n_c = 10L, s = confounded_share,
                               decay = 0.5, rho1 = 0.8, rho2 = 0.6) {
  m <- ncol(B)
  idx_c <- seq_len(n_c)
  sc <- sqrt(nrow(B) / n_c)
  a1 <- sc * matrix(rnorm(n_c * nsims), n_c, nsims)
  a2 <- sc * matrix(rnorm(n_c * nsims), n_c, nsims)
  Ac <- B[, idx_c, drop = FALSE] %*% a1
  U  <- B[, idx_c, drop = FALSE] %*% (rho1 * a1 + rho2 * a2)
  Auc <- draw_basis_field(B, seq.int(n_c + 1L, m), nsims, decay)
  list(Ac = Ac, Auc = rescale_to_share(Ac, Auc, s), U = U)
}

compute_data_TPS_basis <- function(B_tps_full, nsims, ...) {
  compute_data_basis(B_tps_full, nsims, ...)
}

compute_data_GL_basis <- function(B_gl_full, nsims, ...) {
  compute_data_basis(B_gl_full, nsims, ...)
}

# Two-confounder basis mechanisms (3: TPS, 7: GL). Same partition and A^uc as
# compute_data_basis(); the combined A^c share is fixed at s.
#   B1c = B[, 1:n_c] a1,  B2c = B[, 1:n_c] a2,  Ac = B1c + 0.8 B2c
#   U1  = rho1 B1c + rho2 B2c + B[, 1:n_c] e1
#   U2  = rho3 B1c + rho4 B2c + B[, 1:n_c] e2
compute_data_basis_2U <- function(B, nsims, n_c = 10L, s = confounded_share,
                                  decay = 0.5, rho1 = 0.8, rho2 = 0.6,
                                  rho3 = 0.5, rho4 = 0.4) {
  m <- ncol(B)
  Bc <- B[, seq_len(n_c), drop = FALSE]
  sc <- sqrt(nrow(B) / n_c)  # per-unit variance ~1, as in compute_data_basis()
  draw_c <- function() Bc %*% (sc * matrix(rnorm(n_c * nsims), n_c, nsims))
  B1c <- draw_c(); B2c <- draw_c()
  Ac <- 1.0 * B1c + 0.8 * B2c
  U1 <- rho1 * B1c + rho2 * B2c + draw_c()
  U2 <- rho3 * B1c + rho4 * B2c + draw_c()
  Auc <- draw_basis_field(B, seq.int(n_c + 1L, m), nsims, decay)
  list(Ac = Ac, Auc = rescale_to_share(Ac, Auc, s), U1 = U1, U2 = U2)
}

compute_data_TPS_2U <- function(B_tps_full, nsims, ...) {
  compute_data_basis_2U(B_tps_full, nsims, ...)
}

compute_data_GL_2U <- function(B_gl_full, nsims, ...) {
  compute_data_basis_2U(B_gl_full, nsims, ...)
}

# Scale-reversed basis mechanisms (4: TPS, 8: GL): the mirror of mechanism
# 1/2. A^uc lives in the n_uc largest-scale columns; A^c and U live in all
# remaining (smaller-scale) columns, with coefficient sd decaying as j^(-decay)
# from the first confounded column (n_uc + 1). simfunc() reverses the basis
# for these mechanisms, so in search order the confounded block is the first
# m - n_uc columns (true n_ac = m - n_uc). The decay puts most confounded
# variance next to the boundary: with equal weights each confounded column
# carries ~0.06% of Var(A), so moving one across the boundary is undetectable
# and the walk would run deep into the confounded block.
#   Auc = B[, 1:n_uc] au, rescaled to share s
#   Ac  = B[, (n_uc+1):m] a1,  U = B[, (n_uc+1):m] (rho1 a1 + rho2 a2)
compute_data_basis_reversed <- function(B, nsims, n_uc = 20L,
                                        s = confounded_share, decay = 0.5,
                                        rho1 = 0.8, rho2 = 0.6) {
  m <- ncol(B)
  idx_c <- seq.int(n_uc + 1L, m)
  k <- length(idx_c)
  # scaled so A^c and U have average per-unit variance ~1
  sds <- seq_len(k)^(-decay) * sqrt(nrow(B) / sum(seq_len(k)^(-2 * decay)))
  a1 <- matrix(rnorm(k * nsims, sd = sds), k, nsims)
  a2 <- matrix(rnorm(k * nsims, sd = sds), k, nsims)
  Ac <- B[, idx_c, drop = FALSE] %*% a1
  U  <- B[, idx_c, drop = FALSE] %*% (rho1 * a1 + rho2 * a2)
  Auc <- draw_basis_field(B, seq_len(n_uc), nsims)
  list(Ac = Ac, Auc = rescale_to_share(Ac, Auc, s), U = U)
}

compute_data_TPS_reversed <- function(B_tps_full, nsims, ...) {
  compute_data_basis_reversed(B_tps_full, nsims, ...)
}

compute_data_GL_reversed <- function(B_gl_full, nsims, ...) {
  compute_data_basis_reversed(B_gl_full, nsims, ...)
}

rleroux_univariate <- function(W, rho_sp, tau_sp = 1) {
  n  <- nrow(W)
  D  <- Diagonal(n, rowSums(W))
  Qs <- tau_sp * ((1 - rho_sp) * Diagonal(n) + rho_sp * (D - W))
  Qs <- forceSymmetric(Qs)
  cf <- Cholesky(Qs, LDL = FALSE, perm = TRUE, super = TRUE)
  Z  <- rnorm(n)
  Yp <- solve(cf, Z, system = "L")
  Xp <- solve(cf, Yp, system = "Lt")
  Xp[order(cf@perm)]
}

rleroux_bivariate <- function(W, rho_sp, tau_sp, R){
  n  <- nrow(W)
  D  <- Diagonal(n, rowSums(W))
  Qs <- tau_sp * ((1 - rho_sp) * Diagonal(n) + rho_sp * (D - W))
  Qs <- forceSymmetric(Qs)

  Q  <- kronecker(R, Qs)           # joint precision (2n x 2n), sparse
  cf <- Cholesky(Q, LDL=FALSE, perm=TRUE, super=TRUE)

  Z  <- rnorm(2*n)
  Yp <- solve(cf, Z, system="L")
  Xp <- solve(cf, Yp, system="Lt")
  X  <- Xp[order(cf@perm)]
  list(phi1 = X[1:n], phi2 = X[(n+1):(2*n)])
}

# Bivariate Leroux CAR mechanism (mechanism 5; basis misspecification).
#   (U, B1_c) ~ bivariate Leroux(rho_c), cross-correlation 0.7
#   B2_c ~ univariate Leroux(rho_c); B1_uc, B2_uc ~ univariate Leroux(rho_uc)
#   Ac = 1.0 B1_c + 0.8 B2_c;  Auc = 1.1 B1_uc + 0.7 B2_uc, rescaled to share s
# rho_c = 0.99 makes A^c very smooth (half its variance in the first ~3
# TPS/GL columns) and rho_uc = 0.6 keeps A^uc rougher (half by ~40 columns),
# so the reference is ~30 columns. Raising rho_uc instead (with rho_c = 0.6)
# would make A^uc smoother than A^c and reverse the premise.
# Each field is centered and scaled within the replicate: near rho = 1 the
# Leroux model is close to intrinsic, so raw draws carry a large random
# overall level (A's mean ranged from -8 to 15 across replicates), leaving
# almost no units near the cutoff.
compute_data_leroux <- function(W, nsims, rho_c = 0.99, rho_uc = 0.6,
                                s = confounded_share) {
  n <- nrow(W)
  stopifnot(ncol(W) == n)
  R_c <- solve(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  standardize <- function(v) (v - mean(v)) / sd(v)
  out <- list(Auc = matrix(NA_real_, n, nsims),
              Ac  = matrix(NA_real_, n, nsims),
              U   = matrix(NA_real_, n, nsims))
  for (i in seq_len(nsims)) {
    smp <- rleroux_bivariate(W = W, rho_sp = rho_c, tau_sp = 1, R = R_c)
    out$U[, i]  <- standardize(smp$phi1)
    out$Ac[, i] <- 1.0 * standardize(smp$phi2) +
      0.8 * standardize(rleroux_univariate(W = W, rho_sp = rho_c))
    out$Auc[, i] <- 1.1 * standardize(rleroux_univariate(W = W, rho_sp = rho_uc)) +
      0.7 * standardize(rleroux_univariate(W = W, rho_sp = rho_uc))
  }
  out$Auc <- rescale_to_share(out$Ac, out$Auc, s)
  out
}

# Coordinate-based mechanism (mechanism 6).
# U = a1 sin(pi lat) + a2 cos(pi long) + a3 lat long, a_k ~ N(0, 1) per
# replicate (normalized coordinates), centered; Ac = U.
# B1_uc, B2_uc ~ independent Matern GP fields (ranges 0.10, 0.05),
# Auc = 1.1 B1_uc + 0.7 B2_uc, rescaled to share s.
compute_data_spatialcoord <- function(lat, long, nsims, distmat,
                                      s = confounded_share) {
  n      <- length(lat)
  lat_n  <- (lat  - min(lat))  / (max(lat)  - min(lat))
  long_n <- (long - min(long)) / (max(long) - min(long))
  f <- cbind(sin(pi * lat_n), cos(pi * long_n), lat_n * long_n)

  phi_fn <- function(r, kappa = 2) r / (2 * sqrt(kappa))
  gp     <- function(K) {
    draws <- MASS::mvrnorm(nsims, rep(0, n), K)
    if (nsims == 1L) matrix(draws, n, 1L) else t(draws)
  }
  Auc <- 1.1 * gp(geoR::matern(u = distmat, phi = phi_fn(0.10), kappa = 2)) +
    0.7 * gp(geoR::matern(u = distmat, phi = phi_fn(0.05), kappa = 2))

  U <- f %*% matrix(rnorm(3 * nsims), 3, nsims)
  U <- sweep(U, 2L, colMeans(U))
  list(Auc = rescale_to_share(U, Auc, s), Ac = U, U = U)
}

# One entry point for all mechanisms, shared by simfunc() and
# compute_true_psi.R so both generate identical data from the same seed.
generate_mechanism_data <- function(confounding_mechanism, nsims, B_tps_full,
                                    B_gl_full, W, lat, lon, distmat) {
  switch(
    as.character(confounding_mechanism),
    "1" = compute_data_TPS_basis(B_tps_full, nsims),
    "2" = compute_data_GL_basis(B_gl_full, nsims),
    "3" = compute_data_TPS_2U(B_tps_full, nsims),
    "4" = compute_data_TPS_reversed(B_tps_full, nsims),
    "5" = compute_data_leroux(W, nsims),
    "6" = compute_data_spatialcoord(lat, lon, nsims, distmat),
    "7" = compute_data_GL_2U(B_gl_full, nsims),
    "8" = compute_data_GL_reversed(B_gl_full, nsims),
    stop("confounding_mechanism must be between 1 and 8.")
  )
}

# Mechanisms whose basis is reversed (smallest scale first) for selection.
reversed_mechanisms <- c(4L, 8L)
