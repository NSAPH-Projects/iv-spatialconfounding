# Load simulation and estimation functions
source("../funcs.R")
source("../R/basis_selection.R")
source("../R/truncated_effect.R")
source("../R/outer_crossfit.R")

library(dplyr)
library(geosphere)
library(foreach)
library(doParallel)
library(Matrix)

load("sim.RData")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3L || length(args) > 10L) {
  stop(
    paste(
      "Usage: Rscript run_simfunc.R",
      "<nsims> <mechanism> <linear|nonlinear>",
      "[select_basis] [n_cores] [results_dir] [seed] [methods] [core_fraction] [alpha]"
    )
  )
}

parse_positive_integer <- function(value, name) {
  out <- suppressWarnings(as.integer(value))
  if (length(out) != 1L || is.na(out) || out < 1L) {
    stop(name, " must be a positive integer.")
  }
  out
}

parse_flag <- function(value, name) {
  value <- toupper(value)
  if (!value %in% c("TRUE", "FALSE")) {
    stop(name, " must be TRUE or FALSE.")
  }
  identical(value, "TRUE")
}

parse_fraction <- function(value, name) {
  out <- suppressWarnings(as.numeric(value))
  if (length(out) != 1L || is.na(out) || !is.finite(out) ||
      out <= 0 || out >= 1) {
    stop(name, " must be a finite number strictly between 0 and 1.")
  }
  out
}

parse_positive_number <- function(value, name) {
  out <- suppressWarnings(as.numeric(value))
  if (length(out) != 1L || is.na(out) || !is.finite(out) || out <= 0) {
    stop(name, " must be a positive finite number.")
  }
  out
}

nsims <- parse_positive_integer(args[[1]], "nsims")
confounding_mechanism <- parse_positive_integer(
  args[[2]],
  "confounding_mechanism"
)

if (!confounding_mechanism %in% 1:8) {
  stop("confounding_mechanism must be between 1 and 8.")
}

option <- match.arg(args[[3]], c("linear", "nonlinear"))

select_basis <- if (length(args) >= 4L) {
  parse_flag(args[[4]], "select_basis")
} else {
  TRUE
}

n_cores <- if (length(args) >= 5L) {
  parse_positive_integer(args[[5]], "n_cores")
} else {
  1L
}

results_dir <- if (length(args) >= 6L) {
  args[[6]]
} else {
  "results_manuscript"
}

seed <- if (length(args) >= 7L) {
  parse_positive_integer(args[[7]], "seed")
} else {
  123L
}

methods <- if (length(args) >= 8L) {
  trimws(strsplit(args[[8]], ",", fixed = TRUE)[[1]])
} else {
  c("oracle", "IV-TPS", "IV-GraphLaplacian")
}

core_fraction <- if (length(args) >= 9L) {
  parse_fraction(args[[9]], "core_fraction")
} else {
  0.9 # SMW change 0825
}

alpha <- if (length(args) >= 10L) {
  parse_positive_number(args[[10]], "alpha")
} else {
  1
}

allowed_methods <- c(
  "baseline",
  "oracle",
  "spatialcoord",
  "trueIV",
  "trueIV-spatialcoord",
  "IV-TPS",
  "IV-GraphLaplacian",
  "IV-TPS-spatialcoord",
  "IV-GraphLaplacian-spatialcoord"
)

invalid_methods <- setdiff(methods, allowed_methods)
if (length(invalid_methods) > 0L) {
  stop("Unsupported methods: ", paste(invalid_methods, collapse = ", "))
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(results_dir)) {
  stop("Could not create results directory: ", results_dir)
}

# simfunc() expects a trailing path separator.
results_dir <- paste0(
  normalizePath(results_dir, mustWork = TRUE),
  .Platform$file.sep
)

set.seed(seed)

message("Simulation configuration:")
print(list(
  nsims = nsims,
  confounding_mechanism = confounding_mechanism,
  option = option,
  select_basis = select_basis,
  n_cores = n_cores,
  results_dir = results_dir,
  seed = seed,
  methods = methods,
  core_fraction = core_fraction,
  alpha = alpha
))

simfunc(
  nsims = nsims,
  lat = simlist$lat,
  lon = simlist$lon,
  confounding_mechanism = confounding_mechanism,
  option = option,
  methods = methods,
  B_tps_full = simlist$B_tps_full,
  B_gl_full = simlist$B_gl_full,
  statemat = simlist$statemat,
  W = simlist$W,
  select_basis = select_basis,
  n_cores = n_cores,
  results_dir = results_dir,
  iv_control = list(core_fraction = core_fraction, alpha = alpha)
)
