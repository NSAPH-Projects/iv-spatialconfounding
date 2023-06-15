source("funcs.R")

outsim <- sim(
    n = 5, outcome = "linear",
    rhox = c(0.9, 0.001),
    decomposition = "nested",
    quiet = TRUE
)
analnested <- analysis(
    n = 5, A = outsim$A,
    X = outsim$X,
    Y = outsim$Y,
    groups = outsim$groups,
    decomposition = "nested",
    quiet = TRUE,
    return_decomps = TRUE
)
analspectral <- analysis(
    n = 5, A = outsim$A,
    X = outsim$X,
    Y = outsim$Y,
    groups = outsim$groups,
    decomposition = "spectral",
    quiet = TRUE,
    return_decomps = TRUE
)
nest <- analnested$nest
spec <- analspectral$spec
