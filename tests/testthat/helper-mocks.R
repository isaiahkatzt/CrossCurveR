# Scope mocks to the real package/dependency namespaces and restore on exit.
with_temp_bindings <- function(bindings, expr) {
  expr <- substitute(expr)
  external <- list(
    Matrix = "bdiag",
    purrr = "transpose",
    vars = "VARselect",
    tsDyn = c("VECM", "rank.test", "coefA", "coefB")
  )

  for (package in names(external)) {
    selected <- intersect(names(bindings), external[[package]])
    if (length(selected)) {
      do.call(testthat::local_mocked_bindings, c(
        bindings[selected], list(.package = package, .env = environment())
      ))
      bindings <- bindings[setdiff(names(bindings), selected)]
    }
  }

  if (length(bindings)) {
    do.call(testthat::local_mocked_bindings, c(
      bindings, list(.package = "CrossCurveR", .env = environment())
    ))
  }

  eval(expr, envir = parent.frame())
}
