#' @importFrom BiocManager install
#' @importFrom remotes install_github
#'
NULL

# utilitis -----------------------------------------------------------------

# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

# Set a default value if an object is NOT null
#
# @param lhs An object to set if it's NOT null
# @param rhs The value to provide if x is NOT null
#
# @return lhs if lhs is null, else rhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}

# Check to ensure a package is installed
#
# @param package Name of pacakge to check
# @param repository Repository that package is available on;
# choose from 'bioconductor', 'github', or 'cran'
# @param ... Extra parameters passed to BiocManager::install, remotes::install_github, or install.packages, depending on \code{repository}
#
#' @importFrom utils menu install.packages
#
CheckPackage <- function(package, repository, ...) {
  if (!requireNamespace(package = basename(path = package), quietly = TRUE)) {
    if (interactive()) {
      message('Package ', package, ' is not yet installed')
      message('Install now?')
      choice <- menu(choices = c('yes', 'no'))
      if (choice == 1) {
        repository <- match.arg(
          arg = tolower(x = repository),
          choices = c('github', 'bioconductor', 'cran')
        )
        switch(
          EXPR = repository,
          'github' = remotes::install_github(repo = package, ...),
          'bioconductor' = BiocManager::install(pkgs = package, ...),
          'cran' = install.packages(pkgs = package, ...),
          stop('Unknown repository ', repository, call. = FALSE)
        )
        return(invisible(x = NULL))
      }
    }
    stop('Unable to find package ', package, ', please install', call. = FALSE)
  }
}

# Check if a matrix is empty
#
# Takes a matrix and asks if it's empty (either 0x0 or 1x1 with a value of NA)
#
# @param x A matrix
#
# @return Whether or not \code{x} is empty
#
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}

# Cell cycle scoring -----------------------------------------------------------------

#' Get back ground list
#'
GetBackgroundList <- function(features, N, expr.bin) {
  res <- list()
  features.bin.tab <- table(expr.bin[features])
  for (i in 1:N) {
    res[[i]] <- unlist(lapply(names(features.bin.tab), function(b) {
      sel <- which(expr.bin == as.numeric(b) & !(names(expr.bin) %in% features))
      sample(names(expr.bin)[sel], features.bin.tab[b])
    }))
  }
  return(res)
}

#' Calculate enrich z-score
#'
EnrichScore <- function(expr, features, bg.list) {
  features.mean <- apply(expr[features, ], 2, mean)
  bg.mean <- sapply(1:length(bg.list), function(i) apply(expr[bg.list[[i]], ], 2, mean))
  return((features.mean - apply(bg.mean, 1, mean)) / apply(bg.mean, 1, sd))
}

# # AlignCCA -----------------------------------------------------------------
#
# # bicor helper function to standardize the two vectors and perform common
# # calculations.
# #
# # @author Patrick Roelli
# # @param x Vector to prep
# # @param verbose If TRUE, prints a warning when falling back on robust
# # standardization when MAD(x) is 0.
# #
# # @return returns the prepped vector
# #
# BicorPrep <- function(x, verbose = FALSE){
#   if (stats::mad(x) == 0) {
#     if (verbose){
#       warning('mad == 0, using robust standardization')
#     }
#     xat <- x - mean(x = x)
#     xab <- sqrt(x = sum((x - mean(x = x)) ^ 2))
#     result <- xat / xab
#     return (result)
#   } else {
#     ua <- (x - stats::median(x = x)) /
#       (9 * stats::mad(x = x) *
#          stats::qnorm(p = 0.75))
#     i.x <- ifelse(test = ua <= -1 | ua >= 1, yes = 0, no = 1)
#     wax <- ((1 - (ua ^ 2)) ^ 2) * i.x
#     xat <- (x - stats::median(x = x)) * wax
#     xab <- sqrt(x = sum(xat ^ 2))
#     result <- xat / xab
#     return(result)
#   }
# }
#
# # Calculate the biweight midcorrelation (bicor) of two vectors using
# # implementation described in Langfelder, J Stat Sotfw. 2012. If MAD of one of
# # the two vectors is 0, falls back on robust standardization.
# #
# # @author Patrick Roelli
# # @param x First vector
# # @param y Second vector
# #
# # @return returns the biweight midcorrelation of x and y
# #
# BiweightMidcor <- function(x, y){
#   resx <- BicorPrep(x)
#   resy <- BicorPrep(y)
#   result <- sum(resx * resy)
#   return(result)
# }
#
# # Calculate position along a defined reference range for a given vector of numerics. Will range from 0 to 1.
# #
# # @param x      Vector of numeric type
# # @param lower  Lower end of reference range
# # @param upper  Upper end of reference range
# #
# #' @importFrom stats quantile
# #
# # @return       Returns a vector that describes the position of each element in x along the defined reference range
# ReferenceRange <- function(x, lower = 0.025, upper = 0.975) {
#   return((x - quantile(x = x, probs = lower)) /
#            (quantile(x = x, probs = upper) - quantile(x = x, probs = lower)))
# }
