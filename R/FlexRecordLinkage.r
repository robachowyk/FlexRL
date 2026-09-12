
# ============================================================================
# FlexRL — flexible probabilistic record linkage with diagnostics for
# downstream inference on linked data.
# ============================================================================


#' Simulate two linked data sources for record linkage benchmarking
#'
#' Creates two synthetic data sources of given sizes sharing a given number
#' of common entities ("links"), each described by a set of Partially
#' Identifying Variables (PIVs). For every PIV you choose the number of
#' possible values, the proportion of mistakes and of missing values, and
#' whether it is stable over time, flexible (may change but change is not
#' modelled), or structured (expected to change over time, with a survival
#' model for the hazard of change). For structured PIVs, `enforceEstimability`
#' forces half of the linked pairs to have a near-zero time gap, which helps
#' separate "mistake" from "genuine change over time" when fitting the model.
#'
#' @param PIVs_config Named list, one entry per PIV, each a list with:
#'   `dynamics` (`"stable"`, `"flexible"`, or `"structured"`),
#'   `boundMistakes` (length-2 numeric/NA, upper bound on the mistake
#'   probability in file 1 / file 2), `fixMistakes` (length-2 numeric/NA,
#'   mistake probability fixed to this value in file 1 / file 2), and, only
#'   for `dynamics = "structured"`, `condHazardCov` (a list with `cov1` and
#'   `cov2`, the names of covariates in file 1 / file 2 used to model the
#'   hazard of change).
#' @param Nval Integer vector, number of unique values per PIV (same order
#'   as `PIVs_config`).
#' @param NRecords Integer vector of length 2, number of records to generate
#'   in file 1 and file 2 (file 2 must be the larger of the two).
#' @param Nlinks Integer, number of records shared between the two files.
#' @param Pmistake Named list (one entry per PIV) of length-2 numeric vectors,
#'   proportion of mistakes to introduce in file 1 / file 2.
#' @param Pmissing Named list (one entry per PIV) of length-2 numeric vectors,
#'   proportion of missing values to introduce in file 1 / file 2.
#' @param condHazard_params Named list (one entry per PIV) of numeric vectors
#'   with the log-hazard coefficients (baseline hazard, then the coefficients
#'   for the `cov1` covariates, then for the `cov2` covariates); only used for
#'   `dynamics = "structured"` PIVs, and must have length
#'   `1 + length(cov1) + length(cov2)`.
#' @param enforceEstimability Logical; if `TRUE`, half of the linked pairs are
#'   given a near-zero time gap to help estimate the instability parameters.
#'
#' @return A list with: `dataSet1`, `dataSet2` the two simulated (encoded)
#'   data frames, `Nvalues` number of unique values per PIV, `TimeDifference`
#'   time gap between linked records (`NA` if no PIV is structured),
#'   `proba_same_H` matrix (n links x n PIVs) of probabilities that the true
#'   values coincide, `true_pairs` data frame with the true `1`/`2` indices
#'   of the linked records
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V2 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V3 = list(dynamics = "flexible",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V4 = list(dynamics = "structured",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(0.03,0.03),
#'                                 condHazardCov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())
#'                                 )
#' )
#' Nval  <- c(6, 7, 8, 9)
#' Pmistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02))
#' Pmissing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
#' condHazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = c(0.7,0.6,0.5))
#'
#' GenData <- DataCreation(
#'   PIVs_config, Nval, NRecords = c(400, 600), Nlinks = 300,
#'   Pmistake, Pmissing, condHazard_params, enforceEstimability = TRUE
#' )
#' str(GenData, max.level = 1)
DataCreation <- function(PIVs_config, Nval, NRecords, Nlinks, Pmistake, Pmissing,
                         condHazard_params, enforceEstimability) {

  if (NRecords[2] < NRecords[1]) {
    stop("`NRecords[2]` (file B) must be >= `NRecords[1]` (file A).", call. = FALSE)
  }
  if (Nlinks > min(NRecords)) {
    stop("`Nlinks` cannot exceed the size of the smallest file.", call. = FALSE)
  }

  PIVs <- names(PIVs_config)
  PIVs_stable <- sapply(PIVs_config, function(x) x$dynamics != "structured")
  PIVs_conditionalHazardVariables <- sapply(PIVs_config, function(x) {
    if (x$dynamics == "structured") x$condHazardCov else FALSE
  })
  modelDynaPIVs <- which(!PIVs_stable)

  Pmistake1 <- sapply(Pmistake, function(x) x[1])
  Pmistake2 <- sapply(Pmistake, function(x) x[2])
  Pmissing1 <- sapply(Pmissing, function(x) x[1])
  Pmissing2 <- sapply(Pmissing, function(x) x[2])

  # Simulate the true PIV values for both files
  for (i in 1:2) {
    dataSet <- c()
    for (u in seq_along(Nval)) {
      xp    <- exp(0.23 * (0:(Nval[u] - 1)))
      probx <- xp / sum(xp)
      dataSet <- cbind(dataSet, sample(1:Nval[u], NRecords[i], replace = TRUE, prob = probx))
    }
    dataSet <- as.data.frame(dataSet)
    names(dataSet) <- PIVs
    assign(paste0("dataSet", i), dataSet)
  }

  # First `Nlinks` records of file 2 are copies of the first `Nlinks` of file 1
  dataSet2[1:Nlinks, ] <- dataSet1[1:Nlinks, ]

  # Introduce mistakes
  for (x in seq_len(ncol(dataSet1))) {
    biased <- as.logical(stats::rbinom(nrow(dataSet1), 1, Pmistake1[x]))
    if (any(biased)) {
      dataSet1[, x][biased] <- sapply(dataSet1[, x][biased], function(i) sample((1:Nval[x])[-c(i)], 1))
    }
  }
  for (x in seq_len(ncol(dataSet2))) {
    biased <- as.logical(stats::rbinom(nrow(dataSet2), 1, Pmistake2[x]))
    if (any(biased)) {
      dataSet2[, x][biased] <- sapply(dataSet2[, x][biased], function(i) sample((1:Nval[x])[-c(i)], 1))
    }
  }

  # Introduce missing values
  for (x in seq_len(ncol(dataSet1))) {
    biased <- as.logical(stats::rbinom(nrow(dataSet1), 1, Pmissing1[x]))
    if (any(biased)) dataSet1[, x][biased] <- NA
  }
  for (x in seq_len(ncol(dataSet2))) {
    biased <- as.logical(stats::rbinom(nrow(dataSet2), 1, Pmissing2[x]))
    if (any(biased)) dataSet2[, x][biased] <- NA
  }

  dataSet1$change <- FALSE
  dataSet2$change <- FALSE

  if (any(!PIVs_stable)) {

    # Registration times (file 2 always registered after file 1)
    dataSet1$date <- stats::runif(nrow(dataSet1), 0, 3)
    dataSet2$date <- stats::runif(nrow(dataSet2), 3, 6)

    if (enforceEstimability) {
      nullTimeDiff <- as.integer(Nlinks / 2)
      dataSet1[1:nullTimeDiff, "date"] <- stats::runif(nullTimeDiff, 0.00, 0.01)
      dataSet2[1:nullTimeDiff, "date"] <- stats::runif(nullTimeDiff, 0.00, 0.01)
    }

    TimeDifference <- abs(dataSet2[1:Nlinks, "date"] - dataSet1[1:Nlinks, "date"])
    intercept <- rep(1, Nlinks)
    proba_same_H <- matrix(1, Nlinks, length(PIVs))

    for (k in modelDynaPIVs) {

      hasCov <- !is.logical(PIVs_conditionalHazardVariables[[k]]) &&
        any(!sapply(PIVs_conditionalHazardVariables[[k]], is.null))

      if (hasCov) {
        for (covariate in PIVs_config[[k]]$condHazardCov$cov1) {
          dataSet1[, covariate] <- stats::rnorm(nrow(dataSet1), 1, 1)
        }
        for (covariate in PIVs_config[[k]]$condHazardCov$cov2) {
          dataSet2[, covariate] <- stats::rnorm(nrow(dataSet2), 2, 1)
        }
        cov <- cbind(
          intercept,
          dataSet1[1:Nlinks, PIVs_config[[k]]$condHazardCov$cov1, drop = FALSE],
          dataSet2[1:Nlinks, PIVs_config[[k]]$condHazardCov$cov2, drop = FALSE]
        )
      } else {
        cov <- cbind(intercept)
      }

      if (length(condHazard_params[[k]]) != ncol(as.matrix(cov))) {
        stop(sprintf(
          "`condHazard_params[['%s']]` must have length 1 + length(cov1) + length(cov2).",
          PIVs[k]
        ), call. = FALSE)
      }

      proba_same_H[, k] <- exp(-exp(as.matrix(cov) %*% log(condHazard_params[[k]])) * TimeDifference)

      # Generate instability: for each linked pair, decide (via the survival
      # probability above) whether the true value changed between the two
      # registrations
      for (i in seq_len(Nlinks)) {
        is_not_changing <- stats::rbinom(1, 1, proba_same_H[i, k])
        if (!is_not_changing && !is.na(dataSet1[i, x])) {
          dataSet2[i, x] <- sample((1:Nval[x])[-c(dataSet1[i, x])], 1)
          dataSet2[i, "change"] <- TRUE
        }
      }
    }
  } else {
    TimeDifference <- NA
    proba_same_H <- NA
  }

  # Recode the PIVs
  levels_PIVs <- lapply(PIVs, function(x) levels(factor(as.character(c(dataSet1[, x], dataSet2[, x])))))
  for (i in seq_along(PIVs)) {
    dataSet1[, PIVs[i]] <- as.numeric(factor(as.character(dataSet1[, PIVs[i]]), levels = levels_PIVs[[i]]))
    dataSet2[, PIVs[i]] <- as.numeric(factor(as.character(dataSet2[, PIVs[i]]), levels = levels_PIVs[[i]]))
  }
  Nvalues <- sapply(levels_PIVs, length)

  dataSet1$localID <- seq_len(nrow(dataSet1))
  dataSet2$localID <- seq_len(nrow(dataSet2))
  dataSet1$entityID <- c(seq_len(Nlinks), seq.int(Nlinks + 1, nrow(dataSet1)))
  dataSet2$entityID <- c(seq_len(Nlinks), seq.int(nrow(dataSet1) + 1, nrow(dataSet1) + nrow(dataSet2) - Nlinks))
  dataSet1$source <- "1"
  dataSet2$source <- "2"

  true_Delta <- data.frame(`1` = seq_len(Nlinks), `2` = seq_len(Nlinks), check.names = FALSE)

  list(
    dataSet1 = dataSet1, dataSet2 = dataSet2, Nvalues = Nvalues,
    TimeDifference = TimeDifference, proba_same_H = proba_same_H, true_pairs = true_Delta
  )
}

#' Book-keeping data frame for parameters of PIVs dynamics
#'
#' Internal helper used by [stEM()] to accumulate, across Gibbs iterations,
#' the covariates, true-value agreement indicator, and time gaps needed to
#' re-estimate the survival (hazard) parameters of a dynamic structured PIV.
#'
#' @param nCoefUnstable Integer, number of hazard coefficients for this PIV
#'   (1 for the baseline hazard, plus one per covariate from file A and file B).
#' @param stable Logical, whether this PIV is stable
#'   (`dynamics != "structured"`).
#'
#' @return An empty data frame with `nCoefUnstable + 2` columns (the
#'   covariates/intercept, plus `Hequal` and `times`) if `stable` is `FALSE`;
#'   `NULL` if the PIV is stable (nothing to accumulate).
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V2 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V3 = list(dynamics = "flexible",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V4 = list(dynamics = "structured",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(0.03,0.03),
#'                                 condHazardCov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())
#'                                 )
#' )
#' PIVs_stable <- sapply(PIVs_config, function(x) x$dynamics != "structured")
#' nCoefUnstable = c(0,0,0,3)
#' Valpha <- mapply(createDataAlpha, nCoefUnstable = nCoefUnstable,
#'                   stable = PIVs_stable, SIMPLIFY = FALSE)
createDataAlpha <- function(nCoefUnstable, stable) {
  if (stable) return(invisible(NULL))
  nCoef <- nCoefUnstable + 2  # covariates/intercept + Hequal + times
  data.frame(matrix(nrow = 0, ncol = nCoef))
}

#' Log number of possible linkage configurations
#'
#' Computes `log(nB! / (nB - sumD)!)`, i.e. the log number of ways to choose
#' `sumD` ordered links among `Brecords` records of the larger file; used in
#' [loglik()] to normalise the likelihood of the linkage matrix.
#'
#' @param Brecords Integer, number of records in the larger data source (B).
#' @param sumD Integer, number of currently linked records.
#'
#' @return Numeric, `sum(log((Brecords - sumD + 1):Brecords))`, or `0` if `sumD == 0`.
#' @export
#'
#' @examples
#' logPossibleConfig(Brecords = 15, sumD = 5)
logPossibleConfig <- function(Brecords, sumD) {
  if (sumD > 0) sum(log(Brecords:(Brecords - sumD + 1))) else 0
}

#' Log-likelihood of the linkage matrix
#'
#' @param LLL (Sparse) matrix of log-likelihood contributions for linked
#'   records.
#' @param LLA Numeric vector of log-likelihood contributions for non-linked
#'   records from A.
#' @param LLB Numeric vector of log-likelihood contributions for non-linked
#'   records from B.
#' @param links 2-column matrix of indices (A, B) for the currently linked
#'   records.
#' @param sumRowD Logical vector, one entry per record in A: does it have a
#'   link?
#' @param sumColD Logical vector, one entry per record in B: does it have a
#'   link?
#' @param gamma Numeric, proportion of linked records as a fraction of the
#'   smaller file.
#'
#' @return Numeric, the log-likelihood of the linkage matrix.
#' @export
#'
#' @examples
#' LLL <- Matrix::Matrix(0, nrow = 13, ncol = 15, sparse = TRUE)
#' LLA <- runif(13, 0, 2)
#' LLB <- runif(15, 0, 2)
#' links <- as.matrix(data.frame(idxA = c(5, 9, 11, 12, 13),
#'                               idxB = c(5, 9, 11, 13, 15)))
#' LLL[links] <- 0.67
#' sumRowD <- (seq_len(13) %in% links[, 1])
#' sumColD <- (seq_len(15) %in% links[, 2])
#' gamma <- 0.5
#' loglik(LLL, LLA, LLB, links, sumRowD, sumColD, gamma)
loglik <- function(LLL, LLA, LLB, links, sumRowD, sumColD, gamma) {
  if (length(sumColD) - nrow(links) + 1 <= 0) {
    warning(
      "File B has fewer records (", length(sumColD), ") than linked records (",
      nrow(links), "); file A should never be larger than file B.",
      call. = FALSE
    )
  }
  logPossD <- sum(log(gamma) * sumRowD + log(1 - gamma) * (1 - sumRowD)) -
    logPossibleConfig(length(sumColD), nrow(links))
  logPossD + sum(LLA[sumRowD == 0]) + sum(LLB[sumColD == 0]) + sum(LLL[links])
}

#' Simulate the true PIV values underlying the registered records
#'
#' One Gibbs step of the StEM algorithm: draws, for every record, the latent
#' true value of each PIV given the currently registered (possibly mistaken
#' or missing) value, the current linkage status, and the current parameters.
#'
#' @param data List with `encodedA`, `encodedB` (the two encoded data
#'   sources, missing values coded as `0`), `Nvalues`, and `PIVs_config`;
#'   see [stEM()].
#' @param links 2-column matrix of (A, B) indices for the currently linked
#'   records.
#' @param survivalpSameH Matrix (n links x n PIVs); `1` for stable PIVs, and the
#'   survival probability that the true value is unchanged for unstable PIVs.
#' @param sumRowD Logical vector, one entry per record in A: does it have a
#'   link?
#' @param sumColD Logical vector, one entry per record in B: does it have a
#'   link?
#' @param eta List (one per PIV) of the distribution of true values.
#' @param phi List (one per PIV) of length-4 vectors: agreement prob. in A,
#'   agreement prob. in B, missing prob. in A, missing prob. in B.
#'
#' @return List with `truepivsA` and `truepivsB`, matrices (same shape as the
#'   input data) of simulated true PIV values.
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V2 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V3 = list(dynamics = "flexible",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V4 = list(dynamics = "structured",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(0.03,0.03),
#'                                 condHazardCov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())
#'                                 )
#' )
#' Nval  <- c(6, 7, 8, 9)
#' Pmistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02))
#' Pmissing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
#' condHazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = c(0.7,0.6,0.5))
#'
#' GenData <- DataCreation(
#'   PIVs_config, Nval, NRecords = c(400, 600), Nlinks = 300,
#'   Pmistake, Pmissing, condHazard_params, enforceEstimability = TRUE
#' )
#'
#' Data4StEM <- prepare_data(GenData$dataSet1,
#'                           GenData$dataSet2,
#'                           "1",
#'                           "2",
#'                           PIVs_config,
#'                           TRUE,
#'                           "entityID",
#'                           TRUE)
#' PIVs_stable <- sapply(Data4StEM$PIVs_config, function(x)
#'                         x$dynamics != "structured")
#' FlexRL:::initDeltaMap()
#' linksR = base::matrix(0,0,2)
#' linksCpp = linksR
#' sumRowD = rep(0, nrow(Data4StEM$encodedA))
#' sumColD = rep(0, nrow(Data4StEM$encodedB))
#' nlinkrec = 0
#' survivalpSameH = base::matrix(1, nrow(linksR), length(Data4StEM$Nvalues))
#' gamma = 0.5
#' eta = lapply(Data4StEM$Nvalues, function(x) rep(1/x,x))
#' phi = lapply(Data4StEM$Nvalues, function(x)  c(0.9,0.9,0.1,0.1))
#' nCoefUnstable = lapply( seq_along(PIVs_stable), function(idx)
#'  if(PIVs_stable[idx]){ 0 }else{
#'    ncol(Data4StEM$encodedA[, Data4StEM$PIVs_config[[idx]]$condHazardCov$covA,
#'                              drop=FALSE]) +
#'    ncol(Data4StEM$encodedB[, Data4StEM$PIVs_config[[idx]]$condHazardCov$covB,
#'                              drop=FALSE]) + 1 } )
#' alpha = lapply( seq_along(PIVs_stable),
#'                 function(idx) if(PIVs_stable[idx]){ c(-Inf) }else{
#'                   rep(log(0.05), nCoefUnstable[[idx]]) })
#' newTruePivs = simulateH(data=Data4StEM, links=linksCpp,
#'                         survivalpSameH=survivalpSameH,
#'                         sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
#' truepivsA = newTruePivs$truepivsA
#' truepivsB = newTruePivs$truepivsB
simulateH <- function(data, links, survivalpSameH, sumRowD, sumColD, eta, phi) {
  PIVs <- names(data$PIVs_config)
  PIVs_stable <- sapply(data$PIVs_config, function(x) x$dynamics != "structured")
  truePIVs <- sampleH(
    nA = dim(data$encodedA[, PIVs]), nB = dim(data$encodedB[, PIVs]),
    links = links, survivalpSameH = as.matrix(survivalpSameH), pivs_stable = PIVs_stable,
    pivsA = data$encodedA[, PIVs], pivsB = data$encodedB[, PIVs], nvalues = data$Nvalues,
    nonlinkedA = sumRowD == 0, nonlinkedB = sumColD == 0, eta = eta, phi = phi
  )
  list(truepivsA = truePIVs$truepivsA, truepivsB = truePIVs$truepivsB)
}

#' Simulate the linkage matrix D (one Gibbs step)
#'
#' Given the current draw of true PIV values, resamples the linkage matrix D
#' (which pairs, if any, each record from A is linked to in B) from its full
#' conditional distribution, and returns the updated log-likelihood.
#'
#' @param data List with `encodedA`, `encodedB`, `Nvalues`, `PIVs_config`;
#'   see [stEM()].
#' @param linksR 2-column matrix (1-indexed) of the currently linked (A, B)
#'   indices.
#' @param sumRowD Logical vector, one entry per record in A: does it have a
#'   link?
#' @param sumColD Logical vector, one entry per record in B: does it have a
#'   link?
#' @param truepivsA Matrix of true PIV values, as returned by [simulateH()].
#' @param truepivsB Matrix of true PIV values, as returned by [simulateH()].
#' @param gamma Numeric, proportion of linked records as a fraction of the
#'   smaller file.
#' @param eta List (one per PIV) of the distribution of true values.
#' @param alpha List (one per PIV) of hazard coefficients
#'   (see [SurvivalUnstable()]).
#' @param phi List (one per PIV) of registration-error parameters.
#'
#' @return A list (`Dsample`, Rcpp `sampleD()` output) with: `links` updated set
#'   of links, `sumRowD` updated sumRowD, `sumColD` updated sumColD, `loglik`
#'   updated value of the complete log likelihood, `nlinkrec` updated number of
#'   linked records
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V2 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V3 = list(dynamics = "flexible",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V4 = list(dynamics = "structured",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(0.03,0.03),
#'                                 condHazardCov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())
#'                                 )
#' )
#' Nval  <- c(6, 7, 8, 9)
#' Pmistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02))
#' Pmissing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
#' condHazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = c(0.7,0.6,0.5))
#'
#' GenData <- DataCreation(
#'   PIVs_config, Nval, NRecords = c(400, 600), Nlinks = 300,
#'   Pmistake, Pmissing, condHazard_params, enforceEstimability = TRUE
#' )
#'
#' Data4StEM <- prepare_data(GenData$dataSet1,
#'                           GenData$dataSet2,
#'                           "1",
#'                           "2",
#'                           PIVs_config,
#'                           TRUE,
#'                           "entityID",
#'                           TRUE)
#' PIVs_stable <- sapply(Data4StEM$PIVs_config, function(x)
#'                         x$dynamics != "structured")
#' FlexRL:::initDeltaMap()
#' linksR = base::matrix(0,0,2)
#' linksCpp = linksR
#' sumRowD = rep(0, nrow(Data4StEM$encodedA))
#' sumColD = rep(0, nrow(Data4StEM$encodedB))
#' nlinkrec = 0
#' survivalpSameH = base::matrix(1, nrow(linksR), length(Data4StEM$Nvalues))
#' gamma = 0.5
#' eta = lapply(Data4StEM$Nvalues, function(x) rep(1/x,x))
#' phi = lapply(Data4StEM$Nvalues, function(x)  c(0.9,0.9,0.1,0.1))
#' nCoefUnstable = lapply( seq_along(PIVs_stable), function(idx)
#'  if(PIVs_stable[idx]){ 0 }else{
#'    ncol(Data4StEM$encodedA[, Data4StEM$PIVs_config[[idx]]$condHazardCov$covA,
#'                              drop=FALSE]) +
#'    ncol(Data4StEM$encodedB[, Data4StEM$PIVs_config[[idx]]$condHazardCov$covB,
#'                              drop=FALSE]) + 1 } )
#' alpha = lapply( seq_along(PIVs_stable),
#'                 function(idx) if(PIVs_stable[idx]){ c(-Inf) }else{
#'                   rep(log(0.05), nCoefUnstable[[idx]]) })
#' newTruePivs = simulateH(data=Data4StEM, links=linksCpp,
#'                         survivalpSameH=survivalpSameH,
#'                         sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
#' truepivsA = newTruePivs$truepivsA
#' truepivsB = newTruePivs$truepivsB
#' Dsample = simulateD(data=Data4StEM, linksR=linksR, sumRowD=sumRowD,
#'                    sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB,
#'                    gamma=gamma, eta=eta, alpha=alpha, phi=phi)
#' linksCpp = Dsample$links
#' linksR = linksCpp + 1
simulateD <- function(data, linksR, sumRowD, sumColD, truepivsA, truepivsB, gamma, eta, alpha, phi) {
  PIVs <- names(data$PIVs_config)
  PIVs_stable <- sapply(data$PIVs_config, function(x) x$dynamics != "structured")

  # Restrict to pairs that are compatible on the stable PIVs (candidate links)
  UA <- pasteIntoPattern(as.matrix(truepivsA[, PIVs_stable]))
  UB <- pasteIntoPattern(as.matrix(truepivsB[, PIVs_stable]))
  valuesU <- unique(c(UA, UB))
  UA <- as.numeric(factor(UA, levels = valuesU))
  UB <- as.numeric(factor(UB, levels = valuesU))
  tmpA <- indexPatterns(UA, length(valuesU))
  tmpB <- indexPatterns(UB, length(valuesU))
  select <- pairPatterns(tmpA, tmpB, length(tmpA))

  pLink <- rep(gamma, nrow(data$encodedA[, PIVs]))

  # Contribution to the log-likelihood if a record from A is NOT linked
  LLA <- rep(0, nrow(data$encodedA[, PIVs]))
  for (k in seq_along(data$Nvalues)) {
    logpTrue <- log(eta[[k]])[truepivsA[, k]]
    pMissingA <- phi[[k]][3]
    pTypoA <- (1 - pMissingA) * (1 - phi[[k]][1]) / (data$Nvalues[k] - 1)
    pAgreeA <- (1 - pMissingA) * phi[[k]][1]
    contr <- rep(pAgreeA, nrow(data$encodedA[, PIVs]))
    contr[data$encodedA[, PIVs][, k] != truepivsA[, k]] <- pTypoA
    contr[data$encodedA[, PIVs][, k] == 0] <- pMissingA
    LLA <- LLA + logpTrue + log(contr)
  }

  # Contribution to the log-likelihood if a record from B is NOT linked
  LLB <- rep(0, nrow(data$encodedB[, PIVs]))
  for (k in seq_along(data$Nvalues)) {
    logpTrue <- log(eta[[k]])[truepivsB[, k]]
    pMissingB <- phi[[k]][4]
    pTypoB <- (1 - pMissingB) * (1 - phi[[k]][2]) / (data$Nvalues[k] - 1)
    pAgreeB <- (1 - pMissingB) * phi[[k]][2]
    contr <- rep(pAgreeB, nrow(data$encodedB[, PIVs]))
    contr[data$encodedB[, PIVs][, k] != truepivsB[, k]] <- pTypoB
    contr[data$encodedB[, PIVs][, k] == 0] <- pMissingB
    LLB <- LLB + logpTrue + log(contr)
  }

  # Contribution to the log-likelihood if a candidate pair IS linked
  LLL <- Matrix::Matrix(0, nrow = nrow(data$encodedA[, PIVs]), ncol = nrow(data$encodedB[, PIVs]), sparse = TRUE)
  for (k in seq_along(data$Nvalues)) {
    HA <- truepivsA[select[, 1], k]
    HB <- truepivsB[select[, 2], k]
    logpTrue <- log(eta[[k]])[HA]
    pMissingA <- phi[[k]][3]
    pTypoA <- (1 - pMissingA) * (1 - phi[[k]][1]) / (data$Nvalues[k] - 1)
    pAgreeA <- (1 - pMissingA) * phi[[k]][1]
    pMissingB <- phi[[k]][4]
    pTypoB <- (1 - pMissingB) * (1 - phi[[k]][2]) / (data$Nvalues[k] - 1)
    pAgreeB <- (1 - pMissingB) * phi[[k]][2]
    # Contribution to the likelihood of linked observation from A
    helpA <- rep(pAgreeA, length(HA))
    helpA[data$encodedA[, PIVs][select[, 1], k] != HA] <- pTypoA
    helpA[data$encodedA[, PIVs][select[, 1], k] == 0] <- pMissingA
    # Contribution to the likelihood of linked observation from B
    helpB <- rep(pAgreeB, length(HB))
    helpB[data$encodedB[, PIVs][select[, 2], k] != HB] <- pTypoB
    helpB[data$encodedB[, PIVs][select[, 2], k] == 0] <- pMissingB

    LLL[select] <- LLL[select] + logpTrue + log(helpA) + log(helpB)

    # Add unstable part if unstable
    if (!PIVs_stable[k]) {
      times <- abs(data$encodedB[select[, 2], "date"] - data$encodedA[select[, 1], "date"])
      intercept <-  rep(1, nrow(select))
      cov_k <- cbind(
        intercept,
        data$encodedA[select[, 1], data$PIVs_config[[k]]$condHazardCov$covA, drop = FALSE],
        data$encodedB[select[, 2], data$PIVs_config[[k]]$condHazardCov$covB, drop = FALSE]
      )
      pSameH <- SurvivalUnstable(cov_k, alpha[[k]], times)
      helpH <- pSameH^(HA == HB) * ((1 - pSameH) / (data$Nvalues[k] - 1))^(HA != HB)
      LLL[select] <- LLL[select] + log(helpH)
    }
  }

  if (anyNA(LLL[select])) {
    warning("Some entries of the linkage log-likelihood matrix are NA.", call. = FALSE)
  }

  # Complete data likelihood
  LL0 <- loglik(LLL = LLL, LLA = LLA, LLB = LLB, links = linksR, sumRowD = sumRowD, sumColD = sumColD, gamma = gamma)

  Dsample <- sampleD(
    S = as.matrix(select), LLA = LLA, LLB = LLB, LLL = LLL[select],
    gamma = pLink, loglik = LL0, nlinkrec = as.integer(nrow(linksR)),
    sumRowD = sumRowD > 0, sumColD = sumColD > 0
  )
  linksR <- Dsample$links + 1

  # Sanity check: recomputed log-likelihood should match the sampler's own value
  ll_check <- loglik(LLL = LLL, LLA = LLA, LLB = LLB, links = linksR, sumRowD = Dsample$sumRowD, sumColD = Dsample$sumColD, gamma = pLink)
  if (round(Dsample$loglik, 3) != round(ll_check, 3)) {
    stop(
      "Log-likelihood sanity check failed (sampler: ", round(Dsample$loglik, 3),
      ", recomputed: ", round(ll_check, 3), ").",
      call. = FALSE
    )
  }
  Dsample
}

#' Negative log-likelihood of the exponential survival model for instability
#'
#' Used by [stEM()] to re-estimate the hazard coefficients `alpha` at every
#' iteration (via `stats::nlminb`). To use a different survival model, adapt
#' this function together with [SurvivalUnstable()].
#'
#' @param alphas Numeric vector of hazard coefficients (baseline + covariates).
#' @param X Matrix (n linked records x length(alphas)), covariates used to model
#'   the hazard.
#' @param times Numeric vector, time gap between the linked records.
#' @param Hequal Logical vector, whether the true PIV value agrees between the
#'   linked records.
#'
#' @return Numeric, the negative log-likelihood (the StEM algorithm minimises
#'   this).
#' @export
#'
#' @examples
#' alphaInit <- -0.05
#' X <- matrix(runif(5,0,2), nrow = 5, ncol = 1)
#' times <- c(0.001, 0.2, 1.3, 1.5, 2)
#' Hequal <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
#' stats::nlminb(alphaInit, loglikSurvival, X = X, times = times,
#'                 Hequal = Hequal)
loglikSurvival <- function(alphas, X, times, Hequal) {
  ll <- sum( - exp(as.matrix(X) %*% alphas) * times +
               (!Hequal) * log(exp(exp(as.matrix(X) %*% alphas) * times) - 1) )
  -ll
}

#' Exponential survival function for PIV instability
#'
#' For a linked pair `(i, j)` from A and B, models the probability that the
#' true value of an unstable PIV coincides as
#' `P(HAik = HBjk | t_ij) = exp(-exp(X . alpha) * t_ij)`. To use a different
#' survival model, adapt this function together with [loglikSurvival()].
#'
#' @param Xlinksk Matrix (n linked records x (1 + n covariates)), first column:
#'   intercept, remaining columns: covariates used to model instability.
#' @param alphask Numeric vector, hazard coefficients matching the columns
#'   of `Xlinksk`.
#' @param times Numeric vector, time gap between the linked records.
#'
#' @return Numeric vector, probability that the true PIV value is unchanged
#'   for each linked pair.
#' @export
#'
#' @examples
#' cov_k <- cbind(intercept = rep(1, 5))
#' times <- c(0.001, 0.2, 1.3, 1.5, 2)
#' SurvivalUnstable(cov_k, alphask = log(0.28), times = times)
SurvivalUnstable <- function(Xlinksk, alphask, times) {
  exp(-exp(as.matrix(Xlinksk) %*% alphask) * times)
}

#' Stochastic Expectation-Maximisation for record linkage
#'
#' Fits the FlexRL model with a Stochastic EM algorithm: each iteration runs
#' a Gibbs sampler (alternating between simulating the true PIV values and
#' the linkage matrix D) and then updates the model parameters (`gamma`,
#' `eta`, `alpha`, `phi`) from the post-burn-in Gibbs draws. See the
#' methodology paper (\doi{10.1093/jrsssc/qlaf016}) for details.
#'
#' @param data List, typically the output of [prepare_data()], with:
#'   `encodedA` (smaller source, PIVs encoded to natural numbers,
#'    `0` = missing), `encodedB` (larger source, encoded), `Nvalues`,
#'    `sameMistakes` (logical: one mistake parameter shared by A and B?),
#'    `PIVs_config` (named list, one entry per PIV, each a list with: `dynamics`
#'    (`"stable"`, `"flexible"`, or `"structured"`), `boundMistakes` (length-2
#'    numeric/NA, upper bound on the mistake probability in file 1 / file 2),
#'    `fixMistakes` (length-2 numeric/NA, mistake probability fixed to this
#'    value in file 1 / file 2), and, only for `dynamics = "structured"`,
#'    `condHazardCov` (a list with `cov1` and `cov2`, the names of covariates
#'    in file 1 / file 2 used to model the hazard of change).
#' @param StEMIter Integer, total number of StEM iterations (including burn-in).
#' @param StEMBurnin Integer, number of StEM iterations discarded as burn-in.
#' @param GibbsIter Integer, total number of Gibbs iterations per StEM step
#'   (including burn-in).
#' @param GibbsBurnin Integer, number of Gibbs iterations discarded as burn-in
#'   (`0` lets the algorithm auto-detect burn-in from the stabilisation of the
#'   linked count).
#' @param musicOn Logical; if `TRUE`, plays a short tune when the algorithm
#'   inishes.
#' @param newDirectory Path to an existing directory to save progress after
#'   each iteration, or `NULL` to disable.
#' @param saveInfoIter Logical; save the environment at the end of each
#'   iteration (only used if `newDirectory` is not `NULL`).
#' @param gamma0 Optional starting values for `gamma`; at random if `NULL`.
#' @param phiA0 Optional starting values for `phi`; at random if `NULL`.
#' @param phiB0 Optional starting values for `phi`; at random if `NULL`.
#' @param nPostSamp Integer, number of posterior draws used to estimate the
#'   final linkage probabilities `Delta`.
#'
#' @return A list with: `Delta` sparse-matrix summary (`i`, `j`, `x`) of
#'   posterior linkage probabilities; a pair is a valid link candidate once
#'   `x > 0.5` (one-to-one constraint), `gamma`, `eta`, `alpha`, `phi` the StEM
#'   chains for each parameter
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V2 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V3 = list(dynamics = "flexible",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V4 = list(dynamics = "structured",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(0.03,0.03),
#'                                 condHazardCov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())
#'                                 )
#' )
#' Nval  <- c(6, 7, 8, 9)
#' Pmistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02))
#' Pmissing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
#' condHazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = c(0.7,0.6,0.5))
#'
#' GenData <- DataCreation(
#'   PIVs_config, Nval, NRecords = c(400, 600), Nlinks = 300,
#'   Pmistake, Pmissing, condHazard_params, enforceEstimability = TRUE
#' )
#'
#' PrepData <- prepare_data(GenData$dataSet1, GenData$dataSet2, "1", "2",
#'                      PIVs_config, sameMistakes = TRUE, uniqID = "entityID")
#'
#' fit <- stEM(data = PrepData, StEMIter = 10, StEMBurnin = 5,
#'            GibbsIter = 10, GibbsBurnin = 5, musicOn = FALSE)
#' head(fit$Delta[fit$Delta$x > 0.5, ])
stEM <- function(data, StEMIter = 30, StEMBurnin = 15, GibbsIter = 20, GibbsBurnin = 10,
                 musicOn = TRUE, newDirectory = NULL, saveInfoIter = FALSE,
                 gamma0 = NULL, phiA0 = NULL, phiB0 = NULL, nPostSamp = 1000) {

  message("Running FlexRL")

  required <- c("encodedA", "encodedB", "PIVs_config", "Nvalues", "sameMistakes")
  if (!all(required %in% names(data))) {
    stop("`data` must contain: ", paste(required, collapse = ", "), ".", call. = FALSE)
  }

  nGibbsIter <- GibbsIter - GibbsBurnin
  if (StEMIter - StEMBurnin <= 0) stop("`StEMIter` must be greater than `StEMBurnin`.", call. = FALSE)
  if (nGibbsIter <= 0) stop("`GibbsIter` must be greater than `GibbsBurnin`.", call. = FALSE)

  fileA <- data$encodedA
  fileB <- data$encodedB

  PIVs <- names(data$PIVs_config)
  PIVs_stable <- sapply(data$PIVs_config, function(x) x$dynamics != "structured")
  PIVs_boundMistakes <- sapply(data$PIVs_config, function(x) x$boundMistakes)
  PIVs_fixMistakes <- sapply(data$PIVs_config, function(x) x$fixMistakes)
  PIVs_conditionalHazardVariables <- sapply(data$PIVs_config, function(x) {
    if (x$dynamics == "structured") x$condHazardCov else FALSE
  })

  N_PIVs <- length(PIVs)
  modelDynaPIVs <- which(!PIVs_stable)

  if (length(modelDynaPIVs) > 0) {
    if (!"date" %in% colnames(fileA) || !"date" %in% colnames(fileB)) {
      stop("Modelling unstable PIVs requires a `date` column in both file A and file B.", call. = FALSE)
    } else {
      tryCatch(
        { times_test = abs( mean(fileB$date,na.rm=TRUE) - mean(fileA$date,na.rm=TRUE) )
        }, error = function(msg){
          warning("Error in the `date` columns in fileA or fileB, format does not allow to compute absolute time gaps.\n")
        })
    }
    for (k in modelDynaPIVs) {
      covA <- PIVs_conditionalHazardVariables[[k]]$covA
      covB <- PIVs_conditionalHazardVariables[[k]]$covB
      if (length(covA) > 0 && !all(covA %in% colnames(fileA))) {
        stop("Some `condHazardCov$covA` variables for PIV '", PIVs[k], "' are missing from file A.", call. = FALSE)
      }
      if (length(covB) > 0 && !all(covB %in% colnames(fileB))) {
        stop("Some `condHazardCov$covB` variables for PIV '", PIVs[k], "' are missing from file B.", call. = FALSE)
      }
    }
  }

  # --- Initial parameter values -------------------------------------------
  gamma <- if (is.null(gamma0)) stats::runif(1, 0.2, 0.8) else gamma0
  eta <- lapply(data$Nvalues, function(x) rep(1 / x, x))

  nCoefUnstable <- lapply(seq_along(PIVs_stable), function(idx) {
    if (PIVs_stable[idx]) 0 else {
      ncol(fileA[, PIVs_conditionalHazardVariables[[idx]]$covA, drop = FALSE]) +
        ncol(fileB[, PIVs_conditionalHazardVariables[[idx]]$covB, drop = FALSE]) + 1
    }
  })
  alpha <- lapply(seq_along(PIVs_stable), function(idx) {
    if (PIVs_stable[idx]) -Inf else rep(log(0.05), nCoefUnstable[[idx]])
  })

  # (agreement in A, agreement in B, missing in A, missing in B)
  if (!is.null(phiA0) && !is.null(phiB0)) {
    phi <- lapply(data$Nvalues, function(x) c(phiA0, phiB0, 0.1, 0.1))
  } else if (!is.null(phiA0)) {
    phi <- lapply(data$Nvalues, function(x) c(phiA0, stats::runif(1, 0.8, 0.97), 0.1, 0.1))
  } else if (!is.null(phiB0)) {
    phi <- lapply(data$Nvalues, function(x) c(stats::runif(1, 0.8, 0.97), phiB0, 0.1, 0.1))
  } else {
    phi <- lapply(data$Nvalues, function(x) c(stats::runif(1, 0.8, 0.97), stats::runif(1, 0.8, 0.97), 0.1, 0.1))
  }

  NmissingA <- lapply(seq_along(data$Nvalues), function(k) sum(fileA[, PIVs][, k] == 0))
  NmissingB <- lapply(seq_along(data$Nvalues), function(k) sum(fileB[, PIVs][, k] == 0))

  gamma.iter <- array(NA, c(StEMIter, length(gamma)))
  eta.iter   <- lapply(data$Nvalues, function(x) array(NA, c(StEMIter, x)))
  phi.iter   <- lapply(data$Nvalues, function(x) array(NA, c(StEMIter, 4)))
  alpha.iter <- lapply(nCoefUnstable, function(x) array(NA, c(StEMIter, x)))

  boundHitCountA <- stats::setNames(rep(0L, N_PIVs), PIVs)
  boundHitCountB <- stats::setNames(rep(0L, N_PIVs), PIVs)
  boundHitThreshold <- 0.25

  pb_id <- cli::cli_progress_bar(
    format = "Running StEM algorithm {cli::pb_bar} {cli::pb_percent} | iter {cli::pb_current}/{cli::pb_total} [{cli::pb_elapsed}]",
    total = StEMIter, clear = FALSE
  )

  # Runs one Gibbs iteration and refreshes `survivalpSameH` for unstable PIVs
  gibbs_step <- function(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, gamma, eta, alpha, phi) {
    newTruePivs <- simulateH(data = data, links = linksCpp, survivalpSameH = survivalpSameH,
                             sumRowD = sumRowD, sumColD = sumColD, eta = eta, phi = phi)
    Dsample <- simulateD(data = data, linksR = linksR, sumRowD = sumRowD, sumColD = sumColD,
                         truepivsA = newTruePivs$truepivsA, truepivsB = newTruePivs$truepivsB,
                         gamma = gamma, eta = eta, alpha = alpha, phi = phi)
    linksCpp <- Dsample$links
    linksR <- linksCpp + 1
    survivalpSameH <- matrix(1, nrow(linksR), N_PIVs)
    if (length(modelDynaPIVs) > 0 && nrow(linksR) > 0) {
      times <- abs(fileB[linksR[, 2], "date"] - fileA[linksR[, 1], "date"])
      intercept <- rep(1, nrow(linksR))
      for (k in modelDynaPIVs) {
        cov_k <- cbind(intercept,
                       fileA[linksR[, 1], PIVs_conditionalHazardVariables[[k]]$covA, drop = FALSE],
                       fileB[linksR[, 2], PIVs_conditionalHazardVariables[[k]]$covB, drop = FALSE])
        survivalpSameH[, k] <- SurvivalUnstable(cov_k, alpha[[k]], times)
      }
    }
    list(linksCpp = linksCpp, linksR = linksR, sumRowD = Dsample$sumRowD, sumColD = Dsample$sumColD,
         survivalpSameH = survivalpSameH, truepivsA = newTruePivs$truepivsA, truepivsB = newTruePivs$truepivsB,
         loglik = Dsample$loglik, nlinkrec = Dsample$nlinkrec)
  }
  run_auto_burnin <- function(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, nlinkrec, gamma, eta, alpha, phi) {
    burnin_id <- cli::cli_progress_bar(
      format = "  Burn-in (auto) {cli::pb_spin} iter {Burnin_total} | linked records: {nlinkrec}",
      total = 500, current = FALSE
    )
    countBurnin <- 10
    Burnin_total <- 0
    while (countBurnin != 0) {
      Burnin_total <- Burnin_total + 1
      step <- gibbs_step(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, gamma, eta, alpha, phi)
      linksCpp <- step$linksCpp
      linksR <- step$linksR
      sumRowD <- step$sumRowD
      sumColD <- step$sumColD
      survivalpSameH <- step$survivalpSameH
      countBurnin <- if (abs(step$nlinkrec - nlinkrec) < 10) countBurnin - 1 else 15
      nlinkrec <- step$nlinkrec
      cli::cli_progress_update(id = burnin_id)
      if (Burnin_total >= 500) {
        # cli::cli_progress_done(id = burnin_id)
        stop("Auto burn-in exceeded 500 iterations without stabilising; set `GibbsBurnin` explicitly.", call. = FALSE)
      }
    }
    cli::cli_progress_done(id = burnin_id)
    list(linksCpp = linksCpp, linksR = linksR, sumRowD = sumRowD, sumColD = sumColD,
         survivalpSameH = survivalpSameH, nlinkrec = nlinkrec)
  }

  for (iter in seq_len(StEMIter)) {

    initDeltaMap()
    linksR <- matrix(0, 0, 2)
    linksCpp <- linksR
    sumRowD <- rep(0, nrow(fileA))
    sumColD <- rep(0, nrow(fileB))
    nlinkrec <- 0
    survivalpSameH <- matrix(1, nrow(linksR), N_PIVs)

    # --- Burn-in ------------------------------------------------------------
    if (GibbsBurnin == 0) {
      # Auto burn-in: keep sampling until the number of linked records stops
      # moving by more than 10 for 10 consecutive iterations, capped at 500.
      burnin <- run_auto_burnin(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, nlinkrec, gamma, eta, alpha, phi)
      linksCpp <- burnin$linksCpp
      linksR <- burnin$linksR
      sumRowD <- burnin$sumRowD
      sumColD <- burnin$sumColD
      survivalpSameH <- burnin$survivalpSameH
      nlinkrec <- burnin$nlinkrec
    } else {
      for (j in seq_len(GibbsBurnin)) {
        step <- gibbs_step(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, gamma, eta, alpha, phi)
        linksCpp <- step$linksCpp
        linksR <- step$linksR
        sumRowD <- step$sumRowD
        sumColD <- step$sumColD
        survivalpSameH <- step$survivalpSameH
      }
    }

    # --- Post burn-in: accumulate sufficient statistics for the M-step ------
    Vgamma <- c()
    Veta <- lapply(data$Nvalues, function(x) c())
    Valpha <- mapply(createDataAlpha, nCoefUnstable = nCoefUnstable, stable = PIVs_stable, SIMPLIFY = FALSE)
    Vphi <- lapply(data$Nvalues, function(x) c())

    for (j in seq_len(nGibbsIter)) {
      step <- gibbs_step(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, gamma, eta, alpha, phi)
      linksCpp <- step$linksCpp
      linksR <- step$linksR
      sumRowD <- step$sumRowD
      sumColD <- step$sumColD
      survivalpSameH <- step$survivalpSameH
      truepivsA <- step$truepivsA
      truepivsB <- step$truepivsB

      if (length(modelDynaPIVs) > 0 && nrow(linksR) > 0) {
        times <- abs(fileB[linksR[, 2], "date"] - fileA[linksR[, 1], "date"])
        intercept <- rep(1, nrow(linksR))
        for (k in modelDynaPIVs) {
          cov_k <- cbind(intercept,
                         fileA[linksR[, 1], PIVs_conditionalHazardVariables[[k]]$covA, drop = FALSE],
                         fileB[linksR[, 2], PIVs_conditionalHazardVariables[[k]]$covB, drop = FALSE])
          Hequal <- truepivsA[linksR[, 1], k] == truepivsB[linksR[, 2], k]
          Valpha[[k]] <- rbind(Valpha[[k]], cbind(cov_k, times = times, Hequal = Hequal))
        }
      }

      # Update gamma
      Vgamma[j] <- nrow(linksR)
      if (nrow(linksR) == 0) {
        warning("No link was made in a Gibbs iteration (log-lik = ", step$loglik,
                ", gamma so far = ", gamma, "); check for a support mismatch between A and B.", call. = FALSE)
      }

      # Update eta and phi
      for (k in seq_len(N_PIVs)) {
        facpivsA <- factor(truepivsA[, k], levels = 1:data$Nvalues[k])
        facpivsB <- factor(truepivsB[, k], levels = 1:data$Nvalues[k])
        Veta[[k]] <- rbind(Veta[[k]], table(facpivsA[sumRowD == 0]) + table(facpivsB[sumColD == 0]) + table(facpivsA[sumRowD == 1]))
        Vphi[[k]] <- rbind(Vphi[[k]], c(
          sum(truepivsA[, k] == fileA[, PIVs][, k]),
          sum(truepivsB[, k] == fileB[, PIVs][, k]),
          NmissingA[[k]], NmissingB[[k]]
        ))
      }
    }

    # --- M-step -------------------------------------------------------------
    gamma <- sum(Vgamma) / (nGibbsIter * nrow(truepivsA))
    if (gamma == 1) { gamma <- 0.99; warning("`gamma` hit 1 and was reset to 0.99.", call. = FALSE) }
    if (gamma == 0) { gamma <- 0.01; warning("`gamma` hit 0 and was reset to 0.01.", call. = FALSE) }

    for (k in seq_len(N_PIVs)) {
      eta[[k]] <- colSums(Veta[[k]]) / sum(Veta[[k]])
      if (any(eta[[k]] == 0)) {
        warning("PIV '", PIVs[k], "' has one or more values with zero estimated frequency; check for rare categories.", call. = FALSE)
      }
    }

    if (length(modelDynaPIVs) > 0 && nrow(linksR) > 0) {
      for (k in modelDynaPIVs) {
        X <- Valpha[[k]][, 1:nCoefUnstable[[k]]]
        opt <- stats::nlminb(rep(-0.05, nCoefUnstable[[k]]), loglikSurvival,
                             X = X, times = Valpha[[k]]$times, Hequal = Valpha[[k]]$Hequal)
        alpha[[k]] <- opt$par
      }
    }

    for (k in seq_len(N_PIVs)) {
      NtotalA <- nGibbsIter * nrow(fileA)
      NtotalB <- nGibbsIter * nrow(fileB)
      Ntotal  <- NtotalA + NtotalB

      if (data$sameMistakes) {
        phi[[k]][1] <- phi[[k]][2] <-
          (sum(Vphi[[k]][, 1]) + sum(Vphi[[k]][, 2])) / (Ntotal - sum(Vphi[[k]][, 3]) - sum(Vphi[[k]][, 4]))
      } else {
        phi[[k]][1] <- sum(Vphi[[k]][, 1]) / (NtotalA - sum(Vphi[[k]][, 3]))
        phi[[k]][2] <- sum(Vphi[[k]][, 2]) / (NtotalB - sum(Vphi[[k]][, 4]))
      }

      if (!is.na(PIVs_fixMistakes[k][1])) phi[[k]][1] <- 1 - PIVs_fixMistakes[k][1]
      if (!is.na(PIVs_fixMistakes[k][2])) phi[[k]][2] <- 1 - PIVs_fixMistakes[k][2]

      if (!is.na(PIVs_boundMistakes[k][1])) {
        if (phi[[k]][1] < 1 - PIVs_boundMistakes[k][1]) {
          phi[[k]][1] <- 1 - PIVs_boundMistakes[k][1]
          boundHitCountA[k] <- boundHitCountA[k] + 1L
        }
      }
      if (!is.na(PIVs_boundMistakes[k][2])) {
        if (phi[[k]][2] < 1 - PIVs_boundMistakes[k][2]) {
          phi[[k]][2] <- 1 - PIVs_boundMistakes[k][2]
          boundHitCountB[k] <- boundHitCountB[k] + 1L
        }
      }

      phi[[k]][3] <- sum(Vphi[[k]][, 3]) / NtotalA
      phi[[k]][4] <- sum(Vphi[[k]][, 4]) / NtotalB
    }

    gamma.iter[iter, ] <- gamma
    for (k in seq_len(N_PIVs)) {
      eta.iter[[k]][iter, ] <- eta[[k]]
      alpha.iter[[k]][iter, ] <- alpha[[k]]
      phi.iter[[k]][iter, ] <- phi[[k]]
    }

    cli::cli_progress_update(id = pb_id)

    if (!is.null(newDirectory) && saveInfoIter) {
      save.image(file=file.path(newDirectory, "myEnvironment.RData"))
      save(iter, gamma, eta, alpha, phi, gamma.iter, eta.iter, phi.iter, alpha.iter,
           file = file.path(newDirectory, "myEnvironmentLocal.RData"))
    }
  }

  freqA <- boundHitCountA / StEMIter
  freqB <- boundHitCountB / StEMIter
  warning_msg_bound_phi <- c(
    sprintf("PIV '%s', file A hit its mistake bound on %.0f%% of StEM iterations.", PIVs[freqA > boundHitThreshold], 100 * freqA[freqA > boundHitThreshold]),
    sprintf("PIV '%s', file B hit its mistake bound on %.0f%% of StEM iterations.", PIVs[freqB > boundHitThreshold], 100 * freqB[freqB > boundHitThreshold])
  )
  if (length(warning_msg_bound_phi) > 0) {
    warning(paste(warning_msg_bound_phi, collapse = "\n  "), call. = FALSE)
  }

  cli::cli_progress_done(id = pb_id)

  # --- Posterior draws of the linkage matrix --------------------------------
  gamma_avg <- apply(gamma.iter, 2, function(x) mean(x[StEMBurnin:StEMIter], drop = FALSE))
  eta_avg   <- lapply(eta.iter,     function(x) apply(x[StEMBurnin:StEMIter, , drop = FALSE], 2, mean))
  alpha_avg <- lapply(alpha.iter,   function(x) apply(x[StEMBurnin:StEMIter, , drop = FALSE], 2, mean))
  phi_avg   <- lapply(phi.iter,     function(x) apply(x[StEMBurnin:StEMIter, , drop = FALSE], 2, mean))

  Delta <- Matrix::Matrix(0, nrow = nrow(fileA), ncol = nrow(fileB), sparse = TRUE)
  pbfinal_id <- cli::cli_progress_bar(
    format = "Drawing Delta          {cli::pb_bar} {cli::pb_percent} [{cli::pb_elapsed}]",
    total = nPostSamp, clear = FALSE
  )
  for (m in seq_len(nPostSamp)) {
    step <- gibbs_step(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, gamma_avg, eta_avg, alpha_avg, phi_avg)
    # Note: eta/alpha/phi are frozen at their posterior mean for these draws
    linksCpp <- step$linksCpp
    linksR <- step$linksR
    sumRowD <- step$sumRowD
    sumColD <- step$sumColD
    survivalpSameH <- step$survivalpSameH

    if (nrow(linksR) > 0) {
      for (l in seq_len(nrow(linksR))) Delta[linksR[l, 1], linksR[l, 2]] <- Delta[linksR[l, 1], linksR[l, 2]] + 1
    }
    cli::cli_progress_update(id = pbfinal_id)
  }
  Delta <- Delta / nPostSamp
  cli::cli_progress_done(id = pbfinal_id)

  if (musicOn) utils::browseURL("https://www.youtube.com/watch?v=NTa6Xbzfq1U")

  list(
    Delta = as.data.frame(Matrix::summary(Delta)),
    gamma = gamma.iter, eta = eta.iter, alpha = alpha.iter, phi = phi.iter
  )
}

#' Naive (deterministic exact-match) record linkage
#'
#' Links records that agree exactly on every non-missing PIV. Does not
#' enforce the one-to-one assignment constraint, so should be used only to
#' gauge the difficulty of the linkage task (amount of duplication, and
#' discriminative power of the PIVs together).
#'
#' @param PIVs Character vector, names of the PIVs (columns present in both
#'   files).
#' @param encodedA The data source (PIVs encoded to natural numbers).
#' @param encodedB The data source (PIVs encoded to natural numbers).
#' @param na.match Logical; if `TRUE`, a missing PIV value is treated as
#'   matching any value in the other file (default `TRUE`).
#' @param na.is.zero Logical; if `TRUE`, missing values are already coded as
#'   `0` (FlexRL's convention); if `FALSE`, `NA` will be recoded to `0`
#'   (default `TRUE`).
#'
#' @return Data frame with columns `idxA`, `idxB`: pairs of record indices
#'   that match.
#' @export
#'
#' @examples
#' PIVs_config <- list(V1 = list(dynamics = "stable",
#'                               boundMistakes = c(0.1, 0.1),
#'                               fixMistakes = c(NA, NA)),
#'                     V2 = list(dynamics = "stable",
#'                               boundMistakes = c(0.1, 0.1),
#'                               fixMistakes = c(NA, NA)))
#' GenData <- DataCreation(PIVs_config, Nval = c(6,6), NRecords = c(30, 50),
#'                Nlinks = 15, Pmistake = list(V1 = c(0, 0), V2 = c(0, 0)),
#'                Pmissing = list(V1 = c(0, 0), V1 = c(0, 0)),
#'                condHazard_params = list(V1 = c()),
#'                enforceEstimability = FALSE)
#' PIVs <- names(PIVs_config)
#' encodedA <- GenData$dataSet1
#' encodedA[,PIVs][ is.na(encodedA[,PIVs]) ] = 0
#' encodedB <- GenData$dataSet2
#' encodedB[,PIVs][ is.na(encodedB[,PIVs]) ] = 0
#' NaiveLinkage(PIVs, encodedA, encodedB)
NaiveLinkage <- function(PIVs, encodedA, encodedB, na.match = TRUE, na.is.zero = TRUE) {

  if (!na.is.zero) {
    encodedA[PIVs][is.na(encodedA[PIVs])] <- 0
    encodedB[PIVs][is.na(encodedB[PIVs])] <- 0
  }
  if (!isTRUE(na.match) && !isFALSE(na.match)) {
    stop("`na.match` must be TRUE or FALSE.", call. = FALSE)
  }

  rownames(encodedA) <- seq_len(nrow(encodedA))
  rownames(encodedB) <- seq_len(nrow(encodedB))
  DeltaNaiveLinked <- data.frame(idxA = character(0), idxB = character(0))

  match_block <- function(idA, idB) {
    valuesU <- unique(c(idA$U, idB$U))
    a <- as.numeric(factor(idA$U, levels = valuesU))
    b <- as.numeric(factor(idB$U, levels = valuesU))
    tmpA <- indexPatterns(a, length(valuesU))
    tmpB <- indexPatterns(b, length(valuesU))
    select <- pairPatterns(tmpA, tmpB, length(tmpA))
    if (nrow(select) == 0) return(NULL)
    data.frame(idxA = idA$ID[as.integer(select[, 1])], idxB = idB$ID[as.integer(select[, 2])])
  }

  # A (non-missing) vs. B (exact match on all PIVs), and vice versa
  isNotMissingA <- apply(encodedA[, PIVs] != 0, 1, all)
  isNotMissingB <- apply(encodedB[, PIVs] != 0, 1, all)
  DeltaNaiveLinked <- rbind(
    DeltaNaiveLinked,
    match_block(
      list(U = pasteIntoPattern(as.matrix(encodedA[isNotMissingA, PIVs])), ID = rownames(encodedA)[isNotMissingA]),
      list(U = pasteIntoPattern(as.matrix(encodedB[, PIVs])), ID = rownames(encodedB))
    ),
    match_block(
      list(U = pasteIntoPattern(as.matrix(encodedA[, PIVs])), ID = rownames(encodedA)),
      list(U = pasteIntoPattern(as.matrix(encodedB[isNotMissingB, PIVs])), ID = rownames(encodedB)[isNotMissingB])
    )
  )

  # Records with exactly one missing PIV: match on the remaining PIVs
  # (only run when na.match = TRUE; missing values otherwise never count
  # as matches), (there is no match when several values are missing)
  if (isTRUE(na.match)) {
    for (k in seq_along(PIVs)) {
      isMissingA_k <- encodedA[, PIVs][, k] == 0
      isMissingB_k <- encodedB[, PIVs][, k] == 0
      DeltaNaiveLinked <- rbind(
        DeltaNaiveLinked,
        match_block(
          list(U = pasteIntoPattern(as.matrix(encodedA[isMissingA_k, PIVs][, -k])), ID = rownames(encodedA)[isMissingA_k]),
          list(U = pasteIntoPattern(as.matrix(encodedB[, PIVs][, -k])), ID = rownames(encodedB))
        ),
        match_block(
          list(U = pasteIntoPattern(as.matrix(encodedA[, PIVs][, -k])), ID = rownames(encodedA)),
          list(U = pasteIntoPattern(as.matrix(encodedB[isMissingB_k, PIVs][, -k])), ID = rownames(encodedB)[isMissingB_k])
        )
      )
    }
  }

  DeltaNaiveLinked[!duplicated(DeltaNaiveLinked), ]
}

# ============================================================================
# FDP ("false discovery proportion") estimation via synthetic-record
# augmentation, for a range of external record-linkage packages.
#
# Every FDPinRL_* wrapper below has the same contract: take the two data
# sets (`NewA`, `NewB`, as produced by [synthesise()]) and an `arguments`
# list of extra parameters forwarded to the underlying package, run that
# package's linkage pipeline, and return
#   list(result_index_linked_A, result_index_linked_B, result_vec_linked_scores)
# so that [compute_FDP_RLwithSynth()] can treat every method uniformly.
# `result_vec_linked_scores` may be `NULL` for methods that don't return a
# continuous linkage score (in which case only the default threshold is used).
# Details about this FDP estimation method in https://doi.org/10.1002/sim.70292.
# ============================================================================

#' FDP estimation with FlexRL
#'
#' More details about FlexRL on
#' https://cran.r-project.org/web/packages/FlexRL/index.html
#'
#' @param NewA Data set as returned by [synthesise()].
#' @param NewB Data set as returned by [synthesise()].
#' @param arguments List of extra arguments forwarded to [stEM()] (e.g.
#'   `data`, `StEMIter`, `StEMBurnin`, `GibbsIter`, `GibbsBurnin`).
#'
#' @return List with `result_index_linked_A`, `result_index_linked_B`,
#'   `result_vec_linked_scores` (posterior linkage probabilities).
#' @export
#'
FDPinRL_FlexRL <- function(NewA, NewB, arguments) {
  PIVs <- names(arguments$data$PIVs_config)
  NewA[PIVs][is.na(NewA[PIVs])] <- 0 # FlexRL sentinel for missing values
  NewB[PIVs][is.na(NewB[PIVs])] <- 0
  arguments$data[["A"]] <- NewA
  arguments$data[["B"]] <- NewB
  fit <- do.call(stEM, arguments)
  list(fit$Delta$i, fit$Delta$j, fit$Delta$x)
}

#' FDP estimation with the fedmatch package
#'
#' More details about fedmatch on
#' https://cran.r-project.org/web/packages/fedmatch/index.html
#'
#' @param NewA Data set as returned by [synthesise()].
#' @param NewB Data sets as returned by [synthesise()].
#' @param arguments List of extra arguments forwarded to
#'   `fedmatch::merge_plus()` (e.g. `by`, `match_type`, `unique_key_1`,
#'   `unique_key_2`, `multivar_settings`).
#'
#' @return List with `result_index_linked_A`, `result_index_linked_B`,
#'   `result_vec_linked_scores`.
#' @export
#'
FDPinRL_fedmatch <- function(NewA, NewB, arguments) {
  NewA[, arguments$unique_key_1] <- rownames(NewA)
  NewB[, arguments$unique_key_2] <- rownames(NewB)
  arguments$data1 <- NewA
  arguments$data2 <- NewB
  results <- do.call(fedmatch::merge_plus, arguments)
  list(results$matches$unique_key_A, results$matches$unique_key_B, results$matches$multivar_score)
}

#' FDP estimation with the reclin2 package
#'
#' More details about reclin2 on
#' https://cran.r-project.org/web/packages/reclin2/index.html
#'
#' @param NewA Data set as returned by [synthesise()].
#' @param NewB Data sets as returned by [synthesise()].
#' @param arguments List of extra arguments forwarded across the `reclin2`
#'   pipeline (`pair()`, `compare_pairs()`, `problink_em()`, `predict()`,
#'   `select_threshold()`): e.g. `on`, `formula`, `type`, `add`, `variable`,
#'   `score`, `threshold`.
#'
#' @return List with `result_index_linked_A`, `result_index_linked_B`, `result_vec_linked_scores`.
#' @export
#'
FDPinRL_reclin2 <- function(NewA, NewB, arguments) {
  arguments$x <- NewA
  arguments$y <- NewB
  pairs <- do.call(reclin2::pair, arguments[intersect(names(arguments), names(formals(reclin2::pair)))])

  arguments$pairs <- pairs
  pairs <- do.call(reclin2::compare_pairs, arguments[intersect(names(arguments), names(formals(reclin2::compare_pairs)))])

  arguments$data <- pairs
  model <- do.call(reclin2::problink_em, arguments[intersect(names(arguments), names(formals(reclin2::problink_em)))])

  arguments$object <- model
  arguments$pairs <- pairs
  predict_formals <- union(names(formals(stats::predict)), c("object", "pairs", "newdata", "type", "binary", "add", "comparators", "inplace", "new_name"))
  pairs <- do.call(stats::predict, arguments[intersect(names(arguments), predict_formals)])

  arguments$pairs <- pairs
  selected <- do.call(reclin2::select_threshold, arguments[intersect(names(arguments), names(formals(reclin2::select_threshold)))])
  selected <- as.data.frame(selected)
  selected <- selected[selected[, arguments$variable] == TRUE, ]

  list(selected$.x, selected$.y, selected[, arguments$type])
}

#' FDP estimation with the BRL package
#'
#' More details about BRL on
#' https://cran.r-project.org/web/packages/BRL/index.html
#'
#' @param NewA Data set as returned by [synthesise()].
#' @param NewB Data sets as returned by [synthesise()].
#' @param arguments List of extra arguments forwarded to `BRL::compareRecords()`
#'   and `BRL::bipartiteGibbs()` (e.g. `flds`, `types`, `nIter`).
#'
#' @return List with `result_index_linked_A`, `result_index_linked_B`, `result_vec_linked_scores`.
#' @export
#'
FDPinRL_BRL <- function(NewA, NewB, arguments) {
  arguments$df1 <- NewB
  arguments$df2 <- NewA
  myCompData <- do.call(BRL::compareRecords, arguments[intersect(names(arguments), names(formals(BRL::compareRecords)))])

  arguments$cd <- myCompData
  chain <- do.call(BRL::bipartiteGibbs, arguments[intersect(names(arguments), names(formals(BRL::bipartiteGibbs)))])

  if (!"nIter" %in% names(arguments)) arguments$nIter <- 1000
  n1 <- nrow(NewB)
  n2 <- nrow(NewA)
  Zchain <- chain$Z[, setdiff(seq_len(arguments$nIter), seq_len(0.10 * arguments$nIter)), drop = FALSE]
  Zchain[Zchain > n1 + 1] <- n1 + 1

  tableLabels <- apply(Zchain, 1, tabulate, nbins = n1 + 1) / ncol(Zchain)
  probNoLink <- tableLabels[n1 + 1, ]
  maxProbOption <- apply(tableLabels, 2, which.max)
  probMaxProbOption <- apply(tableLabels, 2, max)
  maxProbOptionIsLink <- maxProbOption <= n1

  # Bayes-optimal decision rule under a symmetric loss
  # (see BRL package documentation)
  lFM1 <- 1
  lFM2 <- 2
  lFNM <- 1
  tholdLink <- lFM1 / (lFM1 + lFNM) + (lFM2 - lFM1 - lFNM) * (1 - probNoLink - probMaxProbOption) / (lFM1 + lFNM)
  isLink <- maxProbOptionIsLink & (probMaxProbOption > tholdLink)
  Zhat <- (n1+1):(n1+n2)
  Zhat[isLink] <- maxProbOption[isLink]

  list(which(Zhat <= n1), Zhat[Zhat <= n1], probMaxProbOption[isLink])
}

#' FDP estimation with the fastLink package
#'
#' More details about fastLink on
#' https://cran.r-project.org/web/packages/fastLink/index.html
#'
#' @param NewA Data set as returned by [synthesise()].
#' @param NewB Data sets as returned by [synthesise()].
#' @param arguments List of extra arguments forwarded to `fastLink::fastLink()`
#'   (e.g. `varnames`, `threshold.match`, `tol.em`, `return.all`).
#'
#' @return List with `result_index_linked_A`, `result_index_linked_B`,
#'   `result_vec_linked_scores`.
#' @export
#'
FDPinRL_fastLink <- function(NewA, NewB, arguments) {
  arguments$dfA <- NewA
  arguments$dfB <- NewB
  out <- do.call(fastLink::fastLink, arguments)
  list(out$matches$inds.a, out$matches$inds.b, out$posterior)
}

#' FDP estimation with the multilink package
#'
#' More details about multilink on
#' https://cran.r-project.org/web/packages/multilink/index.html
#'
#' @param NewA Data set as returned by [synthesise()].
#' @param NewB Data sets as returned by [synthesise()].
#' @param arguments List of extra arguments forwarded across the `multilink`
#'   pipeline (`create_comparison_data()`, `reduce_comparison_data()`,
#'   `specify_prior()`, `gibbs_sampler()`, `find_bayes_estimate()`,
#'   `relabel_bayes_estimate()`): e.g. `records`, `types`, `breaks`, `duplicates`, `n_iter`.
#'
#' @return List with `result_index_linked_A`, `result_index_linked_B`,
#'   `result_vec_linked_scores` (`NULL`: multilink does not return a
#'   per-pair score).
#' @export
#'
FDPinRL_multilink <- function(NewA, NewB, arguments) {
  PIVs <- names(arguments$records)
  arguments$records <- as.data.frame(lapply(rbind(NewA[, PIVs], NewB[, PIVs]), as.character))
  arguments$file_sizes <- c(nrow(NewA), nrow(NewB))
  comparison_list <- do.call(multilink::create_comparison_data, arguments[intersect(names(arguments), names(formals(multilink::create_comparison_data)))])

  arguments$comparison_list <- comparison_list
  arguments$pairs_to_keep <- rep(TRUE, nrow(NewA) * nrow(NewB))
  reduced <- do.call(multilink::reduce_comparison_data, arguments[intersect(names(arguments), names(formals(multilink::reduce_comparison_data)))])

  if (!"dup_count_prior_family" %in% names(arguments)) arguments$dup_count_prior_family <- rep(NA, reduced$K)
  if (!"dup_count_prior_pars" %in% names(arguments)) arguments$dup_count_prior_pars <- rep(NA, reduced$K)
  prior <- do.call(multilink::specify_prior, arguments[intersect(names(arguments), names(formals(multilink::specify_prior)))])

  arguments$comparison_list <- reduced
  arguments$prior_list <- prior
  results <- do.call(multilink::gibbs_sampler, arguments[intersect(names(arguments), names(formals(multilink::gibbs_sampler)))])

  arguments$burn_in <- if ("n_iter" %in% names(arguments)) as.integer(arguments$n_iter / 5) else 400
  arguments$partitions <- results$partitions
  estimate <- do.call(multilink::find_bayes_estimate, arguments[intersect(names(arguments), names(formals(multilink::find_bayes_estimate)))])

  arguments$reduced_comparison_list <- reduced
  arguments$bayes_estimate <- estimate
  relabel <- do.call(multilink::relabel_bayes_estimate, arguments[intersect(names(arguments), names(formals(multilink::relabel_bayes_estimate)))])

  linked_pairs <- data.frame(idxA = integer(0), idxB = integer(0))
  for (i in seq_along(relabel$link_id)) {
    l <- relabel$link_id[i]
    if (l > 0 && sum(relabel$link_id == l) > 1) {
      idx <- which(relabel$link_id == l)
      linked_pairs <- rbind(linked_pairs, data.frame(idxA = idx[1], idxB = idx[2] - nrow(NewA)))
    }
  }
  linked_pairs <- unique(linked_pairs)
  list(linked_pairs$idxA, linked_pairs$idxB, NULL)
}

#' FDP estimation with the diyar package
#'
#' More details about diyar on
#' https://cran.r-project.org/web/packages/diyar/index.html
#'
#' @param NewA Data set as returned by [synthesise()].
#' @param NewB Data sets as returned by [synthesise()].
#' @param arguments List of extra arguments forwarded to
#'   `diyar::prob_score_range()` and `diyar::links_wf_probabilistic()` (e.g.
#'   `attribute`, `probabilistic`, `return_weights`).
#'
#' @return List with `result_index_linked_A`, `result_index_linked_B`,
#'   `result_vec_linked_scores` (`NULL`: not exposed per-pair by this wrapper).
#' @export
#'
FDPinRL_diyar <- function(NewA, NewB, arguments) {
  allrecords <- rbind(NewA, NewB)
  PIVs <- names(arguments$attribute)
  arguments$attribute <- as.list(allrecords[, PIVs])
  linkscores <- do.call(diyar::prob_score_range, arguments[intersect(names(arguments), names(formals(diyar::prob_score_range)))])

  arguments$score_threshold <- linkscores$mid_scorce
  arguments$data_source <- allrecords$source
  idlinkage <- do.call(diyar::links_wf_probabilistic, arguments[intersect(names(arguments), names(formals(diyar::links_wf_probabilistic)))])

  linked_pairs <- data.frame(idxA = integer(0), idxB = integer(0))
  for (clusterid in unique(idlinkage$pid)) {
    idx <- which(idlinkage$pid == clusterid)
    if (length(idx) == 2 && idx[1] <= nrow(NewA) && idx[2] > nrow(NewA)) {
      linked_pairs <- rbind(linked_pairs, data.frame(idxA = idx[1], idxB = idx[2] - nrow(NewA)))
    }
  }
  linked_pairs <- unique(linked_pairs)
  list(linked_pairs$idxA, linked_pairs$idxB, NULL)
}

#' Agreement rate of linked pairs on shared variables
#'
#' For a set of linked pairs (and, optionally, the true pairs), computes how
#' often the two records agree on each variable in `common_vars`.
#'
#' @param data1 Data frame containing `common_vars`.
#' @param data2 Data frame containing `common_vars`.
#' @param common_vars Character vector, names of the variables to compare
#'   (must exist in both `data1` and `data2`).
#' @param linked_pairs Data frame/matrix/list with 2 columns of indices
#'   (into `data1`, `data2`) for the pairs to evaluate.
#' @param known_truth Logical; if `TRUE` and `true_pairs` is supplied, also
#'   compute the agreement rate for the true pairs.
#' @param true_pairs Data frame/matrix/list with 2 columns of indices for the
#'   true matches; only used if `known_truth = TRUE`.
#' @param na.rm Logical; if `TRUE` (default), pairs with a missing value on a
#'   variable are excluded from that variable's agreement rate.
#' @param na.match Logical, required if `na.rm = FALSE`: should a missing value
#'   be treated as agreeing (`TRUE`) or disagreeing (`FALSE`) with any value?
#'
#' @return List with `linked_agreements` (named numeric vector, one entry per
#'   variable in `common_vars`) and, if available, `true_agreements`.
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V2 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V3 = list(dynamics = "flexible",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V4 = list(dynamics = "structured",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(0.03,0.03),
#'                                 condHazardCov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())
#'                                 )
#' )
#' Nval  <- c(6, 7, 8, 9)
#' Pmistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02))
#' Pmissing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
#' condHazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = c(0.7,0.6,0.5))
#'
#' GenData <- DataCreation(
#'   PIVs_config, Nval, NRecords = c(400, 600), Nlinks = 300,
#'   Pmistake, Pmissing, condHazard_params, enforceEstimability = TRUE
#' )
#'
#' PrepData <- prepare_data(GenData$dataSet1, GenData$dataSet2, "1", "2",
#'                      PIVs_config, sameMistakes = TRUE, uniqID = "entityID")
#'
#' fit <- stEM(data = PrepData, StEMIter = 10, StEMBurnin = 5,
#'            GibbsIter = 10, GibbsBurnin = 5, musicOn = FALSE)
#' linked_pairs <- fit$Delta[fit$Delta$x > 0.5, ]
#'
#' rl_agreement(PrepData$encodedA, PrepData$encodedB,
#'                names(PIVs_config), linked_pairs)
#' rl_agreement(PrepData$encodedA, PrepData$encodedB,
#'                names(PIVs_config), linked_pairs, TRUE, GenData$true_pairs)
rl_agreement <- function(data1, data2, common_vars, linked_pairs,
                         known_truth = FALSE, true_pairs = NULL, na.rm = TRUE, na.match = NULL) {

  stopifnot(is.data.frame(data1), is.data.frame(data2))
  if (!is.character(common_vars) || length(common_vars) == 0) {
    stop("`common_vars` must be a non-empty character vector.", call. = FALSE)
  }
  missing1 <- setdiff(common_vars, names(data1))
  missing2 <- setdiff(common_vars, names(data2))
  if (length(missing1) || length(missing2)) {
    stop("`common_vars` missing from data1: [", paste(missing1, collapse = ", "),
         "]; from data2: [", paste(missing2, collapse = ", "), "]", call. = FALSE)
  }
  linked_pairs <- as.data.frame(linked_pairs)
  if (ncol(linked_pairs) < 2) stop("`linked_pairs` must have two columns.", call. = FALSE)

  out <- list(linked_agreements = .rl_agreement_rates(data1, data2, common_vars, linked_pairs, na.rm, na.match))
  if (isTRUE(known_truth) && !is.null(true_pairs)) {
    true_pairs <- as.data.frame(true_pairs)
    out$true_agreements <- .rl_agreement_rates(data1, data2, common_vars, true_pairs, na.rm, na.match)
  } else if (isTRUE(known_truth)) {
    message("`known_truth = TRUE` but `true_pairs` is NULL: returning `linked_agreements` only.")
  }
  out
}

# Proportion of agreement, one value per variable, for the pairs given by
# (idx1[k], idx2[k]). Compared as character so that factor level sets differing
# across sources don't produce spurious mismatches.
.rl_agreement_rates <- function(data1, data2, common_vars, pairs, na.rm = TRUE, na.match = NULL) {
  idx1 <- pairs[[1]]
  idx2 <- pairs[[2]]
  vapply(common_vars, function(v) {
    is_match <- data1[idx1, v] == data2[idx2, v]
    if (isTRUE(na.rm)) {
      is_match <- is_match[!is.na(is_match)]
    } else {
      is_match[is.na(is_match)] <- isTRUE(na.match)
    }
    if (length(is_match) == 0) NA_real_ else mean(is_match)
  }, numeric(1))
}

#' Prepare two data sources for [stEM()]
#'
#' Wraps the data-preparation steps needed before calling [stEM()]: tags each
#' source with a `source` column, relabels the larger source as `B` (FlexRL
#' requires B to be the larger file — adjusting `condHazardCov`/`boundMistakes`/
#' `fixMistakes` accordingly if A and B are swapped), drops records whose PIV
#' values fall outside the support shared by both files, warns if two PIVs
#' are strongly associated (Cramer's V > 0.3, which may degrade record linkage
#' performance), encodes every PIV to natural numbers using levels pooled
#' across both sources, and recodes missing values to `0`.
#'
#' @param data1 Data frame, the raw data source (whichever has more rows
#' becomes `B`).
#' @param data2 Data frame, the raw data source (whichever has more rows
#' becomes `B`).
#' @param label1 Character, label recorded in the `source` column
#'   for `data1`.
#' @param label2 Character, label recorded in the `source` column
#'   for `data2`.
#' @param PIVs_config Named list describing each PIV — see [DataCreation()].
#' @param sameMistakes Logical, will A and B share one mistake-probability
#'   parameter per PIV.
#' @param uniqID Optional column name (present in both files) with the true
#'   entity identifier, used to build `true_pairs` for evaluation; `NULL` if
#'   unavailable.
#' @param restrict_support_intersection Logical; if `TRUE` (default), records
#'   with an out-of-common-support PIV value are dropped (and a warning issued);
#'   if `FALSE`, only the warning is issued.
#'
#' @return A list ready to use as the `data` argument of [stEM()]: `encodedA`,
#'   `encodedB`, `Nvalues`, `PIVs_config`, `sameMistakes`, and `true_pairs`
#'   (`NULL` if `uniqID` was not supplied).
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V2 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V3 = list(dynamics = "flexible",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V4 = list(dynamics = "structured",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(0.03,0.03),
#'                                 condHazardCov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())
#'                                 )
#' )
#' Nval  <- c(6, 7, 8, 9)
#' Pmistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02))
#' Pmissing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
#' condHazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = c(0.7,0.6,0.5))
#'
#' GenData <- DataCreation(
#'   PIVs_config, Nval, NRecords = c(400, 600), Nlinks = 300,
#'   Pmistake, Pmissing, condHazard_params, enforceEstimability = TRUE
#' )
#'
#' PrepData <- prepare_data(GenData$dataSet1, GenData$dataSet2, "1", "2",
#'                      PIVs_config, sameMistakes = TRUE, uniqID = "entityID")
#' str(PrepData, max.level = 1)
prepare_data <- function(data1, data2, label1, label2, PIVs_config,
                         sameMistakes = TRUE, uniqID = NULL,
                         restrict_support_intersection = TRUE) {

  rownames(data1) <- seq_len(nrow(data1))
  rownames(data2) <- seq_len(nrow(data2))

  # --- Ground truth (optional) --------------------------------------------
  true_Delta <- NULL
  if (!is.null(uniqID)) {
    true_links <- intersect(data1[[uniqID]], data2[[uniqID]])
    true_Delta <- data.frame(matrix(0, nrow = 0, ncol = 2))
    for (id in true_links) {
      id1 <- which(data1[[uniqID]] == id)
      id2 <- which(data2[[uniqID]] == id)
      true_Delta <- rbind(true_Delta, cbind(rownames(data1[id1, ]), rownames(data2[id2, ])))
    }
    colnames(true_Delta) <- c(label1, label2)
  }

  PIVs <- names(PIVs_config)
  n_pivs <- length(PIVs)

  # --- PIVs must be categorical, restricted to the common support ---------
  for (i in seq_len(n_pivs)) {
    unique1 <- unique(data1[, PIVs[i]])
    unique2 <- unique(data2[, PIVs[i]])
    common <- intersect(unique1, unique2)
    if (length(common) > 100) {
      stop("PIV '", PIVs[i], "' has more than 100 unique values in the common support; PIVs must be categorical.", call. = FALSE)
    }
    missing1 <- setdiff(unique1, common)
    missing1 <- missing1[!is.na(missing1)]
    missing2 <- setdiff(unique2, common)
    missing2 <- missing2[!is.na(missing2)]
    if (length(missing1) > 0 || length(missing2) > 0) {
      warning(sprintf(
        "PIV '%s': values out of the common support%s%s.",
        PIVs[i],
        if (length(missing1) > 0) sprintf(" in data1 [%s]", paste(missing1, collapse = ", ")) else "",
        if (length(missing2) > 0) sprintf(" in data2 [%s]", paste(missing2, collapse = ", ")) else ""
      ), call. = FALSE)
      if (restrict_support_intersection) {
        data1 <- data1[data1[, PIVs[i]] %in% c(NA, common), ]
        data2 <- data2[data2[, PIVs[i]] %in% c(NA, common), ]
      }
    }
  }
  rownames(data1) <- seq_len(nrow(data1))
  rownames(data2) <- seq_len(nrow(data2))

  # --- Warn about strongly associated PIVs
  # (FlexRL assumes conditional independence) ---
  for (i in seq_len(n_pivs)) for (j in seq_len(i - 1)) {
    v_A <- rcompanion::cramerV(table(data1[, PIVs[i]], data1[, PIVs[j]]))
    v_B <- rcompanion::cramerV(table(data2[, PIVs[i]], data2[, PIVs[j]]))
    if (v_A > 0.3 || v_B > 0.3) {
      warning(sprintf(
        "PIVs '%s' and '%s' are associated (Cramer's V = %.2f in data1, %.2f in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.",
        PIVs[i], PIVs[j], v_A, v_B
      ), call. = FALSE)
    }
  }

  # --- Validate PIVs_config -------------------------------------------------
  stopifnot(is.data.frame(data1), is.data.frame(data2))
  if (!is.list(PIVs_config) || length(PIVs_config) == 0 || is.null(names(PIVs_config))) {
    stop("`PIVs_config` must be a non-empty named list, one entry per PIV.", call. = FALSE)
  }
  allowed <- c("dynamics", "boundMistakes", "fixMistakes", "condHazardCov")
  for (k in seq_len(n_pivs)) {
    if (!all(names(PIVs_config[[k]]) %in% allowed)) {
      stop("`PIVs_config[['", PIVs[k], "']]` may only contain: ", paste(allowed, collapse = ", "), ".", call. = FALSE)
    }
  }
  missing1 <- setdiff(PIVs, names(data1))
  missing2 <- setdiff(PIVs, names(data2))
  if (length(missing1) || length(missing2)) {
    stop("PIVs missing from data1: [", paste(missing1, collapse = ", "),
         "]; from data2: [", paste(missing2, collapse = ", "), "]", call. = FALSE)
  }

  bad_dynamics <- !vapply(PIVs_config, function(x) is.character(x$dynamics) && !is.na(x$dynamics), logical(1))
  if (any(bad_dynamics)) {
    stop("`dynamics` must be one of 'stable', 'flexible', 'structured'; problem for: ",
         paste(PIVs[bad_dynamics], collapse = ", "), ".", call. = FALSE)
  }
  is_numeric_or_na <- function(v) is.vector(v) && length(v) == 2 && all(is.na(v) | grepl("^-?[0-9.]+$", v))
  bad_fix <- !vapply(PIVs_config, function(x) is_numeric_or_na(x$fixMistakes), logical(1))
  if (any(bad_fix)) stop("`fixMistakes` must be a length-2 numeric/NA vector; problem for: ", paste(PIVs[bad_fix], collapse = ", "), ".", call. = FALSE)
  bad_bound <- !vapply(PIVs_config, function(x) is_numeric_or_na(x$boundMistakes), logical(1))
  if (any(bad_bound)) stop("`boundMistakes` must be a length-2 numeric/NA vector; problem for: ", paste(PIVs[bad_bound], collapse = ", "), ".", call. = FALSE)

  PIVs_stable <- vapply(PIVs_config, function(x) x$dynamics != "structured", logical(1))
  modelDynaPIVs <- PIVs[!PIVs_stable]

  if (length(modelDynaPIVs) > 0) {
    cov_issues <- character(0)
    for (p in modelDynaPIVs) {
      cov1 <- PIVs_config[[p]]$condHazardCov$cov1
      cov2 <- PIVs_config[[p]]$condHazardCov$cov2
      if (length(cov1) > 0) {
        if (!all(cov1 %in% names(data1))) cov_issues <- c(cov_issues, sprintf("PIV '%s': cov1 not found in data1: %s", p, paste(setdiff(cov1, names(data1)), collapse = ", ")))
        if (!"date" %in% names(data1)) cov_issues <- c(cov_issues, sprintf("PIV '%s' is structured but data1 has no `date` column.", p))
      }
      if (length(cov2) > 0) {
        if (!all(cov2 %in% names(data2))) cov_issues <- c(cov_issues, sprintf("PIV '%s': cov2 not found in data2: %s", p, paste(setdiff(cov2, names(data2)), collapse = ", ")))
        if (!"date" %in% names(data2)) cov_issues <- c(cov_issues, sprintf("PIV '%s' is structured but data2 has no `date` column.", p))
      }
    }
    if (length(cov_issues) > 0) stop("Invalid `condHazardCov` configuration:\n  ", paste(cov_issues, collapse = "\n  "), call. = FALSE)
  }

  if (!is.logical(sameMistakes)) stop("`sameMistakes` must be logical.", call. = FALSE)
  if (sameMistakes) {
    for (k in seq_len(n_pivs)) {
      bm <- PIVs_config[[k]]$boundMistakes
      fm <- PIVs_config[[k]]$fixMistakes
      if (!identical(bm[1], bm[2]) || !identical(fm[1], fm[2])) {
        stop("`sameMistakes = TRUE` requires identical `boundMistakes`/`fixMistakes` for file 1 and file 2 (PIV '", PIVs[k], "').", call. = FALSE)
      }
    }
  }

  # --- Advisory warnings for atypical dynamics/parameter combinations ------
  # For each PIV, `boundMistakes`/`fixMistakes` should match its `dynamics`:
  # stable     -> recommend `boundMistakes`, discourage `fixMistakes`
  # flexible   -> discourage both (nothing to bound/fix: change is not modelled)
  # structured -> discourage `boundMistakes`, recommend `fixMistakes`
  #               (to avoid confounding mistakes with genuine change over time)
  cov_issues <- c()
  for (k in seq_len(n_pivs)) {
    dyn <- PIVs_config[[k]]$dynamics
    bm  <- PIVs_config[[k]]$boundMistakes
    fm  <- PIVs_config[[k]]$fixMistakes

    if (dyn == "stable") {
      if (!all(is.numeric(bm))) {
        cov_issues <- c(cov_issues, sprintf("We recommend bounding the mistakes with `boundMistakes` when `dynamics` is `stable` (PIV: %s).", PIVs[k]))
      }
      if (!all(is.na(fm))) {
        cov_issues <- c(cov_issues, sprintf("We do not recommend fixing the mistakes with `fixMistakes` when `dynamics` is `stable` (PIV: %s).", PIVs[k]))
      }
    } else if (dyn == "flexible") {
      if (!all(is.na(bm))) {
        cov_issues <- c(cov_issues, sprintf("We do not recommend bounding the mistakes with `boundMistakes` when `dynamics` is `flexible` (PIV: %s).", PIVs[k]))
      }
      if (!all(is.na(fm))) {
        cov_issues <- c(cov_issues, sprintf("We do not recommend fixing the mistakes with `fixMistakes` when `dynamics` is `flexible` (PIV: %s).", PIVs[k]))
      }
    } else if (dyn == "structured") {
      if (!all(is.na(bm))) {
        cov_issues <- c(cov_issues, sprintf("We do not recommend bounding the mistakes with `boundMistakes` when `dynamics` is `structured` (PIV: %s).", PIVs[k]))
      }
      if (!all(is.numeric(fm))) {
        cov_issues <- c(cov_issues, sprintf("We recommend fixing the mistakes with `fixMistakes` when `dynamics` is `structured` (PIV: %s).", PIVs[k]))
      }
    }
  }
  if (length(cov_issues) > 0) warning(paste(cov_issues, collapse = "\n  "), call. = FALSE)

  # Configuration entries outside the 4 recognised fields are silently ignored
  # downstream; warn so the user notices a likely typo.
  useless_config <- vapply(PIVs_config, function(x) {
    any(!names(x) %in% allowed)
  }, logical(1))
  if (any(useless_config)) {
    warning(
      "Configuration fields outside of `dynamics`, `boundMistakes`, `fixMistakes`, `condHazardCov` are ignored (PIV: ",
      paste(PIVs[useless_config], collapse = ", "), ").",
      call. = FALSE
    )
  }

  # `condHazardCov` only has an effect when `dynamics = "structured"`; warn if
  # it was supplied for a stable/flexible PIV, since it will be silently unused.
  useless_given_cov <- vapply(PIVs_config, function(x) {
    non_empty <- any(vapply(x[["condHazardCov"]], length, integer(1)) > 0)
    is_stable_or_flexible <- x$dynamics %in% c("stable", "flexible")
    non_empty && is_stable_or_flexible
  }, logical(1))
  if (any(useless_given_cov)) {
    warning(
      sprintf(
        "`condHazardCov` was supplied but `dynamics` is not `structured`, so no dynamics will be modelled (PIV: %s).",
        paste(PIVs[useless_given_cov], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # --- source column + assign the larger file as B --------------------------
  if (!"source" %in% names(data1)) data1$source <- label1
  if (!"source" %in% names(data2)) data2$source <- label2

  swap <- nrow(data1) > nrow(data2)
  encodedA <- if (swap) data2 else data1
  encodedB <- if (swap) data1 else data2
  message(sprintf("'%s' is the larger source, saved as B; '%s' saved as A.", if (swap) label1 else label2, if (swap) label2 else label1))
  if (swap) {
    for (p in modelDynaPIVs) {
      names(PIVs_config[[p]]$condHazardCov) <- gsub("^cov1$", "covB", gsub("^cov2$", "covA", names(PIVs_config[[p]]$condHazardCov)))
    }
    for (k in seq_len(n_pivs)) {
      PIVs_config[[k]]$boundMistakes <- rev(PIVs_config[[k]]$boundMistakes)
      PIVs_config[[k]]$fixMistakes <- rev(PIVs_config[[k]]$fixMistakes)
    }
  } else {
    for (p in modelDynaPIVs) {
      names(PIVs_config[[p]]$condHazardCov) <- gsub("^cov1$", "covA", gsub("^cov2$", "covB", names(PIVs_config[[p]]$condHazardCov)))
    }
  }

  # --- Encode PIVs to natural numbers, pooling levels across A and B --------
  levels_PIVs <- stats::setNames(lapply(PIVs, function(x) levels(factor(as.character(c(encodedA[[x]], encodedB[[x]]))))), PIVs)
  for (x in PIVs) {
    encodedA[[x]] <- as.numeric(factor(as.character(encodedA[[x]]), levels = levels_PIVs[[x]]))
    encodedB[[x]] <- as.numeric(factor(as.character(encodedB[[x]]), levels = levels_PIVs[[x]]))
  }
  nvalues <- stats::setNames(vapply(levels_PIVs, length, integer(1)), PIVs)
  encodedA[PIVs][is.na(encodedA[PIVs])] <- 0  # FlexRL sentinel for missing values
  encodedB[PIVs][is.na(encodedB[PIVs])] <- 0

  list(
    encodedA = encodedA, encodedB = encodedB, Nvalues = nvalues,
    PIVs_config = PIVs_config, sameMistakes = sameMistakes, true_pairs = true_Delta
  )
}

#' Compare distributions of shared variables across data sets
#' (histograms/barplots)
#'
#' Overlays, for each variable in `common_vars`, the empirical distribution
#' in every data set of `data_list` (e.g. a baseline file vs. the linked set).
#'
#' @param data_list Named list of data frames, each containing `common_vars`.
#' @param common_vars Character vector, variables to compare.
#' @param colours_vec Optional vector of colours, one per element of `data_list`.
#'
#' @return `NULL`, invisibly; called for its plotting side effect.
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V2 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V3 = list(dynamics = "flexible",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V4 = list(dynamics = "structured",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(0.03,0.03),
#'                                 condHazardCov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())
#'                                 )
#' )
#' Nval  <- c(6, 7, 8, 9)
#' Pmistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02))
#' Pmissing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
#' condHazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = c(0.7,0.6,0.5))
#'
#' GenData <- DataCreation(
#'   PIVs_config, Nval, NRecords = c(400, 600), Nlinks = 300,
#'   Pmistake, Pmissing, condHazard_params, enforceEstimability = TRUE
#' )
#'
#' PrepData <- prepare_data(GenData$dataSet1, GenData$dataSet2, "1", "2",
#'                      PIVs_config, sameMistakes = TRUE, uniqID = "entityID")
#'
#' fit <- stEM(data = PrepData, StEMIter = 10, StEMBurnin = 5,
#'            GibbsIter = 10, GibbsBurnin = 5, musicOn = FALSE)
#' threshold_strict <- stats::quantile(fit$Delta$x, 0.95)
#' DFLinkedStrict = data.frame( cbind( data.frame(
#'       PrepData$encodedA[fit$Delta[fit$Delta$x > threshold_strict, "i"],]),
#'                                     data.frame(
#'       PrepData$encodedB[fit$Delta[fit$Delta$x > threshold_strict, "j"],]) ) )
# data_list = list(dataBaseline=PrepData$encodedA, dataSelect=DFLinkedStrict)
# common_vars = names(PIVs_config)
# hist_comp(data_list, common_vars)
hist_comp <- function(data_list, common_vars, colours_vec = NULL) {
  if (!is.character(common_vars) || length(common_vars) == 0) {
    stop("`common_vars` must be a non-empty character vector.", call. = FALSE)
  }
  for (data in data_list) {
    stopifnot(is.data.frame(data))
    missing <- setdiff(common_vars, names(data))
    if (length(missing)) stop("`common_vars` not found in a data set: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  N <- length(data_list)
  if (is.null(colours_vec)) colours_vec <- grDevices::rgb(stats::runif(N), stats::runif(N), stats::runif(N), alpha = 0.6)

  for (v in common_vars) {
    if (length(unique(data_list[[1]][, v])) > 20 && is.numeric(data_list[[1]][, v])) {
      b <- diff(range(unique(unlist(data_list[[1]][, v])), na.rm = TRUE)) + 2
      h <- graphics::hist(as.numeric(data_list[[1]][, v]), ylim = c(0, 1), breaks = b, col = colours_vec[1],
                          main = sprintf("Distribution of %s", v), prob = TRUE, xlab = v, ylab = "Density")
      for (d in 2:N) {
        graphics::hist(as.numeric(data_list[[d]][, v]), ylim = c(0, 1), col = colours_vec[d], breaks = h$breaks, prob = TRUE, add = TRUE)
      }
    } else {
      levels_v <- sort(unique(unlist(lapply(data_list, function(d) unique(d[[v]])))))
      dens <- prop.table(table(factor(data_list[[1]][[v]], levels_v)))
      graphics::barplot(dens, ylim = c(0, 1), main = sprintf("Distribution of %s", v), col = colours_vec[1], xlab = v, ylab = "Density")
      for (d in 2:N) {
        graphics::barplot(prop.table(table(factor(data_list[[d]][[v]], levels_v))), ylim = c(0, 1), col = colours_vec[d], add = TRUE)
      }
    }
    graphics::legend("topright", paste(names(data_list), sapply(data_list, nrow), sep = ": obs. "), col = colours_vec, lwd = 10)
  }
  invisible(NULL)
}

#' Maximum Mean Discrepancy between two sets of variables
#'
#' Joint (multivariate) measure of distributional discrepancy between `x`
#' and `y`, using a Gaussian (RBF) kernel with bandwidth set by the median
#' pairwise distance (biased estimator, Gretton et al. 2012).
#'
#' @param x,y Numeric matrices with the same number of columns.
#'
#' @return Numeric, the (biased) MMD estimate.
#' @export
#'
#' @examples
#' x <- matrix(rnorm(50), ncol = 5)
#' y <- matrix(rnorm(30) + 0.5, ncol = 5)
#' mmd(x, y)
mmd <- function(x, y) {
  stopifnot(ncol(x) == ncol(y))
  n <- nrow(x)
  m <- nrow(y)
  dists <- as.matrix(stats::dist(rbind(x, y)))
  sigma <- stats::median(dists) / 2
  k <- exp(-(dists^2) / (2 * sigma^2)) + diag(1e-5, n + m)

  k_x <- k[1:n, 1:n]
  k_y <- k[(n + 1):(n + m), (n + 1):(n + m)]
  k_xy <- k[1:n, (n + 1):(n + m)]

  sum(k_x) / (n * (n - 1)) + sum(k_y) / (m * (m - 1)) - 2 * sum(k_xy) / (n * m)
}


#' Intersection-over-union of two histogram supports
#'
#' @param h1 `histogram` object (as returned by [graphics::hist()]).
#' @param h2 `histogram` object (as returned by [graphics::hist()]).
#'
#' @return Numeric, IoU of the ranges over which `h1` and `h2`
#'   have non-zero counts.
#' @export
#'
#' @examples
#' h1 <- hist(rnorm(200), plot = FALSE)
#' h2 <- hist(rnorm(200) + 1, plot = FALSE)
#' compute_histogram_support_iou(h1, h2)
compute_histogram_support_iou <- function(h1, h2) {
  r1 <- range(h1$breaks[h1$counts > 0])
  r2 <- range(h2$breaks[h2$counts > 0])
  intersection <- max(0, min(r1[2], r2[2]) - max(r1[1], r2[1]))
  union <- max(r1[2], r2[2]) - min(r1[1], r2[1])
  intersection / union
}

#' Proportion of a data set with a given variable at a given level
#'
#' @param df Data frame.
#' @param var Character, column name in `df`.
#' @param level Value to match against `df[[var]]`.
#'
#' @return Numeric, proportion of rows where `df[[var]] == level` (na.rm=TRUE).
#' @export
#'
#' @examples
#' df <- data.frame(colour = sample(c("orange", "purple"), 100, replace = TRUE))
#' compute_proba_control(df, "colour", "purple")
compute_proba_control <- function(df, var, level) {
  mean(df[, var] == level, na.rm=TRUE)
}

#' Augment file B with synthetic records
#'
#' Fits a generative model on `encodedB[, PIVs]` and draws `syntheticSample`
#' new synthetic records from it, appended to `encodedB` with `source =
#' "synthetic"`; used by [compute_FDP_RLwithSynth()] to estimate the false
#' discovery proportion of a record-linkage method without ground truth.
#'
#' @param method One of `"arf"` ([arf::adversarial_rf()]), `"synthpop"`
#'   ([synthpop::syn()]), or `"mice"` ([mice::mice()]).
#' @param encodedA The (already prepared / encoded) data source.
#' @param encodedB The (already prepared / encoded) data source.
#' @param PIVs Character vector, names of the PIVs to synthesise.
#' @param syntheticSample Integer, number of synthetic records to generate.
#' @param restrict_support_intersection Logical; drop synthetic records whose
#'   PIV values fall outside the support shared with `encodedA` (default `TRUE`).
#'
#' @return List with `NewA` (unchanged `encodedA`, support-restricted) and
#'   `NewB` (`encodedB` plus the synthetic records).
#' @export
#'
synthesise <- function(method, encodedA, encodedB, PIVs, syntheticSample, restrict_support_intersection = TRUE) {

  if (!method %in% c("arf", "synthpop", "mice")) {
    stop("`method` must be one of 'arf', 'synthpop', 'mice'.", call. = FALSE)
  }

  if (method == "arf") {
    arf_model <- arf::adversarial_rf(encodedB[, PIVs])
    psi <- arf::forde(arf_model, encodedB[, PIVs])
    syntheticNewB <- arf::forge(psi, syntheticSample)
  } else {
    # synthpop / mice work better with continuous coding for
    # very-high-cardinality PIVs
    for (p in PIVs) {
      if (length(unique(encodedB[, p])) >= 60) encodedB[, p] <- as.numeric(encodedB[, p])
      if (length(unique(encodedA[, p])) >= 60) encodedA[, p] <- as.numeric(encodedA[, p])
    }
    if (method == "synthpop") {
      syntheticNewB <- synthpop::syn(encodedB[, PIVs], k = syntheticSample)$syn
    } else {
      empty <- matrix(NA, nrow = syntheticSample, ncol = length(PIVs), dimnames = list(NULL, PIVs))
      imputed <- mice::complete(mice::mice(rbind(encodedB[, PIVs], empty), m = 1))
      syntheticNewB <- imputed[(nrow(encodedB) + 1):nrow(imputed), ]
    }
  }
  rownames(syntheticNewB) <- seq_len(nrow(syntheticNewB))

  extraCols <- setdiff(names(encodedB), c(PIVs, "localID", "source"))
  for (col in extraCols) syntheticNewB[, col] <- NA

  syntheticNewB$localID <- nrow(encodedB) + seq_len(nrow(syntheticNewB))
  syntheticNewB$source <- "synthetic"
  encodedNewB <- rbind(encodedB, syntheticNewB)

  for (p in PIVs) {
    encodedNewB[, p] <- as.integer(encodedNewB[, p])
    encodedA[, p] <- as.integer(encodedA[, p])
  }

  for (p in PIVs) {
    common <- intersect(unique(encodedA[, p]), unique(encodedNewB[, p]))
    missing1 <- setdiff(unique(encodedA[, p]), common)
    missing2 <- setdiff(unique(encodedNewB[, p]), common)
    if (length(missing1) > 0 || length(missing2) > 0) {
      warning(sprintf("PIV '%s': synthetic records out of the common support were %s.",
                      p, if (restrict_support_intersection) "removed" else "kept (support not restricted)"), call. = FALSE)
      if (restrict_support_intersection) {
        encodedA <- encodedA[encodedA[, p] %in% c(NA, common), ]
        encodedNewB <- encodedNewB[encodedNewB[, p] %in% c(NA, common), ]
      }
    }
  }
  rownames(encodedA) <- seq_len(nrow(encodedA))
  rownames(encodedNewB) <- seq_len(nrow(encodedNewB))

  list(NewA = encodedA, NewB = encodedNewB)
}

#' Standardised mean difference between a selected set and a baseline
#'
#' @param dataSelect Data frame to compare (e.g. the linked
#'   set vs. the original file).
#' @param dataBaseline Data frame to compare (e.g. the linked
#'   set vs. the original file).
#' @param var Character, column name to compare.
#' @param continuous Logical; if `TRUE`, compares means of `var` directly; if
#'   `FALSE`, compares the proportion at each observed level of `var`
#'   (default `TRUE`).
#'
#' @return Named list, one SMD per variable (`continuous = TRUE`) or per level
#'   (`continuous = FALSE`).
#' @export
#'
#' @examples
#' base <- data.frame(age = rnorm(200, 40, 10), sex = sample(c("M", "F"),
#'                     200, replace = TRUE))
#' select <- data.frame(age = rnorm(80, 43, 10), sex = sample(c("M", "F"), 80,
#'                     replace = TRUE, prob = c(0.6, 0.4)))
#' SMD(select, base, "age")
#' SMD(select, base, "sex", continuous = FALSE)
SMD <- function(dataSelect, dataBaseline, var, continuous = TRUE) {
  res <- list()
  if (continuous) {
    res[[var]] <- (mean(dataSelect[, var], na.rm = TRUE) - mean(dataBaseline[, var], na.rm = TRUE)) / stats::sd(dataBaseline[, var], na.rm = TRUE)
  } else {
    for (val in sort(unique(dataSelect[, var]))) {
      res[[paste(var, val, sep = "_")]] <-
        (mean(dataSelect[, var] == val, na.rm = TRUE) - mean(dataBaseline[, var] == val, na.rm = TRUE)) /
        stats::sd(dataBaseline[, var] == val, na.rm = TRUE)
    }
  }
  res
}

#' Plot the distribution of linkage scores
#'
#' @param total_nbr_pairs Integer, total number of candidate pairs considered
#'   (`nrow(A) * nrow(B)`).
#' @param vec_linked_scores Numeric vector, linkage scores of the pairs above
#'   0 (e.g. `Delta$x`).
#'
#' @return `NULL`, invisibly; called for its plotting side effect.
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V2 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V3 = list(dynamics = "flexible",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V4 = list(dynamics = "structured",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(0.03,0.03),
#'                                 condHazardCov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())
#'                                 )
#' )
#' Nval  <- c(6, 7, 8, 9)
#' Pmistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02))
#' Pmissing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
#' condHazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = c(0.7,0.6,0.5))
#'
#' GenData <- DataCreation(
#'   PIVs_config, Nval, NRecords = c(400, 600), Nlinks = 300,
#'   Pmistake, Pmissing, condHazard_params, enforceEstimability = TRUE
#' )
#'
#' PrepData <- prepare_data(GenData$dataSet1, GenData$dataSet2, "1", "2",
#'                      PIVs_config, sameMistakes = TRUE, uniqID = "entityID")
#'
#' fit <- stEM(data = PrepData, StEMIter = 10, StEMBurnin = 5,
#'            GibbsIter = 10, GibbsBurnin = 5, musicOn = FALSE)
#'
#' plot_linkage_score(nrow(PrepData$encodedA)*nrow(PrepData$encodedB),
#'                     fit$Delta$x)
plot_linkage_score <- function(total_nbr_pairs, vec_linked_scores) {
  breaks <- seq(0, 1, by = 0.05)
  h <- graphics::hist(vec_linked_scores, plot = FALSE, breaks = breaks)
  h$counts[1] <- h$counts[1] + (total_nbr_pairs - length(vec_linked_scores))
  h$counts <- pmax(log10(h$counts), 0)
  graphics::plot(h, xlim = c(0, 1), main = "Linkage score distribution",
                 xlab = "Posterior linkage score", ylab = "Log10 frequency")
  invisible(NULL)
}

#' Model-based false discovery proportion at a given threshold
#'
#' @param vec_link_scores Numeric vector, linkage scores of the candidate pairs.
#' @param threshold Numeric, score threshold above which a pair is declared
#'   linked.
#'
#' @return Numeric, `1 - mean(score | score > threshold)`, the estimated FDP at
#'   `threshold`.
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V2 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V3 = list(dynamics = "flexible",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V4 = list(dynamics = "structured",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(0.03,0.03),
#'                                 condHazardCov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())
#'                                 )
#' )
#' Nval  <- c(6, 7, 8, 9)
#' Pmistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02))
#' Pmissing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
#' condHazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = c(0.7,0.6,0.5))
#'
#' GenData <- DataCreation(
#'   PIVs_config, Nval, NRecords = c(400, 600), Nlinks = 300,
#'   Pmistake, Pmissing, condHazard_params, enforceEstimability = TRUE
#' )
#'
#' PrepData <- prepare_data(GenData$dataSet1, GenData$dataSet2, "1", "2",
#'                      PIVs_config, sameMistakes = TRUE, uniqID = "entityID")
#'
#' fit <- stEM(data = PrepData, StEMIter = 10, StEMBurnin = 5,
#'            GibbsIter = 10, GibbsBurnin = 5, musicOn = FALSE)
#'
#' compute_FDP_RLmodelSpecific(fit$Delta$x, 0.5)
compute_FDP_RLmodelSpecific <- function(vec_link_scores, threshold) {
  linked <- vec_link_scores > threshold
  1 - sum(vec_link_scores[linked]) / sum(linked)
}


#' Estimate the false discovery proportion of a record-linkage method via
#'synthetic augmentation
#'
#' Repeatedly augments file B with synthetic records ([synthesise()]), runs
#' the chosen record-linkage method ([FDPinRL_FlexRL()] and siblings), and
#' compares the synthetic-vs-real proportion among linked pairs to estimate
#' the false discovery proportion, for a range of score thresholds
#' (0.50 to 0.95). See https://doi.org/10.1002/sim.70292 for the method.
#'
#' @param SynthMethod Passed to [synthesise()]: `"arf"`, `"synthpop"`, or
#'   `"mice"`.
#' @param fileA The prepared data source (`fileA` must be the
#'   smaller one).
#' @param fileB The prepared data source (`fileB` must be the
#'   larger one).
#' @param PIVs Character vector, names of the PIVs.
#' @param subsample_size Integer, number of synthetic records to generate per
#'   iteration (default: 10% of `nrow(fileB)`).
#' @param restrict_support_intersection Passed to [synthesise()].
#' @param maxIter4CV Integer, max number of retries per iteration if no valid
#'   FDP estimate is obtained.
#' @param NIter Integer, number of augmentation iterations to average over.
#' @param RLMethod One of `"multilink"`, `"fastLink"`, `"BRL"`,
#'   `"reclin2"`, `"diyar"`, `"fedmatch"`, `"FlexRL"`.
#' @param ... Extra arguments forwarded to the chosen `FDPinRL_*()` wrapper
#'   (i.e. to the underlying record-linkage package).
#'
#' @return List with `FDP_scoring_estimator`, `FDP_synth_estimator`,
#'   `Linked_obs_pairs`: data frames (`NIter` rows x 10 thresholds).
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V2 = list(dynamics = "stable",
#'                                 boundMistakes = c(0.10,0.10),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V3 = list(dynamics = "flexible",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(NA,NA)
#'                                 ),
#'                     V4 = list(dynamics = "structured",
#'                                 boundMistakes = c(NA,NA),
#'                                 fixMistakes = c(0.03,0.03),
#'                                 condHazardCov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())
#'                                 )
#' )
#' Nval  <- c(10, 11, 12, 13)
#' Pmistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02))
#' Pmissing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
#' condHazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = c(0.7,0.6,0.5))
#'
#' GenData <- DataCreation(
#'   PIVs_config, Nval, NRecords = c(400, 600), Nlinks = 300,
#'   Pmistake, Pmissing, condHazard_params, enforceEstimability = TRUE
#' )
#'
#' PrepData <- prepare_data(GenData$dataSet1, GenData$dataSet2, "1", "2",
#'                      PIVs_config, sameMistakes = TRUE, uniqID = "entityID")
#'
#' fdp_res_brl <- compute_FDP_RLwithSynth("arf", PrepData$encodedA,
#'       PrepData$encodedB, names(PIVs_config),
#'       subsample_size=NULL, restrict_support_intersection=TRUE,
#'       maxIter4CV=3, NIter=5,
#'       RLMethod = "BRL",
#'       flds = names(PIVs_config),
#'       types = rep("bi",length(PIVs_config))
#'       )
compute_FDP_RLwithSynth <- function(SynthMethod, fileA, fileB, PIVs, subsample_size = NULL,
                                    restrict_support_intersection = TRUE, maxIter4CV = 10, NIter = 10,
                                    RLMethod, ...) {

  nbrRealRecordsA <- nrow(fileA)
  nbrRealRecordsB <- nrow(fileB)
  if (nbrRealRecordsA > nbrRealRecordsB) stop("`fileA` must be smaller than `fileB`.", call. = FALSE)

  supported <- c("multilink", "fastLink", "BRL", "reclin2", "diyar", "fedmatch", "FlexRL")
  if (!RLMethod %in% supported) stop("`RLMethod` must be one of: ", paste(supported, collapse = ", "), ".", call. = FALSE)

  if (is.null(subsample_size)) subsample_size <- as.integer(0.10 * nbrRealRecordsB)

  thresholds <- seq(0.5, 0.95, by = 0.05)
  th_names <- c("thresholds: 0.50", sprintf("%.2f", thresholds[-1]))
  FDP_RLwithSynth_results <- stats::setNames(data.frame(matrix(NA, NIter, 10)), th_names)
  FDP_RLmodelSpecific_results <- stats::setNames(data.frame(matrix(NA, NIter, 10)), th_names)
  Real_linked_results <- stats::setNames(data.frame(matrix(NA, NIter, 10)), th_names)
  arguments <- list(...)

  run_RL <- function(NewA, NewB) {
    switch(RLMethod,
           FlexRL     = FDPinRL_FlexRL(NewA, NewB, arguments),
           fedmatch   = FDPinRL_fedmatch(NewA, NewB, arguments),
           reclin2    = FDPinRL_reclin2(NewA, NewB, arguments),
           BRL        = FDPinRL_BRL(NewA, NewB, arguments),
           fastLink   = FDPinRL_fastLink(NewA, NewB, arguments),
           multilink  = FDPinRL_multilink(NewA, NewB, arguments),
           diyar      = FDPinRL_diyar(NewA, NewB, arguments)
    )
  }

  for (i in seq_len(NIter)) {

    Newdata <- synthesise(SynthMethod, fileA[, c(PIVs, "localID", "source")],
                          fileB[, c(PIVs, "localID", "source")], PIVs, subsample_size, restrict_support_intersection)

    anyValidEstimate <- FALSE
    countTmp <- 0

    while (countTmp < maxIter4CV && !anyValidEstimate) {
      res <- run_RL(Newdata$NewA, Newdata$NewB)
      result_index_linked_A <- res[[1]]
      result_index_linked_B <- res[[2]]
      result_vec_linked_scores <- res[[3]]

      synth_fdp <- function(index_linked_A, index_linked_B) {
        real <- index_linked_B <= nbrRealRecordsB & index_linked_A <= nbrRealRecordsA
        synthfp <- sum(!real)
        N_linked <- length(index_linked_A)
        list(
          fdp_synth = (synthfp * (nbrRealRecordsB / subsample_size)) / (N_linked - synthfp),
          n_real = sum(real)
        )
      }

      if (is.null(result_vec_linked_scores)) {
        N_linked <- length(result_index_linked_A)
        if (N_linked > 0) {
          est <- synth_fdp(result_index_linked_A, result_index_linked_B)
          FDP_RLwithSynth_results[i, 1] <- est$fdp_synth
          Real_linked_results[i, 1] <- est$n_real
        } else {
          warning(sprintf("Nothing linked at iteration %s with default parameters.", i), call. = FALSE)
          FDP_RLwithSynth_results[i, 1] <- 0
          Real_linked_results[i, 1] <- 0
        }
      } else {
        for (j in seq_along(thresholds)) {
          keep <- result_vec_linked_scores > thresholds[j]
          N_linked <- sum(keep, na.rm = TRUE)
          if (N_linked > 0) {
            est <- synth_fdp(result_index_linked_A[keep], result_index_linked_B[keep])
            FDP_RLwithSynth_results[i, j] <- est$fdp_synth
            FDP_RLmodelSpecific_results[i, j] <- compute_FDP_RLmodelSpecific(result_vec_linked_scores, thresholds[j])
            Real_linked_results[i, j] <- est$n_real
          } else {
            if (j == 1) warning(sprintf("Nothing linked at iteration %s at threshold 0.50.", i), call. = FALSE)
            FDP_RLwithSynth_results[i, j:10] <- 0
            FDP_RLmodelSpecific_results[i, j:10] <- 0
            Real_linked_results[i, j:10] <- 0
            break
          }
        }
      }

      countTmp <- countTmp + 1
      anyValidEstimate <- any(!is.na(FDP_RLwithSynth_results[i, ]) & FDP_RLwithSynth_results[i, ] <= 1 & Real_linked_results[i, ] > 0)
    }

    if (countTmp == maxIter4CV && !anyValidEstimate) {
      stop(sprintf(
        "No valid FDP estimate after %s attempts at iteration %s. Increase `maxIter4CV`, or the estimator may be unreliable for this method/data.",
        maxIter4CV, i
      ), call. = FALSE)
    }
  }

  ToShow <- data.frame(
    `FDP model score estimator` = round(colMeans(FDP_RLmodelSpecific_results, na.rm = TRUE), 2),
    `FDP synth data estimator`  = round(colMeans(FDP_RLwithSynth_results, na.rm = TRUE), 2),
    `Linked obs. pairs`         = round(colMeans(Real_linked_results, na.rm = TRUE), 2),
    check.names = FALSE
  )
  message(sprintf("%s results (average over %s iterations):", RLMethod, NIter))
  print(t(ToShow))

  list(
    FDP_scoring_estimator = FDP_RLmodelSpecific_results,
    FDP_synth_estimator = FDP_RLwithSynth_results,
    Linked_obs_pairs = Real_linked_results
  )
}
