#' FlexRL
#'
#' A Flexible Model For Record Linkage
#'
#' The example below aims to link 2 synthetic data sources, with 5 PIVs.
#' PIVs are stable (not changing over time, probability of mistakes could be
#' bounded), flexible (dynamic but no information to model changes over time)
#' or structured (dynamic and information to model changes over time,
#' probability of mistakes could be fixed). We may need to fix the mistake
#' parameter of the 5th dynamic PIV to avoid estimability problems here.
#' We know the true linkage structure in this example so we can compute
#' performances of the method at the end.
#'
#' Methodological paper: \doi{10.1093/jrsssc/qlaf016}.
#' Experiments repository of the methodological paper:
#' https://github.com/robachowyk/FlexRL-experiments.
#' More details in the documentation of the main algorithm ?FlexRL::StEM.
#'
#' @author Kayané Robach
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib FlexRL, .registration=TRUE
#' @name FlexRL
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
#' fit <- stEM(data = PrepData, StEMIter = 20, StEMBurnin = 10,
#'            GibbsIter = 20, GibbsBurnin = 10, musicOn = FALSE)
#' head(fit$Delta[fit$Delta$x > 0.5, ])
#'
#' DeltaResult = fit$Delta
#' colnames(DeltaResult) = c("idxA","idxB","LinkageScores")
#' DeltaResult = DeltaResult[DeltaResult$LinkageScores>0.5,]
#'
#' results = data.frame( matrix(NA, nrow=5, ncol=0) )
#' rownames(results) = c("tp","fp","fn","fdp","sensitivity")
#' if(nrow(DeltaResult)>1){
#'   linked_pairs    = do.call(paste, c(DeltaResult[,c("idxA","idxB")], list(sep="_")))
#'   true_pairs      = do.call(paste, c(PrepData$true_pairs, list(sep="_")))
#'   truepositive    = length( intersect(linked_pairs, true_pairs) )
#'   falsepositive   = length( setdiff(linked_pairs, true_pairs) )
#'   falsenegative   = length( setdiff(true_pairs, linked_pairs) )
#'   fdp             = falsepositive / (truepositive + falsepositive)
#'   sensitivity     = truepositive / (truepositive + falsenegative)
#'   results[,"FlexRL"] = c(truepositive,falsepositive,falsenegative,fdp,sensitivity)
#' }
#'
#' rl_agreement(PrepData$encodedA, PrepData$encodedB, names(PIVs_config),
#'                        DeltaResult[,c("idxA","idxB")], PrepData$true_pairs)
#'
#' diag <- rl_diagnostics(fit, PrepData$encodedA, PrepData$encodedB,
#'                           names(PIVs_config), 0.75, PrepData$true_pairs, 5)
#' diag
#' summary(diag)
#' plot(diag,"scores")
#' plot(diag,"distributions")
#' plot(diag,"smd")
#' plot(diag,"convergence", ask=FALSE)
#'
NULL

.onLoad <- function(...) {
  base::packageStartupMessage("If you are happy with FlexRL, please cite us! Also, if you are unhappy, please cite us anyway.\nHERE ADD\nbibtex format in CITATION.", appendLF = TRUE)
}
