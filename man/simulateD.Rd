% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FlexRecordLinkage.r
\name{simulateD}
\alias{simulateD}
\title{simulateD}
\usage{
simulateD(
  data,
  linksR,
  sumRowD,
  sumColD,
  truepivsA,
  truepivsB,
  gamma,
  eta,
  alpha,
  phi
)
}
\arguments{
\item{data}{A list with elements:
\itemize{
\item A: the smallest data source (encoded: the categorical values of the PIVs have to be mapped to sets of natural numbers and missing values are encoded as 0).
\item B: the largest data source (encoded).
\item Nvalues: A vector (of size number of PIVs) with the number fo unique values per PIVs (in the order of the PIVs defined in PIVs_config).
\item PIVs_config: A list (of size number of PIVs) where element names are the PIVs and element values are lists with elements: stable (boolean for whether the PIV is stable), conditionalHazard (boolean for whether there are external covariates available to model instability, only required if stable is FALSE), pSameH.cov.A and pSameH.covB (vectors with strings corresponding to the names of the covariates to use to model instability from file A and file B, only required if stable is FALSE, empty vectors may be provided if conditionalHazard is FALSE).
\item controlOnMistakes: A vector (of size number of PIVs) of booleans indicating potential bounds on the mistakes probabilities for each PIV. For each PIV, if TRUE there will be control on mistake and the mistake probability will not go above 10\%. If FALSE there is no bound on the probability of mistake. WATCH OUT, if you suspect that a variable is unstable but you do not have data to model its dynamics the boolean value should be set to FALSE to allow the parameter for mistake to adapt for the instability. However if you model this instability, the boolean value should be set to TRUE to help the algorithm differenciate the mistakes from the changes over time.
\item sameMistakes: A boolean value for whether there should be one parameter for the mistakes in A and B or whether each source should have its own parameter. Setting sameMistakes=TRUE is recommended in case of small data sources; the estimation with 2 parameters in that case will fail to capture the mistakes correctly while 1 parameter will be more adapted.
\item phiMistakesAFixed A vector (of size number of PIVs) of booleans indicating whether the parameters for mistakes should be fixed in case of instability. It should be FALSE, except for unstable PIVs for which it may be set to TRUE in order to avoid estimability problems between the parameter for mistake and the parameter for changes across time.
\item phiMistakesBFixed A vector (of size number of PIVs) of booleans indicating whether the parameters for mistakes should be fixed in case of instability. It should be FALSE, except for unstable PIVs for which it may be set to TRUE in order to avoid estimability problems between the parameter for mistake and the parameter for changes across time.
\item phiForMistakesA A vector (of size number of PIVs) of NA or fixed values for the parameters for mistakes. It should be NA, except for unstable PIVs for which one wants to fix the parameter to avoid estimability problem (as indicated with the boolean values in phiMistakesAFixed). In that case it should be set the the expected value for the probability of mistake. If you have no idea: you can put it to 0, the algorithm is quite robust to wrongly fixed parameters.
\item phiForMistakesB A vector (of size number of PIVs) of NA or fixed values for the parameters for mistakes. It should be NA, except for unstable PIVs for which one wants to fix the parameter to avoid estimability problem (as indicated with the boolean values in phiMistakesBFixed). In that case it should be set the the expected value for the probability of mistake. If you have no idea: you can put it to 0, the algorithm is quite robust to wrongly fixed parameters.
}}

\item{linksR}{A matrix of 2 columns with indices of the linked records.}

\item{sumRowD}{A boolean vector indicating, for each row of the linkage matrix, i.e. for each record in the smallest file A, whether the record has a link in B or not.}

\item{sumColD}{A boolean vector indicating, for each column of the linkage matrix, i.e. for each record in the largest file B, whether the record has a link in A or not.}

\item{truepivsA}{A matrix of the shape of data source A, representing the true values of the PIVs underlying the registered values present in A.}

\item{truepivsB}{A matrix of the shape of data source B, representing the true values of the PIVs underlying the registered values present in B.}

\item{gamma}{The proportion of linked records as a fraction of the smallest file.}

\item{eta}{The distribution weights for the PIVs.}

\item{alpha}{The parameter involved in the survival model for the probability of true values to coincide (parameter for instability).}

\item{phi}{The proportion of mistakes and missing for the PIVs.}
}
\value{
A list with:
\itemize{
\item new set of links
\item new sumRowD
\item new sumColD
\item new value of the complete log likelihood
\item new number fo linked records
}
}
\description{
simulateD
}
\examples{
PIVs_config = list( V1 = list(stable = TRUE),
                    V2 = list(stable = TRUE),
                    V3 = list(stable = TRUE),
                    V4 = list(stable = TRUE),
                    V5 = list( stable = FALSE,
                               conditionalHazard = FALSE,
                               pSameH.cov.A = c(),
                               pSameH.cov.B = c()) )
PIVs = names(PIVs_config)
PIVs_stable = sapply(PIVs_config, function(x) x$stable)
Nval = c(6, 7, 8, 9, 15)
NRecords = c(13, 15)
Nlinks = 6
PmistakesA = c(0.02, 0.02, 0.02, 0.02, 0.02)
PmistakesB = c(0.02, 0.02, 0.02, 0.02, 0.02)
PmissingA = c(0.007, 0.007, 0.007, 0.007, 0.007)
PmissingB = c(0.007, 0.007, 0.007, 0.007, 0.007)
moving_params = list(V1=c(),V2=c(),V3=c(),V4=c(),V5=c(0.28))
enforceEstimability = TRUE
DATA = DataCreation( PIVs_config,
                     Nval,
                     NRecords,
                     Nlinks,
                     PmistakesA,
                     PmistakesB,
                     PmissingA,
                     PmissingB,
                     moving_params,
                     enforceEstimability)
A                    = DATA$A
B                    = DATA$B
Nvalues              = DATA$Nvalues

encodedA = A
encodedB = B

encodedA[,PIVs][ is.na(encodedA[,PIVs]) ] = 0
encodedB[,PIVs][ is.na(encodedB[,PIVs]) ] = 0

dataForStEM = list( A                    = encodedA,
                    B                    = encodedB,
                    Nvalues              = Nvalues,
                    PIVs_config          = PIVs_config,
                    controlOnMistakes    = c(TRUE, TRUE, TRUE, TRUE, TRUE),
                    sameMistakes         = TRUE,
                    phiMistakesAFixed    = TRUE,
                    phiMistakesBFixed    = TRUE,
                    phiForMistakesA      = c(NA, NA, NA, NA, 0),
                    phiForMistakesB      = c(NA, NA, NA, NA, 0)
)

initDeltaMap()
linksR = base::matrix(0,0,2)
linksCpp = linksR
sumRowD = rep(0, nrow(dataForStEM$A))
sumColD = rep(0, nrow(dataForStEM$B))
nlinkrec = 0
survivalpSameH = base::matrix(1, nrow(linksR), length(dataForStEM$Nvalues))
gamma = 0.5
eta = lapply(dataForStEM$Nvalues, function(x) rep(1/x,x))
phi = lapply(dataForStEM$Nvalues, function(x)  c(0.9,0.9,0.1,0.1))
nCoefUnstable = lapply( seq_along(PIVs_stable),
                        function(idx)
                        if(PIVs_stable[idx]){ 0 }
                        else{
                        ncol(dataForStEM$A[, dataForStEM$PIVs_config[[idx]]$pSameH.cov.A,
                             drop=FALSE]) +
                        ncol(dataForStEM$B[, dataForStEM$PIVs_config[[idx]]$pSameH.cov.B,
                             drop=FALSE]) +
                        1 } )
alpha = lapply( seq_along(PIVs_stable),
                function(idx) if(PIVs_stable[idx]){ c(-Inf) }else{
                  rep(log(0.05), nCoefUnstable[[idx]]) })
newTruePivs = simulateH(data=dataForStEM, links=linksCpp, survivalpSameH=survivalpSameH,
                        sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
truepivsA = newTruePivs$truepivsA
truepivsB = newTruePivs$truepivsB
Dsample = simulateD(data=dataForStEM, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD,
                    truepivsA=truepivsA, truepivsB=truepivsB,
                    gamma=gamma, eta=eta, alpha=alpha, phi=phi)
linksCpp = Dsample$links
linksR = linksCpp + 1
}
