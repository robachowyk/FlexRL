#' Title
#'
#' @param PIVs_config A list (of size number of PIVs) where element names are the PIVs and element values are lists with elements: stable (boolean for whether the PIV is stable), conditionalHazard (boolean for whether there are external covariates available to model instability, only required if stable is FALSE), pSameH.cov.A and pSameH.covB (vectors with strings corresponding to the names of the covariates to use to model instability from file A and file B, only required if stable is FALSE, empty vectors may be provided if conditionalHazard is FALSE)
#' @param Nval A vector (of size number of PIVs) with the number fo unique values per PIVs (in the order of the PIVs defined in PIVs_config)
#' @param NRecords A vector (of size 2) with the number of records to be generated in file A and in file B
#' @param Nlinks An integer with the number of links (record referring to the same entities) to be generated
#' @param PmistakesA A vector (of size number of PIVs) with the proportion of mistakes to be generated per PIVs in file A (in the order of the PIVs defined in PIVs_config)
#' @param PmistakesB A vector (of size number of PIVs) with the proportion of mistakes to be generated per PIVs in file B (in the order of the PIVs defined in PIVs_config)
#' @param PmissingA A vector (of size number of PIVs) with the proportion of missing to be generated per PIVs in file A (in the order of the PIVs defined in PIVs_config)
#' @param PmissingB A vector (of size number of PIVs) with the proportion of missing to be generated per PIVs in file A (in the order of the PIVs defined in PIVs_config)
#' @param moving_params A list (of size number of PIVs) where element names are the PIVs and element values are vectors (of size: 1 + number of covariates to use from A + number of covariates to use from B) with the log hazards coefficient (1st one: log baseline hazard, then: the coefficients for conditional hazard covariates from A, then: the coefficients for conditional hazard covariates from B)
#' @param enforceEstimability A boolean value for whether half of the links should have near-0 time gaps (useful for modeling instability and avoiding estimability issues as discussed in the paper)
#'
#' @return A list with generated dataframe A, dataframe B, vector of Nvalues (Nval), vector of TimeDifference (for the links, when thewre is instability), matrix proba_same_H (number of links, number fo PIVs) with the proba that true values coincide (e.g. 1 - proba of moving)
#' @export
#'
#' @examples
#' PIVs_config = list( V1 = list(stable = TRUE),
#'                     V2 = list(stable = TRUE),
#'                     V3 = list(stable = TRUE),
#'                     V4 = list(stable = TRUE),
#'                     V5 = list( stable = FALSE,
#'                                conditionalHazard = FALSE,
#'                                pSameH.cov.A = c(),
#'                                pSameH.cov.B = c()) )
#' Nval = c(6, 7, 8, 9, 15)
#' NRecords = c(500, 800)
#' Nlinks = 300
#' PmistakesA = c(0.02, 0.02, 0.02, 0.02, 0.02)
#' PmistakesB = c(0.02, 0.02, 0.02, 0.02, 0.02)
#' PmissingA = c(0.007, 0.007, 0.007, 0.007, 0.007)
#' PmissingB = c(0.007, 0.007, 0.007, 0.007, 0.007)
#' moving_params = list(V1=c(),V2=c(),V3=c(),V4=c(),V5=c(0.28))
#' enforceEstimability = TRUE
#' DataCreation( PIVs_config,
#'               Nval,
#'               NRecords,
#'               Nlinks,
#'               PmistakesA,
#'               PmistakesB,
#'               PmissingA,
#'               PmissingB,
#'               moving_params,
#'               enforceEstimability)
DataCreation = function(PIVs_config, Nval, NRecords, Nlinks, PmistakesA, PmistakesB, PmissingA, PmissingB, moving_params, enforceEstimability){

  PIVs = names(PIVs_config)
  PIVs_stable = sapply(PIVs_config, function(x) x$stable)
  PIVs_conditionalHazard_variables = sapply(PIVs_config, function(x) if(!x$stable){x$conditionalHazard}else{FALSE})

  for(i in 1:2)
  {
    dataSet=c()
    for(u in 1:length(Nval))
    {
      # Simulate different probabilities for the true values (xp = exp(0 *(0:(Nval[u]-1))) for uniform distribution of PIVs)
      xp            = exp(0.23 *(0:(Nval[u]-1)))
      probx         = xp/sum(xp)
      dataSet       = cbind(dataSet, sample(1:Nval[u], NRecords[i], replace=TRUE, prob=probx))
    }
    dataSet         = as.data.frame(dataSet)
    names(dataSet)  = PIVs
    assign( paste("dataSet", i, sep=""), dataSet )
  }

  # Add overlapping units
  dataSet2[1:Nlinks,] = dataSet1[1:Nlinks,]

  # Add typos
  for(x in 1:ncol(dataSet1))
  {
    biased = as.logical(stats::rbinom(nrow(dataSet1),1,PmistakesA[x]))
    if(sum(biased)>0)
      dataSet1[,x][biased] = sapply(dataSet1[,x][biased], function(i) sample((1:Nval[x])[-c(i)], 1))
  }

  for(x in 1:ncol(dataSet2))
  {
    biased = as.logical(stats::rbinom(nrow(dataSet2),1,PmistakesB[x]))
    if(sum(biased)>0)
      dataSet2[,x][biased] = sapply(dataSet2[,x][biased], function(i) sample((1:Nval[x])[-c(i)], 1))
  }

  # Add missings
  for(x in 1:ncol(dataSet1))
  {
    biased = as.logical(stats::rbinom(nrow(dataSet1),1,PmissingA[x]))
    if(sum(biased)>0)
      dataSet1[,x][biased] = NA
  }

  for(x in 1:ncol(dataSet2))
  {
    biased = as.logical(stats::rbinom(nrow(dataSet2),1,PmissingB[x]))
    if(sum(biased)>0)
      dataSet2[,x][biased] = NA
  }

  instability = any(!PIVs_stable)
  dataSet1$change                   = FALSE
  dataSet2$change                   = FALSE
  if(instability)
  {
    # Add registration time
    dataSet1$date                     = stats::runif(nrow(dataSet1), 0, 3)
    dataSet2$date                     = stats::runif(nrow(dataSet2), 3, 6)

    if(enforceEstimability){
      nullTimeDiff                      = as.integer(Nlinks/2)
      dataSet1[1:nullTimeDiff, "date"]  = stats::runif(nullTimeDiff, 0.00, 0.01)
      dataSet2[1:nullTimeDiff, "date"]  = stats::runif(nullTimeDiff, 0.00, 0.01)
    }

    TimeDifference                    = abs( dataSet2[1:Nlinks, "date"] - dataSet1[1:Nlinks, "date"] )

    intercept                         = rep(1, Nlinks)
    proba_same_H                      = base::matrix(1,Nlinks,length(PIVs))
    for(k in 1:length(PIVs)){
      if(!PIVs_stable[[k]]){
        # unstable piv
        if(PIVs_conditionalHazard_variables[[k]]){
          # create covariates
          for(c in 1:length(PIVs_config[[k]]$pSameH.cov.A)){
            covariate = PIVs_config[[k]]$pSameH.cov.A[[c]]
            dataSet1[, covariate] = stats::rnorm(nrow(dataSet1), 1, 1)
          }
          for(c in 1:length(PIVs_config[[k]]$pSameH.cov.B)){
            covariate = PIVs_config[[k]]$pSameH.cov.B[[c]]
            dataSet2[, covariate] = stats::rnorm(nrow(dataSet2), 2, 1)
          }
          cov                         = cbind( intercept,
                                               dataSet1[1:Nlinks, PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                                               dataSet2[1:Nlinks, PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
        }else{
          cov                         = cbind( intercept )
        }

        # generate instability
        for(i in 1:Nlinks)
        {
          is_not_moving = stats::rbinom(1, 1, proba_same_H[i,k])
          if(!is_not_moving)
          {
            if(!is.na(dataSet1[i,x])){
              dataSet2[i,x]         = sample((1:Nval[x])[-c(dataSet1[i,x])], 1)
              dataSet2[i,"change"]  = TRUE
            }
          }
        }
      }
    }
  }

  A = dataSet1
  B = dataSet2

  # Recode the PIVs
  levels_PIVs = lapply(PIVs, function(x) levels(factor(as.character(c(A[,x], B[,x])))))

  for(i in 1:length(PIVs))
  {
    A[,PIVs[i]] = as.numeric(factor(as.character(A[,PIVs[i]]), levels=levels_PIVs[[i]]))
    B[,PIVs[i]] = as.numeric(factor(as.character(B[,PIVs[i]]), levels=levels_PIVs[[i]]))
  }

  Nvalues       = sapply(levels_PIVs, length)

  A$localID     = 1:nrow(A)
  B$localID     = 1:nrow(B)

  A$source      = "A"
  B$source      = "B"

  list(A=A, B=B, Nvalues=Nvalues, TimeDifference=TimeDifference, proba_same_H=proba_same_H)
}

#' Title
#'
#' @param nCoefUnstable .
#' @param stable .
#'
#' @return .
#' @export
#'
#' @examples .
createDataAlpha <- function(nCoefUnstable, stable){
  nCoef = nCoefUnstable + 2 # cov A + cov B + 1 + Hequal + times
  if (!stable){ return ( data.frame(base::matrix(nrow = 0, ncol = nCoef)) ) } }

#' Title
#'
#' @param Brecords .
#' @param sumD .
#'
#' @return .
#' @export
#'
#' @examples .
logPossibleConfig = function(Brecords,sumD)
{
  return = 0
  if(sumD>0)
    return = sum( log(Brecords:(Brecords-sumD+1))  )
  return
}

#' Title
#'
#' @param LLL .
#' @param LLA .
#' @param LLB .
#' @param D .
#' @param links .
#' @param sumRowD .
#' @param sumColD .
#' @param gamma .
#'
#' @return .
#' @export
#'
#' @examples .
loglik = function(LLL, LLA, LLB, D, links, sumRowD, sumColD, gamma)
{
  # Sanity check for using logPossibleConfig(.) below
  if(ncol(D) - sum(D) + 1 <= 0){
    cat("\nThe number of records in B has to be >= to the number of linked records.\nNumber of records in B:", ncol(D), "\nNumber of linked records:", sum(D), "\nThe problem may come from the fact that file A is bigger than file B, which should not happen.\nNumber of records in A:", nrow(D), "\n")
  }
  logPossD = sum(log(gamma) * sumRowD + log(1-gamma) * (1-sumRowD)) - logPossibleConfig(ncol(D),sum(D))
  logPossD + sum(LLA[sumRowD==0]) + sum(LLB[sumColD==0]) + sum(LLL[links])
}

#' Title
#'
#' @param data .
#' @param D .
#' @param links .
#' @param survivalpSameH .
#' @param sumRowD .
#' @param sumColD .
#' @param eta .
#' @param phi .
#'
#' @return .
#' @export
#'
#' @examples .
simulateH = function(data, D, links, survivalpSameH, sumRowD, sumColD, eta, phi)
{
  nonlinkedA = sumRowD==0
  nonlinkedB = sumColD==0
  truePIVs = sampleH(nA=dim(data$A[,PIVs]), nB=dim(data$B[,PIVs]), links=links, survivalpSameH=as.matrix(survivalpSameH), pivs_stable=PIVs_stable, pivsA=data$A[,PIVs], pivsB=data$B[,PIVs], nvalues=data$Nvalues, D=D, nonlinkedA=nonlinkedA, nonlinkedB=nonlinkedB, eta=eta, phi=phi)
  list(truepivsA=truePIVs$truepivsA, truepivsB=truePIVs$truepivsB)
}

#' Title
#'
#' @param data .
#' @param D .
#' @param linksR .
#' @param sumRowD .
#' @param sumColD .
#' @param truepivsA .
#' @param truepivsB .
#' @param gamma .
#' @param eta .
#' @param alpha .
#' @param phi .
#'
#' @return .
#' @export
#'
#' @examples .
simulateD = function(data, D, linksR, sumRowD, sumColD, truepivsA, truepivsB, gamma, eta, alpha, phi)
{
  #Determine which observation pairs can be matches (or not) based on the true values of the PIVS
  UA = sspaste2(as.matrix(truepivsA[,PIVs_stable]))
  UB = sspaste2(as.matrix(truepivsB[,PIVs_stable]))
  tmpUA = UA
  tmpUB = UB
  valuesU = unique(c(UA,UB))
  UA = as.numeric(factor(UA,levels=valuesU))
  UB = as.numeric(factor(UB,levels=valuesU))
  tmpA = F2(UA, length(valuesU))
  tmpB = F2(UB, length(valuesU))
  select = F33(tmpA, tmpB, length(tmpA))
  pLink = rep(gamma, nrow(data$A[,PIVs]))
  # What should we add/substract from the loglikelihood if an observation is a not linked
  LLA = rep(0,nrow(data$A[,PIVs]))
  for(k in 1:length(data$Nvalues))
  {
    logpTrue = log(eta[[k]])[truepivsA[,k]]
    pMissingA = phi[[k]][3]
    pTypoA = (1-pMissingA) * (1-phi[[k]][1]) / (data$Nvalues[k]-1)
    pAgreeA = (1-pMissingA) * phi[[k]][1]
    # Contribution to the likelihood
    contr = rep(pAgreeA, nrow(data$A[,PIVs]))
    contr[data$A[,PIVs][,k] != truepivsA[,k]] = pTypoA
    contr[data$A[,PIVs][,k] == 0] = pMissingA
    LLA = LLA + logpTrue + log(contr)
  }
  # What should we add/substract from the loglikelihood if an observation is a not linked
  LLB = rep(0,nrow(data$B[,PIVs]))
  for(k in 1:length(data$Nvalues))
  {
    logpTrue = log(eta[[k]])[truepivsB[,k]]
    pMissingB = phi[[k]][4]
    pTypoB = (1-pMissingB) * (1-phi[[k]][2]) / (data$Nvalues[k]-1)
    pAgreeB = (1-pMissingB) * phi[[k]][2]
    # Contribution to the likelihood
    contr = rep(pAgreeB, nrow(data$B[,PIVs]))
    contr[data$B[,PIVs][,k] != truepivsB[,k]] = pTypoB
    contr[data$B[,PIVs][,k] == 0] = pMissingB
    LLB = LLB + logpTrue + log(contr)
  }
  # What do we add/substract from the loglikelihood if a pair is linked
  LLL = Matrix::Matrix(0, nrow=nrow(data$A[,PIVs]), ncol=nrow(data$B[,PIVs]), sparse=TRUE)
  for(k in 1:length(data$Nvalues))
  {
    HA = truepivsA[ select[,1],k ]
    HB = truepivsB[ select[,2],k ]
    logpTrue = log(eta[[k]])[HA]
    pMissingA = phi[[k]][3]
    pTypoA = (1-pMissingA) * (1-phi[[k]][1]) / (data$Nvalues[k]-1)
    pAgreeA = (1-pMissingA) * phi[[k]][1]
    pMissingB = phi[[k]][4]
    pTypoB = (1-pMissingB) * (1-phi[[k]][2]) / (data$Nvalues[k]-1)
    pAgreeB = (1-pMissingB) * phi[[k]][2]
    # Contribution to the likelihood of linked observation from A
    helpA = rep(pAgreeA, length(HA))
    helpA[data$A[,PIVs][select[,1],k] != HA] = pTypoA
    helpA[data$A[,PIVs][select[,1],k] == 0] = pMissingA
    # Contribution to the likelihood of linked observation from B
    helpB = rep(pAgreeB, length(HB))
    helpB[data$B[,PIVs][select[,2],k] != HB] = pTypoB
    helpB[data$B[,PIVs][select[,2],k] == 0] = pMissingB
    LLL[select] = LLL[select] + logpTrue + log(helpA) + log(helpB)
    # Add unstable part if unstable
    if(!PIVs_stable[k])
    {
      times = abs(data$B[data$B$source!="synthetic",][select[,2], "date"] - data$A[data$A$source!="synthetic",][select[,1], "date"])
      intercept = rep(1, nrow(select))
      cov_k = cbind( intercept,
                     data$A[select[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                     data$B[select[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
      pSameH = SurvivalUnstable(cov_k, alpha[[k]], times)
      helpH = pSameH^(HA==HB) * ((1-pSameH)/(data$Nvalues[k]-1))^(HA!=HB)
      LLL[select] = LLL[select] + log(helpH)
    }
  }
  LLL[select][is.na(LLL[select])] = Inf
  # Complete data likelihood
  LL0 = loglik(LLL=LLL, LLA=LLA, LLB=LLB, D=D, links=linksR, sumRowD=sumRowD, sumColD=sumColD, gamma=gamma)
  # Single run through D
  Dsample = sampleD(S=as.matrix(select),
                    LLA=LLA,
                    LLB=LLB,
                    LLL=LLL,
                    gamma=pLink,
                    loglik=LL0,
                    D=D,
                    nlinkrec=as.integer(sum(D)),
                    sumRowD=sumRowD>0,
                    sumColD=sumColD>0)
  linksR = Dsample$links+1
  # Sanity check: does it give the same likelihood?
  if (round(Dsample$loglik, digits = 3) != round(loglik( LLL=LLL, LLA=LLA, LLB=LLB, D=Dsample$D, links=linksR, sumRowD=Dsample$sumRowD, sumColD=Dsample$sumColD, gamma=pLink), digits = 3))
  {
    cat("\nSanity check failed.\nLog likelihood associated with new Delta:", Dsample$loglik, "\nLog likelihood computed on the new Delta:", loglik( LLL=LLL, LLA=LLA, LLB=LLB, D=Dsample$D, links=linksR, sumRowD=Dsample$sumRowD, sumColD=Dsample$sumColD, gamma=pLink), "\n")
  }
  Dsample
}

#' Title
#'
#' @param alphas A vector of size 1+cov in A+cov in B with coefficients of the hazard (baseline hazard and conditional hazard)
#' @param X A matrix with number of linked records rows and 1+cov in A+cov in B columns (first column: intercept, following columns: covariates from A and then from B to model instability) (used for optimisation: X concatenate the X obtained in each iteration of the Gibbs sampler)
#' @param times A vector of size number of linked records with the time gaps between the record from each sources (used for optimisation: times concatenate the times vectors obtained in each iteration of the Gibbs sampler)
#' @param Hequal A vector of size number of linked records with boolean values indicating wether the values in A and in B coincide (used for optimisation: times concatenate the times vectors obtained in each iteration of the Gibbs sampler)
#'
#' @return The value of the opposite (-) of the log likelihood associated with the survival function defining the probabilities that true values coincide (as defined in the paper) (the algorithm minimises -loglik i.e. maximises the loglik)
#' @export
#'
#' @examples .
loglikSurvival = function(alphas, X, times, Hequal)
{
  loglikelihood = sum( - exp( as.matrix(X) %*% alphas ) * times + (!Hequal) * log( exp( exp( as.matrix(X) %*% alphas ) * times ) - 1 ) )
  - loglikelihood
}

#' Title
#'
#' @param Xlinksk A matrix with number of linked records rows and 1+cov in A+cov in B columns, with a first column filled with 1 (intercept), and following columns filled with the values of the covariates useful for modelling instability for the linked records
#' @param alphask A vector of size 1+cov in A+cov in B, with as first element the baseline hazard and following elements being the coefficient of the conditional hazard associated with the covariates given in X
#' @param times A vector of size number of linked records with the time gaps between the record from each sources
#'
#' @return A vector (for an unstable PIV) of size number of linked records with the probabilities that true values coincide (e.g. 1 - proba to move if the PIV is postal code) defined according to the survival function with exponential model proposed in the paper
#' @export
#'
#' @examples .
SurvivalUnstable = function(Xlinksk, alphask, times)
{
  # returns the proba that true values of PIV indexed by k coincide
  # P(HAik = HBjk | t_ij) = exp( - exp(X.ALPHA) t_ij )
  # first column in Xlinks is the intercept
  # following columns in Xlinks are the covariates to include
  # number of rows in Xlinks is the number of linked records
  exp( - exp( as.matrix(Xlinksk) %*% alphask ) * times )
}

#' Title
#'
#' @param data A list with elements A, B, Nvalues, PIVs_config, controlOnMistakes,
#' sameMistakes, phiMistakesAFixed, phiMistakesBFixed, phiForMistakesA, phiForMistakesB
#' @param StEMIter An integer with the total number of iterations of the Stochastic
#' EM algorithm (including the period to discard as burn-in)
#' @param StEMBurnin An integer with the number of iterations to discard as burn-in
#' @param GibbsIter An integer with the total number of iterations of the Gibbs sampler
#' (done in each iteration of the StEM) (including the period to discard as burn-in)
#' @param GibbsBurnin An integer with the number of iterations to discard as burn-in
#' @param sparseOutput A boolean value, to handle the final linkage matrix as sparse or
#' not, for large data sets it should be set to TRUE to avoid R crashing or memory errors
#' @param musicOn A boolean value, if TRUE the algorithm will play music at the end of
#' the algorithm, useful if you have to wait for the record linkage to run and to act as
#' an alarm when record linkage is done
#' @param newDirectory A NULL value or: A string with the name of (or path to) the directory
#' (which should already exist) where to save the environment variables at the end of each
#' iteration (useful when record linkage is very long, to not loose everything and not restart
#' from scratch in case your computer shut downs before record linkage is finished)
#' @param saveInfoIter A boolean value to indicate whether you want the environment variables
#' to be saved at the end of each iteration (useful when record linkage is very long, to not
#' loose everything and not restart from scratch in case your computer shut downs before
#' record linkage is finished)
#'
#' @return A list with: Delta, (sparse) matrix with the pairs of records linked and their posterior probabilities to be linked (select the pairs where the proba>0.5 to get a valid set of linked records), gamma, a vector with the chain of the parameter gamma representing the proportion of linked records as a fraction of the smallest file, eta, a vector with the chain of the parameter eta representing the distribution of the PIVs, alpha, a vector with the chain of the parameter alpha representing the hazard coefficient of the model for instability, phi, a vector with the chain of the parameter phi representing the registration errors parameters)
#' @export
#'
#' @examples data = list( A=encodedA,
#'                        B=encodedB,
#'                        Nvalues=nvalues,
#'                        PIVs_config=PIVs_config,
#'                        controlOnMistakes=c(TRUE,TRUE,FALSE,FALSE,FALSE),
#'                        sameMistakes=TRUE,
#'                        phiMistakesAFixed=FALSE,
#'                        phiMistakesBFixed=FALSE,
#'                        phiForMistakesA=c(NA,NA,NA,NA,NA),
#'                        phiForMistakesB=c(NA,NA,NA,NA,NA))
#'            fit = StEM( data,
#'                        100,
#'                        70,
#'                        200,
#'                        100,
#'                        TRUE,
#'                        TRUE,
#'                        NULL,
#'                        FALSE )
stEM = function(data, StEMIter, StEMBurnin, GibbsIter, GibbsBurnin, sparseOutput=TRUE, musicOn=TRUE, newDirectory=NULL, saveInfoIter=FALSE)
{

  PIVs = names(data$PIVs_config)
  PIVs_stable = sapply(data$PIVs_config, function(x) x$stable)
  PIVs_conditionalHazard_variables = sapply(data$PIVs_config, function(x) if(!x$stable){x$conditionalHazard}else{FALSE})

  # data$A, data$B, data$PIVs_stable, PIVs are used globally without necessarily being passed as parameters

  # ANY CONTROL ON MISTAKES? i.e. BOUND ON PHI
  if(any(data$controlOnMistakes)){
    controlMistakes = which(data$controlOnMistakes)
  }else{
    controlMistakes = c()
  }

  # ANY UNSTABLE PIV?
  instability = any(!PIVs_stable)
  if(instability){
    unstablePIVs = which(!PIVs_stable)
    if(any(PIVs_conditionalHazard_variables)){
      for(k in 1:length(unstablePIVs)){
        unstablePIV = unstablePIVs[k]
        if( length(data$PIVs_config[[unstablePIV]]$pSameH.cov.A)>0 )
          testit::assert( "Some variables from A to include in the survival model for unstable PIVs do not exist.", data$PIVs_config[[unstablePIV]]$pSameH.cov.A %in% colnames(A) )
        if( length(data$PIVs_config[[unstablePIV]]$pSameH.cov.B)>0 )
          testit::assert( "Some variables from B to include in the survival model for unstable PIVs do not exist.", data$PIVs_config[[unstablePIV]]$pSameH.cov.B %in% colnames(B) )
      }
    }
  }

  nGibbsIter = GibbsIter - GibbsBurnin
  testit::assert( "Number of iterations for StEM should be positive.", StEMIter - StEMBurnin > 0 )
  testit::assert( "Number of iterations for StEM Gibbs sampler should be positive.", nGibbsIter > 0 )

  # Parameters for PIVs

  # Parameter for the probability for a pair to be linked
  gamma = 0.5

  # Parameters for the distributions of the true values
  eta = lapply(data$Nvalues, function(x) rep(1/x,x))

  # Parameters for the survival model describing pSameH for unstable PIVs (over time and potentially more covariates)
  # Number of coefficients: baseline + covariates from A + covariates from B
  nCoefUnstable = lapply(seq_along(PIVs_stable), function(idx) if(PIVs_stable[idx]){ 0 }else{ ncol(data$A[, data$PIVs_config[[idx]]$pSameH.cov.A, drop=FALSE]) + ncol(data$B[, data$PIVs_config[[idx]]$pSameH.cov.B, drop=FALSE]) + 1 } )
  alpha = lapply(seq_along(PIVs_stable), function(idx) if(PIVs_stable[idx]){ c(-Inf) }else{ rep(log(0.05), nCoefUnstable[[idx]]) })

  # Parameters for the registration errors (agreement in A, agreement in B, missing in A, missing in B)
  phi = lapply(data$Nvalues, function(x)  c(0.9,0.9,0.1,0.1))

  NmissingA = lapply(seq_along(data$Nvalues), function(k) sum(data$A[,PIVs][,k]==0))
  NmissingB = lapply(seq_along(data$Nvalues), function(k) sum(data$B[,PIVs][,k]==0))

  gamma.iter = array(NA, c(StEMIter, length(gamma)))
  eta.iter = lapply(data$Nvalues, function(x) array(NA, c(StEMIter, x)))
  phi.iter = lapply(data$Nvalues, function(x) array(NA, c(StEMIter, 4)))
  alpha.iter = lapply(nCoefUnstable, function(x) array(NA, c(StEMIter, x)))

  time.iter=c()
  pb = progress::progress_bar$new(format = "Running StEM algorithm [:bar] :percent in :elapsed",       total = StEMIter, clear = FALSE, width= 60)

  # Stochastic EM
  for(iter in 1:StEMIter)
  {
    tijdM = Sys.time()

    D = Matrix::Matrix(0, nrow=nrow(data$A), ncol=nrow(data$B), sparse=TRUE)
    linksR = which(D==1, arr.ind=TRUE)
    linksCpp = linksR
    sumRowD = rowSums(D)
    sumColD = colSums(D)
    nlinkrec = 0
    survivalpSameH = base::matrix(1, nrow(linksR), length(data$Nvalues))

    # if burnin value is 0, the algorithm will explore the necessary burnin for the number of linked records to stagnate
    countBurnin = 10
    Burnin_total = 0

    # Burn-in period Gibbs sampler
    if(GibbsBurnin == 0){
      pb_burnin_explore = progress::progress_bar$new(format = "(:spin)     Burn-in period: :whatburnin     Number of linked records: :whatlinkrec", total = 300, clear = FALSE, width= 100)
      while(countBurnin != 0)
      {
        Burnin_total = Burnin_total + 1
        newTruePivs = simulateH(data=data, D=D, links=linksCpp, survivalpSameH=survivalpSameH, sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
        truepivsA = newTruePivs$truepivsA
        truepivsB = newTruePivs$truepivsB
        Dsample = simulateD(data=data, D=D, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma, eta=eta, alpha=alpha, phi=phi)

        D = Dsample$D
        linksCpp = Dsample$links
        linksR = linksCpp + 1
        sumRowD = Dsample$sumRowD
        sumColD = Dsample$sumColD
        loglikelihood = Dsample$loglik
        new_nlinkred = Dsample$nlinkrec

        if( abs(new_nlinkred - nlinkrec) < 10 ){
          countBurnin = countBurnin - 1
          nlinkrec = new_nlinkred
        }else{
          countBurnin = 15
          nlinkrec = new_nlinkred
        }

        survivalpSameH = base::matrix(1, nrow(linksR), length(data$Nvalues))
        if(instability){
          if(nrow(linksR)>0){
            times = abs(data$B[data$B$source!="synthetic",][linksR[,2], "date"] - data$A[data$A$source!="synthetic",][linksR[,1], "date"])
            intercept = rep(1, nrow(linksR))
            for(k in 1:length(PIVs)){
              if(k %in% unstablePIVs){
                cov_k = cbind( intercept,
                               data$A[linksR[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                               data$B[linksR[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
                survivalpSameH[,k] = SurvivalUnstable(cov_k, alpha[[k]], times)
              }
            }
          }
        }

        pb_burnin_explore$tick(tokens = list(whatburnin = Burnin_total, whatlinkrec = nlinkrec))
        if(Burnin_total>=300){
          countBurnin = 0
          cat("\nThe burn-in period exceeded 300 iterations, this is the default maximum burn-in in that case.\nIf you want to increase the burn-in period fix it with the parameter 'GibbsBurnin' of the stEM function.\n")
          pb_burnin_explore$terminate()
          GibbsBurnin = 300
          stop()
        }
      }
    }else{
      for(j in 1:GibbsBurnin)
      {
        newTruePivs = simulateH(data=data, D=D, links=linksCpp, survivalpSameH=survivalpSameH, sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
        truepivsA = newTruePivs$truepivsA
        truepivsB = newTruePivs$truepivsB
        Dsample = simulateD(data=data, D=D, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma, eta=eta, alpha=alpha, phi=phi)

        D = Dsample$D
        linksCpp = Dsample$links
        linksR = linksCpp + 1
        sumRowD = Dsample$sumRowD
        sumColD = Dsample$sumColD
        loglikelihood = Dsample$loglik
        nlinkred = Dsample$nlinkrec

        survivalpSameH = base::matrix(1, nrow(linksR), length(data$Nvalues))
        if(instability){
          if(nrow(linksR)>0){
            times = abs(data$B[data$B$source!="synthetic",][linksR[,2], "date"] - data$A[data$A$source!="synthetic",][linksR[,1], "date"])
            intercept = rep(1, nrow(linksR))
            for(k in 1:length(PIVs)){
              if(k %in% unstablePIVs){
                cov_k = cbind( intercept,
                               data$A[linksR[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                               data$B[linksR[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
                survivalpSameH[,k] = SurvivalUnstable(cov_k, alpha[[k]], times)
              }
            }
          }
        }
      }
    }

    # Administration for the M-step
    Vgamma = c()
    Veta = lapply(data$Nvalues, function(x) c())
    Valpha = mapply(createDataAlpha, nCoefUnstable = nCoefUnstable, stable = PIVs_stable, SIMPLIFY=FALSE)
    Vphi = lapply(data$Nvalues, function(x) c())

    for(j in 1:nGibbsIter)
    {
      newTruePivs = simulateH(data=data, D=D, links=linksCpp, survivalpSameH=survivalpSameH, sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
      truepivsA = newTruePivs$truepivsA
      truepivsB = newTruePivs$truepivsB
      Dsample = simulateD(data=data, D=D, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma, eta=eta, alpha=alpha, phi=phi)

      D = Dsample$D
      linksCpp = Dsample$links
      linksR = linksCpp + 1
      sumRowD = Dsample$sumRowD
      sumColD = Dsample$sumColD
      loglikelihood = Dsample$loglik
      nlinkrec = Dsample$nlinkrec

      survivalpSameH = base::matrix(1, nrow(linksR), length(data$Nvalues))
      if(instability){
        if(nrow(linksR)>0){
          times = abs(data$B[data$B$source!="synthetic",][linksR[,2], "date"] - data$A[data$A$source!="synthetic",][linksR[,1], "date"])
          intercept = rep(1, nrow(linksR))
          for(k in 1:length(PIVs)){
            if(k %in% unstablePIVs){
              cov_k = cbind( intercept,
                             data$A[linksR[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                             data$B[linksR[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
              survivalpSameH[,k] = SurvivalUnstable(cov_k, alpha[[k]], times)
              Hequal = truepivsA[linksR[,1],k] == truepivsB[linksR[,2],k]
              Vtmp_k = cbind( cov_k,
                              times,
                              Hequal )
              Valpha[[k]] = rbind(Valpha[[k]], Vtmp_k)
            }
          }
        }
      }

      # Update gamma
      Vgamma[j] = sum(D)
      if(sum(D)==0){
        cat("\nloglikelihood = ", loglikelihood, "\nSo far, gamma = ", gamma, "\nIn the last iteration no link has been made.\nThis is probably due to a difference in the support of some PIV between file A and file B.\n")
      }

      # Update eta
      for(k in 1:length(data$Nvalues))
      {
        facpivsA = factor(truepivsA[,k], levels=1:data$Nvalues[k])
        facpivsB = factor(truepivsB[,k], levels=1:data$Nvalues[k])
        Veta[[k]] = rbind(Veta[[k]], table(facpivsA[sumRowD==0]) + table(facpivsB[sumColD==0]) + table(facpivsA[sumRowD==1]))
      }

      # Update phi: count agreements / missings
      for(k in 1:length(data$Nvalues)){
        agreements_in_A = sum(truepivsA[,k]==data$A[,PIVs][,k])
        agreements_in_B = sum(truepivsB[,k]==data$B[,PIVs][,k])
        Vphi[[k]] = rbind( Vphi[[k]],
                           c( agreements_in_A,
                              agreements_in_B,
                              NmissingA[[k]],
                              NmissingB[[k]]) )
      }
    }

    # Calculate new parameters gamma/eta/phi/alpha

    # New gamma
    gamma = sum(Vgamma) / (nGibbsIter*nrow(truepivsA))

    # New eta
    for(k in 1:length(data$Nvalues)){
      eta[[k]] = colSums(Veta[[k]])/sum(Veta[[k]])
      if(sum(eta[[k]]==0)>0){
        cat("\nAn error may appear due to some rare values in PIV", k, "\nSo far, eta ", k, ":\n", eta[[k]], "\nCount values of PIV", k, "in A: \n", table(encodedA[,PIVs][,k]), "\nand in B: \n", table(encodedB[,PIVs][,k]), "\n")
      }
    }

    # New alpha
    if(instability){
      if(nrow(linksR)>0){
        for(k in 1:length(PIVs)){
          if(!PIVs_stable[k]){
            alphaInit = rep(-0.05, nCoefUnstable[[k]])
            X = Valpha[[k]][,1:nCoefUnstable[[k]]]
            times = Valpha[[k]]$times
            Hequal = Valpha[[k]]$Hequal
            optim = stats::nlminb(alphaInit, loglikSurvival, control=list(trace=FALSE), X=X, times=times, Hequal=Hequal)
            alpha[[k]] = optim$par
          }
        }
      }
    }

    # New phi
    for(k in 1:length(data$Nvalues))
    {
      NtotalA = nGibbsIter * nrow(data$A)
      NtotalB = nGibbsIter * nrow(data$B)
      Ntotal  = nGibbsIter * (nrow(data$A) + nrow(data$B))
      if(k %in% controlMistakes){
        phi[[k]][1] = max( 0.9, sum(Vphi[[k]][,1]) / (NtotalA - sum(Vphi[[k]][,3])) )
        phi[[k]][2] = max( 0.9, sum(Vphi[[k]][,2]) / (NtotalB - sum(Vphi[[k]][,4])) )
      }else{
        phi[[k]][1] = sum(Vphi[[k]][,1]) / (NtotalA - sum(Vphi[[k]][,3]))
        phi[[k]][2] = sum(Vphi[[k]][,2]) / (NtotalB - sum(Vphi[[k]][,4]))
      }
      if(data$sameMistakes){
        phi[[k]][1] = (sum(Vphi[[k]][,1]) + sum(Vphi[[k]][,2])) / (Ntotal - sum(Vphi[[k]][,3]) - sum(Vphi[[k]][,4]))
        phi[[k]][2] = (sum(Vphi[[k]][,1]) + sum(Vphi[[k]][,2])) / (Ntotal - sum(Vphi[[k]][,3]) - sum(Vphi[[k]][,4]))
      }
      if(data$phiMistakesAFixed | data$phiMistakesBFixed){
        for(i in 1:length(unstablePIVs)){
          if(data$phiMistakesAFixed){
            phi[[ unstablePIVs[i] ]][1] = 1-data$phiForMistakesA[[ unstablePIVs[i] ]]
          }
          if(data$phiMistakesBFixed){
            phi[[ unstablePIVs[i] ]][2] = 1-data$phiForMistakesB[[ unstablePIVs[i] ]]
          }
        }
      }
      phi[[k]][3] = sum(Vphi[[k]][,3]) / (nGibbsIter*nrow(data$A))
      phi[[k]][4] = sum(Vphi[[k]][,4]) / (nGibbsIter*nrow(data$B))
    }

    # Administration

    gamma.iter[iter,] = gamma

    for(k in 1:length(data$Nvalues))
    {
      eta.iter[[k]][iter,] = eta[[k]]
      alpha.iter[[k]][iter,] = alpha[[k]]
      phi.iter[[k]][iter,] = phi[[k]]
    }

    pb$tick()

    # save current environment, global variables and local variables
    if(!is.null(newDirectory) & saveInfoIter)
    {
      save.image(file=file.path(newDirectory, 'myEnvironment.RData'))
      save(iter, gamma, eta, alpha, phi, gamma.iter, eta.iter, phi.iter, alpha.iter, file=file.path(newDirectory, 'myEnvironmentLocal.RData'))
    }
  }

  pb$terminate()

  Delta = Matrix::Matrix(0, nrow=nrow(data$A), ncol=nrow(data$B), sparse=TRUE)

  gamma_avg = apply(gamma.iter, 2, function(x) mean(x[StEMBurnin:StEMIter , drop=FALSE]))
  eta_avg = lapply(eta.iter, function(x) apply(x[StEMBurnin:StEMIter, , drop=FALSE], 2, mean))
  alpha_avg = lapply(alpha.iter, function(x) apply(x[StEMBurnin:StEMIter, , drop=FALSE], 2, mean))
  phi_avg = lapply(phi.iter, function(x) apply(x[StEMBurnin:StEMIter, , drop=FALSE], 2, mean))

  pbfinal = progress::progress_bar$new(format = "Drawing Delta          [:bar] :percent in :elapsed",       total = 1000, clear = FALSE, width= 60)

  for(m in 1:1000)
  {
    newTruePivs = simulateH(data=data, D=D, links=linksCpp, survivalpSameH=survivalpSameH, sumRowD=sumRowD, sumColD=sumColD, eta=eta_avg, phi=phi_avg)
    truepivsA = newTruePivs$truepivsA
    truepivsB = newTruePivs$truepivsB
    Dsample = simulateD(data=data, D=D, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma_avg, eta=eta_avg, alpha=alpha_avg, phi=phi_avg)

    D = Dsample$D
    linksCpp = Dsample$links
    linksR = linksCpp+1
    sumRowD = Dsample$sumRowD
    sumColD = Dsample$sumColD
    loglikelihood = Dsample$loglik
    nlinkrec = Dsample$nlinkrec

    survivalpSameH = base::matrix(1, nrow(linksR), length(data$Nvalues))
    if(instability){
      if(nrow(linksR)>0){
        times = abs(data$B[data$B$source!="synthetic",][linksR[,2], "date"] - data$A[data$A$source!="synthetic",][linksR[,1], "date"])
        intercept = rep(1, nrow(linksR))
        for(k in 1:length(PIVs)){
          if(k %in% unstablePIVs){
            cov_k = cbind( intercept,
                           data$A[linksR[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                           data$B[linksR[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
            survivalpSameH[,k] = SurvivalUnstable(cov_k, alpha[[k]], times)
          }
        }
      }
    }

    Delta = Delta + D
    pbfinal$tick()
  }
  Delta = Delta / 1000

  if(sparseOutput){
    Delta = summary(Delta)
  }

  pbfinal$terminate()

  if(musicOn){
    utils::browseURL('https://www.youtube.com/watch?v=NTa6Xbzfq1U')
  }

  list(Delta=Delta, gamma=gamma.iter, eta=eta.iter, alpha=alpha.iter, phi=phi.iter)
}

#' Title
#'
#' @param PIVs .
#' @param encodedA .
#' @param encodedB .
#'
#' @return .
#' @export
#'
#' @examples .
launchNaive <- function(PIVs, encodedA, encodedB){
  testit::assert( "NA in A should be encoded as 0.", sum(is.na(encodedA))==0 )
  testit::assert( "NA in B should be encoded as 0.", sum(is.na(encodedB))==0 )
  rownames(encodedA) = 1:nrow(encodedA)
  rownames(encodedB) = 1:nrow(encodedB)
  DeltaNaiveLinked = data.frame( Matrix::Matrix(0, nrow=0, ncol=2, sparse=TRUE) )
  colnames(DeltaNaiveLinked) = c("idxA", "idxB")
  # A NOT MISSING THAT MATCH WITH B
  isNotMissingA = apply(encodedA[,PIVs]!=0, 1, all)
  A_PIVs_notMissing = encodedA[isNotMissingA,PIVs]
  A_PIVs_notMissing_ID = rownames(encodedA[isNotMissingA,])
  UA = sspaste2(as.matrix(A_PIVs_notMissing))
  UB = sspaste2(as.matrix(encodedB[,PIVs]))
  valuesU = unique(c(UA,UB))
  UA = as.numeric(factor(UA,levels=valuesU))
  UB = as.numeric(factor(UB,levels=valuesU))
  tmpA = F2(UA, length(valuesU))
  tmpB = F2(UB, length(valuesU))
  select = F33(tmpA, tmpB, length(tmpA))
  if(nrow(select>0)){
    for (l in 1:nrow(select))
    {
      idxA = as.integer(select[l,1])
      idxA = A_PIVs_notMissing_ID[idxA]
      idxB = as.integer(select[l,2])
      idxB = rownames(encodedB[idxB,])
      DeltaNaiveLinked = rbind(DeltaNaiveLinked, cbind(idxA=idxA, idxB=idxB))
    }
  }
  # A MISSING THAT MATCH WITH B
  for(k in 1:length(PIVs)){
    isMissingA_k = encodedA[,k]==0
    A_PIVs_k_Missing = encodedA[isMissingA_k,PIVs]
    A_PIVs_k_Missing_ID = rownames(encodedA[isMissingA_k,])
    UA = sspaste2(as.matrix(A_PIVs_k_Missing[,-k]))
    UB = sspaste2(as.matrix(encodedB[,PIVs][,-k]))
    valuesU = unique(c(UA,UB))
    UA = as.numeric(factor(UA,levels=valuesU))
    UB = as.numeric(factor(UB,levels=valuesU))
    tmpA = F2(UA, length(valuesU))
    tmpB = F2(UB, length(valuesU))
    select = F33(tmpA, tmpB, length(tmpA))
    if(nrow(select>0)){
      for (l in 1:nrow(select))
      {
        idxA = as.integer(select[l,1])
        idxA = A_PIVs_k_Missing_ID[idxA]
        idxB = as.integer(select[l,2])
        idxB = rownames(encodedB[idxB,])
        DeltaNaiveLinked = rbind(DeltaNaiveLinked, cbind(idxA=idxA, idxB=idxB))
      }
    }
  }
  # B NOT MISSING THAT MATCH WITH A
  isNotMissingB = apply(encodedB[,PIVs]!=0, 1, all)
  B_PIVs_notMissing = encodedB[isNotMissingB,PIVs]
  B_PIVs_notMissing_ID = rownames(encodedB[isNotMissingB,])
  UB = sspaste2(as.matrix(B_PIVs_notMissing))
  UA = sspaste2(as.matrix(encodedA[,PIVs]))
  valuesU = unique(c(UA,UB))
  UA = as.numeric(factor(UA,levels=valuesU))
  UB = as.numeric(factor(UB,levels=valuesU))
  tmpA = F2(UA, length(valuesU))
  tmpB = F2(UB, length(valuesU))
  select = F33(tmpA, tmpB, length(tmpA))
  if(nrow(select>0)){
    for (l in 1:nrow(select))
    {
      idxA = as.integer(select[l,1])
      idxA = rownames(encodedA[idxA,])
      idxB = as.integer(select[l,2])
      idxB = B_PIVs_notMissing_ID[idxB]
      DeltaNaiveLinked = rbind(DeltaNaiveLinked, cbind(idxA=idxA, idxB=idxB))
    }
  }
  # B MISSING THAT MATCH WITH A
  for(k in 1:length(PIVs)){
    isMissingB_k = encodedB[,k]==0
    B_PIVs_k_Missing = encodedB[isMissingB_k,PIVs]
    B_PIVs_k_Missing_ID = rownames(encodedB[isMissingB_k,])
    UB = sspaste2(as.matrix(B_PIVs_k_Missing[,-k]))
    UA = sspaste2(as.matrix(encodedA[,PIVs][,-k]))
    valuesU = unique(c(UA,UB))
    UA = as.numeric(factor(UA,levels=valuesU))
    UB = as.numeric(factor(UB,levels=valuesU))
    tmpA = F2(UA, length(valuesU))
    tmpB = F2(UB, length(valuesU))
    select = F33(tmpA, tmpB, length(tmpA))
    if(nrow(select>0)){
      for (l in 1:nrow(select))
      {
        idxA = as.integer(select[l,1])
        idxA = rownames(encodedA[idxA,])
        idxB = as.integer(select[l,2])
        idxB = B_PIVs_k_Missing_ID[idxB]
        DeltaNaiveLinked = rbind(DeltaNaiveLinked, cbind(idxA=idxA, idxB=idxB))
      }
    }
  }
  DeltaNaiveLinked = DeltaNaiveLinked[!duplicated(DeltaNaiveLinked), ]
  DeltaNaiveLinked
}

