% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FlexRecordLinkage.r
\name{SurvivalUnstable}
\alias{SurvivalUnstable}
\title{The survival function with exponential model}
\usage{
SurvivalUnstable(Xlinksk, alphask, times)
}
\arguments{
\item{Xlinksk}{A matrix with number of linked records rows and 1+cov in A+cov in B columns, with a first column filled with 1 (intercept), and following columns filled with the values of the covariates useful for modelling instability for the linked records}

\item{alphask}{A vector of size 1+cov in A+cov in B, with as first element the baseline hazard and following elements being the coefficient of the conditional hazard associated with the covariates given in X}

\item{times}{A vector of size number of linked records with the time gaps between the record from each sources}
}
\value{
A vector (for an unstable PIV) of size number of linked records with the probabilities that true values coincide (e.g. 1 - proba to move if the PIV is postal code) defined according to the survival function with exponential model proposed in the paper
}
\description{
The survival function with exponential model (as proposed in our paper) is representing the probability that true values of a pair of records referring to the same entity coincide. For a linked pair of record (i,j) from sources A, B respectively, P(HAik = HBjk | t_ij, ...) = exp( - exp(X.ALPHA) t_ij ).
This function is only used if the PIV is unstable and evolve over time. If so the true values of a linked pair of records may not coincide.
If you want to use a different survival function to model instability, you can change this function 'SurvivalUnstable' as well as the associated log(likelihood) function 'loglikSurvival'.
Also see ?FlexRL::loglikSurvival.
}
\details{
The simplest model (without covariates) just writes P(HAik = HBjk | t_ij, ...) = exp( - exp(alpha) . t_ij )
The more complex model (with covariates) writes P(HAik = HBjk | t_ij, ...) = exp( - exp(X.ALPHA) . t_ij ) and uses a matrix X (nrow=nbr of linked records, ncol=1 + nbr of cov from A + nbr of cov from B) where the first column is filled with 1 (intercept) and the subsequent columns are the covariates values from source A and/or from source B to be used. The ALPHA in this case is a vector of parameters, the first one being associated with the intercept is the same one than for the simplest model, the subsequent ones are associated with the covariates from A and/or from B.
}
\examples{
\dontrun{
for(iter in 1:StEMIter)
{

  for(j in 1:GibbsBurnin)
  {

    ...

  }

  for(j in 1:nGibbsIter)
  {
    ...

    survivalpSameH = base::matrix(1, nrow(linksR), length(data$Nvalues))
    if(instability){
      if(nrow(linksR)>0){
        times = abs(data$B[data$B$source!="synthetic",][linksR[,2], "date"] -
                    data$A[data$A$source!="synthetic",][linksR[,1], "date"])
        intercept = rep(1, nrow(linksR))
        for(k in 1:length(PIVs)){
          if(k \%in\% unstablePIVs){
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

    ...

  }

  ...

  # Update the alpha chain wth a new parameter:
  if(instability){
    if(nrow(linksR)>0){
      for(k in 1:length(PIVs)){
        if(!PIVs_stable[k]){
          alphaInit = rep(-0.05, nCoefUnstable[[k]])
          X = Valpha[[k]][,1:nCoefUnstable[[k]]]
          times = Valpha[[k]]$times
          Hequal = Valpha[[k]]$Hequal
          optim = stats::nlminb(alphaInit, loglikSurvival, control=list(trace=FALSE),
                                X=X, times=times, Hequal=Hequal)
          alpha[[k]] = optim$par
        }
      }
    }
  }

  ...

}
}
}
