##############################################
#
# General description:
#
#   The following function fits different models for the selection coefficient
#   frequency of an allele under selection (equation after Boissinot at al. 2006 PNAS)

# Input:
#
#     PredictMat: matrix of predictor variables for selection coefficient
#                 (1 row per LINE-1)
#     s: selection coefficient
#     N: population size

# Comment:
#     This function requires the function AlleleFreqTime

##############################################

FitSelectionModels_quadgr <- function(PredictMat,  Freqs, 
                               Counts, PopSize, SampleSize,
                               blnIns, 
                               LogRegCoeff,
                               DetectProb,
                               aBorder = 0.003, 
                               bBorder = 10^(-6), 
                               cBorder = 10^(-6)

){
    

  # Estimate maximum likelihood for a single selection coefficient
  cat("Estimate maximum likelihood for a single selection coefficient\n")
  ML_a <-  constrOptim(theta = c(a = 0),
                          f = function(x) -AlleleFreqLogLik_4Par_quadgr(
                            Freqs = Freqs,
                            Counts = rep(1, length(Freqs)),
                            Predict = PredictMat,
                            a = x[1], b = 0, c = 0, d = 0, N = PopSize,
                            SampleSize = SampleSize,
                            blnIns = blnIns, 
                            LogRegCoeff = LogRegCoeff,
                            DetectProb = DetectProb),
                          grad = NULL,
                          ui = rbind(1, -1),
                          ci = c(a = -aBorder, a = -aBorder),
                          method = "Nelder-Mead")
  cat("done!\n")
  
  
  # Get maximum likelihood estimate for effect of L1 start on selection
  cat("Estimate effect of", colnames(PredictMat)[1], "on selection ...\n")
  ML_ab <-  constrOptim(theta = c(a = ML_a$par, b = 0),
                             f = function(x) -AlleleFreqLogLik_4Par_quadgr(
                               Freqs = Freqs,
                               Counts = rep(1, length(Freqs)),
                               Predict = PredictMat,
                               a = x[1], b = x[2], c = 0, d = 0, N = PopSize, 
                               SampleSize = SampleSize,
                               blnIns = blnIns, 
                               LogRegCoeff = LogRegCoeff,
                               DetectProb = DetectProb),
                             grad = NULL,
                             ui = rbind(c(1, 0),  c(0, 1),   
                                        c(-1, 0), c(0, -1)),
                             ci = c(a = -aBorder, b = -bBorder, 
                                    a = -aBorder, b = -bBorder),
                             method = "Nelder-Mead")
  cat("done!\n")
  
  # Get maximum likelihood estimate for effect of full-length L1 on selection
  cat("Estimate effect of", colnames(PredictMat)[2], "on selection ...\n")
  ML_ac <-  constrOptim(theta = c(a = ML_a$par, c = 0),
                            f = function(x) -AlleleFreqLogLik_4Par_quadgr(
                              Freqs = Freqs,
                              Counts = rep(1, length(Freqs)),
                              Predict = PredictMat,
                              a = x[1], b = 0, c = x[2], d = 0, N = PopSize, 
                              SampleSize = SampleSize,
                              blnIns = blnIns, 
                              LogRegCoeff = LogRegCoeff,
                              DetectProb = DetectProb),
                            grad = NULL,
                            ui = rbind(c(1, 0),  c(0, 1),   
                                       c(-1, 0), c(0, -1)),
                            ci = c(a = -aBorder, c = -cBorder, 
                                   a = -aBorder, c = -cBorder),
                            method = "Nelder-Mead")
  cat("done!\n")
  
  # Determine maximum likelihood with 3 parameters (selection coefficient as 
  # function of L1 start and indicator for full-length)
  cat("Estimate effect of", paste(colnames(PredictMat)[1:2], sep = " and "), 
      "on selection ...\n")
  ML_abc <- constrOptim(theta = c(a = ML_ab$par[1], 
                                            b = ML_ab$par[2], 
                                            c = ML_ac$par[2]),
                                  f = function(x) -AlleleFreqLogLik_4Par_quadgr(
                                    Freqs = Freqs,
                                    Counts = rep(1, length(Freqs)),
                                    Predict = PredictMat,
                                    a = x[1], b = x[2], c = x[3], d = 0, N = PopSize, 
                                    SampleSize = SampleSize,
                                    blnIns = blnIns, 
                                    LogRegCoeff = LogRegCoeff,
                                    DetectProb = DetectProb),
                                  grad = NULL,
                                  ui = rbind(c(1, 0, 0),  c(0, 1, 0),  c(0, 0, 1), 
                                             c(-1, 0, 0), c(0, -1, 0), c(0, 0, -1)),
                                  ci = c(a = -aBorder, b = -bBorder, c = -cBorder, 
                                         a = -aBorder, b = -bBorder, c = -cBorder),
                                  method = "Nelder-Mead")
  cat("done!\n")
  
  ###################################################
  #                                                 #
  #  Summarize results                              #
  #                                                 #
  ###################################################
  
  # Function to extract AIC from optim results
  GetAIC <- function(OptimResults){
    round(2 * (length(OptimResults$par) + OptimResults$value), 2)
  }
  GetParVals <- function(OptimResults){
    Results <- paste(names(OptimResults$par), 
                     format(OptimResults$par, digits = 2), sep = " = ",
                     collapse = ", ")
  }
  GetNPar <- function(OptimResults){
    length(OptimResults$par)
  }
  
  # Get columns of AIC and parameter values
  Cols2Append <- t(sapply(list(ML_a, 
                               ML_ab, 
                               ML_ac, 
                               ML_abc), function(x){
                                 c(AIC = GetAIC(x), Pars = GetParVals(x))
                               }))

    # Combine AIC values into a table
  AICTab <- cbind(data.frame(
         NrParameters = c(1, 2, 2, 3),
         Predictor = c("none", 
                       colnames(PredictMat)[1], 
                       colnames(PredictMat)[2], 
                       paste(colnames(PredictMat)[1:2], collapse = " and ")),
         stringsAsFactors = F),
         Cols2Append)

  # Put all results in a list     
  list(ML_a = ML_a, ML_ab = ML_ab, ML_ac = ML_ac, ML_abc = ML_abc, 
       AICTab = AICTab)
}


