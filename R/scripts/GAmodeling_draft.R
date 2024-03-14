predictCentilesForAgeRange_GAbins <- function(gamModel, ageRange, euler=0, cent=0.5, by.sex = TRUE, df = NULL){
  if (euler == 0){
    newDataM_VPM <- data.frame(logAge=ageRange,
                               sex=c(rep(as.factor("M"),  length(ageRange))), GAbins_recode=c(rep(as.factor("VPM"),  length(ageRange))))
    newDataF_VPM <- data.frame(logAge=ageRange,
                               sex=c(rep(as.factor("F"),  length(ageRange))), GAbins_recode=c(rep(as.factor("VPM"),  length(ageRange))))
    
    newDataM_LPM <- data.frame(logAge=ageRange,
                               sex=c(rep(as.factor("M"), length(ageRange))), GAbins_recode=c(rep(as.factor("LPM"),  length(ageRange))))
    newDataF_LPM <- data.frame(logAge=ageRange,
                               sex=c(rep(as.factor("F"), length(ageRange))), GAbins_recode=c(rep(as.factor("LPM"),  length(ageRange))))
    
    newDataM_Term <- data.frame(logAge=ageRange,
                                sex=c(rep(as.factor("M"),  length(ageRange))), GAbins_recode=c(rep(as.factor("Term"),  length(ageRange))))
    newDataF_Term <- data.frame(logAge=ageRange,
                                sex=c(rep(as.factor("F"),  length(ageRange))), GAbins_recode=c(rep(as.factor("Term"),  length(ageRange))))
  } else  {
    newDataM <- data.frame(logAge=ageRange,
                           SurfaceHoles=c(rep(euler, length(ageRange))),
                           sex=c(rep(as.factor("M"),  length(ageRange))))
    
    newDataF <- data.frame(logAge=ageRange,
                           SurfaceHoles=c(rep(euler, length(ageRange))),
                           sex=c(rep(as.factor("F"),  length(ageRange))))
  } 
  
  # Predict phenotype values for set age range for each sex
  gammModelM_VPM <- predictAll(gamModel, newdata=newDataM_VPM, type="response",data = df)
  gammModelF_VPM <- predictAll(gamModel, newdata=newDataF_VPM, type="response",data = df)
  gammModelM_LPM <- predictAll(gamModel, newdata=newDataM_LPM, type="response",data = df)
  gammModelF_LPM <- predictAll(gamModel, newdata=newDataF_LPM, type="response",data = df)
  gammModelM_Term <- predictAll(gamModel, newdata=newDataM_Term, type="response",data = df)
  gammModelF_Term <- predictAll(gamModel, newdata=newDataF_Term, type="response",data = df)
  
  # Calculate the `cent`th centiles for the sex models
  phenoMedianPredsM_VPM <- qGG(c(cent), 
                               mu=gammModelM_VPM$mu, 
                               sigma=gammModelM_VPM$sigma, 
                               nu=gammModelM_VPM$nu)
  
  phenoMedianPredsF_VPM <- qGG(c(cent), 
                               mu=gammModelF_VPM$mu, 
                               sigma=gammModelF_VPM$sigma, 
                               nu=gammModelF_VPM$nu)
  
  phenoMedianPredsM_LPM <- qGG(c(cent), 
                               mu=gammModelM_LPM$mu, 
                               sigma=gammModelM_LPM$sigma, 
                               nu=gammModelM_LPM$nu)
  
  phenoMedianPredsF_LPM <- qGG(c(cent), 
                               mu=gammModelF_LPM$mu, 
                               sigma=gammModelF_LPM$sigma, 
                               nu=gammModelF_LPM$nu)
  
  phenoMedianPredsM_Term <- qGG(c(cent), 
                                mu=gammModelM_Term$mu, 
                                sigma=gammModelM_Term$sigma, 
                                nu=gammModelM_Term$nu)
  
  phenoMedianPredsF_Term <- qGG(c(cent), 
                                mu=gammModelF_Term$mu, 
                                sigma=gammModelF_Term$sigma, 
                                nu=gammModelF_Term$nu)
  
  # Average the two calculated centile curves 
  # For visualizations and correlations with the LBCC data
  phenoMedianPreds_VPM <- (phenoMedianPredsF_VPM + phenoMedianPredsM_VPM)/2
  phenoMedianPreds_LPM <- (phenoMedianPredsF_LPM + phenoMedianPredsM_LPM)/2
  phenoMedianPreds_Term<- (phenoMedianPredsF_Term + phenoMedianPredsM_Term)/2
  
  centOut <- list(VPM = NULL, LPM = NULL, Term = NULL, logAge = ageRange)
  centOut$VPM <- phenoMedianPreds_VPM
  centOut$LPM <- phenoMedianPreds_LPM
  centOut$Term <- phenoMedianPreds_Term
  # Return the calculated centile curve
  # Return a DF that returns the age with it
  return(centOut)
}