#-------------------------------------------------------------------------------
# This script was assembled to walk a user through building a GAMLSS model
# using their own neuroimaging phenotype and demographic data. If you customize
# it and use it in your research, we politely ask that you cite the following 
# papers (at least the first one if your citations are limited):
#
# [SLIP paper](https://www.medrxiv.org/content/10.1101/2023.01.13.23284533v1)
# 
# [Lifespan Nature paper](https://www.nature.com/articles/s41586-022-04554-y)
#
#
# The models build in this script require the following data:
# - logAge: numeric type with log(post conception age in days, base=10). In (1),
#           we use a conversion factor of 325.25 days/year and a post conception 
#           offset of 280 days if post conception age is not available.
# - sex: factor with M or F values.
# - phenotype: numeric type with the measurements of the phenotype of interest
#
# For FreeSurfer (FS) neuroimaging phenotypes, the measurement SurfaceHoles is  
# included in the GAMLSS model. SynthSeg (SS) does not produce SurfaceHoles.
#-------------------------------------------------------------------------------

# gc()
# library(ggplot2)
# library(gamlss) #to fit model
# library(mgcv) # helps with the gam models
# library(tidymv) # helps with the gam models
# source("lib_mpr_analysis.r")

growthchart_model <- function(p, df, covs = c("sex"), agevar = "logAge", formula = "~fp(logAge, npoly=3) + sex - 1", n.cyc = 200, desiredCentiles = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), predictFUN = predictCentilesForAgeRange) {
  df <- df %>% select(all_of(c(p, covs, agevar)))
  df <- df[order(get(agevar, df)),] %>% na.omit(.) %>% rename(logAge = agevar)
  #The df should have no NAs. So create df for the specific variables needed. Or remove all the variables with NAs.
  formulaSs <- as.formula(paste0(p, formula))
  growthChartModel <-gamlss(formula = formulaSs,
                            sigma.formula = formulaSs,
                            nu.formula = as.formula(paste0(p, "~1")),
                            family = GG,
                            data = df,
                            control = gamlss.control(n.cyc = n.cyc),  # See (2)
                            trace = F)
  
  # Predict a set of centiles for the model
  centileCurvesSs <- c() 
  desiredCentiles <- desiredCentiles
  for (i in c(1:length(desiredCentiles))){
    centileCurvesSs[[i]] <- do.call(predictFUN, list(gamModel = growthChartModel, ageRange = df$logAge, df = df,
                                                       cent=desiredCentiles[[i]]))}
  out <- list(model = NULL, centileCurves = NULL)
  out$model <- growthChartModel
  out$centileCurves <- centileCurvesSs %>% set_names(desiredCentiles)
  return(out)
}