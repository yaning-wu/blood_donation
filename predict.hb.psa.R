# function for predicting haemoglobin levels at specific timepoints during follow-up for probabilistic sensitivity analyses - coded with assistance of ChatGPT
predict.hb.psa <- function(formula_hb, psa.hb.hb, b, newdata, ranef.hb){
  # defining regression parameters for haemoglobin submodel drawn from multivariate normal distribution
  newbeta <- psa.hb.hb$pars.hb[psa.hb.hb$b == b]
  
  # prediction using fixed effects
  model.matrix.hb <- model.matrix(~ ethnicDesc_bin + ABORH + HGB_g_dL_bl + followup + pct_ebv_removed_lagged + hb_defer_lagged, data = newdata)
  fix.ef.pred <- model.matrix.hb %*% newbeta
  
  # prediction using random effects
  ran.ef.dt <- data.table(identifier = as.numeric(rownames(ranef.hb)), ran.ef = ranef.hb[[1]])
  
  setDT(newdata)  
  setkey(ran.ef.dt, identifier)
  setkey(newdata, identifier)
  
  newdata <- ran.ef.dt[newdata, on = "identifier"]
  newdata[is.na(ran.ef), ran.ef := 0]
  
  # summing fixed and random effects
  pred <- fix.ef.pred + newdata$ran.ef
  return(pred)
}