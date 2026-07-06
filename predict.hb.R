# function for predicting haemoglobin levels at specific timepoints during follow-up for main analyses - coded with assistance of ChatGPT
predict.hb <- function(formula_hb, newdata, ranef.hb, beta.hb){
  # prediction using fixed effects
  model.matrix.hb <- model.matrix(~ ethnicDesc_bin + ABORH + HGB_g_dL_bl + followup + pct_ebv_removed_lagged + hb_defer_lagged, data = newdata)
  fix.ef.pred <- model.matrix.hb %*% beta.hb
  
  # prediction using random effects
  ran.ef.dt <- data.table(identifier = as.numeric(rownames(ranef.hb)), ran.ef = ranef.hb[[1]])
  
  newdata <- setDT(newdata)  
  setkey(ran.ef.dt, identifier)
  setkey(newdata, identifier)
  
  newdata <- ran.ef.dt[newdata, on = "identifier"]
  newdata[is.na(ran.ef), ran.ef := 0]
  
  # summing fixed and random effects
  pred <- fix.ef.pred + newdata$ran.ef
  return(pred)
}
