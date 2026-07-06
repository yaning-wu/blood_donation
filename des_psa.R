# function for running DES model in 1K hypothetical donors for probabilistic sensitivity analyses
des_boot_new_sens_bmi <- function(b,
                         N, 
                         donorsex,
                         group,
                         vol,
                         ind_data,
                         attend,
                         timehorizon,
                         ranef.hb,
                         formula_hb,
                         psa.hb.hb,
                         psa.hb.dropout,
                         psa.hb.defer,
                         psa.hb.othdefer,
                         psa.hb.vvr,
                         psa.hb.fpmod,
                         thresh,
                         curr.recall,
                         lowhbdef.recall,
                         vlowhbdef.recall,
                         othdef.recall,
                         streams) 
{
  
  # defining matrix to denote mutually exclusive events occurring during 18-month follow-up for donation attendances (excludes VVRs, which are not mutually exclusive)
  allevents <- matrix(nrow = N, ncol = 9) 
  colnames(allevents) <- c("baseline", "attend", "test", "recall", "donate", "defer lowhb", "bleed under threshold", "dropout", "defer other") 
  eventhistory <- data.frame(event = NULL, time = NULL)
  
  # sampling 10K hypothetical donors from the STRIDES PDT dataset
  sample_data <- ind_data[sample(.N, N, replace = T)]
  
  # applying the DES function to 1K hypothetical donors  
  result <- lapply(1:N, 
                   des.donor_boot_sens_bmi, 
                   donorsex, 
                   group,
                   vol,
                   sample_data, 
                   attend, 
                   timehorizon, 
                   ranef.hb, 
                   formula_hb,
                   psa.hb.hb,
                   psa.hb.dropout,
                   psa.hb.defer,
                   psa.hb.othdefer,
                   psa.hb.vvr,
                   psa.hb.fpmod,
                   b,
                   thresh, 
                   curr.recall, 
                   lowhbdef.recall, 
                   vlowhbdef.recall, 
                   othdef.recall,
                   streams)
  
  # documenting number of events occurring in 1K hypothetical donors during 18-month follow-up
  eventlist <- sapply(result,"[", "eventlist")
  timelist <- sapply(result,"[","timelist")
  donornumlist <- sapply(result,"[","donornumlist")
  collectvollist <- sapply(result, "[", "collectvollist")
  newdonorlist <- sapply(result, "[", "newdonorlist")
  pulsevvrlist <- sapply(result, "[", "pulsevvrlist")
  
  hblist <- sapply(result,"[","hblist")
  hbcurrlist <- sapply(result,"[","hbcurrlist") 
  collectvollist <- sapply(result,"[","collectvollist") 
  
  eventhistory <- data.frame(
    donornum   = as.numeric(unlist(donornumlist)),
    event      = unlist(eventlist),
    time       = unlist(timelist),
    hbcurr     = unlist(hbcurrlist),
    collectvol = unlist(collectvollist),
    newdonor = unlist(newdonorlist),
    pulsevvr    = unlist(pulsevvrlist),
    stringsAsFactors = F
  )
  
  allevents <- t(sapply(eventlist, function(x){
    table(factor(x, 
                 levels = colnames(allevents)))
  }))
  eventhistory$event <- factor(eventhistory$event, 
                               levels = colnames(allevents))
  return(list(eventhistory = eventhistory, 
              allevents = allevents,
              data = data)) 
} 

# function for running DES model for each hypothetical donor
des.donor_boot_sens_bmi <- function(i, 
                           donorsex,
                           group,
                           vol,
                           sample_data, 
                           attend, 
                           timehorizon, 
                           ranef.hb,
                           formula_hb,
                           psa.hb.hb,
                           psa.hb.dropout,
                           psa.hb.defer,
                           psa.hb.othdefer,
                           psa.hb.vvr,
                           psa.hb.fpmod,
                           b,
                           thresh,
                           curr.recall,
                           lowhbdef.recall,
                           vlowhbdef.recall,
                           othdef.recall,
                           streams){
  
  # ensuring that random numbers are generated in the same way for each run of the function to ensure reproducibility - lines 112-113 coded with assistance from ChatGPT 
  oldseed <- .Random.seed
  .Random.seed <<- streams[[i]]
  
  eventlist <- vector()
  timelist <- vector()
  nexttime <- vector()
  hblist <- vector()
  hbcurrlist <- vector()
  collectvollist <- vector()
  newdonorlist <- vector()
  donornumlist <- vector()
  pulsevvrlist <- vector()
  
  # defining initial parameters prior to entry into follow-up
  j <- 0
  currenttime <- 0
  
  last.donation <- 0
  last.attend <- 0
  timesince.donation <- 0
  initial.donation <- 1
  
  nexttime[c("attend", "test", "donate", "defer lowhb", 
             "bleed under threshold", "dropout", "defer other")] <- Inf

  covars.i <- as.data.table(sample_data[i, ])
  
  donornum <- i
  
  nexttime["recall"] <- last.attend + curr.recall
  
  # defining parameters and outcomes immediately prior to baseline attendance during follow-up
  j <- 1
  timelist[j] <- 0
  
  event <- "baseline"
  eventlist[j] <- event
  donornumlist[j] <- i
  hb.curr <- covars.i$HGB_g_dL_bl
  hbcurrlist[j] <- hb.curr
  
  collectvollist[j] <- 0
  
  pulsevvr.curr <- covars.i$vvr
  pulsevvrlist[j] <- pulsevvr.curr
  
  newdonorlist[j] <- covars.i$new_donor
  
  # defining parameters and outcomes at baseline attendance during follow-up
  j <- 2
  timelist[j] <- 0.01
  set(covars.i, j = "agePulse", value = covars.i$agePulse + 0.01)
  
  event <- "donate"
  
  eventlist[j] <- event
  donornumlist[j] <- i
  
  hb.curr <- covars.i$HGB_g_dL_bl
  covars.i[, followup := 0.01]
  
  ## defining proportion of EBV removed at attendance
  covars.i[, attendance_outcome := "Donate_Normal"]
  covars.i[, pct_ebv_removed := vol_strat(
    data = covars.i,
    group = group, 
    vol = vol,
    ebv_var = "ebv"
  )]
  covars.i[, pct_ebv_removed_lagged := covars.i$pct_ebv_removed]
  covars.i[, vol_removed_lagged := covars.i$vol_removed]
  covars.i[, hb_defer_lagged := covars.i$hb_defer]
  
  ## predicting haemoglobin levels immediately after baseline attendance using regression parameters drawn from multivariate normal distribution
  hb.curr <- predict.hb.psa(formula_hb = formula_hb, psa.hb.hb = psa.hb.hb, b = b, newdata = covars.i, ranef.hb = ranef.hb)
  hbcurrlist[j] <- hb.curr
  
  set(covars.i, j = "HGB_g_dL_bl", value = hb.curr)
  
  collectvollist[j] <- pmax((covars.i$ebv*covars.i$pct_ebv_removed_lagged/100) - 0.025, 0)
  
  covars.i[, vvr_lagged := pulsevvr.curr]
  
  ## manually predicting risk of VVR at baseline attendance using regression parameters drawn from multivariate normal distribution
  prob.pulsevvr <- inv.logit(psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][1] + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][2]*covars.i$pct_ebv_removed + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][3]*covars.i$agePulse + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][4]*(as.numeric(covars.i$ethnicDesc_bin) - 1) + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][5]*covars.i$site_type + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][6]*ifelse(covars.i$bmi_cat == "obesity", 1, 0) + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][7]*ifelse(covars.i$bmi_cat == "overweight", 1, 0) + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][8]*ifelse(covars.i$bmi_cat == "underweight", 1, 0))
  pulsevvr.curr <- rbinom(1, 1, prob.pulsevvr)
  
  pulsevvrlist[j] <- pulsevvr.curr 
  
  covars.i.new <- covars.i
  
  covars.i.new[, vvr_lagged := ifelse(pulsevvr.curr >= 0.5, 1, 0)]
  covars.i.new[, VVR_hist_bin := ifelse(pulsevvr.curr >= 0.5, 1, covars.i.new$VVR_hist_bin)]
  
  newdonorlist[j] <- covars.i$new_donor
  
  # simulating individuals' repeated attendances during follow-up
  repeat {
    j <- j + 1
    event <- names(which.min(nexttime))
    time <- nexttime[which.min(nexttime)]
    
    currenttime <- time
    
    if(currenttime > as.numeric(timehorizon)){ 
      break 
    }
    
    covars.i[, agePulse_new := covars.i$agePulse + currenttime/52]
    
    covars.i[, new_donor_new := 0]
    
    # defining individual-level variables, parameters, and subsequent DES flow at successful donation (including over- and under-threshold donations)
    if(event %in% c("donate", "bleed under threshold")){
      last.donation <- currenttime
      initial.donation <- 0
      
      covars.i.new[, VVR_hist_bin := ifelse(pulsevvr.curr >= 0.5, 1, covars.i.new$VVR_hist_bin)]
      covars.i.new[, agePulse := covars.i$agePulse_new]
      covars.i.new[, new_donor := covars.i$new_donor_new]
      covars.i.new[, followup := currenttime]
      covars.i.new[, attendance_outcome := "Donate_Normal"]
      covars.i.new[, pct_ebv_removed := vol_strat(
        data = covars.i.new,
        group = group, 
        vol = vol,
        ebv_var = "ebv"
      )]
      
      ## lines 249-253 optimised with assistance from ChatGPT
      collectvollist[j] <- pmax((covars.i.new$ebv*covars.i.new$pct_ebv_removed/100) - 0.025, 0) 
      
      covars.i.new[, pct_ebv_removed_lagged := covars.i.new$pct_ebv_removed]
      covars.i.new[, vol_removed_lagged := (covars.i.new$pct_ebv_removed/100)*covars.i.new$ebv]
      
      prob.pulsevvr <- inv.logit(psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][1] + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][2]*covars.i.new$pct_ebv_removed + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][3]*covars.i.new$agePulse + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][4]*(as.numeric(covars.i.new$ethnicDesc_bin) - 1) + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][5]*covars.i.new$site_type + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][6]*ifelse(covars.i.new$bmi_cat == "obesity", 1, 0) + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][7]*ifelse(covars.i.new$bmi_cat == "overweight", 1, 0) + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][8]*ifelse(covars.i.new$bmi_cat == "underweight", 1, 0))
      pulsevvr.curr <- rbinom(1, 1, prob.pulsevvr)
      
      pulsevvrlist[j] <- pulsevvr.curr
      
      covars.i.new[, vvr_lagged := ifelse(pulsevvr.curr >= 0.5, 1, 0)]
      covars.i.new[, hb_defer_lagged := 0]
      covars.i.new[, hb_defer := 0]
      
      last.attend <- currenttime
      
      nexttime["recall"] <- last.attend + curr.recall
      
      newdonorlist[j] <- covars.i.new$new_donor
    }
    
    timesince.donation <- currenttime - last.donation
    
    timesince.donation.r <- round(timesince.donation, digits = 0)
    if(timesince.donation.r > as.numeric(timehorizon)){ 
      break
    } 
    
    nexttime <- subset(nexttime, !(names(nexttime) %in% event))
    
    # defining individual-level variables, parameters, and subsequent DES flow at donor recall (i.e., invitation to donate)
    if(event == "recall"){
      covars.i.new[, vvr := ifelse(pulsevvr.curr >= 0.5, 1, 0)] 
      
      prob.dropout <- inv.logit(psa.hb.dropout$pars.dropout[psa.hb.dropout$b == b][1] + psa.hb.dropout$pars.dropout[psa.hb.dropout$b == b][2]*(pulsevvr.curr >= 0.5))
      dropout <- rbinom(1, 1, prob.dropout)
      
      covars.i.new[, agePulse := covars.i$agePulse_new]
      
      coefs <- c(psa.hb.fpmod$pars.fpmod[psa.hb.fpmod$b == b])
      names(coefs) <- names(fpmod$coefficients)
      
      nexttime["attend"] <- attend(currenttime, eventlist, timelist, covi = covars.i.new, fpmod, coefs = coefs, logcumhaz, curr.recall = curr.recall)
      
      if(dropout == 1 | covars.i.new$pct_ebv_removed_lagged == 0){ 
        nexttime["attend"] <- Inf
        nexttime["dropout"] <- currenttime + 0.01
      }
      
      nexttime["recall"] <- Inf
      
      collectvollist[j] <- 0
      newdonorlist[j] <- covars.i.new$new_donor
      
    }
    
    # defining individual-level variables, parameters, and subsequent DES flow at individual attendance at donation appointment
    if(event == "attend"){
      prob.othdefer <- inv.logit(psa.hb.othdefer$pars.othdefer[psa.hb.othdefer$b == b])
      defer.other <- rbinom(1, 1, prob.othdefer)
      
      if(defer.other == 1){
        nexttime["defer other"] <- currenttime + 0.01
      } else {
        nexttime["test"] <- currenttime + 0.01
      }
      
      collectvollist[j] <- 0
      newdonorlist[j] <- covars.i.new$new_donor
    }
    
    # defining individual-level variables, parameters, and subsequent DES flow at haemoglobin testing
    if(event == "test"){ 
      nexttime["test"] <- Inf
      
      covars.i.new[, agePulse := covars.i$agePulse_new]
      
      covars.i.new[, followup := currenttime]
      
      hb.curr <- predict.hb.psa(formula_hb = formula_hb, psa.hb.hb = psa.hb.hb, b = b, newdata = covars.i.new, ranef.hb = ranef.hb)
      
      if(hb.curr < thresh){ 
        pdefer = inv.logit(psa.hb.defer$pars.defer[psa.hb.defer$b == b])
        d <- runif(1, 0, 1)
        if (d <= pdefer){
          nexttime["defer lowhb"] <- currenttime + 0.01
        }
        if (d > pdefer){
          nexttime["bleed under threshold"] <- currenttime + 0.01
        }
      } 
      else {
        nexttime["donate"] <- currenttime + 0.01
      }
      
      collectvollist[j] <- 0
      newdonorlist[j] <- covars.i.new$new_donor
    }
    
    # defining individual-level variables, parameters, and subsequent DES flow at deferral due to low haemoglobin
    if(event == "defer lowhb"){
      covars.i.new[, agePulse := covars.i$agePulse_new]
      
      covars.i.new[, new_donor := covars.i$new_donor_new]
      covars.i.new[, followup := currenttime]
      covars.i.new[, pct_ebv_removed := 0]
      covars.i.new[, pct_ebv_removed_lagged := 0]
      covars.i.new[, vol_removed_lagged := 0]
      
      collectvollist[j] <- 0 
      
      prob.pulsevvr <- inv.logit(psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][1] + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][2]*covars.i.new$pct_ebv_removed + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][3]*covars.i.new$agePulse + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][4]*(as.numeric(covars.i.new$ethnicDesc_bin) - 1) + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][5]*covars.i.new$site_type + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][6]*ifelse(covars.i.new$bmi_cat == "obesity", 1, 0) + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][7]*ifelse(covars.i.new$bmi_cat == "overweight", 1, 0) + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][8]*ifelse(covars.i.new$bmi_cat == "underweight", 1, 0))
      pulsevvr.curr <- rbinom(1, 1, prob.pulsevvr)
      
      pulsevvrlist[j] <- pulsevvr.curr
      
      covars.i.new[, vvr_lagged := ifelse(pulsevvr.curr >= 0.5, 1, 0)]
      covars.i.new[, hb_defer_lagged := 1]
      covars.i.new[, other_defer_lagged := 0]
      covars.i.new[, hb_defer := 1]
      
      last.attend <- currenttime
      
      nexttime["recall"] <- last.attend + lowhbdef.recall
      if(hb.curr < thresh - 1){
        nexttime["recall"] <- last.attend + vlowhbdef.recall
      }
      
      newdonorlist[j] <- covars.i.new$new_donor
    }
    
    # defining individual-level variables, parameters, and subsequent DES flow at deferral due to non-haemoglobin-related reasons (e.g., health or travel history)
    if(event == "defer other"){
      covars.i.new[, agePulse := covars.i$agePulse_new]
      
      covars.i.new[, new_donor := covars.i$new_donor_new]
      covars.i.new[, followup := currenttime]
      covars.i.new[, pct_ebv_removed := 0]
      covars.i.new[, pct_ebv_removed_lagged := 0]
      covars.i.new[, vol_removed_lagged := 0]
      
      collectvollist[j] <- 0
      
      prob.pulsevvr <- inv.logit(psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][1] + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][2]*covars.i.new$pct_ebv_removed + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][3]*covars.i.new$agePulse + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][4]*(as.numeric(covars.i.new$ethnicDesc_bin) - 1) + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][5]*covars.i.new$site_type + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][6]*ifelse(covars.i.new$bmi_cat == "obesity", 1, 0) + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][7]*ifelse(covars.i.new$bmi_cat == "overweight", 1, 0) + psa.hb.vvr$pars.vvr[psa.hb.vvr$b == b][8]*ifelse(covars.i.new$bmi_cat == "underweight", 1, 0))
      pulsevvr.curr <- rbinom(1, 1, prob.pulsevvr)
      
      pulsevvrlist[j] <- pulsevvr.curr
      
      covars.i.new$vvr_lagged <- ifelse(pulsevvr.curr >= 0.5, 1, 0)
      covars.i.new[, hb_defer_lagged := 0]
      covars.i.new[, other_defer_lagged := 1]
      covars.i.new[, hb_defer := 0]
      
      last.attend <- currenttime
      
      nexttime["recall"] <- last.attend + othdef.recall
      
      newdonorlist[j] <- covars.i.new$new_donor
    }
    
    # defining individual-level variables, parameters, and subsequent DES flow at donor dropout (i.e., subsequent non-return for donation attendances during follow-up)
    if(event == "dropout"){
      collectvollist[j] <- pmax((covars.i$ebv*covars.i$pct_ebv_removed_lagged/100) - 0.025, 0) # this is NOT actually the collected volume, just used to determine whether low EBV deferral occurred
    }
    
    timelist[j] <- time
    eventlist[j] <- event
    donornumlist[j] <- i
    hbcurrlist[j] <- hb.curr
    collectvollist[j] <- collectvollist[j]
    newdonorlist[j] <- newdonorlist[j]
    pulsevvrlist[j] <- pulsevvr.curr
  }
  
  ## lines 427-428 coded with assistance from ChatGPT
  streams[[i]] <<- .Random.seed   
  .Random.seed <<- oldseed        
  
  # documenting outcomes incurred by individual hypothetical donor during follow-up
  return(list(donornumlist = donornumlist, 
              eventlist = eventlist, 
              timelist = timelist, 
              currenttime = currenttime, 
              hbcurrlist = hbcurrlist,
              collectvollist = collectvollist,
              newdonorlist = newdonorlist,
              pulsevvrlist = pulsevvrlist
  ))
}

# function to simulate attendance times of individual hypothetical donors using flexible parametric survival model with regression parameters drawn from multivariate normal distrbution
attend.dist_boot <- function(currenttime, eventlist, timelist, covi, fpmod, coefs, logcumhaz, curr.recall){
  
  attenddelay <- tryCatch({
    simsurv(betas = coefs,
            x = covi,             
            knots = fpmod$knots,    
            logcumhazard = logcumhaz,  
            maxt = 100000,               
            interval = c(1e-8, 100001))
    simsurv(betas = coefs,
            x = covi,             
            knots = fpmod$knots,    
            logcumhazard = logcumhaz,  
            maxt = 100000,               
            interval = c(1e-8, 100001))
  }, error = function(e) 0)
  
  attenddelay <- tryCatch({
    ((attenddelay$eventtime/24)/7)
  }, error = function(e) 0)
  
  delay <- (max((curr.recall - 12), attenddelay)) - (max(curr.recall - 12))
  currenttime + delay
}