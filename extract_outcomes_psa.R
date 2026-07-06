# function to extract outcomes from output of DES for probabilistic sensitivity analyses
extract_outcomes_psa <- function(i, outcome_name, psa_outcome_set, donorsex){
  outcome_set <- as.data.frame(psa_outcome_set[[i]][[1]])
  des_results <- outcome_set %>%
    mutate(
      eligible_pulsevvr = ifelse(pulsevvr >= 0.5 & event %in% c("donate", "bleed under threshold", "defer other", "defer lowhb"), 1, 0)) %>%
    group_by(donornum) %>%
    arrange(time, .by_group = T) %>%
    mutate(
      total_vvrs = cumsum(eligible_pulsevvr)
    ) %>%
    ungroup()
  
  outcome <- case_when(
    outcome_name == "donate" ~ sum(des_results$event == "donate" & des_results$time >= ifelse(donorsex == 1, 12, 16)),
    outcome_name == "bleed under threshold" ~ sum(des_results$event == "bleed under threshold" & des_results$time >= ifelse(donorsex == 1, 12, 16)),
    outcome_name == "donate and bleed under threshold" ~ sum(des_results$event == "donate" & des_results$time >= ifelse(donorsex == 1, 12, 16)) + sum(des_results$event == "bleed under threshold" & des_results$time >= ifelse(donorsex == 1, 12, 16)),
    outcome_name == "collectvol" ~ sum(as.numeric(des_results$collectvol[des_results$time >= ifelse(donorsex == 1, 12, 16) & des_results$event %in% c("donate", "bleed under threshold")])),
    outcome_name == "defer lowhb" ~ sum(des_results$event == "defer lowhb" & des_results$time >= ifelse(donorsex == 1, 12, 16)),
    outcome_name == "defer other" ~ sum(des_results$event == "defer other" & des_results$time >= ifelse(donorsex == 1, 12, 16)),
    outcome_name == "under threshold" ~ sum(des_results$event == "bleed under threshold" & des_results$time >= ifelse(donorsex == 1, 12, 16)) + sum(des_results$event == "defer lowhb" & des_results$time >= ifelse(donorsex == 1, 12, 16)),
    outcome_name == "vvr" ~ sum(des_results$eligible_pulsevvr == 1 & des_results$time >= ifelse(donorsex == 1, 12, 16)),
    outcome_name == "collectvol under threshold" ~ sum(as.numeric(des_results$collectvol[des_results$time >= ifelse(donorsex == 1, 12, 16) & des_results$event %in% c("bleed under threshold")])),
    outcome_name == "collectvol vvr" ~ sum(as.numeric(des_results$collectvol[des_results$time >= ifelse(donorsex == 1, 12, 16) & des_results$eligible_pulsevvr == 1])),
    outcome_name == "defer low ebv" ~ sum(des_results$event %in% c("dropout") & des_results$collectvol == 0)
  )
  
  return(outcome)
}
