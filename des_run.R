# code for running DES and PSA for all strategies
# loading libraries
invisible(suppressPackageStartupMessages(lapply(c("devtools", "tidyverse", "fastDummies", "lme4", "arm", "boot", "data.table", "MASS", "flexsurv", "simsurv", "future.apply", "logistf", "parallel"), 
                                                library, 
                                                character.only = T)))

# loading self-defined functions
allfuns <- c("des_main.R",
             "des_psa.R",
             "extract_outcomes_des.R",
             "extract_outcomes_psa.R",
             "predict.hb.R",
             "predict.hb.psa.R",
             "vol_strat.R",
             "vol_strat_usable.R")

invisible(lapply(allfuns, source))

# loading cleaned STRIDES PDT data
data <- as.data.table(read_csv("cleaned_data.csv", show_col_types = F))

# cleaning data 
## defining correct data types
data <- data %>%
  mutate(
    ABORH = factor(ABORH, levels = c("A-", "A+", "O+", "O-", "B+", "B-", "AB+", "AB-", NA)),
    followup = as.numeric(followup),
    timebetween_don_weeks = as.numeric(timebetween_don_weeks),
    ethnicDesc_bin = as.factor(ethnicDesc_bin),
    new_donor = as.numeric(new_donor),
    nLowHB_hx2_grouped = as.factor(nLowHB_hx2_grouped),
    nDonation_hx2_grouped = as.factor(nDonation_hx2_grouped),
    donationDate = as.Date(donationDate)
  ) %>%
  setDT()

## defining new variables
data <- data %>%
  group_by(identifier) %>%
  dplyr::arrange(donationDate, .by_group = T) %>%
  mutate(
    # new variable denoting whether VVR occurred at previous attendance
    vvr_lagged = lag(vvr),
    # new variable denoting % of EBV removed under current strategy (dependent on donation outcome)
    pct_ebv_removed = ifelse(
      attendance_outcome %in% c("Donate_Normal", "Donate_Under"),
      100*0.500/ebv,
      0
    ),
    # new variable denoting % of EBV removed at previous attendance under current strategy
    pct_ebv_removed_lagged = ifelse(
      donationDate == min(donationDate),
      NA, 
      lag(pct_ebv_removed)
    )
  ) %>%
  ungroup() %>%
  mutate(
    # new variable denoting volume of blood collected under current strategy (dependent on donation outcome)
    vol_removed = ifelse(attendance_outcome %in% c("Donate_Normal", "Donate_Under"),
                         0.500,
                         0)) %>%
  group_by(identifier) %>%
  dplyr::arrange(donationDate, .by_group = T) %>%
  mutate(
    # new variable denoting volume of blood collected at previous attendance under current strategy 
    vol_removed_lagged = lag(vol_removed)
  ) %>%
  ungroup() %>%
  group_by(identifier) %>%
  dplyr::arrange(donationDate, .by_group = T) %>% 
  mutate(
    # new variable denoting whether deferral due to low haemoglobin occurred at previous attendance
    hb_defer_lagged = ifelse(
      lag(hb_defer) == 1, 
      1,
      0
    )
  ) %>%
  ungroup()

# extracting function arguments from input into high-performance computing (HPC) environment - lines 22, 24-26 coded with the assistance of ChatGPT
args <- R.utils::commandArgs(trailingOnly = T)

strattype <- as.character(args[2]) 
stratvol <- as.numeric(args[4])
donorsex <- as.numeric(args[6])

# defining inputs into DES model
timehorizon <- 78 # 18-month follow-up
NDES <- 10000 # 10K hypothetical donors randomly drawn for main analyses
NUI <- 1000 # 1K hypothetical donors randomly drawn for probabilistic sensitivity analyses
B <- 400 # regression parameters drawn from multivariate normal distributions 400 times for probabilistic sensitivity analyses

# defining formulae for regression models for haemoglobin and VVR
formula_hb <- "ethnicDesc_bin + ABORH + HGB_g_dL_bl + followup + pct_ebv_removed_lagged + hb_defer_lagged"
formula_vvr <- "pct_ebv_removed + agePulse + ethnicDesc_bin + site_type + bmi_cat"

# setting random seed
set.seed(9)

# filtering dataset to include only male or female donors
data <- data[sexPulse == donorsex]

# haemoglobin submodel (linear mixed model)
hb_mod <- lmer(as.formula(paste0("HGB_g_dL_pd ~ (1|identifier) + ", formula_hb)),
               data = data %>%
                 drop_na(ethnicDesc_bin, nDonation_hx2) %>%
                 filter(dropout == 0) %>%
                 group_by(identifier) %>%
                 dplyr::arrange(donationDate, .by_group = T) %>%
                 filter(donationDate != min(donationDate)) %>%
                 ungroup() %>%
                 setDT())
ranef.hb <- ranef(hb_mod)$identifier
beta.hb <- hb_mod@beta

# VVR submodel (Firth's logistic regression model with FLIC)
vvr_lr_mod <- logistf(as.formula(paste0("vvr ~ ", formula_vvr)),
                      data = data %>%
                        group_by(identifier) %>%
                        dplyr::arrange(donationDate, .by_group = T) %>%
                        filter(donationDate != min(donationDate) & dropout == 0) %>%
                        ungroup() %>%
                        drop_na(ethnicDesc_bin, nDonation_hx2) %>%
                        setDT(),
                      flic = T)

# "imputing" over- vs under-threshold haemoglobin
## dataset without initial donation
post_index_data <- data %>%
  group_by(identifier) %>%
  dplyr::arrange(donationDate, .by_group = T) %>%
  filter(donationDate != min(donationDate) & dropout == 0) %>%
  ungroup() %>%
  drop_na(ethnicDesc_bin, nDonation_hx2)

## mutating variable to denote index donation
data <- data %>%
  group_by(identifier) %>%
  dplyr::arrange(donationDate, .by_group = T) %>%
  mutate(
    index_donation_indicator = ifelse(donationDate == min(donationDate),
                                      1,
                                      0)
  ) %>%
  ungroup()

## "imputing" over- vs under-threshold donation for donations with missing haemoglobin values
prop_under_threshold_don <- nrow(post_index_data[!is.na(post_index_data$hb_under) & post_index_data$hb_under == 1 & post_index_data$attendance_outcome %in% c("Donate_Normal", "Donate_Under") & !is.na(post_index_data$ethnicDesc_bin) & !is.na(post_index_data$nDonation_hx2),])/nrow(post_index_data[!is.na(post_index_data$hb_under) & post_index_data$attendance_outcome %in% c("Donate_Normal", "Donate_Under") & !is.na(post_index_data$ethnicDesc_bin) & !is.na(post_index_data$nDonation_hx2),])
rows_to_impute <- data[data$attendance_outcome %in% c("Donate_Normal", "Donate_Under") & is.na(data$hb_under) & data$index_donation_indicator == 0 & data$dropout == 0 & !is.na(data$ethnicDesc_bin) & !is.na(data$nDonation_hx2),]
data$hb_under[data$attendance_outcome %in% c("Donate_Normal", "Donate_Under") & is.na(data$hb_under) & data$index_donation_indicator == 0 & data$dropout == 0 & !is.na(data$ethnicDesc_bin) & !is.na(data$nDonation_hx2)] <- sample(c(rep(1, round(prop_under_threshold_don*nrow(rows_to_impute), 0)), rep(0, round((1 - prop_under_threshold_don)*nrow(rows_to_impute), 0))))

# low haemoglobin deferral submodel (among attendances where donors presented with under-threshold haemoglobin; logistic regression model)
hb_defer_mod <- glm(hb_defer ~ NULL,
                    family = "binomial",
                    data = data %>%
                      group_by(identifier) %>%
                      dplyr::arrange(donationDate, .by_group = T) %>%
                      filter(hb_under == 1) %>%
                      filter(donationDate != min(donationDate) & dropout == 0) %>%
                      filter(attendance_outcome %in% c("Donate_Normal", "Donate_Under", "Defer_LowHB")) %>%
                      ungroup() %>%
                      drop_na(ethnicDesc_bin, nDonation_hx2) %>%
                      setDT())

# non-haemoglobin-related deferral submodel (logistic regression model)
other_defer_mod <- glm(other_defer ~ NULL,
                       family = "binomial",
                       data = data %>%
                         group_by(identifier) %>%
                         dplyr::arrange(donationDate, .by_group = T) %>%
                         filter(donationDate != min(donationDate) & dropout == 0) %>%
                         ungroup() %>%
                         drop_na(ethnicDesc_bin, nDonation_hx2) %>%
                         setDT())

# dropout (i.e., non-return within 18-month follow-up; logistic regression model)
dropout_mod <- glm(dropout ~ vvr + hb_defer,
                   family = "binomial",
                   data = data %>%
                     group_by(identifier) %>%
                     dplyr::arrange(donationDate, .by_group = T) %>%
                     filter(donationDate == min(donationDate)) %>%
                     ungroup() %>%
                     drop_na(ethnicDesc_bin, nDonation_hx2) %>%
                     setDT())

# defining default parameters for DES
params <- list()

params[["defaults"]] <- c( 
  # haemoglobin threshold for acceptance (sex-specific)
  thresh = ifelse(donorsex == 1,
                  13.5,
                  12.5), 
  # recall times between donations (sex-specific)
  curr.recall = ifelse(donorsex == 1,
                       12,
                       16),
  # recall times after low haemoglobin deferral where haemoglobin levels were <=1 g/dL lower than thresholds 
  lowhbdef.recall = 12, 
  # recall times after low haemoglobin deferral where haemoglobin levels were >1 g/dL lower than thresholds
  vlowhbdef.recall = 52, 
  # assumed recall times after deferral due to non-haemoglobin-related reasons (e.g., health or travel history)
  othdef.recall = 4) 

thresh <- params[["defaults"]]["thresh"]
curr.recall <- params[["defaults"]]["curr.recall"]
lowhbdef.recall <- params[["defaults"]]["lowhbdef.recall"]
vlowhbdef.recall <- params[["defaults"]]["vlowhbdef.recall"]
othdef.recall <- params[["defaults"]]["othdef.recall"]

# Function to simulate time to return within DES model using regression parameters from flexible parametric survival model (below)
logcumhaz <- function(t, x, betas, knots){
  basis <- flexsurv::basis(knots, log(t)) 
  res <-
    betas[["gamma0"]]*basis[[1]] +
    betas[["gamma1"]]*basis[[2]] +
    betas[["gamma2"]]*basis[[3]] +
    betas[["gamma3"]]*basis[[4]] +
    betas[["gamma4"]]*basis[[5]] +
    betas[["gamma5"]]*basis[[6]] +
    betas[["agePulse"]]*x[["agePulse"]] +
    betas[["ABORHA+"]]*x[["ABORH_A+"]] +
    betas[["ABORHAB-"]]*x[["ABORH_AB-"]] +
    betas[["ABORHAB+"]]*x[["ABORH_AB+"]] +
    betas[["ABORHB-"]]*x[["ABORH_B-"]] +
    betas[["ABORHB+"]]*x[["ABORH_AB+"]] +
    betas[["ABORHO-"]]*x[["ABORH_O-"]] +
    betas[["ABORHO+"]]*x[["ABORH_O+"]] +
    betas[["ethnicDesc_bin1"]]*(as.numeric(x[["ethnicDesc_bin"]]) - 1) +
    betas[["nLowHB_hx2_grouped1"]]*ifelse(x[["nLowHB_hx2_grouped"]] == 1, 1, 0) +
    betas[["nLowHB_hx2_grouped2"]]*ifelse(x[["nLowHB_hx2_grouped"]] == 2, 1, 0) +
    betas[["nDonation_hx2_grouped1"]]*ifelse(x[["nDonation_hx2_grouped"]] == 1, 1, 0) +
    betas[["nDonation_hx2_grouped2"]]*ifelse(x[["nDonation_hx2_grouped"]] == 2, 1, 0) +
    betas[["HGB_g_dL_bl"]]*x[["HGB_g_dL_bl"]] +
    betas[["vvr_lagged"]]*x[["vvr_lagged"]] +
    betas[["hb_defer_lagged"]]*x[["hb_defer_lagged"]]
}

# time to return submodel (flexible parametric survival model with four knots)
fpmod <- flexsurv::flexsurvspline(Surv(timebetween_recall_hours, return) ~ agePulse + ABORH + ethnicDesc_bin + nLowHB_hx2_grouped + nDonation_hx2_grouped + HGB_g_dL_bl + vvr_lagged + hb_defer_lagged,
                                  data = data %>%
                                    # eliminating ties between "survival times"
                                    mutate(
                                      timebetween_recall_hours = jitter(timebetween_recall_hours,
                                                                        factor = 1e-30) 
                                    ) %>% 
                                    # adding four weeks to recall times for female donors (in original DES model, some female donors permitted to return every 12 weeks)
                                    mutate(
                                      timebetween_recall_hours = ifelse(sexPulse == 1,
                                                                        timebetween_recall_hours,
                                                                        timebetween_recall_hours + 4*7*24)
                                    ) %>% 
                                    # removing unnecessary rows denoting dummy return attendance among individuals who dropped out after first attendance
                                    filter(!(dropout == 1 & donationDate == as.Date("2022-11-09"))) %>%
                                    # removing individuals with missing values for ethnicity and number of low haemoglobin deferrals in the past two years
                                    drop_na(ethnicDesc_bin, nLowHB_hx2), 
                                  k = 4)

# mutating dataset to include relative EBV removed (and previous instance of EBV removed) according to evaluated strategy
data <- data %>%
  mutate(
    pct_ebv_removed = vol_strat(data = data,
                                group = as.character(strattype),
                                vol = as.numeric(stratvol),
                                ebv_var = "ebv")
  ) %>%
  group_by(identifier) %>%
  dplyr::arrange(donationDate, .by_group = T) %>%
  mutate(
    pct_ebv_removed_lagged = ifelse(
      donationDate == min(donationDate),
      NA, 
      lag(pct_ebv_removed)
    )
  ) %>%
  ungroup()

# mutating dataset to include absolute volume removed (full volume if successful donation; 0 if unsuccessful donation/deferral) according to evaluated strategy
data <- data %>%
  mutate(
    vol_removed = (vol_strat(data = data,
                             group = as.character(strattype),
                             vol = as.numeric(stratvol),
                             ebv_var = "ebv")/100)*ebv) %>%
  group_by(identifier) %>%
  dplyr::arrange(donationDate, .by_group = T) %>%
  mutate(
    vol_removed_lagged = lag(vol_removed)
  ) %>%
  ungroup()

# ensuring that random numbers are generated in the same way for each run of the DES function to ensure reproducibility - lines 301-302 and 304-308 coded with assistance from ChatGPT 
RNGkind("L'Ecuyer-CMRG")
set.seed(9)

streams <- vector("list", NDES)
for (i in seq_len(NDES)) {
  streams[[i]] <- .Random.seed
  .Random.seed <- nextRNGStream(.Random.seed)
}

# running DES models for 10K donors
model.hb <- des_new(N = NDES, 
                    donorsex = donorsex,
                    group = strattype,
                    vol = as.numeric(stratvol),
                    ind_data = data %>%
                      drop_na(ethnicDesc_bin, nLowHB_hx2) %>%
                      group_by(identifier) %>%
                      arrange(donationDate, .by_group = T) %>%
                      filter(donationDate == min(donationDate)) %>%
                      ungroup() %>%
                      setDT(),
                    attend = attend.dist,
                    timehorizon = timehorizon,
                    ranef.hb = ranef.hb,
                    beta.hb = beta.hb,
                    formula_hb = formula_hb,
                    vvr_lr_mod = vvr_lr_mod,
                    hb_defer_mod = hb_defer_mod,
                    other_defer_mod = other_defer_mod,
                    dropout_mod = dropout_mod,
                    thresh = thresh,
                    curr.recall = curr.recall,
                    lowhbdef.recall = lowhbdef.recall,
                    vlowhbdef.recall = vlowhbdef.recall,
                    othdef.recall = othdef.recall,
                    streams = streams)

# collating beta coefficients from above regression submodels
hb_mod_b <- c(cons = hb_mod@beta[1],
              ethnicDesc_bin = hb_mod@beta[2],
              ABORHAplus = hb_mod@beta[3],
              ABORHOplus = hb_mod@beta[4],
              ABORHOminus = hb_mod@beta[5],
              ABORHBplus = hb_mod@beta[6],
              ABORHBminus = hb_mod@beta[7],
              ABORHABplus = hb_mod@beta[8],
              ABORHABminus = hb_mod@beta[9],
              HGB_g_dL_bl = hb_mod@beta[10],
              followup = hb_mod@beta[11],
              pct_ebv_removed_lagged = hb_mod@beta[12],
              hb_defer_lagged = hb_mod@beta[13])

dropout_mod_b <- c(cons = dropout_mod$coefficients["(Intercept)"],
                   vvr = dropout_mod$coefficients["vvr"],
                   hb_defer = dropout_mod$coefficients["hb_defer"])

hb_defer_mod_b <- c(cons = hb_defer_mod$coefficients["(Intercept)"])

other_defer_mod_b <- c(cons = other_defer_mod$coefficients["(Intercept)"])

vvr_mod_b <- c(cons = vvr_lr_mod$coefficients["(Intercept)"],
               pct_ebv_removed = vvr_lr_mod$coefficients["pct_ebv_removed"],
               agePulse = vvr_lr_mod$coefficients["agePulse"],
               ethnicDesc_bin = vvr_lr_mod$coefficients["ethnicDesc_bin1"],
               site_type = vvr_lr_mod$coefficients["site_type"],
               bmi_catobesity = vvr_lr_mod$coefficients["bmi_catobesity"],
               bmi_catoverweight = vvr_lr_mod$coefficients["bmi_catoverweight"],
               bmi_catunderweight = vvr_lr_mod$coefficients["bmi_catunderweight"])

fpmod_mod_b <- c(gamma0 = fpmod$coefficients["gamma0"],
                 gamma1 = fpmod$coefficients["gamma1"],
                 gamma2 = fpmod$coefficients["gamma2"],
                 gamma3 = fpmod$coefficients["gamma3"],
                 gamma4 = fpmod$coefficients["gamma4"],
                 gamma5 = fpmod$coefficients["gamma5"],
                 agePulse = fpmod$coefficients["agePulse"],
                 `ABORHA+` = fpmod$coefficients["ABORHA+"],
                 `ABORHO+` = fpmod$coefficients["ABORHO+"],
                 `ABORHO-` = fpmod$coefficients["ABORHO-"],
                 `ABORHB+` = fpmod$coefficients["ABORHB+"],
                 `ABORHB-` = fpmod$coefficients["ABORHB-"],
                 `ABORHAB+` = fpmod$coefficients["ABORHAB+"],
                 `ABORHAB-` = fpmod$coefficients["ABORHAB-"],
                 ethnicDesc_bin = fpmod$coefficients["ethnicDesc_bin1"],
                 nLowHB_hx2_grouped1 = fpmod$coefficients["nLowHB_hx2_grouped1"],
                 nLowHB_hx2_grouped2 = fpmod$coefficients["nLowHB_hx2_grouped2"],
                 nDonation_hx2_grouped1 = fpmod$coefficients["nDonation_hx2_grouped1"],
                 nDonation_hx2_grouped2 = fpmod$coefficients["nDonation_hx2_grouped2"],
                 HGB_g_dL_bl = fpmod$coefficients["HGB_g_dL_bl"],
                 vvr_lagged = fpmod$coefficients["vvr_lagged"],
                 hb_defer_lagged = fpmod$coefficients["hb_defer_lagged"])

psa.b.params = list("hb" = hb_mod_b,
                    "dropout" = dropout_mod_b,
                    "defer" = hb_defer_mod_b,
                    "othdefer" = other_defer_mod_b,
                    "vvr" = vvr_mod_b,
                    "fpmod" = fpmod_mod_b
)

# collating variances from above regression submodels
hb_mod_v <- vcov(hb_mod)
dropout_mod_v <- vcov(dropout_mod)
hb_defer_mod_v <- vcov(hb_defer_mod)
other_defer_mod_v <- vcov(other_defer_mod)
vvr_mod_v <- vcov(vvr_lr_mod)
fpmod_mod_v <- vcov(fpmod)

psa.v.params = list("hb" = hb_mod_v,
                    "dropout" = dropout_mod_v,
                    "defer" = hb_defer_mod_v,
                    "othdefer" = other_defer_mod_v,
                    "vvr" = vvr_mod_v,
                    "fpmod" = fpmod_mod_v
)

# drawing 400 regression parameters from multivariate normal distributions for probabilistic sensitivity analyses
psa.hb <- psa.des(B = B, 
                  params.b = psa.b.params,
                  params.v = psa.v.params)

psa.hb.hb <- psa.hb[[2]]
psa.hb.dropout <- psa.hb[[3]]
psa.hb.defer <- psa.hb[[4]]
psa.hb.othdefer <- psa.hb[[5]]
psa.hb.vvr <- psa.hb[[6]]
psa.hb.fpmod <- psa.hb[[7]]

# ensuring that random numbers are generated in the same way for each run of the DES function to ensure reproducibility - lines 430-431 and 433-437 coded with assistance from ChatGPT 
RNGkind("L'Ecuyer-CMRG")
set.seed(9)

streams <- vector("list", NUI)
for (i in seq_len(NUI)) {
  streams[[i]] <- .Random.seed
  .Random.seed <- nextRNGStream(.Random.seed)
}

# running DES models for probabilistic sensitivity analyses
psa_outcome_set <- lapply(1:B,
                          des_boot_new_sens_bmi,
                          N = NUI,
                          donorsex = donorsex,
                          group = strattype,
                          vol = stratvol,
                          ind_data = data %>%
                            drop_na(ethnicDesc_bin, nLowHB_hx2) %>%
                            group_by(identifier) %>%
                            arrange(donationDate, .by_group = T) %>%
                            filter(donationDate == min(donationDate)) %>%
                            ungroup() %>%
                            setDT(),
                          attend = attend.dist_boot,
                          timehorizon = timehorizon,
                          ranef.hb = ranef.hb,
                          formula_hb = formula_hb,
                          psa.hb.hb = psa.hb.hb,
                          psa.hb.dropout = psa.hb.dropout,
                          psa.hb.defer = psa.hb.defer,
                          psa.hb.othdefer = psa.hb.othdefer,
                          psa.hb.vvr = psa.hb.vvr,
                          psa.hb.fpmod = psa.hb.fpmod,
                          thresh = thresh,
                          curr.recall = curr.recall,
                          lowhbdef.recall = lowhbdef.recall,
                          vlowhbdef.recall = vlowhbdef.recall,
                          othdef.recall = othdef.recall,
                          streams = streams)

# extracting outcomes from main analyses
total_over_threshold_donations <- extract_outcomes_des(outcome_name = "donate", donorsex = donorsex)
total_under_threshold_donations <- extract_outcomes_des(outcome_name = "bleed under threshold", donorsex = donorsex)
total_donations <- extract_outcomes_des(outcome_name = "donate and bleed under threshold", donorsex = donorsex)
total_vvrs <- extract_outcomes_des(outcome_name = "vvr", donorsex = donorsex)
total_collectvol <- extract_outcomes_des(outcome_name = "collectvol", donorsex = donorsex)
total_hb_deferrals <- extract_outcomes_des(outcome_name = "defer lowhb", donorsex = donorsex)
total_other_deferrals <- extract_outcomes_des(outcome_name = "defer other", donorsex = donorsex)
total_under_threshold <- extract_outcomes_des(outcome_name = "under threshold", donorsex = donorsex)
total_collectvol_under_threshold <- extract_outcomes_des(outcome_name = "collectvol under threshold", donorsex = donorsex)
total_collectvol_vvrs <- extract_outcomes_des(outcome_name = "collectvol vvr", donorsex = donorsex)
total_ebv_deferrals <- extract_outcomes_des(outcome_name = "defer low ebv", donorsex = donorsex)

# extracting outcomes from probabilistic sensitivity analyses
psa_over_threshold_donations <- unlist(lapply(1:B, extract_outcomes_psa, outcome_name = "donate", donorsex = donorsex, psa_outcome_set = psa_outcome_set))
psa_under_threshold_donations <- unlist(lapply(1:B, extract_outcomes_psa, outcome_name = "bleed under threshold", donorsex = donorsex, psa_outcome_set = psa_outcome_set))
psa_donations <- unlist(lapply(1:B, extract_outcomes_psa, outcome_name = "donate and bleed under threshold", donorsex = donorsex, psa_outcome_set = psa_outcome_set))
psa_vvrs <- unlist(lapply(1:B, extract_outcomes_psa, outcome_name = "vvr", donorsex = donorsex, psa_outcome_set = psa_outcome_set))
psa_collectvol <- unlist(lapply(1:B, extract_outcomes_psa, outcome_name = "collectvol", donorsex = donorsex, psa_outcome_set = psa_outcome_set))
psa_hb_deferrals <- unlist(lapply(1:B, extract_outcomes_psa, outcome_name = "defer lowhb", donorsex = donorsex, psa_outcome_set = psa_outcome_set))
psa_other_deferrals <- unlist(lapply(1:B, extract_outcomes_psa, outcome_name = "defer other", donorsex = donorsex, psa_outcome_set = psa_outcome_set))
psa_under_threshold <- unlist(lapply(1:B, extract_outcomes_psa, outcome_name = "under threshold", donorsex = donorsex, psa_outcome_set = psa_outcome_set))
psa_collectvol_under_threshold <- unlist(lapply(1:B, extract_outcomes_psa, outcome_name = "collectvol under threshold", donorsex = donorsex, psa_outcome_set = psa_outcome_set))
psa_collectvol_vvrs <- unlist(lapply(1:B, extract_outcomes_psa, outcome_name = "collectvol vvr", donorsex = donorsex, psa_outcome_set = psa_outcome_set))
psa_ebv_deferrals <- unlist(lapply(1:B, extract_outcomes_psa, outcome_name = "defer low ebv", donorsex = donorsex, psa_outcome_set = psa_outcome_set))

# creating dataframe with all probabilistic sensitivity analysis results
psa_results_df <- data.frame(
  psa_over_threshold_donations = psa_over_threshold_donations,
  psa_under_threshold_donations = psa_under_threshold_donations,
  psa_donations = psa_donations,
  psa_vvrs = psa_vvrs,
  psa_collectvol = psa_collectvol,
  psa_hb_deferrals = psa_hb_deferrals,
  psa_other_deferrals = psa_other_deferrals,
  psa_under_threshold = psa_under_threshold,
  psa_collectvol_under_threshold = psa_collectvol_under_threshold,
  psa_collectvol_vvrs = psa_collectvol_vvrs,
  psa_ebv_deferrals = psa_ebv_deferrals
)

# writing above dataframe to CSV
write.csv(psa_results_df, file = paste0("sens_bmi_psa_results_", strattype, "_", stratvol, "_", ifelse(donorsex == 1, "m", "f"), ".csv"))

# generating uncertainty intervals from probabilistic sensitivity analysis outcomes
extracted_ci_over_threshold_donations <- quantile(psa_over_threshold_donations, probs = c(0.025, 0.975))
extracted_ci_under_threshold_donations <- quantile(psa_under_threshold_donations, probs = c(0.025, 0.975))
extracted_ci_donations <- quantile(psa_donations, probs = c(0.025, 0.975))
extracted_ci_vvrs <- quantile(psa_vvrs, probs = c(0.025, 0.975))
extracted_ci_collectvol <- quantile(psa_collectvol, probs = c(0.025, 0.975))
extracted_ci_hb_deferrals <- quantile(psa_hb_deferrals, probs = c(0.025, 0.975))
extracted_ci_other_deferrals <- quantile(psa_other_deferrals, probs = c(0.025, 0.975))
extracted_ci_under_threshold <- quantile(psa_under_threshold, probs = c(0.025, 0.975))
extracted_ci_collectvol_under_threshold <- quantile(psa_collectvol_under_threshold, probs = c(0.025, 0.975))
extracted_ci_collectvol_vvrs <- quantile(psa_collectvol_vvrs, probs = c(0.025, 0.975))
extracted_ci_ebv_deferrals <- quantile(psa_ebv_deferrals, probs = c(0.025, 0.975))

# extracting uncertainty intervals
lci_over_threshold_donations <- extracted_ci_over_threshold_donations[1]
lci_under_threshold_donations <- extracted_ci_under_threshold_donations[1]
lci_donations <- extracted_ci_donations[1]
lci_vvrs <- extracted_ci_vvrs[1]
lci_hb_deferrals <- extracted_ci_hb_deferrals[1]
lci_other_deferrals <- extracted_ci_other_deferrals[1]
lci_under_threshold <- extracted_ci_under_threshold[1]
lci_collectvol <- extracted_ci_collectvol[1]
lci_collectvol_under_threshold <- extracted_ci_collectvol_under_threshold[1]
lci_collectvol_vvrs <- extracted_ci_collectvol_vvrs[1]
lci_ebv_deferrals <- extracted_ci_ebv_deferrals[1]

uci_over_threshold_donations <- extracted_ci_over_threshold_donations[2]
uci_under_threshold_donations <- extracted_ci_under_threshold_donations[2]
uci_donations <- extracted_ci_donations[2]
uci_vvrs <- extracted_ci_vvrs[2]
uci_hb_deferrals <- extracted_ci_hb_deferrals[2]
uci_other_deferrals <- extracted_ci_other_deferrals[2]
uci_under_threshold <- extracted_ci_under_threshold[2]
uci_collectvol <- extracted_ci_collectvol[2]
uci_collectvol_under_threshold <- extracted_ci_collectvol_under_threshold[2]
uci_collectvol_vvrs <- extracted_ci_collectvol_vvrs[2]
uci_ebv_deferrals <- extracted_ci_ebv_deferrals[2]

# compiling outcomes and outcome uncertainty intervals into single dataframe
results_df <- data.frame(
  strattype = strattype,
  stratvol = stratvol,
  donorsex = donorsex,
  outcome_name = c("total over-threshold donations", "total under-threshold donations", "total donations", "total vvrs", "total hb deferrals", "total other deferrals", "total under threshold", "total volume collected", "total volume collected under-threshold", "total volume collected vvrs", "total ebv deferrals"), 
  outcome = c(total_over_threshold_donations, total_under_threshold_donations, total_donations, total_vvrs, total_hb_deferrals, total_other_deferrals, total_under_threshold, total_collectvol, total_collectvol_under_threshold, total_collectvol_vvrs, total_ebv_deferrals),
  outcome_lci = c(lci_over_threshold_donations, lci_under_threshold_donations, lci_donations, lci_vvrs, lci_hb_deferrals, lci_other_deferrals, lci_under_threshold, lci_collectvol, lci_collectvol_under_threshold, lci_collectvol_vvrs, lci_ebv_deferrals),
  outcome_uci = c(uci_over_threshold_donations, uci_under_threshold_donations, uci_donations, uci_vvrs, uci_hb_deferrals, uci_other_deferrals, uci_under_threshold, uci_collectvol, uci_collectvol_under_threshold, uci_collectvol_vvrs, uci_ebv_deferrals)
)

outdir <- R.utils::commandArgs()[23] 

# writing dataframe as CSV
write.csv(results_df, file = paste0("results_sens_bmi_", strattype, "_", stratvol, "_", ifelse(donorsex == 1, "m", "f"), ".csv"))