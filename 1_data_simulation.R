###############################################################################
# THIS SCRIPT SIMULATES AND SAVES THE DATA FOR TYPE-1 ERROR RATE AND          #
# PERFORMANCE EVALUATION STUDY                                                #
###############################################################################

source("0_util.R")

# Type-1 Error Rate and Bias ------------------------------------------------

n_values <- c(180, 360, 720, 1080)
#n_values <- 1080
p1year_values <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)

system.time({
for (n in n_values) {
  ## assumed parameters
  Evpar_C <- list(med_surv = 24, shape = 1, lambda_cens = 0.004288308)
  Evpar_I <- list(med_surv = 24, shape = 1, lambda_cens = 0.004288308)
  hr0 <- calc_hazard_ratio(p_med_survI = Evpar_I$med_surv, p_med_survC = Evpar_C$med_surv)
  p1yearC0 <- med_surv_to_one_year_event_prob(Evpar_C$med_surv)

  ## Median survival time for a range of 1-year control group event probability values
  medC_values <- medI_values <- c()
  for (p in p1year_values) {
    med_survC <- one_year_event_prob_to_med_surv(p)
    p1yearI <- calc_one_year_event_prob_I(p_hazard_ratio = hr0, p_1yearprob = p)
    med_survI <- one_year_event_prob_to_med_surv(p1yearI)
    medC_values <- c(medC_values, med_survC)
    medI_values <- c(medI_values, med_survI)
  }

  surv_parameters <- list(evProb = data.frame(medC_values, medI_values))

  ## simulation settings
  n_sim <- 10
  t_max_acc <- 60
  acc_rate <- n / 36
  min_fu <- 12
  seed <- 23012025

  ## loop over scenarios
  for (range_param in c("evProb")) {
    for (i_param in 1:length(surv_parameters[[range_param]]$medI_values)) {
      Evpar_C <- list(med_surv = surv_parameters[[range_param]]$medC_values[i_param], shape = 1, lambda_cens = 0.004288308)
      Evpar_I <- list(med_surv = surv_parameters[[range_param]]$medI_values[i_param], shape = 1, lambda_cens = 0.004288308)
      data_name <- paste0(range_param, i_param)

      ## simulate patient data
      data_all <- NULL
      for (i_sim in 1:n_sim) {
        set.seed(seed + i_sim - 1)
        data <- simulate_stage(
          p_acc_rate = acc_rate,
          p_max_acc_time = t_max_acc,
          p_min_fu = min_fu,
          p_EvparI = Evpar_I,
          p_EvparC = Evpar_C
        )
        data$sim_idx <- i_sim
        data_all %<>% rbind(data)

        if (i_sim %% 1000 == 0) {
          print(paste0("n=", n, " ", range_param, " simulation run ", i_sim, " of ", n_sim))
        }
      }

      ## save patient data
      fwrite(data_all, paste0("./Simulation Study/Data/Type-1 Error Rate/n=", n, "/", data_name, ".csv"))

      ## save setting parameters
      saveRDS(list(
        n_sim = n_sim,
        seed = seed + i_sim - 1,
        t_max_acc = t_max_acc,
        acc_rate = acc_rate,
        min_fu = min_fu,
        Evpar_C = Evpar_C,
        Evpar_I = Evpar_I,
        hazard_ratio = calc_hazard_ratio(p_med_survI = Evpar_I$med_surv, p_med_survC = Evpar_C$med_surv),
        one_year_event_prob_C = med_surv_to_one_year_event_prob(Evpar_C$med_surv)
      ), paste0("./Simulation Study/Data/Type-1 Error Rate/n=", n, "/", data_name, "_parameters.RDS"))
    }
  }
}
})

# Performance Study ---------------------------------------------------------

study_setups <- list(
  list(
    name = "Study 1_Pocock",
    Evpar_C_med = 37.28113,
    Evpar_I_med = 78.94791,
    p1year_values = c(0.05, 0.1, 0.2, 0.4, 0.6),
    acc_rate = 208 / 36,
    seed = 23012025
  ),
  list(
    name = "Study 2_Pocock",
    Evpar_C_med = 9.08047,
    Evpar_I_med = 16.28234,
    p1year_values = c(0.1, 0.2, 0.4, 0.6, 0.7, 0.8),
    acc_rate = 138 / 36,
    seed = 23012025
  ),
  list(
    name = "Study 1_OBrien-Fleming",
    Evpar_C_med = 37.28113,
    Evpar_I_med = 78.94791,
    p1year_values = c(0.05, 0.1, 0.2, 0.4, 0.6),
    acc_rate = 184 / 36,
    seed = 13082025
  ),
  list(
    name = "Study 2_OBrien-Fleming",
    Evpar_C_med = 9.08047,
    Evpar_I_med = 16.28234,
    p1year_values = c(0.1, 0.2, 0.4, 0.6, 0.7, 0.8),
    acc_rate = 122 / 36,
    seed = 13082025
  )
)

## Constants
n_sim <- 10
t_max_acc <- 60
min_fu <- 12
lambda_cens <- 0.0008403983 # corresponds t0 1% 1-year censoring probability

## Main loop
for (setup in study_setups) {
  ## assumed parameters
  Evpar_C <- list(med_surv = setup$Evpar_C_med, shape = 1, lambda_cens = lambda_cens)
  Evpar_I <- list(med_surv = setup$Evpar_I_med, shape = 1, lambda_cens = lambda_cens)
  hr0 <- calc_hazard_ratio(p_med_survI = Evpar_I$med_surv, p_med_survC = Evpar_C$med_surv)
  p1yearC0 <- med_surv_to_one_year_event_prob(Evpar_C$med_surv)

  ## --- survival parameters: event probabilities ---
  medC_values <- medI_values <- c()
  for (p in setup$p1year_values) {
    med_survC <- one_year_event_prob_to_med_surv(p)
    p1yearI <- calc_one_year_event_prob_I(p_hazard_ratio = hr0, p_1yearprob = p)
    med_survI <- one_year_event_prob_to_med_surv(p1yearI)
    medC_values <- c(medC_values, med_survC)
    medI_values <- c(medI_values, med_survI)
  }
  surv_parameters <- list(evProb = data.frame(medC_values, medI_values))

  ## --- survival parameters: hazard ratios ---
  hazard_ratio <- c(seq(0.2, 0.5, 0.1), hr0, seq(0.6, 1.1, 0.1))
  medC_values <- medI_values <- c()
  for (hr in hazard_ratio) {
    med_survC <- Evpar_C$med_surv
    p1yearI <- calc_one_year_event_prob_I(p_hazard_ratio = hr, p_1yearprob = p1yearC0)
    med_survI <- one_year_event_prob_to_med_surv(p1yearI)
    medC_values <- c(medC_values, med_survC)
    medI_values <- c(medI_values, med_survI)
  }
  surv_parameters[["hr"]] <- data.frame(medC_values, medI_values)

  ## --- simulation ---
  for (range_param in c("evProb", "hr")) {
    for (i_param in 1:length(surv_parameters[[range_param]]$medI_values)) {
      Evpar_C <- list(med_surv = surv_parameters[[range_param]]$medC_values[i_param], shape = 1, lambda_cens = lambda_cens)
      Evpar_I <- list(med_surv = surv_parameters[[range_param]]$medI_values[i_param], shape = 1, lambda_cens = lambda_cens)
      data_name <- paste0(range_param, i_param)

      ## simulate data
      data_all <- NULL
      for (i_sim in 1:n_sim) {
        set.seed(setup$seed + i_sim - 1)
        data <- simulate_stage(
          p_acc_rate = setup$acc_rate, p_max_acc_time = t_max_acc,
          p_min_fu = min_fu, p_EvparI = Evpar_I, p_EvparC = Evpar_C
        )
        data$sim_idx <- i_sim
        data_all %<>% rbind(data)

        if (i_sim %% 1000 == 0) {
          print(paste0(setup$name, " - ", range_param, " run ", i_sim, " of ", n_sim))
        }
      }

      ## save outputs
      fwrite(data_all, paste0("./Simulation Study/Data/Performance/", setup$name, "/", data_name, ".csv"))
      saveRDS(
        list(
          n_sim = n_sim,
          seed = setup$seed + i_sim - 1,
          t_max_acc = t_max_acc,
          acc_rate = setup$acc_rate,
          min_fu = min_fu,
          Evpar_C = Evpar_C,
          Evpar_I = Evpar_I,
          hazard_ratio = calc_hazard_ratio(p_med_survI = Evpar_I$med_surv, p_med_survC = Evpar_C$med_surv),
          one_year_event_prob_C = med_surv_to_one_year_event_prob(Evpar_C$med_surv)
        ),
        paste0("./Simulation Study/Data/Performance/", setup$name, "/", data_name, "_parameters.RDS")
      )
    }
  }
}
