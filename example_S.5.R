#-----------------------------------------------------------
# Script: Example in Section 5
#-----------------------------------------------------------

source("0_util.R")

# ---- Constants / Planning ------------------------------------------------
pi1 <- 0.41
pi2 <- 0.6
cens_rate <- 0.025

# Information rate -> interim after 12 months
design <- getDesignInverseNormal(
  kMax = 2,
  typeOfDesign = "P",
  informationRates = c(.153, 1),
  futilityBounds = 0,
  bindingFutility = TRUE
)
boundaries <- round(design$criticalValues, 3)
w1 <- sqrt(design$informationRates[1])
alpha1 <- design$alphaSpent[1]

# sample size (survival)
sample_size <- getSampleSizeSurvival(
  design = design,
  accrualTime = c(0, 24),
  dropoutTime = 12,
  dropoutRate1 = cens_rate,
  dropoutRate2 = cens_rate,
  eventTime = 12,
  pi1 = pi1,
  pi2 = pi2,
  followUpTime = 48
)

# ----- rates from quantiles ------------------------------------------------
lambda1 <- uniroot(function(x) pexp(12, x) - pi1, interval = c(.001, 10))$root
lambda2 <- uniroot(function(x) pexp(12, x) - pi2, interval = c(.001, 10))$root
lambda_cens <- uniroot(function(x) pexp(12, x) - cens_rate, interval = c(.00001, 10))$root

# ensure total sample size is even (preserves original behaviour)
n_total_even <- ceiling(sample_size$maxNumberOfSubjects) + (ceiling(sample_size$maxNumberOfSubjects) %% 2)
n_per_group <- n_total_even / 2

# expected events at interim (both arms)
Dl <- expected_events_exponential(
  overall_sample_size = n_per_group,
  time = 12,
  rate_event = lambda1,
  rate_censoring = lambda_cens,
  end_recruitment_time = 24
) + expected_events_exponential(
  overall_sample_size = n_per_group,
  time = 12,
  rate_event = lambda2,
  rate_censoring = lambda_cens,
  end_recruitment_time = 24
)

event_prob_interim <- (
  event_probability_exponential(n_per_group, 12, lambda1, lambda_cens, 24) +
    event_probability_exponential(n_per_group, 12, lambda2, lambda_cens, 24)
) / 2

HR <- lambda1 / lambda2

# ----- Hybrid boundaries ---------------------------------------------------
# Hybrid 1
Delta_h1 <- 0.1
df_h1 <- get_hybrid1(
  Dl = Dl,
  HR = HR,
  alpha1 = alpha1,
  Delta = Delta_h1,
  n_int = ceiling(sample_size$numberOfSubjects[1])
)
h1_l <- df_h1$h1_l
h1_u <- df_h1$h1_u

# Hybrid 2
Delta_h2 <- 0.1
h2_l <- qbinom(Delta_h2, n_per_group * 2, event_prob_interim, lower.tail = TRUE) - 1
h2_u <- qbinom(1 - Delta_h2, n_per_group * 2, event_prob_interim)

# Hybrid 3
Delta_h3 <- 0.2
h3_l <- floor(Dl - (Dl * Delta_h3))
h3_u <- ceiling(Dl + (Dl * Delta_h3))

# Design setup ------------------------------------------------------------

###### fixed parameters
fixed_par <- list(t_int = 12, min_fu = 48, sec_rec_length = 12) ## equal for all designs

###### fixed design ##############################
design_name <- "fixed"
fixed <- create_hybrid_design(
  p_cfix = qnorm(0.975), p_c1 = NULL, p_c2 = NULL,
  p_f1 = 0, p_du = 10^6, p_dl = -1, p_fixPar = fixed_par, p_w1 = 0
)
saveRDS(fixed, paste0("Example/Designs/", design_name, ".RDS"))

###### Hybrid Design 1 ###########################
design_name <- "hybrid1"
hybrid1 <- create_hybrid_design(
  p_cfix = qnorm(0.975), p_c1 = boundaries[1], p_c2 = boundaries[2],
  p_f1 = 0, p_du = h1_u, p_dl = h1_l, p_fixPar = fixed_par, p_w1 = w1
)
saveRDS(hybrid1, paste0("Example/Designs/", design_name, ".RDS"))

###### Hybrid Design 2 ###########################
design_name <- "hybrid2"
hybrid2 <- create_hybrid_design(
  p_cfix = qnorm(0.975), p_c1 = boundaries[1], p_c2 = boundaries[2],
  p_f1 = 0, p_du = h2_u, p_dl = h2_l, p_fixPar = fixed_par, p_w1 = w1
)
saveRDS(hybrid2, paste0("Example/Designs/", design_name, ".RDS"))

###### Hybrid Design 3 ###########################
design_name <- "hybrid3"
hybrid3 <- create_hybrid_design(
  p_cfix = qnorm(0.975), p_c1 = boundaries[1], p_c2 = boundaries[2],
  p_f1 = 0, p_du = h3_u, p_dl = h3_l, p_fixPar = fixed_par, p_w1 = w1
)
saveRDS(hybrid3, paste0("Example/Designs/", design_name, ".RDS"))

###### Group Sequential ###########################
design_name <- "gs_design"
gs_design <- create_hybrid_design(
  p_cfix = qnorm(0.975), p_c1 = boundaries[1], p_c2 = boundaries[2],
  p_f1 = 0, p_du = -1, p_dl = -1, p_fixPar = fixed_par, p_w1 = w1
)
saveRDS(gs_design, paste0("Example/Designs/", design_name, ".RDS"))


# Data simulation ---------------------------------------------------------
## assumed parameters
surv_parameters <- list()
surv_parameters[["hr"]] <- data.frame(medC_values = qexp(.5, lambda2), medI_values = qexp(.5, lambda1))

######## generate data sets ############################
## simulation parameters; Just nsim=1 for an exemplary study
n_sim <- 1
t_max_acc <- 100
acc_rate <- 152 / 24
min_fu <- 48
seed <- 9122024

for (range_param in c("hr")) {
  for (i_param in 1:length(surv_parameters[[range_param]]$medI_values)) {
    Evpar_C <- list(med_surv = surv_parameters[[range_param]]$medC_values[i_param], shape = 1, lambda_cens = 0.0008403983)
    Evpar_I <- list(med_surv = surv_parameters[[range_param]]$medI_values[i_param], shape = 1, lambda_cens = 0.0008403983)
    data_name <- paste0(range_param, i_param)

    ## simulate patient data
    data_all <- NULL
    for (i_sim in 1:n_sim) {
      ## simulate survival data
      set.seed(seed + i_sim - 1)
      data <- simulate_stage(p_acc_rate = acc_rate, p_max_acc_time = t_max_acc, p_min_fu = min_fu, p_EvparI = Evpar_I, p_EvparC = Evpar_C)

      ## add simulation run index
      data$sim_idx <- i_sim
      # add results to previous simulation runs
      data_all %<>% rbind(data)

      ## position in simulation
      if (i_sim %% 1000 == 0) {
        print(paste0(range_param, " simulation run ", i_sim, " of ", n_sim))
      }
    }
    ## save patient data
    fwrite(data_all, paste0("Example/Data/", data_name, ".csv"))

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
    ), paste0("Example/Data/", data_name, "_parameters", ".RDS"))
  }
}


# Evaluation with the different designs -----------------------------------

## simulation parameters
ev_times_names <- c("hr1")
design_names <- c("fixed", "gs_design", "hybrid1", "hybrid2", "hybrid3")

for (i_data in 1:length(ev_times_names)) {
  ## read event data and survival parameters and obtain number of simulation runs
  data_all <- fread(paste0("Example/Data/", ev_times_names[i_data], ".csv"))
  surv_param <- readRDS(paste0("Example/Data/", ev_times_names[i_data], "_parameters.RDS"))
  n_sim <- max(data_all$sim_idx)

  for (i_design in 1:length(design_names)) {
    ## read design
    design_eval <- readRDS(paste0("Example/Designs/", design_names[i_design], ".RDS"))

    ## calculate summary statistics for
    int_stats_all <- NULL
    sec_stats_all <- NULL
    for (i_sim in 1:n_sim) {
      ## cut-off data for simulation run
      data <- data_all[sim_idx == i_sim, ]

      ## interim analysis at fixed calendar time
      data_int <- get_results_at_time(p_data_res = data, p_t = design_eval$t_int, p_end_rec = design_eval$t_int)
      survdiff_int <- survdiff(Surv(t, status) ~ group, data = data_int)
      res_logrank_int <- logrank_func(survdiff_int)
      cox_int <- coxph(Surv(t, status) ~ group, data = data_int)
      int_stats <- data.frame(
        sim_idx = i_sim, n_pat = nrow(data_int), n_ev = sum(data_int$status), perc_evI = mean(data_int[group == "I", ]$status),
        perc_evC = mean(data_int[group == "C", ]$status), lr_stat = res_logrank_int$log_rank_stat, hazard_ratio = exp(cox_int$coefficients),
        trial_dur = unique(data_int$time_eval), unblinded = design_eval$is_unblinded(sum(data_int$status))
      )

      ## second stage
      sec_rec_length <- design_eval$dec_rec_dur(p_nEv = int_stats$n_ev, p_logrank_int = res_logrank_int$log_rank_stat) # length of second stage recruitement period
      data_sec <- get_results_at_time(p_data_res = data, p_t = design_eval$t_int + sec_rec_length + design_eval$min_fu, p_end_rec = design_eval$t_int + sec_rec_length)
      survdiff_sec <- survdiff(Surv(t, status) ~ group, data = data_sec)
      res_logrank_sec <- logrank_func(survdiff_sec)
      res_inv_norm <- indep_inc_func(p_lr_fin = res_logrank_sec, p_lr_int = res_logrank_int, p_w1 = design_eval$w1)
      cox_sec <- coxph(Surv(t, status) ~ group, data = data_sec)
      surv_group <- summary(survfit(Surv(t, status) ~ group, data = data_sec), times = 12)
      sec_stats <- data.frame(
        sim_idx = i_sim, n_pat = nrow(data_sec), n_ev = sum(data_sec$status), perc_evI = mean(data_sec[group == "I", ]$status),
        perc_evC = mean(data_sec[group == "C", ]$status), lr_stat = res_logrank_sec$log_rank_stat, inv_norm = res_inv_norm$invN_stat,
        hazard_ratio = exp(cox_sec$coefficients), trial_dur = unique(data_sec$time_eval), km_surv_1year_I = surv_group$surv[surv_group$strata == "group=I"],
        km_surv_1year_C = surv_group$surv[surv_group$strata == "group=C"]
      )

      # add results to previous simulation runs
      int_stats_all %<>% rbind(int_stats)
      sec_stats_all %<>% rbind(sec_stats)

      ## position in simulation
      if (i_sim %% 1000 == 0) {
        print(paste0("Simulation run ", i_sim, " of ", n_sim))
      }
    }

    ## join first and second stage results
    names(int_stats_all)[which(!names(int_stats_all) %in% c("sim_idx"))] <- paste0("int_", names(int_stats_all)[which(!names(int_stats_all) %in% c("sim_idx"))])
    names(sec_stats_all)[which(!names(sec_stats_all) %in% c("sim_idx"))] <- paste0("sec_", names(sec_stats_all)[which(!names(sec_stats_all) %in% c("sim_idx"))])
    stats_all <- join(int_stats_all, sec_stats_all, by = "sim_idx")
    stats_all %<>% as.data.table

    ## rejection decisions
    stats_all <- design_eval$rej_func(stats_all)

    ## index of extreme cases (no/all patients have events at interim or at the end)
    idx_cap <- which((stats_all$int_unblinded & (stats_all$int_perc_evC %in% c(0, 1) | stats_all$int_perc_evI %in% c(0, 1))) | stats_all$sec_perc_evC %in% c(0, 1) | stats_all$sec_perc_evI %in% c(0, 1))

    ## calculate performance measures;
    performance <- data.frame(
      av_pat_int = mean(stats_all$int_n_pat),
      av_pat_total = mean(stats_all$sec_n_pat),
      av_events_int = mean(stats_all$int_n_ev),
      av_events_total = mean(stats_all$sec_n_ev),
      prob_effic = mean(stats_all$int_rej),
      power = mean(stats_all$int_rej) + mean((stats_all$int_rej == 0 & (stats_all$int_n_pat != stats_all$sec_n_pat)) & stats_all$sec_rej == 1),
      prob_fut = mean(stats_all$int_n_pat == stats_all$sec_n_pat & stats_all$int_rej == 0),
      med_surv_control = surv_param$Evpar_C$med_surv,
      one_year_event_prob_C = surv_param$one_year_event_prob_C,
      one_year_event_prob_I = med_surv_to_one_year_event_prob(surv_param$Evpar_I$med_surv),
      av_trial_dur = mean(stats_all$sec_trial_dur),
      hazard_ratio = surv_param$hazard_ratio,
      av_est_hazard_ratio = mean(stats_all$sec_hazard_ratio),
      av_est_log_hazard_ratio = mean(log(stats_all$sec_hazard_ratio)),
      av_est_hazard_ratio_cap = mean(stats_all[-idx_cap, ]$sec_hazard_ratio),
      av_est_log_hazard_ratio_cap = mean(log(stats_all[-idx_cap, ]$sec_hazard_ratio)),
      prob_cap = length(idx_cap) / nrow(stats_all),
      prob_unblind = mean(stats_all$int_unblinded),
      av_est_km_1year_survI = mean(stats_all$sec_km_surv_1year_I),
      av_est_km_1year_survC = mean(stats_all$sec_km_surv_1year_C)
    )
    ## save results
    path_res <- paste0("Example/Results/", ev_times_names[i_data])
    if (!file.exists(path_res)) {
      dir.create(path_res)
    }
    path_res <- paste0(path_res, "/", design_names[i_design])
    dir.create(path_res)
    fwrite(stats_all, paste0(path_res, "/stats_all.csv"))
    saveRDS(performance, paste0(path_res, "/performance.RDS"))
  }
}


# Summary -----------------------------------------------------------------

ev_times_names <- "hr1"
design_names <- c("fixed", "gs_design", "hybrid1", "hybrid2", "hybrid3")

res <- data.frame()

for (i_design in 1:length(design_names)) {
  ## define design
  design_name <- design_names[i_design]
  path_res <- paste0("Example/Results/", ev_times_names, "/", design_name, "/stats_all.csv")
  data <- fread(path_res)
  data$design <- design_name
  res <- rbind(res, data)
}

print(
  res %>%
    select(design, int_n_pat, int_n_ev, int_unblinded, int_lr_stat, int_hazard_ratio, sec_n_pat, sec_n_ev, sec_lr_stat)
)
