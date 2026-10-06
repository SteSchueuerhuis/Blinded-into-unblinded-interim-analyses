###############################################################################
# REPRODUCING STIMULATION STUDY 1 (POCOCK BOUNDARIES)                         #
###############################################################################

source("0_util.R")

# Data Simulation ---------------------------------------------------------

## assumed parameters
Evpar_C <- list(med_surv = 37.28113, shape = 1, lambda_cens = 0.0008403983)
Evpar_I <- list(med_surv = 78.94791, shape = 1, lambda_cens = 0.0008403983)
hr0 <- calc_hazard_ratio(p_med_survI = Evpar_I$med_surv, p_med_survC = Evpar_C$med_surv)
p1yearC0 <- med_surv_to_one_year_event_prob(Evpar_C$med_surv)

####### Median survival time for a range of 1-year control group event probability values
p1year_values <- c(0.05, 0.1, 0.2, 0.4, 0.6)

medC_values <- c()
medI_values <- c()
for (i in 1:length(p1year_values)) {
  med_survC <- one_year_event_prob_to_med_surv(p1year_values[i])
  p1yearI <- calc_one_year_event_prob_I(p_hazard_ratio = hr0, p_1yearprob = p1year_values[i])
  med_survI <- one_year_event_prob_to_med_surv(p1yearI)
  medC_values <- c(medC_values, med_survC)
  medI_values <- c(medI_values, med_survI)
}
surv_parameters <- list(evProb = data.frame(medC_values, medI_values))

## Median survival time for a range of hazard ratios
hazard_ratio <- c(seq(0.2, 0.5, 0.1), hr0, seq(0.6, 1.1, 0.1))
medC_values <- c()
medI_values <- c()
for (i in 1:length(hazard_ratio)) {
  med_survC <- Evpar_C$med_surv
  p1yearI <- calc_one_year_event_prob_I(p_hazard_ratio = hazard_ratio[i], p_1yearprob = p1yearC0)
  med_survI <- one_year_event_prob_to_med_surv(p1yearI)
  medC_values <- c(medC_values, med_survC)
  medI_values <- c(medI_values, med_survI)
}
surv_parameters[["hr"]] <- data.frame(medC_values, medI_values)

######## generate data sets ############################
## simulation parameters
n_sim <- 10000
t_max_acc <- 60
acc_rate <- 208 / 36
min_fu <- 12 ## minimum follow-up
seed <- 23012025

for (range_param in c("evProb", "hr")) {
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
    fwrite(data_all, paste0("./Reproduction Example/Data/", data_name, ".csv"))

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
    ), paste0("./Reproduction Example/Data/", data_name, "_parameters", ".RDS"))
  }
}

# Create Designs ----------------------------------------------------------

###### fixed parameters
fixed_par <- list(t_int = 24, min_fu = 12, sec_rec_length = 12) ## equal for all designs

# Fixed
design_name <- "fixed"
fixed <- create_hybrid_design(
  p_cfix = qnorm(0.975), p_c1 = NULL, p_c2 = NULL,
  p_f1 = 0, p_du = 10^6, p_dl = -1, p_fixPar = fixed_par, p_w1 = 0
)
saveRDS(fixed, paste0("./Reproduction Example/Designs/", design_name, ".RDS"))

# Hybrid Design 1
design_name <- "hybrid1"
w1 <- sqrt(0.3)
hybrid1 <- create_hybrid_design(
  p_cfix = qnorm(0.975), p_c1 = 2.195, p_c2 = 2.195,
  p_f1 = 0, p_du = 28, p_dl = 13, p_fixPar = fixed_par, p_w1 = w1
)
saveRDS(hybrid1, paste0("./Reproduction Example/Designs/", design_name, ".RDS"))

# Hybrid Design 2
design_name <- "hybrid2"
w1 <- sqrt(0.3)
hybrid2 <- create_hybrid_design(
  p_cfix = qnorm(0.975), p_c1 = 2.195, p_c2 = 2.195,
  p_f1 = 0, p_du = 26, p_dl = 14, p_fixPar = fixed_par, p_w1 = w1
)
saveRDS(hybrid2, paste0("./Reproduction Example/Designs/", design_name, ".RDS"))

# Hybrid Design 3
design_name <- "hybrid3"
w1 <- sqrt(0.3)
hybrid3 <- create_hybrid_design(
  p_cfix = qnorm(0.975), p_c1 = 2.195, p_c2 = 2.195,
  p_f1 = 0, p_du = 25, p_dl = 16, p_fixPar = fixed_par, p_w1 = w1
)
saveRDS(hybrid3, paste0("./Reproduction Example/Designs/", design_name, ".RDS"))

# Group-Sequential
design_name <- "gs_design"
w1 <- sqrt(0.3)
gs_design <- create_hybrid_design(
  p_cfix = qnorm(0.975), p_c1 = 2.195, p_c2 = 2.195,
  p_f1 = 0, p_du = -1, p_dl = -1, p_fixPar = fixed_par, p_w1 = w1
)
saveRDS(gs_design, paste0("./Reproduction Example/Designs/", design_name, ".RDS"))


# Evaluation --------------------------------------------------------------

## simulation parameters
ev_times_names <- c(paste0("evProb", 1:5), paste0("hr", 1:11))
design_names <- c("fixed", "gs_design", "hybrid1", "hybrid2", "hybrid3")

for (i_data in 1:length(ev_times_names)) {
  ## read event data and survival parameters and obtain number of simulation runs
  data_all <- fread(paste0("./Reproduction Example/Data/", ev_times_names[i_data], ".csv"))
  surv_param <- readRDS(paste0("./Reproduction Example/Data/", ev_times_names[i_data], "_parameters.RDS"))
  n_sim <- max(data_all$sim_idx)

  for (i_design in 1:length(design_names)) {
    ## read design
    design_eval <- readRDS(paste0("./Reproduction Example/Designs/", design_names[i_design], ".RDS"))

    ## calculate summary statistics for
    int_stats_all <- NULL
    sec_stats_all <- NULL
    for (i_sim in 1:n_sim) {
      ## cut-off data for simulation run
      data <- data_all[sim_idx == i_sim, ]

      ## interim analysis at fixed calendar time // + design_eval$t_follow_up
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
      data_sec$time_eval <- ifelse(sec_rec_length == 0, design_eval$t_int, data_sec$time_eval)
      survdiff_sec <- survdiff(Surv(t, status) ~ group, data = data_sec)
      res_logrank_sec <- logrank_func(survdiff_sec)
      res_inv_norm <- indep_inc_func(p_lr_fin = res_logrank_sec, p_lr_int = res_logrank_int, p_w1 = design_eval$w1)
      cox_sec <- coxph(Surv(t, status) ~ group, data = data_sec)
      surv_group <- summary(survfit(Surv(t, status) ~ group, data = data_sec), times = 12)
      sec_stats <- data.frame(
        sim_idx = i_sim, n_pat = nrow(data_sec), n_ev = sum(data_sec$status), perc_evI = mean(data_sec[group == "I", ]$status),
        perc_evC = mean(data_sec[group == "C", ]$status), lr_stat = res_logrank_sec$log_rank_stat, inv_norm = res_logrank_sec$log_rank_stat,
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

    ## calculate performance measures
    performance <- data.frame(
      av_pat_int = mean(stats_all$int_n_pat),
      av_pat_total = mean(stats_all$sec_n_pat),
      av_events_int = mean(stats_all$int_n_ev),
      av_events_total = mean(stats_all$sec_n_ev),
      prob_effic = mean(stats_all$int_rej),
      power = mean(stats_all$int_rej) + mean((stats_all$int_rej == 0 & (stats_all$int_n_pat != stats_all$sec_n_pat)) & stats_all$sec_rej == 1),
      prob_fut = mean((stats_all$int_n_pat == stats_all$sec_n_pat) & stats_all$int_rej == 0),
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
    path_res <- paste0("./Reproduction Example/Results/", ev_times_names[i_data])
    if (!file.exists(path_res)) {
      dir.create(path_res)
    }
    path_res <- paste0(path_res, "/", design_names[i_design])
    dir.create(path_res)
    fwrite(stats_all, paste0(path_res, "/stats_all.csv"))
    saveRDS(performance, paste0(path_res, "/performance.RDS"))
  }
}


# Visualization -----------------------------------------------------------

ev_times_names <- c(paste0("evProb", 1:5), paste0("hr", 1:11))
design_names <- c("fixed", "hybrid1", "hybrid2", "hybrid3", "gs_design")

for (i_design in 1:length(design_names)) {
  ## define design
  design_name <- design_names[i_design]

  ## performance parameters of interest
  one_year_event_prob_C <- c()
  power <- c()
  av_trial_dur <- c()
  hazard_ratio <- c()
  perc_bias <- c()
  prob_unblind <- c()
  av_patients <- c()
  for (i_evt_set in 1:length(ev_times_names)) {
    ## read performance file
    path_res <- paste0("./Reproduction Example/Results/", ev_times_names[i_evt_set], "/", design_name)
    performance <- readRDS(paste0(path_res, "/performance.RDS"))

    ## modify entries for caped hazard ratio, when no extreme case happened (maybe not important)
    if (performance$prob_cap == 0) {
      performance$av_est_hazard_ratio_cap <- performance$av_est_hazard_ratio
      performance$av_est_log_hazard_ratio_cap <- performance$av_est_log_hazard_ratio
    }

    ## get performance parameters
    one_year_event_prob_C <- c(one_year_event_prob_C, performance$one_year_event_prob_C)
    power <- c(power, performance$power)
    av_trial_dur <- c(av_trial_dur, performance$av_trial_dur)
    hazard_ratio <- c(hazard_ratio, performance$hazard_ratio)
    perc_bias <- c(perc_bias, (performance$av_est_hazard_ratio - performance$hazard_ratio) / performance$hazard_ratio * 100)
    surv_diff_true <- (1 - performance$one_year_event_prob_I) - (1 - performance$one_year_event_prob_C)
    surv_diff_est <- performance$av_est_km_1year_survI - performance$av_est_km_1year_survC
    prob_unblind <- c(prob_unblind, performance$prob_unblind)
    av_patients <- c(av_patients, performance$av_pat_total)
  }

  assign(x = paste0(design_name, "_perf"), value = data.table(
    one_year_event_prob_C = one_year_event_prob_C,
    power = power, av_trial_dur = av_trial_dur,
    hazard_ratio = hazard_ratio, perc_bias = perc_bias,
    prob_unblind = prob_unblind, av_patients = av_patients,
    par_range = substr(ev_times_names, 1, 2)
  ))
}

################## Plot for parameters of primary interest #####################

pdf(file = "./Reproduction Example/Plots/Figure2.pdf", width = 12, height = 8)

layout(matrix(c(1, 1, 1, 2, 3, 4, 5, 5, 5, 6, 7, 8, 9, 9, 9), ncol = 3, byrow = TRUE), heights = c(0.06, 0.4, 0.06, 0.4, 0.06))
plot_colors <- c("black", "darkorange", "forestgreen", "dodgerblue2", "firebrick4")

## title
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, paste0("Performance for HR=", round(fixed_perf[par_range == "ev", ]$hazard_ratio[1], 2), " and varying values of control event probability"), cex = 2, font = 1.5)

## power plot
par(mar = c(3, 3, 2, 2))
fixed_perf[par_range == "ev", ] %$% plot(one_year_event_prob_C, power,
  type = "l", ylim = c(0, 1),
  lwd = 3, main = "Power", mgp = c(1.5, 0.5, 0), xlab = "1-year event probability in control group", col = plot_colors[1]
)
fixed_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, rep(0.8, length(one_year_event_prob_C)), col = "red", lty = 3)
hybrid1_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, power, col = plot_colors[2], lty = 1, lwd = 3)
hybrid2_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, power, col = plot_colors[3], lty = 1, lwd = 3)
hybrid3_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, power, col = plot_colors[4], lty = 1, lwd = 3)
gs_design_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, power, col = plot_colors[5], lty = 1, lwd = 3)


## trial duration plot
par(mar = c(3, 3, 2, 2))
fixed_perf[par_range == "ev", ] %$% plot(one_year_event_prob_C, av_trial_dur,
  type = "l", lwd = 3, main = "Trial duration",
  ylab = "average trial duration", mgp = c(1.5, 0.5, 0), xlab = "1-year event probability in control group"
)
# fixed_exp1_perf %$% lines(one_year_event_prob_C, rep(0.8, length(one_year_event_prob_C)), col = 'red', lty = 3)
hybrid1_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, av_trial_dur, col = plot_colors[2], lty = 1, lwd = 3)
hybrid2_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, av_trial_dur, col = plot_colors[3], lty = 1, lwd = 3)
hybrid3_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, av_trial_dur, col = plot_colors[4], lty = 1, lwd = 3)
gs_design_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, av_trial_dur, col = plot_colors[5], lty = 1, lwd = 3)


## HR bias plot
par(mar = c(3, 3, 2, 2))
fixed_perf[par_range == "ev", ] %$% plot(one_year_event_prob_C, perc_bias,
  type = "l", lwd = 3, main = "HR bias",
  ylab = "% Bias", mgp = c(1.5, 0.5, 0), xlab = "1-year event probability in control group", ylim = c(-13, 13), col = plot_colors[1]
)
fixed_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, rep(0, length(one_year_event_prob_C)), col = "red", lty = 3)
hybrid1_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, perc_bias, col = plot_colors[2], lty = 1, lwd = 3)
hybrid2_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, perc_bias, col = plot_colors[3], lty = 1, lwd = 3)
hybrid3_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, perc_bias, col = plot_colors[4], lty = 1, lwd = 3)
gs_design_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, perc_bias, col = plot_colors[5], lty = 1, lwd = 3)


## title
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, paste0("Performance for 1-year control event probability=", round(fixed_perf[par_range == "hr", ]$one_year_event_prob_C[1], 3), " and varying values of HR"), cex = 2, font = 1.5)

## power plot
par(mar = c(3, 3, 2, 2))
fixed_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% plot(hazard_ratio, power,
  type = "l", ylim = c(0, 1),
  lwd = 3, main = "Power", mgp = c(1.5, 0.5, 0), xlab = "Hazard ratio"
)
fixed_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, rep(0.8, length(one_year_event_prob_C)), col = "red", lty = 3)
hybrid1_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, power, col = plot_colors[2], lty = 1, lwd = 3)
hybrid2_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, power, col = plot_colors[3], lty = 1, lwd = 3)
hybrid3_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, power, col = plot_colors[4], lty = 1, lwd = 3)
gs_design_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, power, col = plot_colors[5], lty = 1, lwd = 3)


## trial duration plot
par(mar = c(3, 3, 2, 2))
fixed_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% plot(hazard_ratio, av_trial_dur,
  type = "l", lwd = 3, main = "Trial duration",
  ylab = "average trial duration", mgp = c(1.5, 0.5, 0), xlab = "Hazard ratio"
)
hybrid1_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, av_trial_dur, col = plot_colors[2], lty = 1, lwd = 3)
hybrid2_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, av_trial_dur, col = plot_colors[3], lty = 1, lwd = 3)
hybrid3_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, av_trial_dur, col = plot_colors[4], lty = 1, lwd = 3)
gs_design_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, av_trial_dur, col = plot_colors[5], lty = 1, lwd = 3)


## HR bias plot
par(mar = c(3, 3, 2, 2))
fixed_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% plot(hazard_ratio, perc_bias,
  type = "l", lwd = 3, main = "HR bias",
  ylab = "% Bias", mgp = c(1.5, 0.5, 0), xlab = "Hazard ratio", ylim = c(-13, 13)
)
hybrid1_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, rep(0, length(hazard_ratio)), col = "red", lty = 3, lwd = 1)
hybrid1_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, perc_bias, col = plot_colors[2], lty = 1, lwd = 3)
hybrid2_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, perc_bias, col = plot_colors[3], lty = 1, lwd = 3)
hybrid3_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, perc_bias, col = plot_colors[4], lty = 1, lwd = 3)
gs_design_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, perc_bias, col = plot_colors[5], lty = 1, lwd = 3)


## legend
par(mar = c(1, 1, 0.5, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
plot_colors <- c("black", plot_colors[2], plot_colors[3], plot_colors[4], plot_colors[5])
legend(
  x = "bottom", inset = 0,
  legend = c("Single-stage design", "Hybrid design I", "Hybrid design II", "Hybrid design III", "GSD"),
  col = plot_colors, lty = 1, lwd = 3, cex = 1.5, horiz = TRUE, bty = "n"
)

dev.off()


############### Plot for parameters of secondary interest ############################

pdf(file = "./Reproduction Example/Plots/Figure3.pdf", width = 12, height = 8)


layout(matrix(c(1, 1, 2, 3, 4, 4, 5, 6, 7, 7), ncol = 2, byrow = TRUE), heights = c(0.06, 0.4, 0.06, 0.4, 0.06))

## title
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, paste0("Performance for HR=", round(fixed_perf[par_range == "ev", ]$hazard_ratio[1], 2), " and varying values of control event probability"), cex = 2, font = 1.5)

## number of patients
par(mar = c(3, 3, 2, 2))
fixed_perf[par_range == "ev", ] %$% plot(one_year_event_prob_C, av_patients,
  type = "l",
  lwd = 3, main = "Average number of patients", mgp = c(1.5, 0.5, 0),
  xlab = "1-year event probability in control group", ylim = c(0, 200), ylab = "average number of patients"
)
hybrid1_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, av_patients, col = plot_colors[2], lty = 1, lwd = 3)
hybrid2_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, av_patients, col = plot_colors[3], lty = 1, lwd = 3)
hybrid3_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, av_patients, col = plot_colors[4], lty = 1, lwd = 3)
gs_design_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, av_patients, col = plot_colors[5], lty = 1, lwd = 3)

## probability of unblinding
par(mar = c(3, 3, 2, 2))
fixed_perf[par_range == "ev", ] %$% plot(one_year_event_prob_C, prob_unblind,
  type = "l",
  lwd = 3, main = "Probability of unblinding", mgp = c(1.5, 0.5, 0),
  xlab = "1-year event probability in control group", ylim = c(0, 1), ylab = "prob. unblinding"
)
hybrid1_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, prob_unblind, col = plot_colors[2], lty = 1, lwd = 3)
hybrid2_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, prob_unblind, col = plot_colors[3], lty = 1, lwd = 3)
hybrid3_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, prob_unblind, col = plot_colors[4], lty = 1, lwd = 3)
gs_design_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, prob_unblind, col = plot_colors[5], lty = 1, lwd = 3)

## title
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, paste0("Performance for 1-year control event probability=", round(fixed_perf[par_range == "hr", ]$one_year_event_prob_C[1], 3), " and varying values of HR"), cex = 2, font = 1.5)

## number of patients
par(mar = c(3, 3, 2, 2))
fixed_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% plot(hazard_ratio, av_patients,
  type = "l", lwd = 3, main = "Average number of patients",
  ylab = "average number of patients", mgp = c(1.5, 0.5, 0), xlab = "Hazard ratio", ylim = c(0, 200)
)
hybrid1_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, av_patients, col = plot_colors[2], lty = 1, lwd = 3)
hybrid2_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, av_patients, col = plot_colors[3], lty = 1, lwd = 3)
hybrid3_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, av_patients, col = plot_colors[4], lty = 1, lwd = 3)
gs_design_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, av_patients, col = plot_colors[5], lty = 1, lwd = 3)

## probability of unblinding
par(mar = c(3, 3, 2, 2))
fixed_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% plot(hazard_ratio, prob_unblind,
  type = "l", lwd = 3, main = "Probability of unblinding",
  ylab = "prob. unblinding", mgp = c(1.5, 0.5, 0), xlab = "Hazard ratio", ylim = c(0, 1)
)
hybrid1_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, prob_unblind, col = plot_colors[2], lty = 1, lwd = 3)
hybrid2_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, prob_unblind, col = plot_colors[3], lty = 1, lwd = 3)
hybrid3_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, prob_unblind, col = plot_colors[4], lty = 1, lwd = 3)
gs_design_perf[par_range == "hr", ] %>% arrange(hazard_ratio) %$% lines(hazard_ratio, prob_unblind, col = plot_colors[5], lty = 1, lwd = 3)

## legend
par(mar = c(1, 1, 0.5, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(
  x = "bottom", inset = 0,
  legend = c("Single-stage design", "Hybrid design I", "Hybrid design II", "Hybrid design III", "GSD"),
  col = plot_colors, lty = 1, lwd = 3, cex = 1.5, horiz = TRUE, bty = "n"
)

dev.off()
