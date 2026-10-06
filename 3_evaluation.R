###############################################################################
# THIS SCRIPT IMPORTS THE DATA AND DESIGNS AND COMPUTES THE PERFORMANCE       #
# CHARACTERISTICS ACCORDINGLY                                                 #
###############################################################################

source("0_util.R")

# Type-1 Error Rate and Bias ------------------------------------------------

## Common simulation parameters
ev_times_names <- paste0("evProb", 1:8)
design_names <- c("fixed", "gs_design", "hybrid1", "hybrid2", "hybrid3")
n_values <- c(180, 360, 720, 1080)
system.time({
for (n in n_values) {
  for (i_data in seq_along(ev_times_names)) {
    ## read event data and survival parameters and obtain number of simulation runs
    data_all <- fread(paste0("./Simulation Study/Data/Type-1 Error Rate/n=", n, "/", ev_times_names[i_data], ".csv"))
    surv_param <- readRDS(paste0("./Simulation Study/Data/Type-1 Error Rate/n=", n, "/", ev_times_names[i_data], "_parameters.RDS"))
    n_sim <- max(data_all$sim_idx)

    for (i_design in seq_along(design_names)) {
      ## read design
      design_eval <- readRDS(paste0("./Simulation Study/Designs/Type-1 Error Rate/n=", n, "/", design_names[i_design], "_", i_data, ".RDS"))

      ## collect results
      int_stats_all <- NULL
      sec_stats_all <- NULL

      for (i_sim in 1:n_sim) {
        ## extract simulation iteration
        data <- data_all[sim_idx == i_sim, ]

        ## interim analysis
        data_int <- get_results_at_time(p_data_res = data, p_t = design_eval$t_int, p_end_rec = design_eval$t_int)
        survdiff_int <- survdiff(Surv(t, status) ~ group, data = data_int)
        res_logrank_int <- logrank_func(survdiff_int)
        cox_int <- coxph(Surv(t, status) ~ group, data = data_int)
        int_stats <- data.frame(
          sim_idx = i_sim,
          n_pat = nrow(data_int),
          n_ev = sum(data_int$status),
          perc_evI = mean(data_int[group == "I", ]$status),
          perc_evC = mean(data_int[group == "C", ]$status),
          lr_stat = res_logrank_int$log_rank_stat,
          hazard_ratio = exp(cox_int$coefficients),
          trial_dur = unique(data_int$time_eval),
          unblinded = design_eval$is_unblinded(sum(data_int$status))
        )

        ## second stage
        sec_rec_length <- design_eval$dec_rec_dur(
          p_nEv = int_stats$n_ev,
          p_logrank_int = res_logrank_int$log_rank_stat
        )
        data_sec <- get_results_at_time(
          p_data_res = data,
          p_t = design_eval$t_int + sec_rec_length + design_eval$min_fu,
          p_end_rec = design_eval$t_int + sec_rec_length
        )
        survdiff_sec <- survdiff(Surv(t, status) ~ group, data = data_sec)
        res_logrank_sec <- logrank_func(survdiff_sec)
        res_inv_norm <- indep_inc_func(p_lr_fin = res_logrank_sec, p_lr_int = res_logrank_int, p_w1 = design_eval$w1)
        cox_sec <- coxph(Surv(t, status) ~ group, data = data_sec)
        surv_group <- summary(survfit(Surv(t, status) ~ group, data = data_sec), times = 12)
        sec_stats <- data.frame(
          sim_idx = i_sim,
          n_pat = nrow(data_sec),
          n_ev = sum(data_sec$status),
          perc_evI = mean(data_sec[group == "I", ]$status),
          perc_evC = mean(data_sec[group == "C", ]$status),
          lr_stat = res_logrank_sec$log_rank_stat,
          inv_norm = res_inv_norm$invN_stat,
          hazard_ratio = exp(cox_sec$coefficients),
          trial_dur = unique(data_sec$time_eval),
          km_surv_1year_I = surv_group$surv[surv_group$strata == "group=I"],
          km_surv_1year_C = surv_group$surv[surv_group$strata == "group=C"]
        )

        int_stats_all %<>% rbind(int_stats)
        sec_stats_all %<>% rbind(sec_stats)

        if (i_sim %% 1000 == 0) {
          print(paste0("n=", n, " - Simulation run ", i_sim, " of ", n_sim))
        }
      }

      ## join stages
      names(int_stats_all)[!names(int_stats_all) %in% "sim_idx"] <- paste0("int_", names(int_stats_all)[!names(int_stats_all) %in% "sim_idx"])
      names(sec_stats_all)[!names(sec_stats_all) %in% "sim_idx"] <- paste0("sec_", names(sec_stats_all)[!names(sec_stats_all) %in% "sim_idx"])
      stats_all <- join(int_stats_all, sec_stats_all, by = "sim_idx") %>% as.data.table()

      ## rejection
      stats_all <- design_eval$rej_func(stats_all)

      ## extreme cases
      idx_cap <- which(
        (stats_all$int_unblinded & (stats_all$int_perc_evC %in% c(0, 1) | stats_all$int_perc_evI %in% c(0, 1))) |
          stats_all$sec_perc_evC %in% c(0, 1) | stats_all$sec_perc_evI %in% c(0, 1)
      )

      ## performance
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
      path_res <- paste0("./Simulation Study/Results/Type-1 Error Rate/n=", n, "/", ev_times_names[i_data])
      if (!file.exists(path_res)) dir.create(path_res, recursive = TRUE)
      path_res <- paste0(path_res, "/", design_names[i_design])
      dir.create(path_res)
      fwrite(stats_all, paste0(path_res, "/stats_all.csv"))
      saveRDS(performance, paste0(path_res, "/performance.RDS"))
    }
  }
}
})

# Performance Study ---------------------------------------------------------

design_names <- c("fixed", "gs_design", "hybrid1", "hybrid2", "hybrid3")

studies <- list(
  list(name = "Study 1_Pocock", ev_times = c(paste0("evProb", 1:5), paste0("hr", 1:11))),
  list(name = "Study 2_Pocock", ev_times = c(paste0("evProb", 1:6), paste0("hr", 1:11))),
  list(name = "Study 1_OBrien-Fleming", ev_times = c(paste0("evProb", 1:5), paste0("hr", 1:11))),
  list(name = "Study 2_OBrien-Fleming", ev_times = c(paste0("evProb", 1:6), paste0("hr", 1:11)))
)

for (study in studies) {
  ev_times_names <- study$ev_times

  for (i_data in seq_along(ev_times_names)) {
    ## read event data + survival parameters
    data_all <- fread(paste0("./Simulation Study/Data/Performance/", study$name, "/", ev_times_names[i_data], ".csv"))
    surv_param <- readRDS(paste0("./Simulation Study/Data/Performance/", study$name, "/", ev_times_names[i_data], "_parameters.RDS"))
    n_sim <- max(data_all$sim_idx)

    for (i_design in seq_along(design_names)) {
      ## read design
      design_eval <- readRDS(paste0("./Simulation Study/Designs/Performance/", study$name, "/", design_names[i_design], ".RDS"))

      int_stats_all <- NULL
      sec_stats_all <- NULL

      for (i_sim in 1:n_sim) {
        ## extract simulation iteration
        data <- data_all[sim_idx == i_sim, ]

        ## interim analysis
        data_int <- get_results_at_time(p_data_res = data, p_t = design_eval$t_int, p_end_rec = design_eval$t_int)
        survdiff_int <- survdiff(Surv(t, status) ~ group, data = data_int)
        res_logrank_int <- logrank_func(survdiff_int)
        cox_int <- coxph(Surv(t, status) ~ group, data = data_int)
        int_stats <- data.frame(
          sim_idx = i_sim,
          n_pat = nrow(data_int),
          n_ev = sum(data_int$status),
          perc_evI = mean(data_int[group == "I", ]$status),
          perc_evC = mean(data_int[group == "C", ]$status),
          lr_stat = res_logrank_int$log_rank_stat,
          hazard_ratio = exp(cox_int$coefficients),
          trial_dur = unique(data_int$time_eval),
          unblinded = design_eval$is_unblinded(sum(data_int$status))
        )

        ## second stage
        sec_rec_length <- design_eval$dec_rec_dur(
          p_nEv = int_stats$n_ev,
          p_logrank_int = res_logrank_int$log_rank_stat
        )
        data_sec <- get_results_at_time(
          p_data_res = data,
          p_t = design_eval$t_int + sec_rec_length + design_eval$min_fu,
          p_end_rec = design_eval$t_int + sec_rec_length
        )
        data_sec$time_eval <- ifelse(sec_rec_length == 0, design_eval$t_int, data_sec$time_eval)
        survdiff_sec <- survdiff(Surv(t, status) ~ group, data = data_sec)
        res_logrank_sec <- logrank_func(survdiff_sec)
        res_inv_norm <- indep_inc_func(p_lr_fin = res_logrank_sec, p_lr_int = res_logrank_int, p_w1 = design_eval$w1)
        cox_sec <- coxph(Surv(t, status) ~ group, data = data_sec)
        surv_group <- summary(survfit(Surv(t, status) ~ group, data = data_sec), times = 12)
        sec_stats <- data.frame(
          sim_idx = i_sim,
          n_pat = nrow(data_sec),
          n_ev = sum(data_sec$status),
          perc_evI = mean(data_sec[group == "I", ]$status),
          perc_evC = mean(data_sec[group == "C", ]$status),
          lr_stat = res_logrank_sec$log_rank_stat,
          inv_norm = res_inv_norm$invN_stat,
          hazard_ratio = exp(cox_sec$coefficients),
          trial_dur = unique(data_sec$time_eval),
          km_surv_1year_I = surv_group$surv[surv_group$strata == "group=I"],
          km_surv_1year_C = surv_group$surv[surv_group$strata == "group=C"]
        )

        int_stats_all %<>% rbind(int_stats)
        sec_stats_all %<>% rbind(sec_stats)

        if (i_sim %% 1000 == 0) {
          print(paste0(study$name, " - Simulation run ", i_sim, " of ", n_sim))
        }
      }

      ## join stages
      names(int_stats_all)[!names(int_stats_all) %in% "sim_idx"] <- paste0("int_", names(int_stats_all)[!names(int_stats_all) %in% "sim_idx"])
      names(sec_stats_all)[!names(sec_stats_all) %in% "sim_idx"] <- paste0("sec_", names(sec_stats_all)[!names(sec_stats_all) %in% "sim_idx"])
      stats_all <- join(int_stats_all, sec_stats_all, by = "sim_idx") %>% as.data.table()

      ## rejection
      stats_all <- design_eval$rej_func(stats_all)

      ## extreme cases
      idx_cap <- which(
        (stats_all$int_unblinded &
          (stats_all$int_perc_evC %in% c(0, 1) | stats_all$int_perc_evI %in% c(0, 1))) |
          stats_all$sec_perc_evC %in% c(0, 1) |
          stats_all$sec_perc_evI %in% c(0, 1)
      )

      ## performance
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
      path_res <- paste0("./Simulation Study/Results/Performance/", study$name, "/", ev_times_names[i_data])
      if (!file.exists(path_res)) dir.create(path_res, recursive = TRUE)
      path_res <- paste0(path_res, "/", design_names[i_design])
      dir.create(path_res)
      fwrite(stats_all, paste0(path_res, "/stats_all.csv"))
      saveRDS(performance, paste0(path_res, "/performance.RDS"))
    }
  }
}