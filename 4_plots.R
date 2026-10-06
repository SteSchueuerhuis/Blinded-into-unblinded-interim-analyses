###############################################################################
# THIS SCRIPT CREATES THE PLOTS SHOWN IN THE PAPER                            #
###############################################################################

source("0_util.R")

# Type-1 Error Rate and Bias ------------------------------------------------

ev_times_names <- c(paste0("evProb", 1:8))
design_names <- c("fixed", "hybrid1", "hybrid2", "hybrid3", "gs_design")
plot_colors <- c("black", "darkorange", "forestgreen", "dodgerblue2", "firebrick4")

# sample sizes to iterate over
n_values <- c(180, 360, 720, 1080)

# bias y-limits (matched to your original code)
bias_ylims <- list(
  `180` = c(-1, 20), # none plotted for n=180
  `360` = c(-1, 20),
  `720` = c(-1, 10),
  `1080` = c(-1, 10)
)

# bias figure numbers (for consistency with original numbering)
bias_figs <- list(
  `180` = "Figure9-0.pdf",
  `360` = "Figure9-1.pdf",
  `720` = "Figure9-2.pdf",
  `1080` = "Figure9-3.pdf"
)

for (n in n_values) {
  # collect performance results per design
  for (i_design in 1:length(design_names)) {
    design_name <- design_names[i_design]

    one_year_event_prob_C <- c()
    power <- c()
    av_trial_dur <- c()
    hazard_ratio <- c()
    perc_bias <- c()
    prob_unblind <- c()
    av_patients <- c()

    for (i_evt_set in 1:length(ev_times_names)) {
      path_res <- paste0("./Simulation Study/Results/Type-1 Error Rate/n=", n, "/", ev_times_names[i_evt_set], "/", design_name)
      performance <- readRDS(paste0(path_res, "/performance.RDS"))

      if (performance$prob_cap == 0) {
        performance$av_est_hazard_ratio_cap <- performance$av_est_hazard_ratio
        performance$av_est_log_hazard_ratio_cap <- performance$av_est_log_hazard_ratio
      }

      one_year_event_prob_C <- c(one_year_event_prob_C, performance$one_year_event_prob_C)
      power <- c(power, performance$power)
      av_trial_dur <- c(av_trial_dur, performance$av_trial_dur)
      hazard_ratio <- c(hazard_ratio, performance$hazard_ratio)
      perc_bias <- c(perc_bias, (performance$av_est_hazard_ratio - performance$hazard_ratio) / performance$hazard_ratio * 100)
      prob_unblind <- c(prob_unblind, performance$prob_unblind)
      av_patients <- c(av_patients, performance$av_pat_total)
    }

    assign(
      x = paste0(design_name, "_perf"),
      value = data.table(
        one_year_event_prob_C = one_year_event_prob_C,
        power = power, av_trial_dur = av_trial_dur,
        hazard_ratio = hazard_ratio, perc_bias = perc_bias,
        prob_unblind = prob_unblind, av_patients = av_patients,
        par_range = substr(ev_times_names, 1, 2)
      )
    )
  }

  ##  -------- Type-1 error plots ----------
  pdf(
    file = paste0(
      "./Simulation Study/Plots/Type-1 Error Rate/n=", n, "/Figure",
      ifelse(n == 180, 5, ifelse(n == 360, 6, ifelse(n == 720, 7, 8))), ".pdf"
    ),
    width = 8, height = 5.5
  )

  fixed_perf[par_range == "ev", ] %$% plot(one_year_event_prob_C, power,
    type = "l", ylim = c(0.015, 0.037),
    ylab = "Simulated type-1 error rate", lwd = 3, mgp = c(1.5, 0.5, 0),
    xlab = "1-year event probability", col = plot_colors[1]
  )

  fixed_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C,
    rep(
      0.025 - qnorm(1 - 0.025) * sqrt((0.025 * (1 - 0.025)) / 10000),
      length(one_year_event_prob_C)
    ),
    col = "grey", lty = 3, lwd = 2
  )
  fixed_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C,
    rep(
      0.025 + qnorm(1 - 0.025) * sqrt((0.025 * (1 - 0.025)) / 10000),
      length(one_year_event_prob_C)
    ),
    col = "grey", lty = 3, lwd = 2
  )
  fixed_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, rep(0.025, length(one_year_event_prob_C)),
    col = "red", lty = 3, lwd = 2
  )

  hybrid1_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, power, col = plot_colors[2], lty = 1, lwd = 2)
  hybrid2_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, power, col = plot_colors[3], lty = 1, lwd = 2)
  hybrid3_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, power, col = plot_colors[4], lty = 1, lwd = 2)
  gs_design_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, power, col = plot_colors[5], lty = 1, lwd = 2)

  legend(0.05, 0.0375,
    legend = c("Target significance level", "Single-stage design", "Hybrid design I", "Hybrid design II", "Hybrid design III", "GSD"),
    col = c("red", plot_colors), box.lty = 0, lty = c(3, rep(1, 5)), lwd = 3, cex = 1
  )

  dev.off()

  # -------- Bias plots ----------
  pdf(
    file = paste0("./Simulation Study/Plots/Type-1 Error Rate/n=", n, "/", bias_figs[[as.character(n)]]),
    width = 8, height = 5.5
  )

  fixed_perf[par_range == "ev", ] %$% plot(one_year_event_prob_C, perc_bias,
    type = "l",
    ylim = bias_ylims[[as.character(n)]],
    ylab = "% Bias", lwd = 3, mgp = c(1.5, 0.5, 0),
    xlab = "1-year event probability", col = plot_colors[1]
  )
  fixed_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, rep(0.0, length(one_year_event_prob_C)),
    col = "red", lty = 3, lwd = 2
  )
  hybrid1_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, perc_bias, col = plot_colors[2], lty = 1, lwd = 2)
  hybrid2_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, perc_bias, col = plot_colors[3], lty = 1, lwd = 2)
  hybrid3_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, perc_bias, col = plot_colors[4], lty = 1, lwd = 2)
  gs_design_perf[par_range == "ev", ] %$% lines(one_year_event_prob_C, perc_bias, col = plot_colors[5], lty = 1, lwd = 2)

  legend(0.27, max(bias_ylims[[as.character(n)]]) - 0.5,
    legend = c("Single-stage design", "Hybrid design I", "Hybrid design II", "Hybrid design III", "GSD"),
    col = plot_colors, box.lty = 0, lty = rep(1, 5), lwd = 3, cex = 1
  )

  dev.off()
}



# Performance Study ---------------------------------------------------------

## Pocock -------------------------------------------------------------------

### Study 1 -----------------------------------------------------------------

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
    path_res <- paste0("./Simulation Study/Results/Performance/Study 1_Pocock/", ev_times_names[i_evt_set], "/", design_name)
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

################### Figure 2 #####################

pdf(file = "./Simulation Study/Plots/Performance/Study 1_Pocock/Figure2.pdf", width = 12, height = 8)

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


############### Figure 3 ############################

pdf(file = "./Simulation Study/Plots/Performance/Study 1_Pocock/Figure3.pdf", width = 12, height = 8)


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

### Study 2 -----------------------------------------------------------------

ev_times_names <- c(paste0("evProb", 1:6), paste0("hr", 1:11))
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
    path_res <- paste0("./Simulation Study/Results/Performance/Study 2_Pocock/", ev_times_names[i_evt_set], "/", design_name)
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

################## Figure 11 #####################


pdf(file = "./Simulation Study/Plots/Performance/Study 2_Pocock/Figure11.pdf", width = 12, height = 8)

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


############### Figure 12 ############################

pdf(file = "./Simulation Study/Plots/Performance/Study 2_Pocock/Figure12.pdf", width = 12, height = 8)


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


## O'Brien-Fleming ----------------------------------------------------------

### Study 1 -----------------------------------------------------------------

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
    path_res <- paste0("./Simulation Study/Results/Performance/Study 1_OBrien-Fleming/", ev_times_names[i_evt_set], "/", design_name)
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

################## Supplement Figure 1 #####################

pdf(file = "./Simulation Study/Plots/Performance/Study 1_OBrien-Fleming/Supplement_Figure1.pdf", width = 12, height = 8)

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


################## Supplement Figure 2 #####################

pdf(file = "./Simulation Study/Plots/Performance/Study 1_OBrien-Fleming/Supplement_Figure2.pdf", width = 12, height = 8)


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


### Study 2 -----------------------------------------------------------------

ev_times_names <- c(paste0("evProb", 1:6), paste0("hr", 1:11))
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
    path_res <- paste0("./Simulation Study/Results/Performance/Study 2_OBrien-Fleming/", ev_times_names[i_evt_set], "/", design_name)
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

################## Supplement Figure 3 #####################

pdf(file = "./Simulation Study/Plots/Performance/Study 2_OBrien-Fleming/Supplement_Figure3.pdf", width = 12, height = 8)

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


################## Supplement Figure 4 #####################

pdf(file = "./Simulation Study/Plots/Performance/Study 2_OBrien-Fleming/Supplement_Figure4.pdf", width = 12, height = 8)


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
