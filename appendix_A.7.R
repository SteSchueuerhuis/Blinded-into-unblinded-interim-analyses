#-----------------------------------------------------------
# Script: Creating Table in Appendix A.7
#-----------------------------------------------------------

ev_times_names <- c(paste0("evProb", 1:5), paste0("hr", 1:11))
design_names <- c("hybrid1", "hybrid2", "hybrid3", "gs_design")

# Empty list to collect all results
results_list <- list()

for (i_design in design_names) {
  for (i_evt_set in ev_times_names) {
    # define file path
    path_res <- paste0(
      "./Simulation Study/Results/Performance/Study 1_Pocock/",
      i_evt_set, "/", i_design, "/"
    )

    # read file
    stats_all <- fread(paste0(path_res, "stats_all.csv"))

    # select only desired columns
    df_sub <- stats_all[, .(int_n_ev, int_lr_stat, int_unblinded, sec_trial_dur)]

    # name of element in the list
    list_name <- paste(i_evt_set, i_design, sep = "_")

    # store into list
    results_list[[list_name]] <- df_sub
  }
}

# Now create merged data.frames per event set
final_datasets <- list()

for (ev in ev_times_names) {
  df <- data.frame(
    int_n_ev = results_list[[paste0(ev, "_hybrid1")]]$int_n_ev,
    int_lr_stat = results_list[[paste0(ev, "_hybrid1")]]$int_lr_stat,
    unblinding_hybrid1 = results_list[[paste0(ev, "_hybrid1")]]$int_unblinded,
    unblinding_hybrid2 = results_list[[paste0(ev, "_hybrid2")]]$int_unblinded,
    unblinding_hybrid3 = results_list[[paste0(ev, "_hybrid3")]]$int_unblinded,
    sec_trial_dur = results_list[[paste0(ev, "_gs_design")]]$sec_trial_dur,
    int_rej = results_list[[paste0(ev, "_gs_design")]]$sec_trial_dur == 24
  )

  final_datasets[[ev]] <- df
}

# Survival parameters for study 1
p1year_values <- c(0.05, 0.1, 0.2, 0.4, 0.6)
hazard_ratio <- c(seq(0.2, 0.5, 0.1), 0.47, seq(0.6, 1.1, 0.1))

# assign to evProb datasets
for (i in seq_along(p1year_values)) {
  final_datasets[[paste0("evProb", i)]]$par <- p1year_values[i]
}

# assign to hr datasets
for (i in seq_along(hazard_ratio)) {
  final_datasets[[paste0("hr", i)]]$par <- hazard_ratio[i]
}

# Initialize results collector
summary_results <- data.frame(
  ev_times_name = character(),
  correlation   = numeric()
)

for (nm in names(final_datasets)) {
  df <- final_datasets[[nm]]

  # compute correlation between LR and d_1
  par <- df$par[1]
  cor_val <- cor(df$int_n_ev, df$int_lr_stat, use = "complete.obs")

  # # Stay blinded when GSD did not stop
  p_blind_H1 <- 1 - (df %>% filter(!unblinding_hybrid1) %>% select(int_rej) %>% pull() %>% mean())
  p_blind_H2 <- 1 - (df %>% filter(!unblinding_hybrid2) %>% select(int_rej) %>% pull() %>% mean())
  p_blind_H3 <- 1 - (df %>% filter(!unblinding_hybrid3) %>% select(int_rej) %>% pull() %>% mean())

  # Unblinding when GSD rejects
  p_unblind_H1 <- df %>%
    filter(unblinding_hybrid1) %>%
    select(int_rej) %>%
    pull() %>%
    mean()
  p_unblind_H2 <- df %>%
    filter(unblinding_hybrid2) %>%
    select(int_rej) %>%
    pull() %>%
    mean()
  p_unblind_H3 <- df %>%
    filter(unblinding_hybrid3) %>%
    select(int_rej) %>%
    pull() %>%
    mean()

  # Monitoring the event numbers
  d_stop <- df %>%
    filter(int_rej) %>%
    select(int_n_ev) %>%
    pull() %>%
    mean(na.rm = T)
  d_stop_f <- df %>%
    filter(int_rej, int_lr_stat <= 0) %>%
    select(int_n_ev) %>%
    pull() %>%
    mean(na.rm = T)
  d_stop_e <- df %>%
    filter(int_rej, int_lr_stat > 0) %>%
    select(int_n_ev) %>%
    pull() %>%
    mean(na.rm = T)
  d_cont <- df %>%
    filter(!int_rej) %>%
    select(int_n_ev) %>%
    pull() %>%
    mean(na.rm = T)

  # store results
  summary_results <- rbind(
    summary_results,
    data.frame(
      ev_times_name = nm,
      par = par,
      correlation = cor_val,
      lr_stat = mean(df$int_lr_stat),
      n_unblind_H1 = df %>% filter(unblinding_hybrid1) %>% nrow(),
      n_unblind_H2 = df %>% filter(unblinding_hybrid2) %>% nrow(),
      n_unblind_H3 = df %>% filter(unblinding_hybrid3) %>% nrow(),
      n_rej_GSD = df %>% filter(int_rej) %>% nrow(),
      n_fut = sum(df$int_lr_stat <= 0),
      d = mean(df$int_n_ev),
      p_unblind_H1 = ifelse(is.na(p_unblind_H1), 0, p_unblind_H1),
      p_unblind_H2 = ifelse(is.na(p_unblind_H2), 0, p_unblind_H2),
      p_unblind_H3 = ifelse(is.na(p_unblind_H3), 0, p_unblind_H3),
      p_blind_H1 = ifelse(is.na(p_blind_H1), 0, p_blind_H1),
      p_blind_H2 = ifelse(is.na(p_blind_H2), 0, p_blind_H2),
      p_blind_H3 = ifelse(is.na(p_blind_H3), 0, p_blind_H3)
    )
  )
}


# Summarize performance characteristics for study 1
summary_results_long_1 <- summary_results %>%
  pivot_longer(
    cols = c(
      n_unblind_H1, n_unblind_H2, n_unblind_H3,
      p_unblind_H1, p_unblind_H2, p_unblind_H3,
      p_blind_H1, p_blind_H2, p_blind_H3
    ),
    names_to = c(".value", "method"),
    names_pattern = "(.*)_(H.)"
  ) %>%
  mutate(
    Study = "Study 1",
    ev_times_name = ifelse(grepl("evProb", ev_times_name), "evProb", "HR")
  ) %>%
  select(Study, ev_times_name, par, method, d, lr_stat, correlation, n_unblind, n_rej_GSD, n_fut, p_unblind, p_blind, method) %>%
  filter((ev_times_name == "evProb" & par %in% c(0.05, 0.2, 0.6)) |
    (ev_times_name == "HR" & par %in% c(0.2, 0.47, 1.1)))


ev_times_names <- c(paste0("evProb", 1:6), paste0("hr", 1:11))
design_names <- c("hybrid1", "hybrid2", "hybrid3", "gs_design")

# Empty list to collect all results
results_list <- list()

for (i_design in design_names) {
  for (i_evt_set in ev_times_names) {
    # define file path
    path_res <- paste0(
      "./Simulation Study/Results/Performance/Study 2_Pocock/",
      i_evt_set, "/", i_design, "/"
    )

    # read file
    stats_all <- fread(paste0(path_res, "stats_all.csv"))

    # select only desired columns
    df_sub <- stats_all[, .(int_n_ev, int_lr_stat, int_unblinded, sec_trial_dur)]

    # name of element in the list
    list_name <- paste(i_evt_set, i_design, sep = "_")

    # store into list
    results_list[[list_name]] <- df_sub
  }
}

# Now create merged data.frames per event set
final_datasets <- list()

for (ev in ev_times_names) {
  df <- data.frame(
    int_n_ev = results_list[[paste0(ev, "_hybrid1")]]$int_n_ev,
    int_lr_stat = results_list[[paste0(ev, "_hybrid1")]]$int_lr_stat,
    unblinding_hybrid1 = results_list[[paste0(ev, "_hybrid1")]]$int_unblinded,
    unblinding_hybrid2 = results_list[[paste0(ev, "_hybrid2")]]$int_unblinded,
    unblinding_hybrid3 = results_list[[paste0(ev, "_hybrid3")]]$int_unblinded,
    sec_trial_dur = results_list[[paste0(ev, "_gs_design")]]$sec_trial_dur,
    int_rej = results_list[[paste0(ev, "_gs_design")]]$sec_trial_dur == 24
  )

  final_datasets[[ev]] <- df
}


# Survival parameters for study 2
p1year_values <- c(0.1, 0.2, 0.4, 0.6, 0.7, 0.8)
hazard_ratio <- c(seq(0.2, 0.5, 0.1), 0.56, seq(0.6, 1.1, 0.1))

# assign to evProb datasets
for (i in seq_along(p1year_values)) {
  final_datasets[[paste0("evProb", i)]]$par <- p1year_values[i]
}

# assign to hr datasets
for (i in seq_along(hazard_ratio)) {
  final_datasets[[paste0("hr", i)]]$par <- hazard_ratio[i]
}

# Initialize results collector
summary_results <- data.frame(
  ev_times_name = character(),
  correlation   = numeric()
)

for (nm in names(final_datasets)) {
  df <- final_datasets[[nm]]

  # compute correlation between LR and d_1
  par <- df$par[1]
  cor_val <- cor(df$int_n_ev, df$int_lr_stat, use = "complete.obs")

  # # Stay blinded when GSD did not stop
  p_blind_H1 <- 1 - (df %>% filter(!unblinding_hybrid1) %>% select(int_rej) %>% pull() %>% mean())
  p_blind_H2 <- 1 - (df %>% filter(!unblinding_hybrid2) %>% select(int_rej) %>% pull() %>% mean())
  p_blind_H3 <- 1 - (df %>% filter(!unblinding_hybrid3) %>% select(int_rej) %>% pull() %>% mean())

  # Unblinding when GSD rejects
  p_unblind_H1 <- df %>%
    filter(unblinding_hybrid1) %>%
    select(int_rej) %>%
    pull() %>%
    mean()
  p_unblind_H2 <- df %>%
    filter(unblinding_hybrid2) %>%
    select(int_rej) %>%
    pull() %>%
    mean()
  p_unblind_H3 <- df %>%
    filter(unblinding_hybrid3) %>%
    select(int_rej) %>%
    pull() %>%
    mean()

  # Monitoring the event numbers
  d_stop <- df %>%
    filter(int_rej) %>%
    select(int_n_ev) %>%
    pull() %>%
    mean(na.rm = T)
  d_stop_f <- df %>%
    filter(int_rej, int_lr_stat <= 0) %>%
    select(int_n_ev) %>%
    pull() %>%
    mean(na.rm = T)
  d_stop_e <- df %>%
    filter(int_rej, int_lr_stat > 0) %>%
    select(int_n_ev) %>%
    pull() %>%
    mean(na.rm = T)
  d_cont <- df %>%
    filter(!int_rej) %>%
    select(int_n_ev) %>%
    pull() %>%
    mean(na.rm = T)

  # store results
  summary_results <- rbind(
    summary_results,
    data.frame(
      ev_times_name = nm,
      par = par,
      correlation = cor_val,
      lr_stat = mean(df$int_lr_stat),
      n_unblind_H1 = df %>% filter(unblinding_hybrid1) %>% nrow(),
      n_unblind_H2 = df %>% filter(unblinding_hybrid2) %>% nrow(),
      n_unblind_H3 = df %>% filter(unblinding_hybrid3) %>% nrow(),
      n_rej_GSD = df %>% filter(int_rej) %>% nrow(),
      n_fut = sum(df$int_lr_stat <= 0),
      d = mean(df$int_n_ev),
      p_unblind_H1 = ifelse(is.na(p_unblind_H1), 0, p_unblind_H1),
      p_unblind_H2 = ifelse(is.na(p_unblind_H2), 0, p_unblind_H2),
      p_unblind_H3 = ifelse(is.na(p_unblind_H3), 0, p_unblind_H3),
      p_blind_H1 = ifelse(is.na(p_blind_H1), 0, p_blind_H1),
      p_blind_H2 = ifelse(is.na(p_blind_H2), 0, p_blind_H2),
      p_blind_H3 = ifelse(is.na(p_blind_H3), 0, p_blind_H3)
    )
  )
}

# Summarize performance characteristics for study 2
summary_results_long_2 <- summary_results %>%
  pivot_longer(
    cols = c(
      n_unblind_H1, n_unblind_H2, n_unblind_H3,
      p_unblind_H1, p_unblind_H2, p_unblind_H3,
      p_blind_H1, p_blind_H2, p_blind_H3
    ),
    names_to = c(".value", "method"),
    names_pattern = "(.*)_(H.)"
  ) %>%
  mutate(
    Study = "Study 2",
    ev_times_name = ifelse(grepl("evProb", ev_times_name), "evProb", "HR")
  ) %>%
  select(Study, ev_times_name, par, method, d, lr_stat, correlation, n_unblind, n_rej_GSD, n_fut, p_unblind, p_blind, method) %>%
  filter((ev_times_name == "evProb" & par %in% c(0.1, 0.6, 0.8)) |
    (ev_times_name == "HR" & par %in% c(0.2, 0.56, 1.1)))

# Merge everything
table6_A7 <- rbind(summary_results_long_1, summary_results_long_2) %>% select(-c(d, lr_stat, n_fut))

# Save as .csv
write.csv(table6_A7, file = "Appendix/A.7_table6.csv")
