#-----------------------------------------------------------
# Script: Simulated Survival Curves with Interim Analysis in Appendix A.5
#-----------------------------------------------------------

# Load utility functions (if needed for helper functions)
source("0_util.R")

#-----------------------------------------------------------
# Study 1: Simulation of survival data and KM plot (Pocock)
#-----------------------------------------------------------
set.seed(16082025) # For reproducibility

# Study parameters
n <- 208 / 2 # Number of patients per group
lambda_event_trt <- 0.00878 # Event rate (treatment group)
lambda_event_ctl <- 0.01860 # Event rate (control group)
lambda_cens <- 0.00084 # Dropout/censoring rate
admin_cens <- 48 # Administrative censoring at 48 months

# Simulate exponential survival times for treatment and control groups
T_event_trt <- rexp(n, rate = lambda_event_trt)
T_event_ctl <- rexp(n, rate = lambda_event_ctl)

# Simulate censoring times
T_cens_trt <- rexp(n, rate = lambda_cens)
T_cens_ctl <- rexp(n, rate = lambda_cens)

# Observed times: take the minimum of event, censoring, and administrative censoring
time_trt <- pmin(T_event_trt, T_cens_trt, admin_cens)
time_ctl <- pmin(T_event_ctl, T_cens_ctl, admin_cens)

# Event indicator: 1 if event occurs before any censoring
event_trt <- as.integer((T_event_trt <= T_cens_trt) & (T_event_trt <= admin_cens))
event_ctl <- as.integer((T_event_ctl <= T_cens_ctl) & (T_event_ctl <= admin_cens))

# Combine into a single data frame
df <- data.frame(
  time = c(time_trt, time_ctl),
  event = c(event_trt, event_ctl),
  group = factor(rep(c("Treatment", "Control"), each = n))
)

# Fit Kaplan-Meier survival curves
fit <- survfit(Surv(time, event) ~ group, data = df)

# Plot survival curves with risk table, adding interim analysis line
p1 <- ggsurvplot(
  fit,
  data = df,
  conf.int = FALSE,
  palette = c("red", "green"),
  risk.table = TRUE, # Show number at risk below plot
  risk.table.height = 0.25, # Height of risk table
  ggtheme = theme_minimal(),
  legend.title = "Group",
  legend.labs = c("Control", "Treatment"),
  title = "Exemplary Simulation Run for Study 1 (Pocock)",
  xlab = "Time (months)",
  ylab = "Survival Probability",
  xlim = c(0, 48),
  break.time.by = 12,
  censor.size = 8
)

# Add vertical line at interim analysis (24 months) and annotation
p1$plot <- p1$plot +
  geom_vline(xintercept = 24, linetype = "dashed", color = "black", size = 1) +
  annotate("text",
    x = 24, y = 0.1, label = "Interim analysis",
    color = "black", angle = 90, vjust = -0.5, hjust = 0
  )

#-----------------------------------------------------------
# Study 2: Same structure, different parameters (Pocock)
#-----------------------------------------------------------
set.seed(16082025)

n <- 136 / 2 # Patients per group
lambda_event_trt <- 0.0426 # Event rate (treatment)
lambda_event_ctl <- 0.0764 # Event rate (control)
lambda_cens <- 0.00084
admin_cens <- 48

# Simulate event and censoring times
T_event_trt <- rexp(n, rate = lambda_event_trt)
T_event_ctl <- rexp(n, rate = lambda_event_ctl)
T_cens_trt <- rexp(n, rate = lambda_cens)
T_cens_ctl <- rexp(n, rate = lambda_cens)

# Observed times and event indicators
time_trt <- pmin(T_event_trt, T_cens_trt, admin_cens)
event_trt <- as.integer((T_event_trt <= T_cens_trt) & (T_event_trt <= admin_cens))

time_ctl <- pmin(T_event_ctl, T_cens_ctl, admin_cens)
event_ctl <- as.integer((T_event_ctl <= T_cens_ctl) & (T_event_ctl <= admin_cens))

# Combine into data frame
df <- data.frame(
  time = c(time_trt, time_ctl),
  event = c(event_trt, event_ctl),
  group = factor(rep(c("Treatment", "Control"), each = n))
)

# Kaplan-Meier fit
fit <- survfit(Surv(time, event) ~ group, data = df)

# Plot survival curves with risk table
p2 <- ggsurvplot(
  fit,
  data = df,
  conf.int = FALSE,
  palette = c("red", "green"),
  risk.table = TRUE,
  risk.table.height = 0.25,
  ggtheme = theme_minimal(),
  legend.title = "Group",
  legend.labs = c("Control", "Treatment"),
  title = "Exemplary Simulation Run for Study 2 (Pocock)",
  xlab = "Time (months)",
  ylab = "Survival Probability",
  xlim = c(0, 48),
  break.time.by = 12,
  censor.size = 8
)

# Add vertical line at interim (24 months)
p2$plot <- p2$plot +
  geom_vline(xintercept = 24, linetype = "dashed", color = "black", size = 1) +
  annotate("text",
    x = 24, y = 0.6, label = "Interim analysis",
    color = "black", angle = 90, vjust = -0.5, hjust = 0
  )

#-----------------------------------------------------------
# Arrange both KM plots side by side with risk tables
#-----------------------------------------------------------
res <- arrange_ggsurvplots(
  list(p1, p2),
  ncol = 2, nrow = 1,
  risk.table.height = 0.25,
  print = FALSE
)

# Print combined plot
print(res)

# Save as PDF
ggsave("Appendix/Figure10.pdf", res, width = 12, height = 6)
