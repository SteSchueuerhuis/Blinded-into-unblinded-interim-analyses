#-----------------------------------------------------------
# Script: Figure 4 in Appendix A.1
#-----------------------------------------------------------

# Load helper functions
source("0_util.R")

#-------------------------------
# Set study parameters
#-------------------------------
pi1 <- 0.1 # Event probability per 12 months in experimental group
pi2 <- 0.2 # Event probability per 12 months in control group
cens <- 0.01 # Dropout probability per 12 months

#-------------------------------
# Define group-sequential design using rpact
#-------------------------------
# This design roughly corresponds to interim analyses after 24 and 48 months
# typeOfDesign = "P": Pocock design
# futilityBounds = 0, bindingFutility = TRUE
design <- getDesignInverseNormal(
  kMax = 2,
  typeOfDesign = "P",
  informationRates = c(.3, 1),
  futilityBounds = 0,
  bindingFutility = TRUE
)

#-------------------------------
# Calculate required sample size for survival study
#-------------------------------
# Assumes exponential distribution, 12-month follow-up, accrual over 0-36 months
sampleSize <- getSampleSizeSurvival(
  design = design,
  accrualTime = c(0, 36),
  dropoutTime = 12,
  dropoutRate1 = cens,
  dropoutRate2 = cens,
  eventTime = 12,
  pi1 = pi1,
  pi2 = pi2,
  followUpTime = 12
)

#-------------------------------
# Convert event probabilities to exponential lambda parameters
#-------------------------------
# Solve for lambda such that P(T <= 12) = pi1 or pi2
lambda1 <- uniroot(\(x) pexp(12, x) - pi1, interval = c(0.001, 10))$root
lambda2 <- uniroot(\(x) pexp(12, x) - pi2, interval = c(0.001, 10))$root

# Solve for dropout rate lambda
lambda_cens <- uniroot(\(x) pexp(12, x) - cens, interval = c(0.00001, 10))$root

# Number of patients per group at interim
n_per_group <- round(sampleSize$maxNumberOfSubjects) / 2

#-------------------------------
# Expected number of events at interim (Dl)
#-------------------------------
Dl <- expected_events_exponential(
  overall_sample_size = n_per_group,
  time = 24,
  rate_event = lambda1,
  rate_censoring = lambda_cens,
  end_recruitment_time = 36
) +
  expected_events_exponential(
    overall_sample_size = n_per_group,
    time = 24,
    rate_event = lambda2,
    rate_censoring = lambda_cens,
    end_recruitment_time = 36
  )

#-------------------------------
# Hazard ratio and alpha spent at interim
#-------------------------------
HR <- lambda1 / lambda2
alpha1 <- design$alphaSpent[1]

#-------------------------------
# Hybrid Design 1
#-------------------------------
Delta <- seq(0.05, 0.95, 0.05) # Grid of Delta design parameters

df <- get_hybrid1(
  Dl = Dl,
  HR = HR,
  alpha1 = alpha1,
  Delta = Delta,
  n_int = ceiling(sampleSize$numberOfSubjects[1])
)

# Event probability at interim (average across groups)
event_prob_interim <- (
  event_probability_exponential(n_per_group, 24, lambda1, lambda_cens, 36) +
    event_probability_exponential(n_per_group, 24, lambda2, lambda_cens, 36)
) / 2

#-------------------------------
# Hybrid Design 2
#-------------------------------
# Compute lower and upper thresholds using binomial approximation
h2_l <- qbinom(Delta, n_per_group * 2, event_prob_interim, lower.tail = TRUE) - 1
h2_u <- qbinom(1 - Delta, n_per_group * 2, event_prob_interim)
df$h2_l <- h2_l
df$h2_u <- h2_u

#-------------------------------
# Hybrid Design 3
#-------------------------------
# Thresholds based on Dl ± Dl*Delta
h3_l <- floor(Dl - (Dl * Delta))
h3_u <- ceiling(Dl + (Dl * Delta))
df$h3_l <- h3_l
df$h3_u <- h3_u

#-------------------------------
# Plot thresholds for all hybrid designs
#-------------------------------
pdf(file = "./Appendix/Figure4.pdf", width = 8, height = 5.5)

# Hybrid Design 1
plot(df$Delta, df$h1_l,
  xaxt = "n", xlim = c(0.05, 0.95),
  lty = 1, col = "darkorange", type = "l",
  ylim = c(0, 138), lwd = 2.2,
  xlab = expression(Delta[i]),
  ylab = expression(paste(d[l], " and ", d[u])),
  cex.lab = 1.3, cex.axis = 1.3, mgp = c(2, 0.5, 0)
)
axis(1, at = seq(0, 1, 0.1), cex.axis = 1.3, cex.lab = 1.3, lwd = 0.9)
lines(df$Delta, df$h1_u, col = "darkorange", lwd = 2.2, lty = 1)

# Hybrid Design 2
lines(df$Delta, df$h2_l, col = "forestgreen", lwd = 2.2, lty = 2)
lines(df$Delta, df$h2_u, col = "forestgreen", lwd = 2.2, lty = 2)

# Hybrid Design 3
lines(df$Delta, df$h3_l, col = "dodgerblue2", lwd = 2.2, lty = 3)
lines(df$Delta, df$h3_u, col = "dodgerblue2", lwd = 2.2, lty = 3)

# Legend
legend(0.1, 135,
  legend = c("Hybrid design I", "Hybrid design II", "Hybrid design III"),
  col = c("darkorange", "forestgreen", "dodgerblue2"),
  box.lty = 0, lty = 1:3, lwd = 3, cex = 1.4
)

dev.off()
