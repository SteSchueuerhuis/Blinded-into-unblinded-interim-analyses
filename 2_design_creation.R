###############################################################################
# THIS SCRIPT CREATES AND SAVES THE STUDY DESIGNS THAT ARE TO BE EVALUATED    #
###############################################################################

source("0_util.R")

# Type-1 Error Rate Rate and Bias ------------------------------------------------

pi <- seq(0.05, 0.4, 0.05)
cens <- 0.05
fixed_par <- list(t_int = 24, min_fu = 12, sec_rec_length = 12)

# Censoring rate
lambda_cens <- uniroot(\(x) pexp(12, x) - cens, interval = c(.00001, 10))$root

# Define scenarios
scenarios <- list(
  list(n_total = 180, HR = 0.5),
  list(n_total = 360, HR = 0.5),
  list(n_total = 720, HR = 0.6),
  list(n_total = 1080, HR = 0.66)
)

for (scenario in scenarios) {
  n_total <- scenario$n_total
  HR <- scenario$HR
  n_per_group <- n_total / 2
  output_dir <- paste0("./Simulation Study/Designs/Type-1 Error Rate/n=", n_total)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Initialize dataframe
  df <- data.frame(
    pi = pi,
    n = n_total,
    n_int = NA, int_t = NA, fin_t = NA,
    Dlint = NA, Dlfin = NA, inf = NA,
    c1 = NA, c2 = NA, l = 0, alpha1 = NA,
    h1_u = NA, h1_l = NA,
    h2_u = NA, h2_l = NA,
    h3_u = NA, h3_l = NA
  )

  # Compute events and design parameters
  for (i in seq_along(pi)) {
    lambda <- uniroot(\(x) pexp(12, x) - pi[i], interval = c(.001, 10))$root

    Dlint <- expected_events_exponential(n_per_group, 24, lambda, lambda_cens, 36) * 2
    Dlfin <- expected_events_exponential(n_per_group, 48, lambda, lambda_cens, 36) * 2

    df$Dlint[i] <- Dlint
    df$Dlfin[i] <- Dlfin
    df$inf[i] <- Dlint / Dlfin

    design <- getDesignInverseNormal(
      kMax = 2,
      typeOfDesign = "P",
      informationRates = c(Dlint / Dlfin, 1),
      futilityBounds = 0,
      bindingFutility = TRUE
    )

    df$c1[i] <- design$criticalValues[1]
    df$c2[i] <- design$criticalValues[2]
    df$alpha1[i] <- design$alphaSpent[1]

    pwr <- getPowerSurvival(
      design = design,
      accrualTime = c(0, 36),
      dropoutTime = 12,
      dropoutRate1 = cens,
      dropoutRate2 = cens,
      eventTime = 12,
      pi1 = pi[i],
      pi2 = pi[i],
      maxNumberOfEvents = Dlfin,
      maxNumberOfSubjects = n_total
    )

    df$n_int[i] <- pwr$numberOfSubjects[1]
    df$int_t[i] <- pwr$analysisTime[1]
    df$fin_t[i] <- pwr$analysisTime[2]

    # Hybrid 1 boundaries
    boundaries_h1 <- get_hybrid1(Dl = Dlint, HR = HR, alpha1 = design$alphaSpent[1], n_int = df$n_int[i], Delta = 0.1)
    df$h1_u[i] <- boundaries_h1$h1_u
    df$h1_l[i] <- boundaries_h1$h1_l

    # Hybrid 2 boundaries
    event_prob <- (event_probability_exponential(n_per_group, 24, lambda, lambda_cens, 36) * 2) / 2
    df$h2_l[i] <- qbinom(0.1, n_total, event_prob, lower.tail = TRUE) - 1
    df$h2_u[i] <- qbinom(0.9, n_total, event_prob)

    # Hybrid 3 boundaries
    df$h3_l[i] <- floor(Dlint * 0.8)
    df$h3_u[i] <- ceiling(Dlint * 1.2)
  }

  # Save designs
  for (i in 1:nrow(df)) {
    w1 <- sqrt(df$inf[i])

    saveRDS(
      create_hybrid_design(qnorm(0.975), NULL, NULL, 0, 1e6, -1, fixed_par, 0),
      paste0(output_dir, "/fixed_", i, ".RDS")
    )

    saveRDS(
      create_hybrid_design(qnorm(0.975), df$c1[i], df$c2[i], 0, df$h1_u[i], df$h1_l[i], fixed_par, w1),
      paste0(output_dir, "/hybrid1_", i, ".RDS")
    )

    saveRDS(
      create_hybrid_design(qnorm(0.975), df$c1[i], df$c2[i], 0, df$h2_u[i], df$h2_l[i], fixed_par, w1),
      paste0(output_dir, "/hybrid2_", i, ".RDS")
    )

    saveRDS(
      create_hybrid_design(qnorm(0.975), df$c1[i], df$c2[i], 0, df$h3_u[i], df$h3_l[i], fixed_par, w1),
      paste0(output_dir, "/hybrid3_", i, ".RDS")
    )

    saveRDS(
      create_hybrid_design(qnorm(0.975), df$c1[i], df$c2[i], 0, -1, -1, fixed_par, w1),
      paste0(output_dir, "/gs_design_", i, ".RDS")
    )
  }
}

# Performance Study ---------------------------------------------------------

studies <- list(
  list(name = "Study 1", pi1 = 0.1, pi2 = 0.2, infoP = c(0.3, 1), infoOF = c(0.3, 1)),
  list(name = "Study 2", pi1 = 0.4, pi2 = 0.6, infoP = c(0.393, 1), infoOF = c(0.393, 1))
)

cens <- 0.01
fixed_par <- list(t_int = 24, min_fu = 12, sec_rec_length = 12)

design_types <- c("P", "OF") # Pocock and O'Brien-Fleming

for (study in studies) {
  for (dtype in design_types) {
    infoRates <- if (dtype == "P") study$infoP else study$infoOF
    design <- getDesignInverseNormal(
      kMax = 2,
      typeOfDesign = dtype,
      informationRates = infoRates,
      futilityBounds = 0,
      bindingFutility = TRUE
    )

    boundaries <- round(design$criticalValues, 3)
    w1 <- sqrt(design$informationRates[1])

    sample_size <- getSampleSizeSurvival(
      design = design,
      accrualTime = c(0, 36),
      dropoutTime = 12,
      dropoutRate1 = cens,
      dropoutRate2 = cens,
      eventTime = 12,
      pi1 = study$pi1,
      pi2 = study$pi2,
      followUpTime = 12
    )

    lambda1 <- uniroot(\(x) pexp(12, x) - study$pi1, interval = c(.001, 10))$root
    lambda2 <- uniroot(\(x) pexp(12, x) - study$pi2, interval = c(.001, 10))$root
    lambda_cens <- uniroot(\(x) pexp(12, x) - cens, interval = c(.00001, 10))$root

    n_per_group <- (ceiling(sample_size$maxNumberOfSubjects) +
      ceiling(sample_size$maxNumberOfSubjects) %% 2) / 2

    Dl <- expected_events_exponential(n_per_group, 24, lambda1, lambda_cens, 36) +
      expected_events_exponential(n_per_group, 24, lambda2, lambda_cens, 36)

    event_prob_interim <- (
      event_probability_exponential(n_per_group, 24, lambda1, lambda_cens, 36) +
        event_probability_exponential(n_per_group, 24, lambda2, lambda_cens, 36)
    ) / 2

    HR <- lambda1 / lambda2
    alpha1 <- design$alphaSpent[1]

    # Hybrid 1
    Delta <- if (dtype == "P") 0.1 else 0.05
    df <- get_hybrid1(
      Dl = Dl, HR = HR, alpha1 = alpha1, Delta = Delta,
      n_int = ceiling(sample_size$numberOfSubjects[1])
    )
    h1_l <- df$h1_l
    h1_u <- df$h1_u

    # Hybrid 2
    Delta <- 0.1
    h2_l <- qbinom(Delta, n_per_group * 2, event_prob_interim, lower.tail = TRUE) - 1
    h2_u <- qbinom(1 - Delta, n_per_group * 2, event_prob_interim)

    # Hybrid 3
    Delta <- 0.2
    h3_l <- floor(Dl - (Dl * Delta))
    h3_u <- ceiling(Dl + (Dl * Delta))

    # Save designs
    design_names <- c("fixed", "hybrid1", "hybrid2", "hybrid3", "gs_design")
    du_list <- list(10^6, h1_u, h2_u, h3_u, -1)
    dl_list <- list(-1, h1_l, h2_l, h3_l, -1)
    w1_list <- c(0, w1, w1, w1, w1)

    for (i in seq_along(design_names)) {
      des <- create_hybrid_design(
        p_cfix = qnorm(0.975),
        p_c1 = if (design_names[i] == "fixed") NULL else boundaries[1],
        p_c2 = if (design_names[i] == "fixed") NULL else boundaries[2],
        p_f1 = 0,
        p_du = du_list[[i]],
        p_dl = dl_list[[i]],
        p_fixPar = fixed_par,
        p_w1 = w1_list[i]
      )

      saveRDS(
        des,
        paste0("./Simulation Study/Designs/Performance/", study$name, "_", ifelse(dtype == "P", "Pocock", "OBrien-Fleming"), "/", design_names[i], ".RDS")
      )
    }
  }
}
