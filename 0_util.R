###############################################################################
# HELPER FUNCTIONS FOR THE SIMULATION CODE                                    #
###############################################################################

# Packages
install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

packages <- c(
  "data.table",  
  "magrittr",    
  "survival",    
  "rpact",
  "tidyr", 
  "plyr",        
  "dplyr",       
  "survminer"    
)

lapply(packages, install_and_load)

# TRIAL SIMULATION --------------------------------------------------------

#-----------------------------------------------------------
# Function: get_results_of_group
#-----------------------------------------------------------
# Purpose:
#   Simulate patient-level data for one treatment group 
#   given Weibull-distributed survival times, random censoring, 
#   and administrative censoring.
#
# Arguments:
#   p_n            - Number of patients in this group
#   p_Evpar        - List of event parameters:
#                      $shape       - Weibull shape parameter
#                      $med_surv    - Median survival time
#                      $lambda_cens - Rate of exponential censoring
#   p_acc_times    - Vector of patient-specific accrual times
#   p_ad_cens_time - Administrative censoring time (calendar time)
#
# Returns:
#   data.table with columns:
#     t_cal_ev   - Calendar time of event
#     t_cal_cens - Calendar time of censoring
#     t_acc      - Accrual time
#     t          - Observed follow-up time
#     status     - Event indicator (1 = event, 0 = censored)
#-----------------------------------------------------------
get_results_of_group <- function(p_n, p_Evpar, p_acc_times, p_ad_cens_time){
  
  # Simulate survival times (Weibull, scaled to median survival)
  surv_times <- rweibull(
    p_n,
    shape = p_Evpar[['shape']], 
    scale = p_Evpar[['med_surv']] / ((log(2))^(1 / p_Evpar[['shape']]))
  )
  
  # Censoring: exponential dropout + administrative censoring
  lfup_time <- rexp(p_n, rate = p_Evpar[['lambda_cens']]) + p_acc_times
  t_cal_cens <- pmin(lfup_time, p_ad_cens_time)
  
  # Build patient-level dataset
  data <- data.table(
    t_cal_ev   = surv_times + p_acc_times, 
    t_cal_cens = t_cal_cens, 
    t_acc      = p_acc_times,
    t          = pmin(surv_times, t_cal_cens - p_acc_times),
    status     = ifelse(surv_times + p_acc_times <= t_cal_cens, 1, 0)
  )
  
  return(data)
}


#-----------------------------------------------------------
# Function: simulate_stage
#-----------------------------------------------------------
# Purpose:
#   Simulate one trial stage including:
#     - Accrual process (Poisson arrivals, uniform accrual times)
#     - Randomized group allocation (coin flip)
#     - Patient outcomes for both groups
#
# Arguments:
#   p_acc_rate     - Expected number of patients recruited per unit time
#   p_max_acc_time - Length of accrual period
#   p_min_fu       - Minimum follow-up time after accrual ends
#   p_EvparI       - Event/censoring parameters for intervention group (see get_results_of_group)
#   p_EvparC       - Event/censoring parameters for control group
#
# Returns:
#   data.table with patient-level data including:
#     group   - Group assignment ("I" or "C")
#     pat_id  - Patient ID
#     (plus all columns returned by get_results_of_group)
#-----------------------------------------------------------
simulate_stage <- function(p_acc_rate, p_max_acc_time, p_min_fu, p_EvparI, p_EvparC){
  
  #--- Simulate accrual ---
  num_pat <- sum(rpois(p_max_acc_time, lambda = p_acc_rate))  # total patients
  acc_times <- sort(runif(n = num_pat, min = 0, max = p_max_acc_time))
  
  if(length(acc_times) < 2){
    print('Less than 2 patients recruited during accrual time.')
    return(NULL)
  }
  
  # Assign patients alternately to groups (coin flip for starting group)
  idx_acc_time_sec_group <- 2 * (1:(length(acc_times) %/% 2))
  if(runif(1) < 0.5){
    acc_timeI <- acc_times[idx_acc_time_sec_group]
    acc_timeC <- acc_times[-idx_acc_time_sec_group]
  } else {
    acc_timeI <- acc_times[-idx_acc_time_sec_group]
    acc_timeC <- acc_times[idx_acc_time_sec_group]
  }
  
  #--- Simulate outcomes for intervention & control ---
  dataI <- get_results_of_group(
    p_n = length(acc_timeI), p_Evpar = p_EvparI,
    p_acc_times = acc_timeI, p_ad_cens_time = p_max_acc_time + p_min_fu
  )
  dataC <- get_results_of_group(
    p_n = length(acc_timeC), p_Evpar = p_EvparC,
    p_acc_times = acc_timeC, p_ad_cens_time = p_max_acc_time + p_min_fu
  )
  
  # Add group labels
  dataI$group <- "I"
  dataC$group <- "C"
  
  # Merge and order by accrual time
  data <- rbind(dataI, dataC)
  data <- data[order(t_acc), ]
  data$pat_id <- 1:nrow(data)
  
  return(data)
}


#-----------------------------------------------------------
# Function: get_results_at_time
#-----------------------------------------------------------
# Purpose:
#   Evaluate patient outcomes at a given calendar time,
#   applying interim/final censoring as needed.
#
# Arguments:
#   p_data_res - Patient-level data from simulate_stage()
#   p_t        - Calendar time at which results are evaluated
#   p_end_rec  - End of accrual period
#
# Returns:
#   data.table with updated patient outcomes:
#     t_cal_cens - Updated censoring time (truncated at p_t)
#     t          - Updated follow-up time
#     status     - Updated event indicator
#     end_rec    - End of accrual
#     time_eval  - Evaluation time (p_t)
#-----------------------------------------------------------
get_results_at_time <- function(p_data_res, p_t, p_end_rec){
  
  data_t <- copy(p_data_res)
  
  # Only include patients recruited before end of accrual
  data_t <- data_t[t_acc < p_end_rec, ]
  
  # Apply administrative censoring at time p_t
  data_t[, t_cal_cens := pmin(t_cal_cens, p_t)] 
  
  # Update follow-up and event status
  data_t[, t := pmin(t_cal_ev, t_cal_cens) - t_acc]
  data_t[, status := ifelse(t_cal_ev <= t_cal_cens, 1, 0)]
  data_t[, end_rec := p_end_rec]
  data_t[, time_eval := p_t]
  
  return(data_t)
}


# TEST STATISTICS ---------------------------------------------------------

#-----------------------------------------------------------
# Function: logrank_func
#-----------------------------------------------------------
# Purpose:
#   Extract a signed log-rank test statistic from the output 
#   of the survival::survdiff() function.
#
# Background:
#   survdiff() returns a chi-squared statistic (Z^2).
#   This function recovers the signed Z-statistic.
#
# Arguments:
#   p_survdiff - Object returned by survival::survdiff().
#
# Returns:
#   A list with:
#     $log_rank_stat - Signed log-rank test statistic (Z)
#     $p_value       - Corresponding p-value
#     $n_event       - Total number of observed events
#-----------------------------------------------------------
logrank_func <- function(p_survdiff) {
  
  # Determine sign based on expected vs. observed counts
  sgns <- sign(p_survdiff$exp - p_survdiff$obs)
  
  # Sanity check: groups should have opposite signs
  if (sgns[1] != (-sgns[2])) {
    stop("Cannot determine direction of Z-statistic.")
  }
  
  # Recover signed log-rank statistic
  log_rank_stat <- sgns[2] * sqrt(p_survdiff$chisq)
  
  # Extract p-value and total number of events
  p_value <- p_survdiff$pvalue
  n_event <- sum(p_survdiff$obs)
  
  results <- list(
    log_rank_stat = log_rank_stat,
    p_value = p_value,
    n_event = n_event
  )
  
  return(results)
}


#-----------------------------------------------------------
# Function: indep_inc_func
#-----------------------------------------------------------
# Purpose:
#   Compute a combined log-rank test statistic using 
#   the independent increments method and inverse-normal method.
#
# Arguments:
#   p_lr_fin - Log-rank test statistic at the final analysis
#   p_lr_int - Log-rank test statistic at the  interim analysis
#   p_w1     - Inverse-normal combination test weight
#
# Returns:
#   A list with:
#     $invN_stat - Inverse-normal combination test statistic
#-----------------------------------------------------------
indep_inc_func <- function(p_lr_fin, p_lr_int, p_w1){
  
  # Compute independent increments statistic
  comb_stat <- (
    sqrt(p_lr_fin$n_event) * p_lr_fin$log_rank_stat -
      sqrt(p_lr_int$n_event) * p_lr_int$log_rank_stat
  ) / sqrt(p_lr_fin$n_event - p_lr_int$n_event)
  
  # Inverse-normal combination of interim + increment
  invN_stat <- p_w1 * p_lr_int$log_rank_stat +
    sqrt(1 - p_w1^2) * comb_stat
  
  results <- list(invN_stat = invN_stat)
  
  return(results)
}

# SURVIVAL PARAMETERS -----------------------------------------------------

#-----------------------------------------------------------
# Function: calc_hazard_ratio
#-----------------------------------------------------------
# Purpose:
#   Calculate the hazard ratio between intervention and control arms,
#   assuming both follow exponential survival distributions.
#
# Arguments:
#   p_med_survI - Median survival time in the intervention group
#   p_med_survC - Median survival time in the control group
#
# Returns:
#   Hazard ratio, which simplifies to
#   HR = median_survival_C / median_survival_I.
#-----------------------------------------------------------
calc_hazard_ratio <- function(p_med_survI, p_med_survC){
  return(p_med_survC / p_med_survI)
}

#-----------------------------------------------------------
# Function: med_surv_to_one_year_event_prob
#-----------------------------------------------------------
# Purpose:
#   Convert a median survival time into the 1-year event probability,
#   assuming exponential survival distribution.
#
# Arguments:
#   p_med_surv - Median survival time (in months)
#
# Returns:
#   Probability of experiencing an event within 12 months.
#
# Formula:
#   S(t) = 0.5^(t / median_survival)
#   => 1 - S(12) gives 1-year event probability.
#-----------------------------------------------------------
med_surv_to_one_year_event_prob <- function(p_med_surv){
  return(1 - 0.5^(12 / p_med_surv))
}


#-----------------------------------------------------------
# Function: one_year_event_prob_to_med_surv
#-----------------------------------------------------------
# Purpose:
#   Convert a 1-year event probability into the corresponding
#   median survival time, assuming exponential distribution.
#
# Arguments:
#   p_one_year_ev_prob - Probability of an event within 12 months
#
# Returns:
#   Median survival time (in months).
#
# Formula:
#   median_survival = 12 * log(0.5) / log(S(12))
#-----------------------------------------------------------
one_year_event_prob_to_med_surv <- function(p_one_year_ev_prob){
  return(12 * log(0.5) / log(1 - p_one_year_ev_prob))
}

#-----------------------------------------------------------
# Function: calc_one_year_event_prob_I
#-----------------------------------------------------------
# Purpose:
#   Compute the 1-year event probability for the intervention group,
#   given a hazard ratio and the control group’s 1-year event probability.
#
# Arguments:
#   p_hazard_ratio - Hazard ratio (λI / λC)
#   p_1yearprob    - 1-year event probability in the control group
#
# Returns:
#   1-year event probability in the intervention group.
#
# Formula:
#   Let xc = monthly survival prob. in control = (1 - p_1yearprob)^(1/12)
#   Intervention monthly survival = xi = 1 - HR * (1 - xc)
#   Intervention 1-year survival = xi^12
#   Event probability = 1 - xi^12
#-----------------------------------------------------------
calc_one_year_event_prob_I <- function(p_hazard_ratio, p_1yearprob){
  xc <- (1 - p_1yearprob)^(1/12)      # monthly survival in control group
  xi <- 1 - p_hazard_ratio * (1 - xc) # monthly survival in intervention group
  return(1 - xi^12)                   # 1-year event probability
}


# STUDY DESIGN ------------------------------------------------------------

#-----------------------------------------------------------
# Function: create_hybrid_design
#-----------------------------------------------------------
# Purpose:
#   Creates a hybrid group-sequential design object that incorporates
#   rules for unblinding, rejection decisions, and second-stage 
#   recruitment length. 
#
# Arguments:
#   p_cfix   - Critical value for blinded analysis (fixed boundary)
#   p_c1     - Critical value for interim (first-stage) log-rank test
#   p_c2     - Critical value for second-stage log-rank test (after unblinding)
#   p_f1     - Futility boundary for interim analysis
#   p_du     - Upper boundary for number of events triggering unblinding
#   p_dl     - Lower boundary for number of events triggering unblinding
#   p_fixPar - List of fixed parameters, containing:
#                $t_int          - Calendar time of interim analysis
#                $min_fu         - Minimum follow-up time
#                $sec_rec_length - Default second-stage recruitment length
#   p_w1     - Weight for combination test statistic
#
# Returns:
#   A list (the "hybrid design object") containing:
#     $env           - The environment storing design variables
#     $t_int         - Interim calendar time
#     $min_fu        - Minimum follow-up time
#     $cfix, $c1, 
#     $c2, $f1       - Critical and futility boundaries
#     $du, $dl       - Upper/lower event thresholds for unblinding
#     $w1            - Weight for test combination
#     $sec_rec_length- Default recruitment length after interim
#     $dec_rec_dur   - Function for deciding second-stage recruitment duration
#     $rej_func      - Function for applying rejection rules
#     $is_unblinded  - Function that determines if unblinding occurs
#
#-----------------------------------------------------------
create_hybrid_design <- function(p_cfix, p_c1, p_c2, p_f1, p_du, p_dl, p_fixPar, p_w1){
  
  # Store design-specific environment (keeps functions self-contained)
  design_env <- environment()
  
  # --- Fixed parameters (do not change within simulations) ---
  t_int <- p_fixPar$t_int          # Calendar time for interim analysis
  min_fu <- p_fixPar$min_fu        # Minimum follow-up time
  sec_rec_length <- p_fixPar$sec_rec_length # Default second-stage recruitment length
  
  # --- Flexible parameters (design-specific critical values, boundaries, etc.) ---
  cfix <- p_cfix   # Blinded critical value
  c1   <- p_c1     # Interim boundary
  c2   <- p_c2     # Final boundary (after unblinding)
  f1   <- p_f1     # Futility boundary
  du   <- p_du     # Upper threshold for number of events (unblinding trigger)
  dl   <- p_dl     # Lower threshold for number of events (unblinding trigger)
  w1   <- p_w1     # Weight for test combination
  
  #-----------------------------------------------------------
  # Function: is_unblinded
  # Purpose: Decide whether trial is unblinded based on observed events
  # Input: p_nEv - number of observed events at interim
  # Output: TRUE if trial is unblinded, FALSE otherwise
  #-----------------------------------------------------------
  is_unblinded <- function(p_nEv){
    if(p_nEv > p_du | p_nEv < p_dl){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  #-----------------------------------------------------------
  # Function: rej_func
  # Purpose: Apply rejection rules to interim and final statistics
  # Input: p_stats_all - data.frame containing interim/final statistics
  #        (must include columns: int_n_ev, sec_lr_stat, int_lr_stat)
  # Output: p_stats_all with added rejection indicators:
  #         $sec_rej (secondary rejection) and $int_rej (interim rejection)
  #-----------------------------------------------------------
  rej_func <- function(p_stats_all){
    
    p_stats_all$sec_rej <- 0
    p_stats_all$int_rej <- 0
    
    # --- Case 1: no unblinding ---
    idx_blind <- which(p_stats_all$int_n_ev >= p_dl & p_stats_all$int_n_ev <= p_du)
    p_stats_all[idx_blind, ]$sec_rej <- ifelse(
      p_stats_all[idx_blind, ]$sec_lr_stat > p_cfix, 1, 0
    )
    
    # --- Case 2: unblinding ---
    idx_unbl <- setdiff(1:nrow(p_stats_all), idx_blind)
    p_stats_all[idx_unbl, ]$sec_rej <- ifelse(
      p_stats_all[idx_unbl, ]$sec_lr_stat > p_c2, 1, 0
    )
    p_stats_all[idx_unbl, ]$int_rej <- ifelse(
      p_stats_all[idx_unbl, ]$int_lr_stat > p_c1, 1, 0
    )
    
    return(p_stats_all)
  }
  
  #-----------------------------------------------------------
  # Function: rec_func
  # Purpose: Determine second-stage recruitment length
  # Input: 
  #   p_nEv        - number of observed events at interim
  #   p_logrank_int- interim log-rank statistic
  #   p_env        - design environment (default: current environment)
  # Output:
  #   Recruitment length (0 if trial stops, sec_rec_length if continues)
  #-----------------------------------------------------------
  rec_func <- function(p_nEv, p_logrank_int, p_env = design_env){
    
    # --- Case 1: early unblinding (too few events) ---
    if(p_nEv < p_dl){ 
      if(p_logrank_int < p_f1){      # Futility
        return(0)
      } else if(p_logrank_int > p_c1){ # Efficacy
        return(0)
      } else {                       # Continue
        return(p_fixPar$sec_rec_length)
      }
    }
    
    # --- Case 2: blinded scenario ---
    if(p_nEv >= p_dl & p_nEv <= p_du){ 
      return(p_fixPar$sec_rec_length)
    }
    
    # --- Case 3: late unblinding (too many events) ---
    if(p_nEv > p_du){ 
      if(p_logrank_int < p_f1){      # Futility
        return(0)
      }
      if(p_logrank_int > p_c1){      # Efficacy
        return(0)
      }
      if(p_logrank_int < p_c1){      # Continue
        return(p_fixPar$sec_rec_length)
      }
    }
  }
  
  #-----------------------------------------------------------
  # Return hybrid design object
  #-----------------------------------------------------------
  hybrid <- list(
    env = design_env,
    t_int = t_int,
    min_fu = min_fu,
    cfix = cfix,
    c1 = c1,
    c2 = c2,
    f1 = f1,
    du = du,
    dl = dl,
    w1 = w1,
    sec_rec_length = sec_rec_length,
    dec_rec_dur = rec_func,     # Function for second-stage recruitment duration
    rej_func = rej_func,        # Function for rejection rules
    is_unblinded = is_unblinded # Function for unblinding decision
  )
  
  return(hybrid)
}

#-----------------------------------------------------------
# Function: hybrid1
#-----------------------------------------------------------
# Purpose:
#   Compute the lower (h1_l) and upper (h1_u) boundaries for 
#   the hybrid design 1, based on design parameters.
#
# Arguments:
#   Dl     - Expected number of events at interim
#   HR     - Hazard ratio (treatment vs. control)
#   alpha1 - One-sided alpha spent at interim
#   n_int  - Planned number of patients at interim
#   Delta  - Hybrid design 1 design parameter
#
# Returns:
#   A data frame with:
#     Delta - The design parameter Delta
#     h1_l  - Lower boundary for the number of events
#     h1_u  - Upper boundary for the number of events
#
# Notes:
#   - Boundaries define the region where the study continues 
#     without unblinding. Outside of this range, unblinding occurs.
#   - If calculation fails (NaN), defaults are used:
#       h1_u → n_int (planned interim patients)
#       h1_l → 0
#-----------------------------------------------------------

get_hybrid1 <- function(Dl, HR, alpha1, n_int, Delta) {
  
  tmp <- sqrt(0.25 * Dl * log(HR)^2) - qnorm(1 - alpha1)
  
  # Upper boundary calculation
  u_calc <- (qnorm(pnorm(tmp) + Delta) + qnorm(1 - alpha1))^2 * (4 / log(HR)^2)
  h1_u <- ifelse(is.na(ceiling(u_calc)), n_int, ceiling(u_calc))
  
  # Lower boundary calculation
  l_calc <- (qnorm(pnorm(tmp) - Delta) + qnorm(1 - alpha1))^2 * (4 / log(HR)^2)
  h1_l <- ifelse(is.na(ceiling(l_calc)), 0, floor(l_calc))
  
  # Return tidy results
  return(data.frame(Delta = Delta, h1_l = h1_l, h1_u = h1_u))
}


# EVENT ESTIMATION --------------------------------------------------------

#-----------------------------------------------------------
# Compute the probability of an event occurring at a given time
#-----------------------------------------------------------
# Arguments:
# - event_distribution: function for the event time distribution (e.g., density function of survival times)
# - censoring_survival_function: survival function for censoring times (P(not censored up to time x))
# - accrual_weighting_function: weighting function describing the probability of a subject being recruited 
#                               at a certain time relative to the accrual period
# - time: calendar time at which the probability is evaluated
#
# Returns:
# - Probability that an event occurs at time "time" for a randomly recruited subject,
#   taking into account censoring and accrual.
event_probability <- function(event_distribution,
                              censoring_survival_function,
                              accrual_weighting_function,
                              time){
  integrate(
    \(x) event_distribution(x) *
      censoring_survival_function(x) * 
      accrual_weighting_function(time - x),
    0,
    Inf
  )$value
}

#-----------------------------------------------------------
# Compute the expected number of events in the study at a given time
#-----------------------------------------------------------
# Arguments:
# - event_distribution: see above
# - censoring_survival_function: see above
# - accrual_weighting_function: see above
# - overall_sample_size: total number of recruited subjects
# - time: calendar time at which expected events are calculated
#
# Returns:
# - Expected number of events at time "time" across all subjects.
expected_events <- function(event_distribution,
                            censoring_survival_function,
                            accrual_weighting_function,
                            overall_sample_size,
                            time){
  overall_sample_size * event_probability(event_distribution,
                                          censoring_survival_function,
                                          accrual_weighting_function,
                                          time)
}

#-----------------------------------------------------------
# Event probability under exponential survival and censoring distributions
#-----------------------------------------------------------
# Arguments:
# - overall_sample_size: total number of recruited subjects (not directly used here)
# - time: calendar time at which probability is evaluated
# - rate_event: event rate parameter of the exponential distribution
# - rate_censoring: censoring rate parameter of the exponential distribution
# - end_recruitment_time: maximum accrual time (subjects recruited uniformly over [0, end_recruitment_time])
#
# Returns:
# - Probability of an event at time "time" under exponential assumptions.
event_probability_exponential <- function(overall_sample_size,
                                          time,
                                          rate_event,
                                          rate_censoring,
                                          end_recruitment_time){
  event_probability(event_distribution = \(x) dexp(x, rate = rate_event),
                    censoring_survival_function = \(x) (1 - pexp(x, rate = rate_censoring)),
                    accrual_weighting_function = \(x) punif(x, 0, end_recruitment_time),
                    time
  )
}

#-----------------------------------------------------------
# Expected number of events under exponential assumptions
#-----------------------------------------------------------
# Arguments:
# - overall_sample_size: total number of recruited subjects
# - time: calendar time at which expected events are calculated
# - rate_event: event rate parameter of the exponential distribution
# - rate_censoring: censoring rate parameter of the exponential distribution
# - end_recruitment_time: maximum accrual time (subjects recruited uniformly over [0, end_recruitment_time])
#
# Returns:
# - Expected number of events at time "time" under exponential assumptions.
expected_events_exponential <- function(overall_sample_size,
                                        time,
                                        rate_event,
                                        rate_censoring,
                                        end_recruitment_time){
  expected_events(event_distribution = \(x) dexp(x, rate = rate_event),
                  censoring_survival_function = \(x) (1 - pexp(x, rate = rate_censoring)),
                  accrual_weighting_function = \(x) punif(x, 0, end_recruitment_time),
                  overall_sample_size,
                  time
  )
}