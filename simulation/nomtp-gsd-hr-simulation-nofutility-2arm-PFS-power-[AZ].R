# ============================================================
# Code title: Two-arm PFS simulation without futility
#
# mtp:        none
# gsd:        calendar-time group-sequential efficacy monitoring
# coc:        no
# cop:        no
# simulation: yes
# futility:   no
# endpoints:  PFS under a piecewise-exponential event-time model
# calculation: power
#
# Key inputs: enrollment, dropout, arm allocation, event-time distributions,
#   treatment effect, calendar-time looks, efficacy spending, futility rule,
#   number of simulations, and random-number seed
#
# Key outputs: efficacy and futility operating characteristics, planned
#   analyses, timeline and enrollment, and observed-effect distribution
#
# Engine: base R, survival, gsDesign, gsDesign2, graphicalMCP
# Input:  jsonlite - fromJSON() (trial1, mtp1, gsd1)
# ============================================================

# ==== READING THE DESIGN INPUTS =============================================

design_input <- jsonlite::fromJSON(
  txt = "/Users/oliviazhang/Desktop/AZ_templates/simulation/nomtp-gsd-hr-simulation-nofutility-2arm-PFS-power-[AZ].json",
  simplifyVector = TRUE,
  simplifyDataFrame = FALSE,
  simplifyMatrix = TRUE)
trial1 <- design_input$trial
mtp1   <- design_input$mtp
gsd1   <- design_input$gsd

# ==== RUNNING THE CALCULATION ================================================

# ==== trial: enrollment, dropout, arm allocation, and PFS assumptions =========

# Resolve the hypothesis, endpoint, and GSD configurations
hypothesis_ids <- names(mtp1$nodes)
hypothesis_id <- hypothesis_ids[[1]]
hypothesis_config <- mtp1$nodes[[hypothesis_id]]
endpoint_id <- hypothesis_config$endpoint
endpoint_config <- trial1$endpoints[[endpoint_id]]
gsd_config <- gsd1[[hypothesis_id]]
arm_ids <- names(trial1$arms)
outcome_models <- trial1$sim$outcome_models

# Read simulation controls and convert allocation weights to subject counts
number_of_simulations <- as.integer(trial1$sim$n_sim)
simulation_seed <- as.integer(trial1$sim$seed)

# Use one design-parameter row; calendar-time values define within-trial looks
param_grid <- data.frame(
  design = "specified")

arm_randomization_weights <- vapply(
  X = trial1$arms,
  FUN = function(arm_config) arm_config$randomization_weight,
  FUN.VALUE = numeric(1))
total_subjects <- as.integer(trial1$enrollment$n_accrual_max)
target_subjects_by_arm <-
  total_subjects * arm_randomization_weights / sum(arm_randomization_weights)
subjects_by_arm <- floor(target_subjects_by_arm)
unassigned_subjects <- total_subjects - sum(subjects_by_arm)
if (unassigned_subjects > 0) {
  arms_receiving_extra_subject <- order(
    x = target_subjects_by_arm - subjects_by_arm,
    decreasing = TRUE)[seq_len(unassigned_subjects)]
  subjects_by_arm[arms_receiving_extra_subject] <-
    subjects_by_arm[arms_receiving_extra_subject] + 1L
}
subjects_by_arm <- as.integer(subjects_by_arm)
names(subjects_by_arm) <- arm_ids

# Derive the intervention-to-comparator allocation ratio
allocation_ratio <-
  trial1$arms[[hypothesis_config$intervention]]$randomization_weight /
  trial1$arms[[hypothesis_config$comparator]]$randomization_weight
# Convert cumulative dropout probability to an exponential rate
dropout_rate <- -log(1 - trial1$dropout$probability) /
  trial1$dropout$window

# Read the endpoint's arm-specific simulation model
outcome_model <- outcome_models[[endpoint_id]]
outcome_distributions <- outcome_model$distributions
outcome_effects <- outcome_model$effects
comparator_distribution <-
  outcome_distributions[[hypothesis_config$comparator]]
comparator_cuts <- as.numeric(comparator_distribution$cut$at)
comparator_rates <- if (!is.null(comparator_distribution$rate)) {
  as.numeric(comparator_distribution$rate)
} else {
  log(2) / as.numeric(comparator_distribution$median)
}
intervention_effect <-
  outcome_effects[[hypothesis_config$intervention]]
intervention_cuts <- sort(unique(c(
  comparator_cuts,
  as.numeric(intervention_effect$cut$at))))
intervention_rates <-
  comparator_rates[findInterval(
    x = intervention_cuts,
    vec = comparator_cuts)] *
  as.numeric(intervention_effect$value)[findInterval(
    x = intervention_cuts,
    vec = as.numeric(intervention_effect$cut$at))]
failure_models <- list()
failure_models[[hypothesis_config$comparator]] <- list(
  cuts = comparator_cuts,
  rates = comparator_rates)
failure_models[[hypothesis_config$intervention]] <- list(
  cuts = intervention_cuts,
  rates = intervention_rates)

# ==== mtp: one-node graph without multiplicity adjustment ====================

# Construct the one-node graph using the same representation as a full MTP
familywise_alpha <- trial1$alpha$value
initial_local_alphas <- vapply(
  X = mtp1$nodes,
  FUN = function(node_config) node_config$initial_alpha,
  FUN.VALUE = numeric(1))
names(initial_local_alphas) <- hypothesis_ids
alpha_transfer_matrix <- matrix(
  data = 0,
  nrow = length(hypothesis_ids),
  ncol = length(hypothesis_ids),
  dimnames = list(
    hypothesis_ids,
    hypothesis_ids))
mtp_graph <- graphicalMCP::graph_create(
  hypotheses = unname(initial_local_alphas / familywise_alpha),
  transitions = unname(alpha_transfer_matrix),
  hyp_names = hypothesis_ids)

# Enumerate the local-alpha values reachable through graph transitions
graph_weights <- graphicalMCP::graph_generate_weights(
  graph = mtp_graph)
alpha_states <- setNames(
  object = lapply(
    X = seq_along(hypothesis_ids),
    FUN = function(hypothesis_index) {
      hypothesis_is_active <- graph_weights[, hypothesis_index] == 1
      reachable_alpha <- familywise_alpha *
        graph_weights[hypothesis_is_active,
                      length(hypothesis_ids) + hypothesis_index]
      sort(unique(reachable_alpha[reachable_alpha > 0]))
    }),
  nm = hypothesis_ids)

# ==== gsd: calendar-time group-sequential design =============================

# Read any simulation-only calendar times
additional_calendar_times <- NULL
if (!is.null(trial1$sim$param_grid)) {
  grid_parameters <- vapply(
    X = trial1$sim$param_grid,
    FUN = function(grid_config) grid_config$parameter,
    FUN.VALUE = character(1))
  calendar_grid_index <- which(grid_parameters == "calendar-time")
  if (length(calendar_grid_index) > 0) {
    additional_calendar_times <- as.numeric(
      trial1$sim$param_grid[[calendar_grid_index[[1]]]]$values)
  }
}

# Combine planned and simulation-only analyses in chronological order
analysis_times <- sort(unique(c(
  as.numeric(gsd_config$looks$schedule$values),
  additional_calendar_times)))
n_looks <- length(analysis_times)

# Map planned efficacy and futility looks to the combined analysis grid
planned_looks <- match(
  x = as.numeric(gsd_config$looks$schedule$values),
  table = analysis_times)
efficacy_looks <- planned_looks[
  as.integer(gsd_config$efficacy$active_at)]
n_efficacy_looks <- length(efficacy_looks)
efficacy_times <-
  analysis_times[efficacy_looks]
futility_enabled <- isTRUE(gsd_config$futility$enabled)
futility_looks <- if (futility_enabled) {
  planned_looks[as.integer(gsd_config$futility$active_at)]
} else {
  integer(0)
}

# Resolve the configured gsDesign spending function
efficacy_spending_function <- getExportedValue(
  ns = "gsDesign",
  name = gsd_config$efficacy$family$fn)

# Convert enrollment assumptions to the gsDesign2 interval format
enrollment_month <- seq_len(trial1$enrollment$period)
planned_cumulative_enrollment <-
  (enrollment_month / trial1$enrollment$period)^trial1$enrollment$k *
  total_subjects
planned_enrollment_model <- data.frame(
  stratum = "All",
  duration = rep(
    x = 1,
    times = length(enrollment_month)),
  rate = planned_cumulative_enrollment - c(
    0,
    head(
      x = planned_cumulative_enrollment,
      n = -1)))

# Align comparator and intervention hazards on common intervals
planned_failure_cuts <- sort(unique(c(
  failure_models[[hypothesis_config$comparator]]$cuts,
  failure_models[[hypothesis_config$intervention]]$cuts)))
planned_comparator_rates <-
  failure_models[[hypothesis_config$comparator]]$rates[findInterval(
    x = planned_failure_cuts,
    vec = failure_models[[hypothesis_config$comparator]]$cuts)]
planned_intervention_rates <-
  failure_models[[hypothesis_config$intervention]]$rates[findInterval(
    x = planned_failure_cuts,
    vec = failure_models[[hypothesis_config$intervention]]$cuts)]
planned_failure_model <- data.frame(
  stratum = "All",
  duration = c(
    diff(planned_failure_cuts),
    max(analysis_times) + 1),
  fail_rate = planned_comparator_rates,
  hr = planned_intervention_rates / planned_comparator_rates,
  dropout_rate = dropout_rate)
# Calculate planned information at each efficacy look
planned_information <- gsDesign2::gs_info_ahr(
  enroll_rate = planned_enrollment_model,
  fail_rate = planned_failure_model,
  ratio = allocation_ratio,
  analysis_time = efficacy_times)
planned_information_fraction <-
  planned_information$info0 /
  tail(
    x = planned_information$info0,
    n = 1)

# Construct fixed GSD boundaries for every reachable local alpha
gsd_bounds <- setNames(
  object = vector(
    mode = "list",
    length = length(hypothesis_ids)),
  nm = hypothesis_ids)
gsd_bounds[[hypothesis_id]] <- do.call(
  what = rbind,
  args = lapply(
    X = alpha_states[[hypothesis_id]],
    FUN = function(local_alpha) {
      gsd_arguments <- list(
        k = n_efficacy_looks,
        test.type = 1,
        alpha = local_alpha / hypothesis_config$test_sides,
        sfu = efficacy_spending_function)
      if (n_efficacy_looks > 1) {
        gsd_arguments$timing <- head(
          x = planned_information_fraction,
          n = -1)
      }
      if (!is.null(gsd_config$efficacy$family$param)) {
        gsd_arguments$sfupar <- gsd_config$efficacy$family$param
      }

      local_gsd <- do.call(
        what = gsDesign::gsDesign,
        args = gsd_arguments)
      z_boundary <- as.numeric(local_gsd$upper$bound)
      data.frame(
        local_alpha = local_alpha,
        look = seq_len(n_efficacy_looks),
        analysis_index = efficacy_looks,
        z = z_boundary,
        critical_effect = exp(
          -z_boundary / sqrt(planned_information$info0)),
        nominal_alpha =
          hypothesis_config$test_sides * stats::pnorm(-z_boundary),
        cumulative_alpha =
          hypothesis_config$test_sides *
            cumsum(as.numeric(local_gsd$upper$spend)))
    }))

# Select the initially allocated boundaries for reporting and simulation
initial_gsd_bounds <- gsd_bounds[[hypothesis_id]][
  gsd_bounds[[hypothesis_id]]$local_alpha ==
    initial_local_alphas[[hypothesis_id]], ]

# ==== simulation =============================================================

# ---- Data generation --------------------------------------------------------
draw_piecewise_exponential <- function(n, cuts, rates) {
  # Map exponential draws through the piecewise cumulative hazard
  cuts <- as.numeric(cuts)
  rates <- as.numeric(rates)
  cumulative_hazard <- c(
    0,
    cumsum(head(
      x = rates,
      n = -1) * diff(cuts)))
  exponential_draw <- stats::rexp(
    n = n)
  interval <- findInterval(
    x = exponential_draw,
    vec = cumulative_hazard)
  cuts[interval] +
    (exponential_draw - cumulative_hazard[interval]) / rates[interval]
}

# simulate_PFS: generate one patient-level PFS dataset
simulate_PFS <- function() {
  # Assign treatment and enrollment times
  assigned_arm <- sample(
    x = rep(
      x = arm_ids,
      times = subjects_by_arm))
  enrollment_time <- (
    stats::runif(
      n = total_subjects) *
      trial1$enrollment$period^trial1$enrollment$k) ^
    (1 / trial1$enrollment$k)

  # Generate arm-specific event times
  failure_time <- numeric(total_subjects)
  for (arm_id in arm_ids) {
    selected_subjects <- which(assigned_arm == arm_id)
    failure_model <- failure_models[[arm_id]]
    failure_time[selected_subjects] <- draw_piecewise_exponential(
      n = length(selected_subjects),
      cuts = failure_model$cuts,
      rates = failure_model$rates)
  }

  # Apply dropout and derive observed follow-up
  dropout_time <- stats::rexp(
    n = total_subjects,
    rate = dropout_rate)
  event_observed <- failure_time <= dropout_time
  follow_up_time <- pmin(
    failure_time,
    dropout_time)

  sim_data <- list(
    assigned_arm = assigned_arm,
    enrollment_time = enrollment_time,
    observed_calendar_time = enrollment_time + follow_up_time,
    event_observed = event_observed)
  sim_data
}

# ---- Endpoint test ----------------------------------------------------------

# test_PFS: calculate the logrank p-value and observed hazard ratio at one look
test_PFS <- function(
    sim_data, hypothesis_id, hypothesis_config, look_index, analysis_time) {
  # Build the dataset available at this analysis time
  enrolled_subjects <-
    sim_data$enrollment_time < analysis_time
  assigned_arm <- sim_data$assigned_arm[enrolled_subjects]
  observed_follow_up_time <- pmax(
    0,
    pmin(
      sim_data$observed_calendar_time[enrolled_subjects],
      analysis_time) -
      sim_data$enrollment_time[enrolled_subjects])
  event_observed <- sim_data$event_observed[enrolled_subjects] &
    sim_data$observed_calendar_time[enrolled_subjects] <=
      analysis_time
  intervention_indicator <-
    as.integer(assigned_arm == hypothesis_config$intervention)

  test_results <- data.frame(
    hypothesis = hypothesis_id,
    look = look_index,
    n = length(intervention_indicator),
    events = sum(event_observed),
    estimate = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    row.names = hypothesis_id)

  # Skip testing when the available information is insufficient
  if (
    length(unique(intervention_indicator)) < 2 ||
      sum(event_observed) < 2
  ) {
    return(test_results)
  }

  # Calculate the logrank p-value and Cox-model hazard ratio
  data_at_analysis <- data.frame(
    observed_follow_up_time = observed_follow_up_time,
    event_observed = as.integer(event_observed),
    intervention_indicator = intervention_indicator)

  logrank_test <- tryCatch(
    expr = survival::survdiff(
      formula = survival::Surv(
        time = observed_follow_up_time,
        event = event_observed) ~ intervention_indicator,
      data = data_at_analysis),
    error = function(e) NULL)
  cox_model <- tryCatch(
    expr = suppressWarnings(survival::coxph(
      formula = survival::Surv(
        time = observed_follow_up_time,
        event = event_observed) ~ intervention_indicator,
      data = data_at_analysis)),
    error = function(e) NULL)

  test_statistic <- if (is.null(logrank_test)) NA_real_ else
    unname(logrank_test$chisq)
  p_value <- if (is.na(test_statistic)) NA_real_ else
    stats::pchisq(
      q = logrank_test$chisq,
      df = 1,
      lower.tail = FALSE)
  effect_estimate <- if (
    is.null(cox_model) || length(stats::coef(cox_model)) == 0
  ) NA_real_ else exp(unname(stats::coef(cox_model)[[1]]))

  test_results$estimate <- effect_estimate
  test_results$statistic <- test_statistic
  test_results$p_value <- p_value
  test_results
}

# ---- Initialize simulation results ------------------------------------------

# Create common simulation-by-look-by-hypothesis storage templates
result_dimensions <- c(
  number_of_simulations,
  n_looks,
  length(hypothesis_ids))
result_dimnames <- list(
  simulation = NULL,
  look = seq_len(n_looks),
  hypothesis = hypothesis_ids)
empty_numeric_results <- array(
  data = NA_real_,
  dim = result_dimensions,
  dimnames = result_dimnames)
empty_logical_results <- array(
  data = FALSE,
  dim = result_dimensions,
  dimnames = result_dimnames)

# Allocate storage for every parameter-grid value
sim_store <- vector(
  mode = "list",
  length = nrow(param_grid))

# ---- Run simulated trials, apply the design, and save results ---------------

# Run data generation, testing, and design application in that order
for (grid_index in seq_len(nrow(param_grid))) {
  # Initialize results for this parameter-grid value
  sim_store[[grid_index]] <- list(
    test_results = list(
      n = empty_numeric_results,
      events = empty_numeric_results,
      estimate = empty_numeric_results,
      statistic = empty_numeric_results,
      p_value = empty_numeric_results),
    trial_decisions = list(
      local_alpha = empty_numeric_results,
      adjusted_p = empty_numeric_results,
      efficacy = empty_logical_results,
      futility = empty_logical_results,
      rejected = empty_logical_results))

  # Set the random-number seed once for this parameter-grid value
  set.seed(
    seed = simulation_seed)

  for (simulation_index in seq_len(number_of_simulations)) {
    # Generate all correlated endpoint data for one trial
    sim_data <- simulate_PFS()

    # Initialize the graph and hypothesis status for this trial
    mtp_state <- mtp_graph
    hypothesis_stopped <- setNames(
      object = rep(
        x = FALSE,
        times = length(hypothesis_ids)),
      nm = hypothesis_ids)

    for (look_index in seq_len(n_looks)) {
      analysis_time <- analysis_times[[look_index]]

      # Calculate all local tests scheduled at this analysis
      test_results <- do.call(
        what = rbind,
        args = lapply(
          X = hypothesis_ids,
          FUN = function(hypothesis_id) {
            hypothesis_config <- mtp1$nodes[[hypothesis_id]]
            test_PFS(
              sim_data = sim_data,
              hypothesis_id = hypothesis_id,
              hypothesis_config = hypothesis_config,
              look_index = look_index,
              analysis_time = analysis_time)
          }))

      # Initialize decisions for this analysis
      trial_decisions <- data.frame(
        hypothesis = hypothesis_ids,
        local_alpha = NA_real_,
        adjusted_p = NA_real_,
        efficacy = FALSE,
        futility = FALSE,
        rejected = FALSE,
        row.names = hypothesis_ids)

      # Apply GSD efficacy boundaries and recycle rejected alpha
      repeat {
        current_local_alphas <-
          familywise_alpha * mtp_state$hypotheses[hypothesis_ids]
        efficacy_reached <- setNames(
          object = rep(
            x = FALSE,
            times = length(hypothesis_ids)),
          nm = hypothesis_ids)

        for (hypothesis_id in hypothesis_ids) {
          trial_decisions[hypothesis_id, "local_alpha"] <-
            current_local_alphas[[hypothesis_id]]
          if (hypothesis_stopped[[hypothesis_id]]) next

          boundary_rows <- gsd_bounds[[hypothesis_id]]
          boundary_row <- boundary_rows[
            abs(boundary_rows$local_alpha -
                  current_local_alphas[[hypothesis_id]]) < 1e-12 &
              boundary_rows$analysis_index == look_index, ]
          if (nrow(boundary_row) == 0) next

          p_value <- test_results[hypothesis_id, "p_value"]
          efficacy_reached[[hypothesis_id]] <-
            !is.na(p_value) && p_value <= boundary_row$nominal_alpha[[1]]
        }

        if (!any(efficacy_reached)) break
        trial_decisions[efficacy_reached, "efficacy"] <- TRUE
        trial_decisions[efficacy_reached, "rejected"] <- TRUE
        hypothesis_stopped[efficacy_reached] <- TRUE
        mtp_state <- graphicalMCP::graph_update(
          graph = mtp_state,
          delete = efficacy_reached)$updated_graph
      }

      # Apply futility rules to hypotheses not stopped for efficacy
      for (hypothesis_id in hypothesis_ids) {
        hypothesis_config <- mtp1$nodes[[hypothesis_id]]
        gsd_config <- gsd1[[hypothesis_id]]
        if (hypothesis_stopped[[hypothesis_id]] ||
            !isTRUE(gsd_config$futility$enabled)) next
        hypothesis_futility_looks <- planned_looks[
          as.integer(gsd_config$futility$active_at)]
        if (!look_index %in% hypothesis_futility_looks) next

        effect_estimate <- test_results[hypothesis_id, "estimate"]
        if (is.na(effect_estimate)) next
        futility_threshold <-
          gsd_config$futility$threshold_rule$threshold
        futility_reached <- switch(
          EXPR = gsd_config$futility$threshold_rule$direction,
          ">" = effect_estimate > futility_threshold,
          ">=" = effect_estimate >= futility_threshold,
          "<" = effect_estimate < futility_threshold,
          "<=" = effect_estimate <= futility_threshold,
          stop(
            "Unsupported threshold direction: ",
            gsd_config$futility$threshold_rule$direction))
        if (futility_reached) {
          trial_decisions[hypothesis_id, "futility"] <- TRUE
          hypothesis_stopped[[hypothesis_id]] <- TRUE
        }
      }

      # Save statistical evidence and trial decisions separately
      sim_store[[grid_index]]$test_results$n[
        simulation_index, look_index, ] <- test_results$n
      sim_store[[grid_index]]$test_results$events[
        simulation_index, look_index, ] <- test_results$events
      sim_store[[grid_index]]$test_results$estimate[
        simulation_index, look_index, ] <- test_results$estimate
      sim_store[[grid_index]]$test_results$statistic[
        simulation_index, look_index, ] <- test_results$statistic
      sim_store[[grid_index]]$test_results$p_value[
        simulation_index, look_index, ] <- test_results$p_value
      sim_store[[grid_index]]$trial_decisions$local_alpha[
        simulation_index, look_index, ] <- trial_decisions$local_alpha
      sim_store[[grid_index]]$trial_decisions$adjusted_p[
        simulation_index, look_index, ] <- trial_decisions$adjusted_p
      sim_store[[grid_index]]$trial_decisions$efficacy[
        simulation_index, look_index, ] <- trial_decisions$efficacy
      sim_store[[grid_index]]$trial_decisions$futility[
        simulation_index, look_index, ] <- trial_decisions$futility
      sim_store[[grid_index]]$trial_decisions$rejected[
        simulation_index, look_index, ] <- trial_decisions$rejected
    }
  }
}

# ---- Calculate operating characteristics -----------------------------------

# Index the requested empirical summaries by label
empirical_quantities <- trial1$sim$empirical_quantities
empirical_labels <- vapply(
  X = empirical_quantities,
  FUN = function(empirical_quantity) empirical_quantity$label,
  FUN.VALUE = character(1))
names(empirical_quantities) <- empirical_labels

# Allocate operating-characteristic storage for every parameter-grid value
oc <- vector(
  mode = "list",
  length = nrow(param_grid))

# Summarize every requested quantity from the saved common arrays
for (grid_index in seq_len(nrow(param_grid))) {
  test_results <- sim_store[[grid_index]]$test_results
  trial_decisions <- sim_store[[grid_index]]$trial_decisions
  hypothesis_index <- match(
    x = hypothesis_id,
    table = hypothesis_ids)

  participant_counts <- test_results$n[, , hypothesis_index]
  event_counts <- test_results$events[, , hypothesis_index]
  observed_effects <- test_results$estimate[, , hypothesis_index]
  efficacy_crossings <- trial_decisions$efficacy[, , hypothesis_index]
  futility_stops <- trial_decisions$futility[, , hypothesis_index]

  oc[[grid_index]] <- list()

  # mean-participant-count: mean enrolled participants at each analysis
  oc[[grid_index]]$mean_participant_count <- colMeans(
    x = participant_counts,
    na.rm = TRUE)

  # mean-event-count: mean observed events at each analysis
  oc[[grid_index]]$mean_event_count <- colMeans(
    x = event_counts,
    na.rm = TRUE)

  # mean-observed-effect: mean observed hazard ratio at each analysis
  oc[[grid_index]]$mean_observed_effect <- colMeans(
    x = observed_effects,
    na.rm = TRUE)

  # observed-effect-interval: requested observed hazard-ratio quantiles
  observed_effect_probabilities <- empirical_quantities[[
    "observed-effect-interval"]]$summary$probabilities
  oc[[grid_index]]$observed_effect_interval <- apply(
    X = observed_effects,
    MARGIN = 2,
    FUN = stats::quantile,
    probs = observed_effect_probabilities,
    na.rm = TRUE)

  if (futility_enabled) {
    # futility-stop-probability: probability of first stopping for futility
    oc[[grid_index]]$futility_probability <- colMeans(
      x = futility_stops)
  }

  # power: marginal and cumulative efficacy-crossing probabilities
  oc[[grid_index]]$marginal_power <- colMeans(
    x = efficacy_crossings)
  oc[[grid_index]]$cumulative_power <- cumsum(
    x = oc[[grid_index]]$marginal_power)
}

# Use the single specified design row for the common output tables
oc <- oc[[1]]

# ==== FORMATTING THE RESULTS =================================================

# ---- Display helpers --------------------------------------------------------

# Round numeric values and render missing values as dashes
format_number <- function(x, digits = 4) {
  ifelse(
    test = is.na(x),
    yes = "-",
    no = format(
      x = round(
        x = x,
        digits = digits),
      trim = TRUE,
      scientific = FALSE))
}

display_table <- function(x, digits = 4) {
  formatted_table <- x
  for (column_index in seq_along(formatted_table)) {
    if (is.numeric(formatted_table[[column_index]])) {
      formatted_table[[column_index]] <-
        format_number(
          x = formatted_table[[column_index]],
          digits = digits)
    }
    formatted_table[[column_index]][
      is.na(formatted_table[[column_index]])] <- "-"
  }
  formatted_table
}

# ---- Build result tables ----------------------------------------------------

# Label planned interim and final analyses
analysis_labels <- rep(
  x = NA_character_,
  times = n_looks)
planned_analysis_labels <- if (length(planned_looks) == 1) "FA" else
  c(
    paste0(
      "IA",
      seq_len(length(planned_looks) - 1)),
    "FA")
analysis_labels[planned_looks] <- planned_analysis_labels

information_fraction <- rep(
  x = NA_real_,
  times = n_looks)
information_fraction[efficacy_looks] <-
  planned_information_fraction

# Build the efficacy-boundary table
type1_error_table <- data.frame(
  Analysis = analysis_labels[efficacy_looks],
  Boundary = "efficacy",
  Z = initial_gsd_bounds$z,
  CV_HR = initial_gsd_bounds$critical_effect,
  `Nominal α` = initial_gsd_bounds$nominal_alpha,
  `Cumulative α` = initial_gsd_bounds$cumulative_alpha,
  check.names = FALSE)

# Retain the common alpha-transfer table for the one-node graph
alpha_transfer_table <- data.frame(
  From = NA_character_,
  To = NA_character_,
  Weight = NA_real_,
  On = NA_character_,
  Condition = NA_character_)

# Build the futility table or its common unavailable row
if (futility_enabled) {
  type2_error_table <- data.frame(
    Analysis = analysis_labels[futility_looks],
    Boundary = "futility",
    Z = NA_real_,
    CV_HR = rep(
      x = gsd_config$futility$threshold_rule$threshold,
      times = length(futility_looks)),
    `Nominal β` = NA_real_,
    `Cumulative β` = NA_real_,
    `Stop Probability` = oc$futility_probability[futility_looks],
    check.names = FALSE)
} else {
  type2_error_table <- data.frame(
    Analysis = NA_character_,
    Boundary = NA_character_,
    Z = NA_real_,
    CV_HR = NA_real_,
    `Nominal β` = NA_real_,
    `Cumulative β` = NA_real_,
    `Stop Probability` = NA_real_,
    check.names = FALSE)
}

# Add analysis timing, information, and power estimates
planned_analyses_table <- data.frame(
  Analysis = analysis_labels,
  `Calendar Month` = analysis_times,
  `Information Fraction` = information_fraction,
  `N Events` = oc$mean_event_count,
  `Median N Events` = NA_real_,
  check.names = FALSE)
planned_analyses_table[[sprintf(
  fmt = "Marginal Power [%s]",
  hypothesis_id)]] <-
  ifelse(
    test = seq_len(n_looks) %in% efficacy_looks,
    yes = oc$marginal_power,
    no = NA_real_)
planned_analyses_table[[sprintf(
  fmt = "Cumulative Power [%s]",
  hypothesis_id)]] <-
  ifelse(
    test = seq_len(n_looks) %in% planned_looks,
    yes = oc$cumulative_power,
    no = NA_real_)

# Build enrollment and timing summaries
timeline_table <- data.frame(
  Analysis = analysis_labels,
  `Calendar Month` = analysis_times,
  `N Subjects` = oc$mean_participant_count,
  AHR = NA_real_,
  check.names = FALSE)

# Build the requested observed-effect summaries
observed_effect_table <- data.frame(
  Analysis = analysis_labels,
  `Calendar Month` = analysis_times,
  `Effect Measure` =
    if (endpoint_config$effect_measure == "hazard-ratio") "HR" else
      endpoint_config$effect_measure,
  Mean = oc$mean_observed_effect,
  Median = NA_real_,
  `2.5%` = oc$observed_effect_interval[1, ],
  `97.5%` = oc$observed_effect_interval[2, ],
  check.names = FALSE)

# Calculate Monte Carlo standard errors for power
empirical_precision_table <- data.frame(
  Analysis = analysis_labels,
  `Calendar Month` = analysis_times,
  `Marginal Power MC SE` = sqrt(
    oc$marginal_power *
      (1 - oc$marginal_power) /
      number_of_simulations),
  `Cumulative Power MC SE` = sqrt(
    oc$cumulative_power *
      (1 - oc$cumulative_power) /
      number_of_simulations),
  check.names = FALSE)

# Collect every formatted table in one result object
results <- list(
  type1_error_table = type1_error_table,
  alpha_transfer_table = alpha_transfer_table,
  type2_error_table = type2_error_table,
  planned_analyses_table = planned_analyses_table,
  timeline_table = timeline_table,
  observed_effect_table = observed_effect_table,
  empirical_precision_table = empirical_precision_table)

# ---- Print results ----------------------------------------------------------

# Temporarily widen the console for the result tables
previous_output_width <- getOption("width")
options(
  width = max(
    previous_output_width,
    240))

# Print trial metadata
cat(
  sprintf(
    fmt = "\n================ %s ================\n",
    trial1$label))
cat(
  sprintf(
    fmt = "  allocation ratio (intervention:comparator) %s:1 | stratification: %s | simulation: yes\n",
    format(
      x = allocation_ratio,
      trim = TRUE),
    if (isTRUE(trial1$strata$enabled)) "yes" else "none"))
cat(
  sprintf(
    fmt = "  simulations: %s | seed: %s | enrollment cap: %s | dropout: %.1f%% by month %s\n",
    format(
      x = number_of_simulations,
      big.mark = ","),
    simulation_seed,
    trial1$enrollment$n_accrual_max,
    100 * trial1$dropout$probability,
    trial1$dropout$window))

# Print the hypothesis label and tables in the common output order
calculation_label <- paste(
  hypothesis_id,
  "simulation")
cat(
  sprintf(
    fmt = "\nstate: %s\n%s\n",
    calculation_label,
    strrep(
      x = "-",
      times = nchar(calculation_label) + 7)))
cat(
  sprintf(
    fmt = "\n  %s  (local alpha %s, %s-sided, %s vs %s; looks by %s; simulation approximation; futility: %s; calculation: power)\n",
    hypothesis_config$label,
    format(
      x = hypothesis_config$initial_alpha),
    hypothesis_config$test_sides,
    hypothesis_config$intervention,
    hypothesis_config$comparator,
    gsd_config$looks$trigger,
    if (futility_enabled) "nonbinding threshold rule" else "none"))

cat("    -- type I error control --\n")
print(
  x = display_table(
    x = type1_error_table,
    digits = 5),
  row.names = FALSE)

cat("       alpha-transfer rules\n")
print(
  x = display_table(
    x = alpha_transfer_table,
    digits = 5),
  row.names = FALSE)

cat("    -- type II error / futility control --\n")
print(
  x = display_table(
    x = type2_error_table,
    digits = 5),
  row.names = FALSE)

cat("    -- planned analyses --\n")
print(
  x = display_table(
    x = planned_analyses_table,
    digits = 4),
  row.names = FALSE)

cat("    -- timeline and enrollment --\n")
print(
  x = display_table(
    x = timeline_table,
    digits = 3),
  row.names = FALSE)

cat("    -- observed-effect distribution --\n")
print(
  x = display_table(
    x = observed_effect_table,
    digits = 4),
  row.names = FALSE)
cat("\n")
options(
  width = previous_output_width)
