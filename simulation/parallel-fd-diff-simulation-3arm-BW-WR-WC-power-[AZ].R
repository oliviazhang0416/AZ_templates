# ============================================================
# Code title: Three-arm non-oncology parallel MTP simulation
#
# mtp:        parallel graphical procedure supplied by the JSON input
# gsd:        none - one fixed final analysis
# coc:        no
# cop:        no
# simulation: yes
# futility:   no
# endpoints:  body weight, at least 5% and 10% weight-loss response, and waist
#   circumference
# calculation: power over an enrollment-per-arm grid
#
# Key inputs: arm allocation, strata, endpoint tests, arm-specific outcome
#   models, MTP nodes and alpha-transfer edges, enrollment grid, empirical
#   decision expressions, number of simulations, and random-number seed
#
# Key outputs: initial graph alpha and transfer rules, per-hypothesis rejection
#   probability, user-defined power rules, Monte Carlo precision,
#   planned analyses, timeline and enrollment, and observed-effect distribution
#
# Engine: base R, stats, graphicalMCP - graphicalMCP::graph_create(),
#   graphicalMCP::graph_test_shortcut()
# Input:  jsonlite - fromJSON() (trial1, mtp1, gsd1)
# ============================================================

# ==== READING THE DESIGN INPUTS =============================================

design_input <- jsonlite::fromJSON(
  txt = "/Users/oliviazhang/Desktop/AZ_templates/simulation/parallel-fd-diff-simulation-3arm-BW-WR-WC-power-[AZ].json",
  simplifyVector = TRUE,
  simplifyDataFrame = FALSE,
  simplifyMatrix = TRUE)
trial1 <- design_input$trial
mtp1   <- design_input$mtp
gsd1   <- design_input$gsd

# ==== RUNNING THE CALCULATION ================================================

# ==== trial: arms, strata, endpoint assumptions, and enrollment grid ==========

# Read identifiers and simulation inputs
arm_ids <- names(trial1$arms)
hypothesis_ids <- names(mtp1$nodes)
outcome_models <- trial1$sim$outcome_models
endpoint_ids <- names(outcome_models)
number_of_simulations <- as.integer(trial1$sim$n_sim)
simulation_seed <- as.integer(trial1$sim$seed)

# Convert the stratum definitions to a named proportion vector
stratum_ids <- if (isTRUE(trial1$strata$enabled))
  names(trial1$strata$definitions) else "All"
stratum_proportions <- if (isTRUE(trial1$strata$enabled)) {
  vapply(
    X = trial1$strata$definitions,
    FUN = function(stratum_config) stratum_config$proportion,
    FUN.VALUE = numeric(1))
} else {
  1
}
names(stratum_proportions) <- stratum_ids

# Expand the enrollment parameter grid
grid_parameters <- vapply(
  X = trial1$sim$param_grid,
  FUN = function(grid_config) grid_config$parameter,
  FUN.VALUE = character(1))
enrollment_grid_index <- which(grid_parameters == "enrollment-per-arm")[[1]]
param_grid <- data.frame(
  enrollment_per_arm = as.integer(
    trial1$sim$param_grid[[enrollment_grid_index]]$values))

# ==== mtp: graph construction and alpha transfer =============================

# Read the initial local alpha for every hypothesis
familywise_alpha <- trial1$alpha$value
initial_local_alphas <- vapply(
  X = mtp1$nodes,
  FUN = function(node_config) node_config$initial_alpha,
  FUN.VALUE = numeric(1))
names(initial_local_alphas) <- hypothesis_ids

# Convert the JSON edges to an alpha-transfer matrix
alpha_transfer_matrix <- matrix(
  data = 0,
  nrow = length(hypothesis_ids),
  ncol = length(hypothesis_ids),
  dimnames = list(
    hypothesis_ids,
    hypothesis_ids))
if (length(mtp1$edges) > 0) {
  for (edge in mtp1$edges) {
    if (!identical(
      x = edge$on,
      y = "reject")) next
    alpha_transfer_matrix[edge$from, edge$to] <- edge$weight
  }
}

# Build the graph applied to every simulated trial
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

# ==== gsd: one fixed final analysis ==========================================

# Represent every fixed local test on the same boundary scale used by a GSD
analysis_times <- 1
gsd_bounds <- setNames(
  object = lapply(
    X = hypothesis_ids,
    FUN = function(hypothesis_id) {
      data.frame(
        local_alpha = alpha_states[[hypothesis_id]],
        look = 1L,
        nominal_alpha = alpha_states[[hypothesis_id]])
    }),
  nm = hypothesis_ids)

# ==== simulation =============================================================

# ---- Data generation --------------------------------------------------------
# simulate_BW_WR_WC: generate body-weight, responder, and waist outcomes
simulate_BW_WR_WC <- function(enrollment_per_arm) {
  # Initialize endpoint-by-arm outcome storage
  sim_data <- list()

  # Allocate each arm total across the configured strata
  target_subjects_by_stratum <-
    enrollment_per_arm * stratum_proportions / sum(stratum_proportions)
  subjects_by_stratum <- floor(target_subjects_by_stratum)
  unassigned_subjects <- enrollment_per_arm - sum(subjects_by_stratum)
  if (unassigned_subjects > 0) {
    strata_receiving_extra_subject <- order(
      x = target_subjects_by_stratum - subjects_by_stratum,
      decreasing = TRUE)[seq_len(unassigned_subjects)]
    subjects_by_stratum[strata_receiving_extra_subject] <-
      subjects_by_stratum[strata_receiving_extra_subject] + 1L
  }
  subjects_by_stratum <- as.integer(subjects_by_stratum)

  # Simulate directly modeled continuous endpoints
  for (endpoint_id in endpoint_ids) {
    # Read the endpoint's distributions and intervention effects
    outcome_model <- outcome_models[[endpoint_id]]
    outcome_distributions <- outcome_model$distributions
    outcome_effects <- outcome_model$effects
    first_arm_distribution <- outcome_distributions[[arm_ids[[1]]]]
    distribution_family <- first_arm_distribution$family
    if (is.null(distribution_family)) {
      distribution_family <-
        first_arm_distribution[[stratum_ids[[1]]]]$family
    }
    if (!identical(
      x = distribution_family,
      y = "normal")) next

    simulated_endpoint_outcomes <- setNames(
      object = lapply(
        X = arm_ids,
        FUN = function(arm_id) {
          # Resolve this arm's distribution, using its effect when needed
          arm_distribution <- outcome_distributions[[arm_id]]
          arm_effect <- outcome_effects[[arm_id]]
          uses_common_distribution <- if (!is.null(arm_distribution)) {
            !is.null(arm_distribution$family)
          } else if (!is.null(arm_effect$vs)) {
            !is.null(
              outcome_distributions[[arm_effect$vs]]$family)
          } else {
            FALSE
          }

          simulated_arm_outcomes <- if (uses_common_distribution) {
            # Simulate from one distribution shared across strata
            simulation_distribution <- if (!is.null(arm_distribution)) {
              arm_distribution
            } else {
              comparator_distribution <-
                outcome_distributions[[arm_effect$vs]]
              comparator_distribution$mean <-
                comparator_distribution$mean + arm_effect$value
              comparator_distribution
            }
            stats::rnorm(
              n = enrollment_per_arm,
              mean = simulation_distribution$mean,
              sd = simulation_distribution$sd)
          } else {
            # Simulate each stratum separately and combine its participants
            simulated_strata <- lapply(
              X = seq_along(stratum_ids),
              FUN = function(stratum_index) {
                stratum_id <- stratum_ids[[stratum_index]]
                simulation_distribution <- if (!is.null(arm_distribution)) {
                  arm_distribution[[stratum_id]]
                } else {
                  intervention_effect <- if (!is.null(arm_effect$vs)) {
                    arm_effect
                  } else {
                    arm_effect[[stratum_id]]
                  }
                  comparator_distribution <-
                    outcome_distributions[[intervention_effect$vs]]
                  if (is.null(comparator_distribution$family)) {
                    comparator_distribution <-
                      comparator_distribution[[stratum_id]]
                  }
                  comparator_distribution$mean <-
                    comparator_distribution$mean + intervention_effect$value
                  comparator_distribution
                }
                stats::rnorm(
                  n = subjects_by_stratum[[stratum_index]],
                  mean = simulation_distribution$mean,
                  sd = simulation_distribution$sd)
              })
            unlist(
              x = simulated_strata,
              use.names = FALSE)
          }
          simulated_arm_outcomes
        }),
      nm = arm_ids)
    sim_data[[endpoint_id]] <- simulated_endpoint_outcomes
  }

  # Derive binary outcomes from the simulated continuous outcomes
  for (endpoint_id in endpoint_ids) {
    if (!is.null(sim_data[[endpoint_id]])) next
    outcome_model <- outcome_models[[endpoint_id]]
    outcome_distributions <- outcome_model$distributions
    simulated_endpoint_outcomes <- setNames(
      object = lapply(
        X = arm_ids,
        FUN = function(arm_id) {
          # Apply the configured threshold to the source endpoint
          arm_distribution <- outcome_distributions[[arm_id]]
          source_outcome <-
            sim_data[[arm_distribution$derived_from]][[arm_id]]
          threshold <- arm_distribution$threshold_rule$value
          simulated_arm_outcomes <- switch(
            EXPR = arm_distribution$threshold_rule$operator,
            "<" = source_outcome < threshold,
            "<=" = source_outcome <= threshold,
            ">" = source_outcome > threshold,
            ">=" = source_outcome >= threshold,
            stop(
              "Unsupported derived-outcome threshold operator: ",
              arm_distribution$threshold_rule$operator))
          simulated_arm_outcomes
        }),
      nm = arm_ids)
    sim_data[[endpoint_id]] <- simulated_endpoint_outcomes
  }

  sim_data
}

# ---- Endpoint test ----------------------------------------------------------

# test_BW_WR_WC: calculate one local test for every MTP hypothesis
test_BW_WR_WC <- function(sim_data) {
  # Initialize the common local-test fields
  test_results <- data.frame(
    hypothesis = hypothesis_ids,
    look = 1L,
    n = NA_integer_,
    events = NA_integer_,
    estimate = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    row.names = hypothesis_ids)

  # Compare the intervention and comparator for every hypothesis
  for (hypothesis_id in hypothesis_ids) {
    # Select the endpoint and arm-level outcomes for this hypothesis
    hypothesis_config <- mtp1$nodes[[hypothesis_id]]
    endpoint_config <- trial1$endpoints[[hypothesis_config$endpoint]]
    intervention_outcome <-
      sim_data[[hypothesis_config$endpoint]][[
        hypothesis_config$intervention]]
    comparator_outcome <-
      sim_data[[hypothesis_config$endpoint]][[
        hypothesis_config$comparator]]
    intervention_sample_size <- length(intervention_outcome)
    comparator_sample_size <- length(comparator_outcome)
    test_results[hypothesis_id, "n"] <-
      intervention_sample_size + comparator_sample_size

    # Apply the configured continuous or binary endpoint test
    if (endpoint_config$type == "continuous") {
      # Calculate means, variances, and the configured t statistic
      intervention_mean <- mean(intervention_outcome)
      comparator_mean <- mean(comparator_outcome)
      mean_difference <- intervention_mean - comparator_mean
      intervention_variance <- stats::var(intervention_outcome)
      comparator_variance <- stats::var(comparator_outcome)

      if (isTRUE(endpoint_config$test$equal_variance)) {
        degrees_freedom <-
          intervention_sample_size + comparator_sample_size - 2
        pooled_variance <- (
          (intervention_sample_size - 1) * intervention_variance +
            (comparator_sample_size - 1) * comparator_variance
        ) / degrees_freedom
        standard_error <- sqrt(pooled_variance * (
          1 / intervention_sample_size + 1 / comparator_sample_size))
      } else {
        standard_error_squared <-
          intervention_variance / intervention_sample_size +
          comparator_variance / comparator_sample_size
        degrees_freedom <- standard_error_squared^2 / (
          (intervention_variance / intervention_sample_size)^2 /
            (intervention_sample_size - 1) +
          (comparator_variance / comparator_sample_size)^2 /
            (comparator_sample_size - 1))
        standard_error <- sqrt(standard_error_squared)
      }

      # Convert the t statistic to a two-sided p-value
      test_statistic <- mean_difference / standard_error
      p_value <-
        2 * stats::pt(
          q = -abs(test_statistic),
          df = degrees_freedom)
      effect_estimate <- mean_difference
    } else {
      # Calculate response counts, proportions, and the chi-square statistic
      intervention_events <- sum(intervention_outcome)
      comparator_events <- sum(comparator_outcome)
      intervention_probability <-
        intervention_events / intervention_sample_size
      comparator_probability <- comparator_events / comparator_sample_size
      pooled_probability <-
        (intervention_events + comparator_events) /
        (intervention_sample_size + comparator_sample_size)
      standard_error <- sqrt(
        pooled_probability * (1 - pooled_probability) *
          (1 / intervention_sample_size + 1 / comparator_sample_size))

      # Convert the chi-square statistic to a two-sided p-value
      test_statistic <-
        (intervention_probability - comparator_probability) / standard_error
      p_value <- stats::pchisq(
        q = test_statistic^2,
        df = 1,
        lower.tail = FALSE)
      if (!is.finite(p_value)) p_value <- 1
      effect_estimate <- intervention_probability - comparator_probability
    }

    # Save the common local-test results
    test_results[hypothesis_id, "estimate"] <- effect_estimate
    test_results[hypothesis_id, "statistic"] <- test_statistic
    test_results[hypothesis_id, "p_value"] <- p_value
  }
  test_results
}

# ---- Initialize simulation results ------------------------------------------

# Create common simulation-by-look-by-hypothesis storage templates
result_dimensions <- c(
  number_of_simulations,
  length(analysis_times),
  length(hypothesis_ids))
result_dimnames <- list(
  simulation = NULL,
  look = "FA",
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
  enrollment_per_arm <- param_grid$enrollment_per_arm[[grid_index]]

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

  # Use common random numbers across parameter-grid values
  set.seed(
    seed = simulation_seed)

  for (simulation_index in seq_len(number_of_simulations)) {
    # Generate all correlated endpoint data for one trial
    sim_data <- simulate_BW_WR_WC(
      enrollment_per_arm = enrollment_per_arm)

    # Calculate the unadjusted local tests
    test_results <- test_BW_WR_WC(
      sim_data = sim_data)

    # Apply the fixed local tests and graphical alpha transfers
    mtp_state <- graphicalMCP::graph_test_shortcut(
      graph = mtp_graph,
      p = setNames(
        object = test_results$p_value,
        nm = hypothesis_ids),
      alpha = familywise_alpha,
      test_values = TRUE)
    test_trace <- mtp_state$test_values$results
    local_alpha <- initial_local_alphas
    if (nrow(test_trace) > 0) {
      for (trace_index in seq_len(nrow(test_trace))) {
        local_alpha[test_trace$Hypothesis[[trace_index]]] <-
          test_trace$Alpha[[trace_index]] * test_trace$Weight[[trace_index]]
      }
    }
    trial_decisions <- data.frame(
      hypothesis = hypothesis_ids,
      local_alpha = unname(local_alpha),
      adjusted_p = unname(mtp_state$outputs$adjusted_p[hypothesis_ids]),
      efficacy = unname(mtp_state$outputs$rejected[hypothesis_ids]),
      futility = FALSE,
      rejected = unname(mtp_state$outputs$rejected[hypothesis_ids]),
      row.names = hypothesis_ids)

    # Save statistical evidence and trial decisions separately
    sim_store[[grid_index]]$test_results$n[simulation_index, 1, ] <-
      test_results$n
    sim_store[[grid_index]]$test_results$events[simulation_index, 1, ] <-
      test_results$events
    sim_store[[grid_index]]$test_results$estimate[simulation_index, 1, ] <-
      test_results$estimate
    sim_store[[grid_index]]$test_results$statistic[simulation_index, 1, ] <-
      test_results$statistic
    sim_store[[grid_index]]$test_results$p_value[simulation_index, 1, ] <-
      test_results$p_value
    sim_store[[grid_index]]$trial_decisions$local_alpha[
      simulation_index, 1, ] <-
      trial_decisions$local_alpha
    sim_store[[grid_index]]$trial_decisions$adjusted_p[
      simulation_index, 1, ] <-
      trial_decisions$adjusted_p
    sim_store[[grid_index]]$trial_decisions$efficacy[
      simulation_index, 1, ] <-
      trial_decisions$efficacy
    sim_store[[grid_index]]$trial_decisions$futility[
      simulation_index, 1, ] <-
      trial_decisions$futility
    sim_store[[grid_index]]$trial_decisions$rejected[
      simulation_index, 1, ] <-
      trial_decisions$rejected
  }
}

# ---- Calculate operating characteristics -----------------------------------

# Evaluate nested any-of and all-of rejection rules
evaluate_expression <- function(expression, rejections) {
  if (!is.null(expression$hypothesis_id)) {
    return(rejections[, expression$hypothesis_id])
  }
  if (!is.null(expression$any_of)) {
    values <- lapply(
      X = expression$any_of,
      FUN = evaluate_expression,
      rejections = rejections)
    return(Reduce(
      f = `|`,
      x = values))
  }
  if (!is.null(expression$all_of)) {
    values <- lapply(
      X = expression$all_of,
      FUN = evaluate_expression,
      rejections = rejections)
    return(Reduce(
      f = `&`,
      x = values))
  }
  stop("Unsupported empirical expression")
}

# Index the requested empirical summaries by label
empirical_quantities <- trial1$sim$empirical_quantities
empirical_labels <- vapply(
  X = empirical_quantities,
  FUN = function(empirical_quantity) empirical_quantity$label,
  FUN.VALUE = character(1))
names(empirical_quantities) <- empirical_labels

# Initialize operating-characteristic storage
oc <- list()

# Calculate one requested empirical probability
calculate_empirical_probability <- function(label, rejections) {
  mean(
    x = evaluate_expression(
      expression = empirical_quantities[[label]]$summary$expression,
      rejections = rejections))
}

oc$empirical_estimates <- matrix(
  data = NA_real_,
  nrow = nrow(param_grid),
  ncol = length(empirical_quantities),
  dimnames = list(
    NULL,
    empirical_labels))
oc$hypothesis_estimates <- matrix(
  data = NA_real_,
  nrow = nrow(param_grid),
  ncol = length(hypothesis_ids),
  dimnames = list(
    NULL,
    hypothesis_ids))

# Summarize every parameter-grid value
for (grid_index in seq_len(nrow(param_grid))) {
  rejections <- sim_store[[grid_index]]$trial_decisions$rejected[, 1, ]
  oc$hypothesis_estimates[grid_index, ] <- colMeans(
    x = rejections)

  # power-high-body-weight: reject the high-dose body-weight hypothesis
  oc$empirical_estimates[grid_index, "power-high-body-weight"] <-
    calculate_empirical_probability(
      label = "power-high-body-weight",
      rejections = rejections)

  # power-high-5pct: reject the high-dose 5% responder hypothesis
  oc$empirical_estimates[grid_index, "power-high-5pct"] <-
    calculate_empirical_probability(
      label = "power-high-5pct",
      rejections = rejections)

  # power-high-any-primary: reject either high-dose primary hypothesis
  oc$empirical_estimates[grid_index, "power-high-any-primary"] <-
    calculate_empirical_probability(
      label = "power-high-any-primary",
      rejections = rejections)

  # power-high-10pct: reject the high-dose 10% responder hypothesis
  oc$empirical_estimates[grid_index, "power-high-10pct"] <-
    calculate_empirical_probability(
      label = "power-high-10pct",
      rejections = rejections)

  # power-high-waist: reject the high-dose waist hypothesis
  oc$empirical_estimates[grid_index, "power-high-waist"] <-
    calculate_empirical_probability(
      label = "power-high-waist",
      rejections = rejections)

  # power-low-body-weight: reject the low-dose body-weight hypothesis
  oc$empirical_estimates[grid_index, "power-low-body-weight"] <-
    calculate_empirical_probability(
      label = "power-low-body-weight",
      rejections = rejections)

  # power-low-5pct: reject the low-dose 5% responder hypothesis
  oc$empirical_estimates[grid_index, "power-low-5pct"] <-
    calculate_empirical_probability(
      label = "power-low-5pct",
      rejections = rejections)

  # power-low-any-primary: reject either low-dose primary hypothesis
  oc$empirical_estimates[grid_index, "power-low-any-primary"] <-
    calculate_empirical_probability(
      label = "power-low-any-primary",
      rejections = rejections)

  # power-low-10pct: reject the low-dose 10% responder hypothesis
  oc$empirical_estimates[grid_index, "power-low-10pct"] <-
    calculate_empirical_probability(
      label = "power-low-10pct",
      rejections = rejections)

  # power-low-waist: reject the low-dose waist hypothesis
  oc$empirical_estimates[grid_index, "power-low-waist"] <-
    calculate_empirical_probability(
      label = "power-low-waist",
      rejections = rejections)

  # overall-high-dose-power: satisfy the complete high-dose success rule
  oc$empirical_estimates[grid_index, "overall-high-dose-power"] <-
    calculate_empirical_probability(
      label = "overall-high-dose-power",
      rejections = rejections)

  # overall-low-dose-power: satisfy the complete low-dose success rule
  oc$empirical_estimates[grid_index, "overall-low-dose-power"] <-
    calculate_empirical_probability(
      label = "overall-low-dose-power",
      rejections = rejections)

  # overall-both-doses-power: satisfy both dose-level success rules
  oc$empirical_estimates[grid_index, "overall-both-doses-power"] <-
    calculate_empirical_probability(
      label = "overall-both-doses-power",
      rejections = rejections)
}

# No observed-effect distribution summaries are requested for this design.

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

# Build the initial graphical-MTP table
initial_available_hypotheses <- paste(
  hypothesis_ids,
  collapse = ", ")
type1_error_table <- do.call(
  what = rbind,
  args = lapply(
    X = hypothesis_ids,
    FUN = function(hypothesis_id) {
      hypothesis_config <- mtp1$nodes[[hypothesis_id]]
      endpoint_config <- trial1$endpoints[[hypothesis_config$endpoint]]
      local_alpha <- initial_local_alphas[[hypothesis_id]]
      data.frame(
        State = "Initial",
        Analysis = "FA",
        Rejected = "none",
        Available = initial_available_hypotheses,
        Hypothesis = hypothesis_id,
        `Graph α` = local_alpha,
        `Effective α` = local_alpha,
        `Local Test` = endpoint_config$test$method,
        Status = if (local_alpha > 1e-14) "initially testable" else "gated",
        Boundary = "efficacy",
        Z = NA_real_,
        CV_HR = NA_real_,
        `Nominal α` = local_alpha,
        `Cumulative α` = NA_real_,
        check.names = FALSE)
  }))

# Build the alpha-transfer table
alpha_transfer_table <- if (length(mtp1$edges) == 0) {
  data.frame(
    From = NA_character_,
    To = NA_character_,
    Weight = NA_real_,
    On = NA_character_,
    Condition = NA_character_)
} else {
  do.call(
    what = rbind,
    args = lapply(
      X = mtp1$edges,
      FUN = function(edge) data.frame(
        From = edge$from,
        To = edge$to,
        Weight = edge$weight,
        On = edge$on,
        Condition = if (is.null(edge$condition)) NA_character_ else
          paste(
            unlist(
              x = edge$condition),
            collapse = " "),
        check.names = FALSE)))
}

# Build the unavailable futility table with the common columns
type2_error_table <- data.frame(
  Analysis = NA_character_,
  Boundary = NA_character_,
  Z = NA_real_,
  CV_HR = NA_real_,
  `Nominal β` = NA_real_,
  `Cumulative β` = NA_real_,
  `Stop Probability` = NA_real_,
  check.names = FALSE)

# Add marginal, cumulative, and requested power estimates
planned_analyses_table <- data.frame(
  Analysis = rep(
    x = "FA",
    times = nrow(param_grid)),
  `Enrollment per Arm` = param_grid$enrollment_per_arm,
  `Information Fraction` = 1,
  `N Events` = NA_real_,
  check.names = FALSE)
for (hypothesis_id in hypothesis_ids) {
  planned_analyses_table[[sprintf(
    fmt = "Marginal Power [%s]",
    hypothesis_id)]] <-
    oc$hypothesis_estimates[, hypothesis_id]
  planned_analyses_table[[sprintf(
    fmt = "Cumulative Power [%s]",
    hypothesis_id)]] <-
    oc$hypothesis_estimates[, hypothesis_id]
}
for (quantity_index in seq_along(empirical_quantities)) {
  expression <- empirical_quantities[[quantity_index]]$summary$expression
  if (is.null(expression$hypothesis_id)) {
    label <- empirical_quantities[[quantity_index]]$label
    planned_analyses_table[[sprintf(
      fmt = "Power [%s]",
      label)]] <-
      oc$empirical_estimates[, quantity_index]
  }
}

# Build enrollment and timing summaries
timeline_table <- data.frame(
  Analysis = rep(
    x = "FA",
    times = nrow(param_grid)),
  `Calendar Month` = NA_real_,
  `N Subjects` = length(arm_ids) * param_grid$enrollment_per_arm,
  AHR = NA_real_,
  check.names = FALSE)

# Retain the common observed-effect table when no summary is requested
observed_effect_table <- data.frame(
  Analysis = NA_character_,
  `Calendar Month` = NA_real_,
  `Effect Measure` = NA_character_,
  Mean = NA_real_,
  Median = NA_real_,
  `2.5%` = NA_real_,
  `97.5%` = NA_real_,
  check.names = FALSE)

# Calculate Monte Carlo standard errors for every empirical estimate
precision_rows <- lapply(
  X = seq_len(nrow(param_grid)),
  FUN = function(grid_index) {
    estimate <- oc$empirical_estimates[grid_index, ]
    data.frame(
      `Enrollment per Arm` = param_grid$enrollment_per_arm[[grid_index]],
      Label = names(estimate),
      Estimate = as.numeric(estimate),
      `Monte Carlo SE` = sqrt(
        estimate * (1 - estimate) / number_of_simulations),
      check.names = FALSE)
  })
empirical_precision_table <- do.call(
  what = rbind,
  args = precision_rows)

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
    1200))

# Print trial metadata
cat(
  sprintf(
    fmt = "\n================ %s ================\n",
    trial1$label))
allocation_text <- paste(
  sprintf(
    fmt = "%s=%s",
    arm_ids,
    vapply(
      X = trial1$arms,
      FUN = function(arm_config) arm_config$randomization_weight,
      FUN.VALUE = numeric(1))),
  collapse = ":")
cat(
  sprintf(
    fmt = "  allocation weights %s | stratification: %s | simulation: yes\n",
    allocation_text,
    if (isTRUE(trial1$strata$enabled))
      paste(
        stratum_ids,
        collapse = "/") else "none"))
cat(
  sprintf(
    fmt = "  simulations per enrollment value: %s | seed: %s | calculation: power\n",
    format(
      x = number_of_simulations,
      big.mark = ","),
    simulation_seed))

# Print the MTP label and tables in the common output order
calculation_label <- paste(
  mtp1$procedure$strategy,
  "MTP simulation")
cat(
  sprintf(
    fmt = "\nstate: %s\n%s\n",
    calculation_label,
    strrep(
      x = "-",
      times = nchar(calculation_label) + 7)))

cat("\n    -- type I error control --\n")
cat("       initial graph\n")
print(
  x = display_table(
    x = type1_error_table,
    digits = 5),
  row.names = FALSE)
cat("\n       alpha-transfer rules\n")
print(
  x = display_table(
    x = alpha_transfer_table,
    digits = 5),
  row.names = FALSE)

cat("\n    -- type II error / futility control --\n")
print(
  x = display_table(
    x = type2_error_table,
    digits = 5),
  row.names = FALSE)

cat("\n    -- planned analyses --\n")
print(
  x = display_table(
    x = planned_analyses_table,
    digits = 4),
  row.names = FALSE)

cat("\n    -- timeline and enrollment --\n")
print(
  x = display_table(
    x = timeline_table,
    digits = 3),
  row.names = FALSE)

cat("\n    -- observed-effect distribution --\n")
print(
  x = display_table(
    x = observed_effect_table,
    digits = 4),
  row.names = FALSE)
cat("\n")
options(
  width = previous_output_width)
