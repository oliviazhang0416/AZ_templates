# ============================================================
# Code title: Three-arm OS event-size calculation under non-proportional hazards with a CoC comparison
#
# mtp:        separate Arm A vs C, Arm B vs C, and loose-alpha Arm A vs B calculations
# gsd:        H1/H2 use event-driven sfLDOF boundaries; H3 uses fixed nominal alpha at calendar times
# coc:        yes - Arm B is a component of Arm A, and H3 compares Arm A with Arm B
# cop:        no
# simulation: no
# futility:   no
# endpoints:  OS with delayed effects
#
# Key inputs: target power by hypothesis, 1:1:1 randomization, control median
#   OS 19.2 months, 6-month effect delay, HR 0.79 for B vs C, HR 0.81 for
#   A vs B, 1% annual dropout, reference enrollment shape, and fixed looks
#
# Key outputs: required final events for every hypothesis, expected timing/AHR,
#   efficacy boundaries, and achieved power
#
# Engine: gsDesign2 - gsDesign2::gs_power_ahr(), gsDesign2::gs_design_ahr()
# Input:  jsonlite - fromJSON() (trial1, mtp1, gsd1)
# ============================================================

# ==== READING THE DESIGN INPUTS =============================================

design_input <- jsonlite::fromJSON(
  "/Users/oliviazhang/Desktop/AZ_templates/nomtp-gsd-ahr-gsdesign2-3arm-OS-CoC-samplesize-[AZ].json",
  simplifyVector = TRUE,
  simplifyDataFrame = FALSE,
  simplifyMatrix = TRUE)
trial1 <- design_input$trial
mtp1   <- design_input$mtp
gsd1   <- design_input$gsd
rm(design_input)

# ==== RUNNING THE CALCULATION ================================================

# ==== trial: pairwise enrollment and linked OS distributions =================

allocation_ratio_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  trial1$arms[[nd$intervention]]$randomization_weight /
    trial1$arms[[nd$comparator]]$randomization_weight
}

reference_enrollment_total <- tail(trial1$enrollment$cumulative, 1)

pairwise_enrollment_for <- function(id, n_total = reference_enrollment_total) {
  nd <- mtp1$nodes[[id]]
  arm_weight <- vapply(trial1$arms, function(x) x$randomization_weight, numeric(1))
  pair_fraction <- sum(arm_weight[c(nd$intervention, nd$comparator)]) / sum(arm_weight)
  enrollment_shape <- trial1$enrollment$cumulative / reference_enrollment_total
  cumulative <- round(enrollment_shape * n_total * pair_fraction)

  tibble::tibble(
    stratum = "All",
    duration = rep(1, length(cumulative)),
    rate = cumulative - c(0, head(cumulative, -1)))
}

step_value <- function(cut, value, at) {
  value[findInterval(at, cut)]
}

hazard_curve_for <- function(endpoint, arm) {
  distribution <- endpoint$distributions[[arm]]
  if (!is.null(distribution)) {
    if (distribution$family == "exponential") {
      rate <- if (is.null(distribution$rate)) log(2) / distribution$median else
        distribution$rate
      return(list(cut = 0, rate = rate))
    }
    return(list(cut = distribution$cut$at, rate = distribution$rate))
  }

  effect <- endpoint$effects[[arm]]
  reference <- hazard_curve_for(endpoint, effect$vs)
  effect_cut <- if (is.null(effect$cut$at)) 0 else effect$cut$at
  effect_value <- if (is.null(effect$cut$at)) rep(effect$value, 1) else effect$value
  cut <- sort(unique(c(reference$cut, effect_cut)))

  list(
    cut = cut,
    rate = step_value(reference$cut, reference$rate, cut) *
      step_value(effect_cut, effect_value, cut))
}

failure_table_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  ep <- trial1$endpoints[[nd$endpoint]]
  reference <- hazard_curve_for(ep, nd$comparator)
  effect <- ep$effects[[nd$intervention]]
  effect_cut <- if (is.null(effect$cut$at)) 0 else effect$cut$at
  effect_value <- if (is.null(effect$cut$at)) rep(effect$value, 1) else effect$value
  cut <- sort(unique(c(reference$cut, effect_cut)))

  tibble::tibble(
    stratum = "All",
    duration = c(diff(cut), 100),
    fail_rate = step_value(reference$cut, reference$rate, cut),
    hr = step_value(effect_cut, effect_value, cut),
    dropout_rate = -log(1 - trial1$dropout$probability) / trial1$dropout$window)
}

# ==== mtp: hypothesis-specific standard and loose alpha ======================

ids <- names(mtp1$nodes)
labels <- vapply(ids, function(id) {
  nd <- mtp1$nodes[[id]]
  ep <- trial1$endpoints[[nd$endpoint]]
  sprintf("%s (%s [%s, %s], %s)",
          id, nd$endpoint, ep$type, ep$effect_measure, nd$label)
}, character(1))

alpha_1s_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  nd$initial_alpha / nd$test_sides
}

# ==== gsd: hypothesis-specific event-size calculations =======================

run_power <- function(id, max_events = NULL,
                      n_total = reference_enrollment_total) {
  g      <- gsd1[[id]]
  family <- g$efficacy$family
  schedule <- g$looks$schedule
  event <- if (schedule$scale == "fraction-of-final-events")
    schedule$values * max_events else NULL
  time <- if (schedule$scale == "calendar-time") schedule$values else NULL

  if (g$efficacy$type == "customized") {
    upper <- gsDesign2::gs_b
    upar <- qnorm(1 - family$nominal_values / mtp1$nodes[[id]]$test_sides)
  } else {
    upper <- gsDesign2::gs_spending_bound
    upar <- list(
      sf = getExportedValue("gsDesign", family$fn),
      total_spend = alpha_1s_for(id),
      param = family$param)
  }

  gsDesign2::gs_power_ahr(
    enroll_rate = pairwise_enrollment_for(id, n_total),
    fail_rate = failure_table_for(id),
    ratio = allocation_ratio_for(id),
    analysis_time = time,
    event = event,
    upper = upper,
    upar = upar,
    test_lower = FALSE,
    interval = c(0.01, 100))
}

final_power <- function(x) {
  upper <- x$bound[x$bound$bound == "upper", , drop = FALSE]
  tail(upper$probability, 1)
}

bracket_root <- function(fn, target) {
  lower <- 100
  upper <- 600
  while (fn(upper) < target) {
    upper <- 2 * upper
    if (upper > 1e7) stop("Could not bracket the required event count")
  }
  uniroot(function(x) fn(x) - target,
          interval = c(lower, upper), tol = 0.01)$root
}

run_calendar_event_design <- function(id) {
  g <- gsd1[[id]]
  family <- g$efficacy$family

  if (g$efficacy$type == "customized") {
    upper <- gsDesign2::gs_b
    upar <- qnorm(1 - family$nominal_values / mtp1$nodes[[id]]$test_sides)
  } else {
    upper <- gsDesign2::gs_spending_bound
    upar <- list(
      sf = getExportedValue("gsDesign", family$fn),
      total_spend = alpha_1s_for(id),
      param = family$param)
  }

  gsDesign2::gs_design_ahr(
    enroll_rate = pairwise_enrollment_for(id),
    fail_rate = failure_table_for(id),
    ratio = allocation_ratio_for(id),
    alpha = alpha_1s_for(id),
    beta = g$beta,
    analysis_time = g$looks$schedule$values,
    upper = upper,
    upar = upar,
    test_upper = seq_along(g$looks$schedule$values) %in% g$efficacy$active_at,
    test_lower = FALSE,
    interval = c(0.01, 100))
}

run_event_size <- function(id) {
  g <- gsd1[[id]]
  target <- 1 - g$beta

  if (g$looks$trigger == "event-count") {
    power_at <- function(max_events)
      final_power(run_power(id, max_events = max_events))
    required <- ceiling(bracket_root(power_at, target))
    return(list(
      required_events = required,
      result = run_power(id, max_events = required)))
  }

  design <- run_calendar_event_design(id)
  list(
    required_events = ceiling(tail(design$analysis$event, 1)),
    result = design)
}

results <- setNames(lapply(ids, run_event_size), ids)

# ==== FORMATTING THE RESULTS =================================================

format_optional_output <- function(x, digits = 4) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return("-")
  format(round(as.numeric(x)[1], digits), trim = TRUE, scientific = FALSE)
}

beta_futility_summary <- function(g, look_labels, calculation) {
  futility <- g$futility
  enabled <- isTRUE(futility$enabled)
  beta_spending <- enabled && identical(futility$type, "beta-spending")
  beta_purpose <- "-"
  if (!is.null(g$beta)) {
    if (identical(calculation, "sample-size")) {
      beta_purpose <- if (beta_spending)
        "sample-size target + futility spending" else "sample-size target"
    } else if (beta_spending) {
      beta_purpose <- "futility spending"
    }
  }
  target_power <- if (identical(calculation, "sample-size") &&
                      !is.null(g$beta)) 1 - g$beta else NULL
  active_looks <- if (enabled && length(futility$active_at) > 0)
    paste(look_labels[as.integer(futility$active_at)], collapse = ", ") else "-"
  data.frame(
    `β Input` = format_optional_output(g$beta),
    `β Purpose` = beta_purpose,
    `Target Power` = format_optional_output(target_power),
    Futility = if (enabled) "yes" else "no",
    Binding = if (!enabled || is.null(futility$binding)) "-" else
      if (isTRUE(futility$binding)) "yes" else "no",
    Type = if (enabled && !is.null(futility$type)) futility$type else "-",
    Family = if (enabled && !is.null(futility$family$fn))
      futility$family$fn else "-",
    Parameter = if (enabled)
      format_optional_output(futility$family$param) else "-",
    `Active Looks` = active_looks,
    check.names = FALSE)
}


cat(sprintf("\n================ %s ================\n", trial1$label))
ratio_text <- paste(sprintf("%s=%s:1", ids,
                            vapply(ids, allocation_ratio_for, numeric(1))),
                    collapse = ", ")
cat(sprintf(
  "  allocation ratios (intervention:comparator) %s | stratification: %s | simulation: %s\n",
  ratio_text, if (isTRUE(trial1$strata$enabled)) "yes" else "none",
  if (isTRUE(trial1$sim$enabled)) "yes" else "none (closed form)"))
cat(sprintf("  reference accrual-shape total: %s | dropout: %.1f%% by month %s\n",
            reference_enrollment_total,
            100 * trial1$dropout$probability, trial1$dropout$window))

for (id in ids) {
  nd <- mtp1$nodes[[id]]
  g  <- gsd1[[id]]
  sizing <- results[[id]]
  x  <- sizing$result
  a  <- x$analysis
  b  <- x$bound[x$bound$bound == "upper", , drop = FALSE]
  looks <- if (nrow(a) == 1) "FA" else c(paste0("IA", seq_len(nrow(a) - 1)), "FA")

  boundary <- data.frame(
    Analysis       = looks,
    Boundary       = "efficacy",
    Z              = round(b$z, 3),
    CV_HR          = round(b$`~hr at bound`, 3),
    `Nominal α`    = round(nd$test_sides * b$`nominal p`, 5),
    `Cumulative α` = round(nd$test_sides * b$probability0, 5),
    check.names = FALSE)

  analyses <- data.frame(
    Analysis               = looks,
    `Information Fraction` = round(a$info_frac, 3),
    `N Events`              = round(a$event, 1),
    `Marginal Power`        = round(pnorm(a$theta * sqrt(a$info) - b$z), 4),
    `Cumulative Power`      = round(b$probability, 4),
    check.names = FALSE)

  timeline <- data.frame(
    Analysis                = looks,
    `Calendar Month`        = round(a$time, 1),
    AHR                     = round(a$ahr, 4),
    check.names = FALSE)

  state <- paste(id, "calculation")
  cat(sprintf("\nstate: %s\n%s\n", state, strrep("-", nchar(state) + 7)))
  cat(sprintf("\n  %s  (local alpha %s, target power %.1f%%, required final events %s, %s-sided, %s vs %s; looks by %s; %s approximation; futility: none; calculation: sample size [events])\n",
              labels[id], format(nd$initial_alpha), 100 * (1 - g$beta),
              sizing$required_events, nd$test_sides,
              nd$intervention, nd$comparator, g$looks$trigger, g$approximation))
  cat("    -- type I error control --\n")
  print(boundary, row.names = FALSE)
  cat("    -- type II error / futility control --\n")
  print(beta_futility_summary(gsd1[[id]], looks, "sample-size"), row.names = FALSE)
  cat("    -- planned analyses --\n")
  print(analyses, row.names = FALSE)
  cat("    -- NPH timing and effect characterization --\n")
  print(timeline, row.names = FALSE)
}
cat("\n")
