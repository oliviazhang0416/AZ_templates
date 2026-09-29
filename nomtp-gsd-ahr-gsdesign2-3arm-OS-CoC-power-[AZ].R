# ============================================================
# Code title: Three-arm OS power under non-proportional hazards with a CoC comparison
#
# mtp:        separate Arm A vs C, Arm B vs C, and loose-alpha Arm A vs B calculations
# gsd:        H1/H2 use event-driven sfLDOF boundaries; H3 uses fixed nominal alpha at calendar times
# coc:        yes - Arm B is a component of Arm A, and H3 compares Arm A with Arm B
# cop:        no
# simulation: no
# futility:   no
# endpoints:  OS with delayed effects
#
# Key inputs: 1:1:1 randomization, control median OS 19.2 months, 6-month
#   effect delay, HR 0.79 for B vs C, HR 0.81 for A vs B, 1% annual
#   dropout, reference accrual shape, events 426/531, and a planned 503.6
#   expected events at the month-49.1 CoC analysis
#
# Key outputs: achieved power, implied enrollment for calendar-time power,
#   event trajectory, average HR, efficacy boundary, and nominal alpha
#
# Engine: gsDesign2 - gsDesign2::gs_power_ahr()
# Input:  jsonlite - fromJSON() (trial1, mtp1, gsd1)
# ============================================================

# ==== READING THE DESIGN INPUTS =============================================

design_input <- jsonlite::fromJSON(
  "/Users/oliviazhang/Desktop/AZ_templates/nomtp-gsd-ahr-gsdesign2-3arm-OS-CoC-power-[AZ].json",
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

pairwise_enrollment_for <- function(id, scale = 1) {
  nd <- mtp1$nodes[[id]]
  arm_weight <- vapply(trial1$arms, function(x) x$randomization_weight, numeric(1))
  pair_fraction <- sum(arm_weight[c(nd$intervention, nd$comparator)]) / sum(arm_weight)
  cumulative <- trial1$enrollment$cumulative * pair_fraction * scale

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

# ==== gsd: event- or calendar-driven NPH power ===============================

run_power <- function(id, enrollment_scale = 1) {
  g      <- gsd1[[id]]
  family <- g$efficacy$family
  schedule <- g$looks$schedule
  event <- if (schedule$scale == "count") schedule$values else NULL
  time  <- if (schedule$scale == "calendar-time") schedule$values else NULL

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
    enroll_rate = pairwise_enrollment_for(id, enrollment_scale),
    fail_rate = failure_table_for(id),
    ratio = allocation_ratio_for(id),
    analysis_time = time,
    event = event,
    upper = upper,
    upar = upar,
    test_lower = FALSE,
    interval = c(0.01, 100))
}

final_events <- function(x) tail(x$analysis$event, 1)

calibrate_to_events <- function(id) {
  g <- gsd1[[id]]
  if (g$looks$trigger == "event-count")
    return(list(enrollment_scale = 1, result = run_power(id)))

  target <- g$max_events
  event_at <- function(scale) final_events(run_power(id, scale))
  lower <- 0.01
  upper <- 1
  while (event_at(upper) < target) {
    upper <- 2 * upper
    if (upper > 1e5) stop("Could not match the planned final event count")
  }
  scale <- uniroot(
    function(value) event_at(value) - target,
    interval = c(lower, upper), tol = 1e-06)$root
  list(enrollment_scale = scale, result = run_power(id, scale))
}

results <- setNames(lapply(ids, calibrate_to_events), ids)

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
    `N Subjects`            = round(a$n, 1),
    AHR                     = round(a$ahr, 4),
    check.names = FALSE)

  state <- paste(id, "calculation")
  cat(sprintf("\nstate: %s\n%s\n", state, strrep("-", nchar(state) + 7)))
  implied_total_n <- reference_enrollment_total * sizing$enrollment_scale
  cat(sprintf("\n  %s  (local alpha %s, planned final %sevents %s, implied/reference total enrollment %.1f, %s-sided, %s vs %s; looks by %s; %s approximation; futility: none; calculation: power)\n",
              labels[id], format(nd$initial_alpha),
              if (g$looks$trigger == "calendar-time") "expected " else "",
              g$max_events, implied_total_n, nd$test_sides,
              nd$intervention, nd$comparator, g$looks$trigger, g$approximation))
  cat("    -- type I error control --\n")
  print(boundary, row.names = FALSE)
  cat("    -- type II error / futility control --\n")
  print(beta_futility_summary(gsd1[[id]], looks, "power"), row.names = FALSE)
  cat("    -- planned analyses --\n")
  print(analyses, row.names = FALSE)
  cat("    -- timeline and implied/reference enrollment --\n")
  print(timeline, row.names = FALSE)
}
cat("\n")
