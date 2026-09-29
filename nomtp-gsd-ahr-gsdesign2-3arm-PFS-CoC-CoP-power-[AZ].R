# ============================================================
# Code title: Peri-operative primary, CoC, and CoP NPH calculations
#
# mtp:        none - primary, CoC, and CoP calculations are evaluated separately
# gsd:        primary calculation has four calendar-time looks with nonbinding
#   beta-spending futility at look 1 and efficacy at looks 2-4; CoC/CoP
#   calculations have three efficacy looks
# coc:        yes - Arm 3 is a component of Arm 2
# cop:        yes
# simulation: no
# futility:   nonbinding sfHSD beta spending for H1 only
# endpoints:  PFS under primary, CoC, and CoP NPH assumptions
#
# Key inputs: hypothesis-specific planned final expected events, reference
#   piecewise accrual shape, piecewise control hazards, randomized allocation,
#   piecewise hazard ratios, calendar-time looks, alpha, and dropout
#
# Key outputs: achieved power, implied enrollment, events and AHR, information
#   fractions, efficacy/futility boundaries, and nominal levels
#
# Engine: gsDesign2 - gsDesign2::gs_power_ahr()
# Input:  jsonlite - fromJSON() (trial1, mtp1, gsd1)
# ============================================================

# ==== READING THE DESIGN INPUTS =============================================

design_input <- jsonlite::fromJSON(
  "/Users/oliviazhang/Desktop/AZ_templates/nomtp-gsd-ahr-gsdesign2-3arm-PFS-CoC-CoP-power-[AZ].json",
  simplifyVector = TRUE,
  simplifyDataFrame = FALSE,
  simplifyMatrix = TRUE)
trial1 <- design_input$trial
mtp1   <- design_input$mtp
gsd1   <- design_input$gsd
rm(design_input)

# ==== RUNNING THE CALCULATION ================================================

# ==== trial: enrollment, dropout, allocation, and NPH endpoint assumptions ====

allocation_ratio_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  trial1$arms[[nd$intervention]]$randomization_weight /
    trial1$arms[[nd$comparator]]$randomization_weight
}

reference_enrollment_total <- sum(
  trial1$enrollment$duration * trial1$enrollment$rate)

enroll_rate_for <- function(scale) {
  tibble::tibble(
    stratum = "All",
    duration = trial1$enrollment$duration,
    rate = scale * trial1$enrollment$rate)
}

dropout_rate <- -log(1 - trial1$dropout$probability) / trial1$dropout$window

step_value <- function(cut, value, at) {
  value[findInterval(at, cut)]
}

failure_table_for <- function(id) {
  nd     <- mtp1$nodes[[id]]
  ep     <- trial1$endpoints[[nd$endpoint]]
  ref    <- ep$distributions[[nd$comparator]]
  effect <- ep$effects[[nd$intervention]]
  effect_cut <- if (is.null(effect$cut$at)) 0 else effect$cut$at
  effect_value <- if (is.null(effect$cut$at)) rep(effect$value, 1) else effect$value
  cut <- sort(unique(c(ref$cut$at, effect_cut)))

  tibble::tibble(
    stratum = "All",
    duration = c(diff(cut), Inf),
    fail_rate = step_value(ref$cut$at, ref$rate, cut),
    hr = step_value(effect_cut, effect_value, cut),
    dropout_rate = dropout_rate)
}

# ==== mtp: independent primary, CoC, and CoP calculations ===================

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

# ==== gsd: event-parameterized calendar-time efficacy and power ===============

spending_arguments <- function(family, total) {
  ans <- list(
    sf = getExportedValue("gsDesign", family$fn),
    total_spend = total)
  if (!is.null(family$param)) ans$param <- family$param
  ans
}

run_power <- function(id, enrollment_scale) {
  g <- gsd1[[id]]
  k <- length(g$looks$schedule$values)

  args <- list(
    enroll_rate = enroll_rate_for(enrollment_scale),
    fail_rate = failure_table_for(id),
    ratio = allocation_ratio_for(id),
    analysis_time = g$looks$schedule$values,
    event = NULL,
    upper = gsDesign2::gs_spending_bound,
    upar = spending_arguments(g$efficacy$family, alpha_1s_for(id)),
    test_upper = seq_len(k) %in% g$efficacy$active_at,
    test_lower = rep(FALSE, k),
    r = 80,
    tol = 1e-10)

  if (isTRUE(g$futility$enabled)) {
    args$binding <- g$futility$binding
    args$lower <- gsDesign2::gs_spending_bound
    args$lpar <- spending_arguments(g$futility$family, g$beta)
    args$test_lower <- seq_len(k) %in% g$futility$active_at
  }

  do.call(gsDesign2::gs_power_ahr, args)
}

final_events <- function(x) tail(x$analysis$event, 1)

calibrate_to_events <- function(id) {
  target <- gsd1[[id]]$max_events
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
allocation_text_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  paste(trial1$arms[[nd$intervention]]$randomization_weight,
        trial1$arms[[nd$comparator]]$randomization_weight, sep = ":")
}
ratio_text <- paste(sprintf("%s=%s", ids,
                            vapply(ids, allocation_text_for, character(1))),
                    collapse = ", ")
cat(sprintf(
  "  allocation ratios (intervention:comparator) %s | stratification: %s | simulation: %s\n",
  ratio_text, if (isTRUE(trial1$strata$enabled)) "yes" else "none",
  if (isTRUE(trial1$sim$enabled)) "yes" else "none (closed form)"))
cat(sprintf("  reference accrual-shape total: %.1f | dropout: %.1f%% by month %s\n",
            reference_enrollment_total,
            100 * trial1$dropout$probability, trial1$dropout$window))

for (id in ids) {
  nd <- mtp1$nodes[[id]]
  g  <- gsd1[[id]]
  sizing <- results[[id]]
  x  <- sizing$result
  a  <- x$analysis
  b  <- x$bound
  k  <- nrow(a)
  looks <- if (k == 1) "FA" else c(paste0("IA", seq_len(k - 1)), "FA")
  upper <- b[b$bound == "upper", , drop = FALSE]
  upper_index <- upper$analysis

  boundary <- data.frame(
    Analysis       = looks[upper_index],
    Boundary       = "efficacy",
    Z              = round(upper$z, 3),
    CV_HR          = round(upper$`~hr at bound`, 3),
    `Nominal α`    = round(nd$test_sides * upper$`nominal p`, 5),
    `Cumulative α` = round(nd$test_sides * upper$probability0, 5),
    check.names = FALSE)

  marginal_power <- rep(NA_real_, k)
  cumulative_power <- rep(0, k)
  marginal_power[upper_index] <- pnorm(
    a$theta[upper_index] * sqrt(a$info[upper_index]) - upper$z)
  cumulative_power[upper_index] <- upper$probability

  analyses <- data.frame(
    Analysis               = looks,
    `Information Fraction` = round(a$info_frac, 3),
    `N Events`              = round(a$event, 1),
    `Marginal Power`        = round(marginal_power, 4),
    `Cumulative Power`      = round(cumulative_power, 4),
    check.names = FALSE)

  timeline <- data.frame(
    Analysis                = looks,
    `Calendar Month`        = round(a$time, 1),
    `N Subjects`            = round(a$n, 1),
    AHR                     = round(a$ahr, 4),
    check.names = FALSE)

  state <- paste(id, "calculation")
  cat(sprintf("\nstate: %s\n%s\n", state, strrep("-", nchar(state) + 7)))
  cat(sprintf("\n  %s  (local alpha %s, planned final expected events %s, implied total enrollment %.1f, %s-sided, %s vs %s; looks by %s; %s approximation; futility: %s; calculation: power)\n",
              labels[id], format(nd$initial_alpha), g$max_events,
              reference_enrollment_total * sizing$enrollment_scale,
              nd$test_sides,
              nd$intervention, nd$comparator, g$looks$trigger, g$approximation,
              if (isTRUE(g$futility$enabled)) "nonbinding" else "none"))
  cat("    -- type I error control --\n")
  print(boundary, row.names = FALSE)
  cat("    -- type II error / futility control --\n")
  print(beta_futility_summary(gsd1[[id]], looks, "power"), row.names = FALSE)
  if (isTRUE(g$futility$enabled)) {
    lower <- b[b$bound == "lower", , drop = FALSE]
    lower_index <- lower$analysis
    futility <- data.frame(
      Analysis       = looks[lower_index],
      Boundary       = "futility",
      Z              = round(lower$z, 3),
      CV_HR          = round(lower$`~hr at bound`, 3),
      `Nominal β`    = round(pnorm(
        lower$z - a$theta[lower_index] * sqrt(a$info[lower_index])), 5),
      `Cumulative β` = round(lower$probability, 5),
      check.names = FALSE)
    print(futility, row.names = FALSE)
  }

  cat("    -- planned analyses --\n")
  print(analyses, row.names = FALSE)
  cat("    -- timeline and implied enrollment --\n")
  print(timeline, row.names = FALSE)
}
cat("\n")
