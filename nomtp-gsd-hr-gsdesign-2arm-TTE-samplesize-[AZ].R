# ============================================================
# Code title: Two-arm TTE event-size calculation under a proportional-hazards design assumption
#
# mtp:        none - one TTE calculation
# gsd:        two information-fraction analyses, sfLDOF efficacy boundary
# coc:        no
# cop:        no
# simulation: no
# futility:   no
# endpoints:  one generic time-to-event endpoint
#
# Key inputs: target power 94.4%, 70% interim information, 1:1 randomization,
#   design hazard ratio 0.645, and two-sided alpha 0.05
#
# Key outputs: required final events, efficacy critical values, nominal and
#   cumulative alpha, planned event counts, and achieved power
#
# Engine: gsDesign - gsDesign::gsDesign(), gsDesign::gsProbability()
# Input:  jsonlite - fromJSON() (trial1, mtp1, gsd1)
# ============================================================

# ==== READING THE DESIGN INPUTS =============================================

design_input <- jsonlite::fromJSON(
  "/Users/oliviazhang/Desktop/AZ_templates/nomtp-gsd-hr-gsdesign-2arm-TTE-samplesize-[AZ].json",
  simplifyVector = TRUE,
  simplifyDataFrame = FALSE,
  simplifyMatrix = TRUE)
trial1 <- design_input$trial
mtp1   <- design_input$mtp
gsd1   <- design_input$gsd
rm(design_input)

# ==== RUNNING THE CALCULATION ================================================

# ==== trial: arm allocation and proportional-hazards assumption =============

allocation_ratio_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  trial1$arms[[nd$intervention]]$randomization_weight /
    trial1$arms[[nd$comparator]]$randomization_weight
}

effect_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  trial1$endpoints[[nd$endpoint]]$effects[[nd$intervention]]$value
}

# ==== mtp: single calculation without multiplicity adjustment ================

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

# ==== gsd: event-size solution under gsDesign ================================

run_power <- function(id, max_events) {
  nd     <- mtp1$nodes[[id]]
  g      <- gsd1[[id]]
  events <- g$looks$schedule$values * max_events
  ratio  <- allocation_ratio_for(id)
  hr     <- effect_for(id)
  alpha  <- alpha_1s_for(id)
  family <- g$efficacy$family

  design_arguments <- list(
    k = length(events),
    test.type = 1,
    alpha = alpha,
    beta = g$beta,
    timing = if (length(events) > 1) head(events / max(events), -1) else NULL,
    sfu = getExportedValue("gsDesign", family$fn))
  if (!is.null(family$param)) design_arguments$sfupar <- family$param
  design <- do.call(gsDesign::gsDesign, design_arguments)

  z <- as.numeric(design$upper$bound)
  alpha_nominal <- 1 - pnorm(z)
  alpha_cumulative <- cumsum(as.numeric(design$upper$spend))
  information <- events * ratio / (1 + ratio)^2
  crossing <- gsDesign::gsProbability(
    k = length(events),
    theta = -log(hr),
    n.I = information,
    a = rep(-20, length(events)),
    b = z)

  list(
    events = events,
    z = z,
    critical_hr = exp(-z / sqrt(information)),
    alpha_nominal = alpha_nominal,
    alpha_cumulative = alpha_cumulative,
    power_marginal = pnorm(-log(hr) * sqrt(information) - z),
    power_cumulative = cumsum(as.numeric(crossing$upper$prob)))
}

solve_event_size <- function(id) {
  target_power <- 1 - gsd1[[id]]$beta
  final_power <- function(max_events)
    tail(run_power(id, max_events)$power_cumulative, 1)

  lower <- 2
  upper <- 100
  while (final_power(upper) < target_power) upper <- 2 * upper
  required_events <- ceiling(uniroot(
    function(max_events) final_power(max_events) - target_power,
    interval = c(lower, upper), tol = 0.001)$root)

  list(
    required_events = required_events,
    result = run_power(id, required_events))
}

results <- setNames(lapply(ids, solve_event_size), ids)

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
  ratio_text,
  if (isTRUE(trial1$strata$enabled)) "yes" else "none",
  if (isTRUE(trial1$sim$enabled)) "yes" else "none (closed form)"))

for (id in ids) {
  nd <- mtp1$nodes[[id]]
  g  <- gsd1[[id]]
  sizing <- results[[id]]
  x <- sizing$result
  looks <- if (length(x$events) == 1) "FA" else
    c(paste0("IA", seq_len(length(x$events) - 1)), "FA")

  boundary <- data.frame(
    Analysis = looks,
    Boundary = "efficacy",
    Z = round(x$z, 3),
    CV_HR = round(x$critical_hr, 3),
    `Nominal α` = round(nd$test_sides * x$alpha_nominal, 5),
    `Cumulative α` = round(nd$test_sides * x$alpha_cumulative, 5),
    check.names = FALSE)

  analyses <- data.frame(
    Analysis = looks,
    `Information Fraction` = round(x$events / max(x$events), 3),
    `N Events` = round(x$events, 1),
    `Marginal Power` = round(x$power_marginal, 4),
    `Cumulative Power` = round(x$power_cumulative, 4),
    check.names = FALSE)

  state <- paste(id, "calculation")
  cat(sprintf("\nstate: %s\n%s\n", state, strrep("-", nchar(state) + 7)))
  cat(sprintf(
    "\n  %s  (local alpha %s, target power %.1f%%, required final events %s, design HR %s, %s-sided, %s vs %s; looks by %s; %s approximation; futility: none; calculation: sample size [events])\n",
    labels[id], format(nd$initial_alpha), 100 * (1 - g$beta),
    sizing$required_events, format(effect_for(id)), nd$test_sides,
    nd$intervention, nd$comparator, g$looks$trigger, g$approximation))
  cat("    -- type I error control --\n")
  print(boundary, row.names = FALSE)
  cat("    -- type II error / futility control --\n")
  print(beta_futility_summary(gsd1[[id]], looks, "sample-size"), row.names = FALSE)
  cat("    -- planned analyses --\n")
  print(analyses, row.names = FALSE)
}
cat("\n")
