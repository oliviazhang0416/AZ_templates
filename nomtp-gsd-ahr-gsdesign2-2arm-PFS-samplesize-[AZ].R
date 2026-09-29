# ============================================================
# Code title: Two-arm PFS event-size calculation under a piecewise control hazard and delayed treatment effect
#
# mtp:        none - one PFS calculation
# gsd:        one calendar-time analysis at month 50, sfLDOF efficacy boundary
# coc:        no
# cop:        no
# simulation: no
# futility:   no
# endpoints:  PFS under a delayed treatment effect
#
# Key inputs: target power 85%, 36-month power-model enrollment with shape 1.5,
#   1:1 randomization, piecewise control hazards, 10% annual dropout,
#   one-sided alpha 0.025, and month-50 data cutoff
#
# Key outputs: required final events, average hazard ratio, efficacy critical
#   value, nominal alpha, and achieved power
#
# Engine: gsDesign2 - gsDesign2::gs_design_ahr()
# Input:  jsonlite - fromJSON() (trial1, mtp1, gsd1)
# ============================================================

# ==== READING THE DESIGN INPUTS =============================================

design_input <- jsonlite::fromJSON(
  "/Users/oliviazhang/Desktop/AZ_templates/nomtp-gsd-ahr-gsdesign2-2arm-PFS-samplesize-[AZ].json",
  simplifyVector = TRUE,
  simplifyDataFrame = FALSE,
  simplifyMatrix = TRUE)
trial1 <- design_input$trial
mtp1   <- design_input$mtp
gsd1   <- design_input$gsd
rm(design_input)

# ==== RUNNING THE CALCULATION ================================================

# ==== trial: enrollment, dropout, arm allocation, and PFS assumptions =========

allocation_ratio_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  trial1$arms[[nd$intervention]]$randomization_weight /
    trial1$arms[[nd$comparator]]$randomization_weight
}

enrollment_month <- seq_len(trial1$enrollment$period)
enroll_rate_for <- function() {
  cumulative_enrollment <-
    (enrollment_month / trial1$enrollment$period)^trial1$enrollment$k
  tibble::tibble(
    stratum = "All",
    duration = rep(1, length(enrollment_month)),
    rate = cumulative_enrollment - c(0, head(cumulative_enrollment, -1)))
}

dropout_rate <- -log(1 - trial1$dropout$probability) / trial1$dropout$window

effect_at <- function(effect, at) {
  if (is.null(effect$cut$at)) return(rep(effect$value, length(at)))
  effect$value[findInterval(at, effect$cut$at)]
}

failure_table_for <- function(id) {
  nd     <- mtp1$nodes[[id]]
  ep     <- trial1$endpoints[[nd$endpoint]]
  ref    <- ep$distributions[[nd$comparator]]
  effect <- ep$effects[[nd$intervention]]
  cut    <- ref$cut$at

  tibble::tibble(
    stratum = "All",
    duration = c(diff(cut), 100),
    fail_rate = ref$rate,
    hr = effect_at(effect, cut),
    dropout_rate = dropout_rate)
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

# ==== gsd: calendar-time NPH model with event-size solution ==================

run_event_size <- function(id) {
  nd     <- mtp1$nodes[[id]]
  g      <- gsd1[[id]]
  family <- g$efficacy$family

  design <- gsDesign2::gs_design_ahr(
    enroll_rate = enroll_rate_for(),
    fail_rate = failure_table_for(id),
    ratio = allocation_ratio_for(id),
    alpha = alpha_1s_for(id),
    beta = g$beta,
    analysis_time = g$looks$schedule$values,
    upper = gsDesign2::gs_spending_bound,
    upar = list(
      sf = getExportedValue("gsDesign", family$fn),
      total_spend = alpha_1s_for(id),
      param = family$param),
    test_upper = seq_along(g$looks$schedule$values) %in% g$efficacy$active_at,
    test_lower = FALSE,
    interval = c(0.01, 100))
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
cat(sprintf(
  "  enrollment shape: power model over %s months | dropout: %.1f%% by month %s\n",
  trial1$enrollment$period,
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
