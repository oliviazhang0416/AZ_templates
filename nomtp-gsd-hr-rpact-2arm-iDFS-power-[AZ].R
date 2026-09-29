# ============================================================
# Code title: Two-arm iDFS power calculation
#
# mtp:        none - one iDFS sample-size calculation
# gsd:        two-look sfLDOF efficacy design at information fractions 0.70 and 1
# coc:        no
# cop:        no
# simulation: no
# futility:   no
# endpoints:  iDFS
#
# Key inputs: one-sided alpha 0.025, HR 0.645, 860 subjects, 266 final
#   events, landmark iDFS rates, piecewise accrual, and 5% dropout by month 12
#
# Key outputs: cumulative events, analysis timing, efficacy boundaries,
#   nominal alpha, marginal power, and cumulative power
#
# Engine: rpact - rpact::getDesignGroupSequential(), rpact::getPowerSurvival()
# Input:  jsonlite - fromJSON() (trial1, mtp1, gsd1)
# ============================================================

# ==== READING THE DESIGN INPUTS =============================================

design_input <- jsonlite::fromJSON(
  "/Users/oliviazhang/Desktop/AZ_templates/nomtp-gsd-hr-rpact-2arm-iDFS-power-[AZ].json",
  simplifyVector = TRUE,
  simplifyDataFrame = FALSE,
  simplifyMatrix = TRUE)
trial1 <- design_input$trial
mtp1   <- design_input$mtp
gsd1   <- design_input$gsd
rm(design_input)

# ==== RUNNING THE CALCULATION ================================================

# ==== trial: allocation, accrual, dropout, and landmark iDFS assumptions =====

allocation_ratio_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  trial1$arms[[nd$intervention]]$randomization_weight /
    trial1$arms[[nd$comparator]]$randomization_weight
}

accrual_time <- c(0, cumsum(head(trial1$enrollment$duration, -1)))
accrual_rate <- trial1$enrollment$rate

landmark_hazards <- function(distribution) {
  rate <- -diff(log(distribution$landmark_value)) /
    diff(distribution$landmark_time)
  c(rate, tail(rate, 1))
}

reference_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  trial1$endpoints[[nd$endpoint]]$distributions[[nd$comparator]]
}

effect_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  trial1$endpoints[[nd$endpoint]]$effects[[nd$intervention]]
}

# ==== mtp: single calculation without multiplicity adjustment ================

ids <- names(mtp1$nodes)
labels <- vapply(ids, function(id) {
  nd <- mtp1$nodes[[id]]
  ep <- trial1$endpoints[[nd$endpoint]]
  sprintf("%s (%s [%s, %s], %s)",
          id, nd$endpoint, ep$type, ep$effect_measure, nd$label)
}, character(1))

# ==== gsd: iDFS power calculation ============================================

rpact_design_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  g  <- gsd1[[id]]

  rpact::getDesignGroupSequential(
    kMax = length(g$looks$schedule$values),
    sided = nd$test_sides,
    alpha = nd$initial_alpha,
    # Internal rpact placeholder only; achieved power is solved from the fixed
    # event count and beta remains null in the power JSON.
    beta = if (is.null(g$beta)) 0.1 else g$beta,
    informationRates = g$looks$schedule$values,
    typeOfDesign = "asOF")
}

run_power <- function(id) {
  nd     <- mtp1$nodes[[id]]
  g      <- gsd1[[id]]
  ref    <- reference_for(id)
  design <- rpact_design_for(id)

  args <- list(
    design = design,
    lambda2 = landmark_hazards(ref),
    piecewiseSurvivalTime = ref$landmark_time,
    hazardRatio = effect_for(id)$value,
    dropoutRate1 = trial1$dropout$probability,
    dropoutRate2 = trial1$dropout$probability,
    dropoutTime = trial1$dropout$window,
    allocationRatioPlanned = allocation_ratio_for(id),
    directionUpper = FALSE,
    accrualTime = accrual_time,
    accrualIntensity = accrual_rate,
    maxNumberOfSubjects = trial1$enrollment$n_accrual_max,
    maxNumberOfEvents = g$max_events,
    typeOfComputation = g$approximation)

  args$eventTime <- trial1$dropout$window
  result <- do.call(rpact::getPowerSurvival, args)

  list(design = design, result = result)
}

results <- setNames(lapply(ids, run_power), ids)

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
cat(sprintf("  enrollment cap: %s | dropout: %.1f%% by month %s\n",
            trial1$enrollment$n_accrual_max,
            100 * trial1$dropout$probability, trial1$dropout$window))

for (id in ids) {
  nd <- mtp1$nodes[[id]]
  x  <- results[[id]]
  d  <- x$design
  r  <- x$result
  k  <- length(d$informationRates)
  critical_hr <- if (!is.null(r$criticalValuesEffectScale))
    r$criticalValuesEffectScale else r$criticalValuesEffectScaleLower
  looks <- if (k == 1) "FA" else c(paste0("IA", seq_len(k - 1)), "FA")
  events <- as.numeric(r$cumulativeEventsPerStage)
  ratio <- allocation_ratio_for(id)
  info <- events * ratio / (1 + ratio)^2
  theta <- -log(effect_for(id)$value)

  boundary <- data.frame(
    Analysis       = looks,
    Boundary       = "efficacy",
    Z              = round(d$criticalValues, 3),
    CV_HR          = round(as.numeric(critical_hr), 3),
    `Nominal α`    = round(as.numeric(r$criticalValuesPValueScale), 5),
    `Cumulative α` = round(d$alphaSpent, 5),
    check.names = FALSE)

  analyses <- data.frame(
    Analysis               = looks,
    `Information Fraction` = round(d$informationRates, 3),
    `N Events`              = round(events, 1),
    `Marginal Power`        = round(pnorm(theta * sqrt(info) - d$criticalValues), 4),
    `Cumulative Power`      = round(cumsum(as.numeric(r$rejectPerStage)), 4),
    check.names = FALSE)

  timeline <- data.frame(
    Analysis                = looks,
    `Calendar Month`        = round(as.numeric(r$analysisTime), 1),
    `N Subjects`            = round(as.numeric(r$numberOfSubjects), 1),
    check.names = FALSE)

  state <- paste(id, "calculation")
  cat(sprintf("\nstate: %s\n%s\n", state, strrep("-", nchar(state) + 7)))
  cat(sprintf("\n  %s  (local alpha %s, %s-sided, %s vs %s; looks by %s; %s approximation; futility: none; calculation: power)\n",
              labels[id], format(nd$initial_alpha), nd$test_sides,
              nd$intervention, nd$comparator, gsd1[[id]]$looks$trigger,
              gsd1[[id]]$approximation))
  cat("    -- type I error control --\n")
  print(boundary, row.names = FALSE)
  cat("    -- type II error / futility control --\n")
  print(beta_futility_summary(gsd1[[id]], looks, "power"), row.names = FALSE)
  cat("    -- planned analyses --\n")
  print(analyses, row.names = FALSE)
  cat("    -- timeline and enrollment --\n")
  print(timeline, row.names = FALSE)
}
cat("\n")
