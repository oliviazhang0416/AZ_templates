# ============================================================
# Code title: Three-arm EFS and OS group-sequential power with a hybrid AND gate
#
# mtp:        hybrid - Dunnett for H1/H2, weighted Bonferroni elsewhere, and
#   epsilon edges from EFS to H3, followed by the ordered H3-to-H4 transfer
# gsd:        EFS uses 2-look sfLDOF designs; OS uses 4-look sfHSD designs
# coc:        no
# cop:        no
# simulation: no
# futility:   no
# endpoints:  EFS and OS, each comparing two experimental arms with a common control
#
# Key inputs: two-sided family alpha 0.05, 1:1:1 randomization, Dunnett
#   correlation 0.5, EFS piecewise control hazards, OS control median 40 months,
#   hypothesis-specific HRs, information fractions, fixed event targets,
#   epsilon 0.0001, and 400 subjects
#
# Key outputs: every graph-reachable state, graph and effective nominal alpha
#   for every remaining hypothesis, and one set of group-sequential boundary,
#   analysis-time, and achieved-power tables per unique hypothesis-by-alpha design
#
# Engine: graphicalMCP - graphicalMCP::graph_create(),
#         graphicalMCP::graph_generate_weights(); mvtnorm - mvtnorm::qmvt();
#         rpact - rpact::getDesignGroupSequential(), rpact::getPowerSurvival()
# Input:  jsonlite - fromJSON() (trial1, mtp1, gsd1)
# ============================================================

# ==== READING THE DESIGN INPUTS =============================================

design_input <- jsonlite::fromJSON(
  "/Users/oliviazhang/Desktop/AZ_templates/hybrid-gsd-hr-rpact-3arm-EFS-OS-power-[AZ].json",
  simplifyVector = TRUE,
  simplifyDataFrame = FALSE,
  simplifyMatrix = TRUE)
trial1 <- design_input$trial
mtp1   <- design_input$mtp
gsd1   <- design_input$gsd
rm(design_input)

# ==== RUNNING THE CALCULATION ================================================

# ==== trial: arm allocation and endpoint assumptions =========================

allocation_ratio_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  trial1$arms[[nd$intervention]]$randomization_weight /
    trial1$arms[[nd$comparator]]$randomization_weight
}

effect_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  trial1$endpoints[[nd$endpoint]]$effects[[nd$intervention]]
}

reference_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  trial1$endpoints[[nd$endpoint]]$distributions[[nd$comparator]]
}

# ==== mtp: Dunnett EFS test and epsilon-edge AND gate to OS ===================

ids <- names(mtp1$nodes)
dunnett_test <- mtp1$intersection_tests$H12
dunnett_alpha_1s <- trial1$alpha$value / trial1$alpha$sidedness
dunnett_critical <- mvtnorm::qmvt(
  p = 1 - dunnett_alpha_1s,
  sigma = dunnett_test$correlation,
  tail = "lower.tail",
  df = Inf)$quantile
dunnett_alpha_2s <- 2 * pnorm(-dunnett_critical)

node_alpha <- vapply(mtp1$nodes, function(nd) nd$initial_alpha, numeric(1))
G <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
for (ed in mtp1$edges) G[ed$from, ed$to] <- ed$weight

graph <- suppressWarnings(graphicalMCP::graph_create(
  node_alpha / trial1$alpha$value, G, hyp_names = ids))
graph_weights <- graphicalMCP::graph_generate_weights(graph)
intersections <- graph_weights[, seq_along(ids), drop = FALSE]
local_weights <- graph_weights[, length(ids) + seq_along(ids), drop = FALSE]
colnames(local_weights) <- ids

efs_ids <- dunnett_test$hypotheses
local_alpha <- trial1$alpha$value * local_weights

# Traverse every nonterminal state reachable by rejecting a hypothesis with
# positive graph alpha. This includes mathematically reachable epsilon branches.
alpha_tolerance <- 1e-15
state_key <- function(x) paste(as.integer(x), collapse = "")
row_for_key <- setNames(
  seq_len(nrow(intersections)), apply(intersections, 1, state_key))
queue <- state_key(rep(1L, length(ids)))
seen <- character()
reachable <- list()

while (length(queue)) {
  key <- queue[1]
  queue <- queue[-1]
  if (key %in% seen) next

  seen <- c(seen, key)
  row <- unname(row_for_key[[key]])
  avail <- ids[intersections[row, ] == 1]
  reachable[[length(reachable) + 1L]] <- list(
    key = key, row = row, avail = avail)

  rejectable <- avail[local_alpha[row, avail] > alpha_tolerance]
  for (id in rejectable) {
    child <- intersections[row, ]
    child[match(id, ids)] <- 0
    if (!any(child == 1)) next
    child_key <- state_key(child)
    if (!(child_key %in% c(seen, queue))) queue <- c(queue, child_key)
  }
}

# ==== gsd: rpact group-sequential designs and survival power ==================

rpact_design_for <- function(id, effective_alpha) {
  nd     <- mtp1$nodes[[id]]
  g      <- gsd1[[id]]
  family <- g$efficacy$family

  args <- list(
    kMax = length(g$looks$schedule$values),
    sided = nd$test_sides,
    alpha = effective_alpha,
    # rpact requires beta to instantiate the design object; for an
    # efficacy-only power calculation this internal value does not define the
    # achieved power or efficacy spending boundaries.
    beta = if (is.null(g$beta)) 0.2 else g$beta,
    informationRates = g$looks$schedule$values,
    typeOfDesign = switch(family$fn, sfLDOF = "asOF", sfHSD = "asHSD"))
  if (family$fn == "sfHSD") args$gammaA <- family$param
  do.call(rpact::getDesignGroupSequential, args)
}

power_for <- function(id, effective_alpha) {
  nd     <- mtp1$nodes[[id]]
  ref    <- reference_for(id)
  effect <- effect_for(id)
  design <- rpact_design_for(id, effective_alpha)

  args <- list(
    design = design,
    hazardRatio = effect$value,
    maxNumberOfSubjects = trial1$enrollment$n_accrual_max,
    maxNumberOfEvents = gsd1[[id]]$max_events,
    allocationRatioPlanned = allocation_ratio_for(id))

  if (ref$family == "piecewise-exponential") {
    args$lambda2 <- ref$rate
    args$piecewiseSurvivalTime <- ref$cut$at
  } else {
    args$lambda2 <- log(2) / ref$median
  }

  list(design = design, result = do.call(rpact::getPowerSurvival, args))
}

labels <- vapply(ids, function(id) {
  nd <- mtp1$nodes[[id]]
  ep <- trial1$endpoints[[nd$endpoint]]
  sprintf("%s (%s [%s, %s], %s)",
          id, nd$endpoint, ep$type, ep$effect_measure, nd$label)
}, character(1))

# ==== mtp/gsd: reachable states and unique calculations ======================

local_test_for <- function(avail, id) {
  if (id %in% efs_ids && all(efs_ids %in% avail)) return("dunnett (H12)")
  if (length(avail) == 1) return("singleton")
  "weighted-bonferroni"
}

effective_alpha_for <- function(avail, id, graph_alpha) {
  if (id %in% efs_ids && all(efs_ids %in% avail))
    return(dunnett_alpha_2s)
  graph_alpha
}

# Catalogue all remaining hypotheses at every reachable state, then assign a
# reusable calculation ID to each distinct hypothesis-by-effective-alpha design.
state_records <- do.call(rbind, lapply(seq_along(reachable), function(i) {
  r <- reachable[[i]]
  graph_alpha <- as.numeric(local_alpha[r$row, r$avail])
  local_test <- vapply(r$avail, function(id)
    local_test_for(r$avail, id), character(1))
  effective_alpha <- vapply(seq_along(r$avail), function(j)
    effective_alpha_for(r$avail, r$avail[j], graph_alpha[j]), numeric(1))
  epsilon_level <- graph_alpha > alpha_tolerance &
    graph_alpha <= trial1$alpha$value * mtp1$procedure$epsilon

  data.frame(
    State = paste0("S", i),
    Available = paste(r$avail, collapse = ", "),
    Hypothesis = r$avail,
    GraphAlpha = graph_alpha,
    EffectiveAlpha = effective_alpha,
    LocalTest = local_test,
    Status = ifelse(
      graph_alpha <= alpha_tolerance, "not currently testable",
      ifelse(epsilon_level, "testable (epsilon-level)", "testable")),
    stringsAsFactors = FALSE)
}))

state_records$Calculation <- "-"
configuration_keys <- character()
configurations <- list()
for (i in seq_len(nrow(state_records))) {
  if (state_records$GraphAlpha[i] <= alpha_tolerance) next
  key <- paste(
    state_records$Hypothesis[i],
    format(state_records$EffectiveAlpha[i], digits = 16, scientific = TRUE),
    sep = "|")
  config_index <- match(key, configuration_keys)
  if (is.na(config_index)) {
    configuration_keys <- c(configuration_keys, key)
    config_index <- length(configuration_keys)
    calc_id <- paste0("C", config_index)
    configurations[[calc_id]] <- list(
      id = state_records$Hypothesis[i],
      effective_alpha = state_records$EffectiveAlpha[i])
  }
  state_records$Calculation[i] <- paste0("C", config_index)
}

for (calc_id in names(configurations)) {
  used <- state_records$Calculation == calc_id
  configurations[[calc_id]]$states <- unique(state_records$State[used])
  configurations[[calc_id]]$graph_alphas <- unique(state_records$GraphAlpha[used])
  configurations[[calc_id]]$local_tests <- unique(state_records$LocalTest[used])
}
calculation_results <- lapply(configurations, function(x)
  power_for(x$id, x$effective_alpha))

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


format_alpha <- function(x) format(x, digits = 8, scientific = FALSE, trim = TRUE)
format_table_alpha <- function(x) {
  vapply(x, function(value) {
    if (value == 0) return("0")
    if (abs(value) < 1e-4)
      return(format(value, digits = 5, scientific = TRUE, trim = TRUE))
    formatC(value, format = "f", digits = 5)
  }, character(1))
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
  "  Dunnett H12 correlation: %.3f | critical Z: %.4f | marginal two-sided alpha: %.6f | AND-gate epsilon: %s\n",
  dunnett_test$correlation[1, 2], dunnett_critical, dunnett_alpha_2s,
  format(mtp1$procedure$epsilon)))

cat("\n-- reachable-state catalogue --\n")
state_table <- data.frame(
  State = state_records$State,
  Available = state_records$Available,
  Hypothesis = state_records$Hypothesis,
  `Graph alpha` = vapply(state_records$GraphAlpha, format_alpha, character(1)),
  `Effective alpha` = vapply(
    state_records$EffectiveAlpha, format_alpha, character(1)),
  `Local test` = state_records$LocalTest,
  Status = state_records$Status,
  Calculation = state_records$Calculation,
  check.names = FALSE)
old_width <- getOption("width")
options(width = max(old_width, 180))
print(state_table, row.names = FALSE, right = FALSE)
options(width = old_width)
cat("  terminal state: no hypotheses available\n")

cat("\n-- unique detailed calculations --\n")
for (calc_id in names(configurations)) {
  config <- configurations[[calc_id]]
  id <- config$id
  nd <- mtp1$nodes[[id]]
  x  <- calculation_results[[calc_id]]
  d  <- x$design
  ss <- x$result
  k  <- length(d$informationRates)
  looks <- if (k == 1) "FA" else c(paste0("IA", seq_len(k - 1)), "FA")
  events <- as.numeric(ss$cumulativeEventsPerStage)
  ratio <- allocation_ratio_for(id)
  info <- events * ratio / (1 + ratio)^2
  theta <- -log(effect_for(id)$value)

  boundary <- data.frame(
    Analysis       = looks,
    Boundary       = "efficacy",
    Z              = round(d$criticalValues, 3),
    CV_HR          = round(ss$criticalValuesEffectScaleLower, 3),
    `Nominal α`    = format_table_alpha(ss$criticalValuesPValueScale),
    `Cumulative α` = format_table_alpha(d$alphaSpent),
    check.names = FALSE)

  analyses <- data.frame(
    Analysis               = looks,
    `Information Fraction` = round(d$informationRates, 3),
    `N Events`              = round(events, 1),
    `Marginal Power`        = round(pnorm(theta * sqrt(info) - d$criticalValues), 4),
    `Cumulative Power`      = round(cumsum(as.numeric(ss$rejectPerStage)), 4),
    check.names = FALSE)

  timeline <- data.frame(
    Analysis                = looks,
    `Calendar Month`        = round(as.numeric(ss$analysisTime), 1),
    `N Subjects`            = round(as.numeric(ss$numberOfSubjects), 1),
    check.names = FALSE)

  cat(sprintf("\ncalculation %s: %s\n", calc_id, labels[id]))
  cat(sprintf(
    "  used by: %s | graph alpha(s) %s | effective nominal alpha %s | local test context(s): %s\n",
    paste(config$states, collapse = ", "),
    paste(vapply(config$graph_alphas, format_alpha, character(1)), collapse = ", "),
    format_alpha(config$effective_alpha),
    paste(config$local_tests, collapse = "; ")))
  cat(sprintf("  %s-sided, %s vs %s; looks by %s; %s approximation; futility: none\n",
              nd$test_sides,
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
