# ============================================================
# Code title: Group-sequential event-size calculation for a fixed-sequence PFS/OS graph with a delayed-effect CoC
#
# mtp:        graph - H1 (PFS) -> H2 (OS) -> H3 (CoC vs PFS), weights 4/5, 1/5, 0
# gsd:        H1 2 looks (sfLDOF), H2 3 looks (sfLDOF), H3 3 looks (customized nominal alpha)
# coc:        yes - H3, a non-randomised "component of arm 2" with a delayed
#   treatment effect stated piecewise BY EVENT COUNT
# cop:        no
# simulation: no
# futility:   no
# endpoints:  PFS (hazard-ratio), OS (hazard-ratio)
#
# Key inputs: alpha (0.05), target beta by hypothesis, PFS hazard ratio (0.65), OS hazard ratio (0.70),
#   CoC delayed effect (no effect for 97 events, then a constant hazard ratio
#   calibrated to an average of 0.82 at 380 events), and information fractions
#
# Key outputs: every graph-reachable state, the local alpha and test status of
#   every remaining hypothesis, and required final event counts with boundary
#   and power tables for each unique hypothesis-by-alpha calculation
#
# Engine: gsDesign - gsDesign::gsDesign(), gsDesign::gsProbability() (H1, H2)
#         gsDesign2 - gsDesign2::gs_power_npe() (H3, average HR by event count)
#         graphicalMCP - graphicalMCP::graph_create(),
#         graphicalMCP::graph_generate_weights() (alpha propagation)
# Input:  jsonlite - fromJSON() (trial1, mtp1, gsd1)
#
# Code status: AZ provided, validated for project use
# ============================================================

# ==== READING THE DESIGN INPUTS =============================================

design_input <- jsonlite::fromJSON(
  "/Users/oliviazhang/Desktop/AZ_templates/fixedseq-gsd-hr-gsdesign-3arm-PFS-OS-samplesize-[AZ].json",
  simplifyVector = TRUE,
  simplifyDataFrame = FALSE,
  simplifyMatrix = TRUE)
trial1 <- design_input$trial
mtp1   <- design_input$mtp
gsd1   <- design_input$gsd
rm(design_input)

# ==== RUNNING THE CALCULATION ================================================

# ==== trial: arm allocation and endpoint metadata ============================

randomization_weight_for <- function(arm_id) {
  arm <- trial1$arms[[arm_id]]
  if (is.null(arm$randomization_weight))
    return(trial1$arms[[arm$component_of]]$randomization_weight)
  arm$randomization_weight
}

allocation_ratio_for <- function(id) {
  nd <- mtp1$nodes[[id]]
  randomization_weight_for(nd$intervention) /
    randomization_weight_for(nd$comparator)
}

ids    <- names(mtp1$nodes)
# labels: "H1 (PFS [time-to-event, hazard-ratio], Arm 2 vs 1)" headers; nd$endpoint is a string key into trial1$endpoints
labels <- vapply(ids, function(id) {
  nd <- mtp1$nodes[[id]]
  ep <- trial1$endpoints[[nd$endpoint]]
  sprintf("%s (%s [%s, %s], %s)", id, nd$endpoint, ep$type, ep$effect_measure, nd$label)
}, character(1))

# ==== mtp: alpha propagation across the fixed-sequence graph H1 -> H2 -> H3 ==

# each node's one-sided alpha, normalised to graph_create()'s starting weights - gated on the
# node initial_alpha values actually summing to trial1$alpha$value
sides      <- vapply(mtp1$nodes, function(nd) nd$test_sides, numeric(1))
node_alpha <- vapply(mtp1$nodes, function(nd) nd$initial_alpha, numeric(1))
a_1s       <- ifelse(sides == 2, node_alpha / 2, node_alpha)
alpha_1s   <- if (isTRUE(all.equal(sum(node_alpha), trial1$alpha$value))) sum(a_1s)

if (is.null(alpha_1s))
  stop("Initial node alpha must sum to trial1$alpha$value")
if (any(sides != trial1$alpha$sidedness))
  stop("All nodes must use the sidedness declared in trial1$alpha")

# G: the transition matrix, H1 -> H2 -> H3, from mtp1's edges
m <- length(ids)
G <- matrix(0, m, m, dimnames = list(ids, ids))
for (ed in mtp1$edges) G[ed$from, ed$to] <- G[ed$from, ed$to] + ed$weight

# the closure's weights at every intersection the graph can reach
graph <- suppressWarnings(graphicalMCP::graph_create(
  unname(a_1s / alpha_1s), unname(G), hyp_names = ids))
ws <- graphicalMCP::graph_generate_weights(graph)
inter <- ws[, seq_len(m), drop = FALSE]

# local_alpha: one-sided alpha per hypothesis at each state, rows labelled by who's still available
local_alpha <- alpha_1s * ws[, m + seq_len(m), drop = FALSE]
colnames(local_alpha) <- ids
rownames(local_alpha) <- apply(inter, 1, function(r)
  paste(paste(ids[r == 1], collapse = ", "), "available"))

# Traverse every nonterminal state that can be reached by rejecting a
# hypothesis holding positive local alpha. Breadth-first traversal makes the
# result independent of a hard-coded rejection order.
alpha_tolerance <- 1e-15
state_key <- function(x) paste(as.integer(x), collapse = "")
row_for_key <- setNames(seq_len(nrow(inter)), apply(inter, 1, state_key))
queue <- state_key(rep(1L, m))
seen <- character()
reachable <- list()

while (length(queue)) {
  key <- queue[1]
  queue <- queue[-1]
  if (key %in% seen) next

  seen <- c(seen, key)
  row <- unname(row_for_key[[key]])
  avail <- ids[inter[row, ] == 1]
  reachable[[length(reachable) + 1L]] <- list(
    key = key, row = row, avail = avail)

  rejectable <- avail[local_alpha[row, avail] > alpha_tolerance]
  for (id in rejectable) {
    child <- inter[row, ]
    child[match(id, ids)] <- 0
    if (!any(child == 1)) next
    child_key <- state_key(child)
    if (!(child_key %in% c(seen, queue))) queue <- c(queue, child_key)
  }
}

# ==== gsd: boundary and power per hypothesis, under gsDesign / gsDesign2 =====

has_futility <- function(x) isTRUE(x$enabled)

# gsdesign_surv_power: constant-HR Schoenfeld power under gsDesign (H1, H2)
gsdesign_surv_power <- function(events, alpha_1s, hr, ratio, beta, family) {
  k    <- length(events)
  frac <- events / max(events)
  # gsDesign() requires beta even though fixed-event power does not solve for it.
  beta_gsd <- if (is.null(beta)) 0.1 else beta
  des <- gsDesign::gsDesign(
    k = k, test.type = 1, alpha = alpha_1s, beta = beta_gsd,
                  timing = frac,
                  sfu = getExportedValue("gsDesign", family$fn),
                  sfupar = family$param)
  z         <- des$upper$bound
  alpha_nom <- 1 - pnorm(z)
  spend_cum <- cumsum(as.numeric(des$upper$spend))

  info_i  <- events * ratio / (1 + ratio)^2
  a_bound <- rep(-20, k)
  xp <- gsDesign::gsProbability(
    k = k, theta = -log(hr), n.I = info_i, a = a_bound, b = z)
  power          <- cumsum(as.numeric(xp$upper$prob))
  power_marginal <- pnorm(-log(hr) * sqrt(info_i) - z)
  cv_hr <- exp(qnorm(alpha_nom) * sqrt((1 + ratio)^2 / (ratio * events)))

  list(events = events, z = z, alpha_nom = alpha_nom, spend_cum = spend_cum,
       cv_hr = cv_hr, power = power, power_marginal = power_marginal)
}

# gsdesign2_npe_power: custom nominal-alpha boundary and NPH power for H3
gsdesign2_npe_power <- function(events, ahr, ratio, nominal_alpha) {
  k     <- length(events)
  theta <- -log(ahr)
  info  <- events * ratio / (1 + ratio)^2

  x <- gsDesign2::gs_power_npe(
    theta = theta, info = info, upper = gsDesign2::gs_b,
                    upar = qnorm(1 - nominal_alpha), lower = gsDesign2::gs_b,
                    lpar = rep(-Inf, k), test_lower = FALSE)
  u <- x[x$bound == "upper", , drop = FALSE]
  z     <- as.numeric(u$z)[order(u$analysis)]
  power <- as.numeric(u$probability)[order(u$analysis)]

  x0 <- gsDesign2::gs_power_npe(
    theta = rep(0, k), info = info, upper = gsDesign2::gs_b,
                     upar = qnorm(1 - nominal_alpha), lower = gsDesign2::gs_b,
                     lpar = rep(-Inf, k), test_lower = FALSE)
  u0 <- x0[x0$bound == "upper", , drop = FALSE]
  spend_cum <- as.numeric(u0$probability)[order(u0$analysis)]

  alpha_nom      <- 1 - pnorm(z)
  power_marginal <- pnorm(theta * sqrt(info) - z)
  cv_hr <- exp(qnorm(alpha_nom) * sqrt((1 + ratio)^2 / (ratio * events)))

  list(events = events, z = z, alpha_nom = alpha_nom, spend_cum = spend_cum,
       cv_hr = cv_hr, power = power, power_marginal = power_marginal, ahr = ahr)
}

# effect_for: this hypothesis's endpoint effect
effect_for <- function(id) {
  nd  <- mtp1$nodes[[id]]
  trial1$endpoints[[nd$endpoint]]$effects[[nd$intervention]]
}

eff_H1 <- effect_for("H1")
eff_H2 <- effect_for("H2")
hr_pfs <- eff_H1$value
hr_os  <- eff_H2$value
coc    <- effect_for("H3")

# reference (comparator) distribution per endpoint - every distribution is NULL here
ref_curves <- lapply(ids, function(id) {
  nd <- mtp1$nodes[[id]]
  r  <- trial1$endpoints[[nd$endpoint]]$distributions[[nd$comparator]]
  if (!is.null(r)) r
})
names(ref_curves) <- ids

# ahr_by_events: H3's average HR by event count
ahr_by_events <- function(events) {
  upper <- c(coc$cut$at[-1], Inf)
  vapply(events, function(d) {
    width <- pmax(0, pmin(upper, d) - coc$cut$at)
    exp(sum(width * log(coc$value)) / d)
  }, numeric(1))
}

events_for <- function(id, max_events) {
  g <- gsd1[[id]]
  g$looks$schedule$values * max_events
}

custom_alpha_for <- function(id, alpha_1s) {
  g   <- gsd1[[id]]
  eff <- g$efficacy
  a <- eff$family$nominal_values / mtp1$nodes[[id]]$test_sides
  replace(a, is.na(a), alpha_1s - sum(a, na.rm = TRUE))
}

power_for <- function(id, alpha_1s, max_events) {
  g     <- gsd1[[id]]
  ratio <- allocation_ratio_for(id)
  events <- events_for(id, max_events)
  switch(id,
         H1 = gsdesign_surv_power(events, alpha_1s, hr_pfs,
                                   ratio, g$beta, g$efficacy$family),
         H2 = gsdesign_surv_power(events, alpha_1s, hr_os,
                                   ratio, g$beta, g$efficacy$family),
         H3 = gsdesign2_npe_power(events, ahr_by_events(events), ratio,
                                   custom_alpha_for(id, alpha_1s)))
}

solve_event_size <- function(id, alpha_1s) {
  target <- 1 - gsd1[[id]]$beta
  power_at <- function(max_events)
    tail(power_for(id, alpha_1s, max_events)$power, 1)
  lower <- if (id == "H3") max(10, max(coc$cut$at) + 1) else 10
  upper <- max(100, 2 * lower)
  while (power_at(upper) < target) {
    upper <- 2 * upper
    if (upper > 1e7) stop("Could not bracket the required event count")
  }
  required_events <- ceiling(uniroot(
    function(max_events) power_at(max_events) - target,
    interval = c(lower, upper), tol = 0.01)$root)
  list(
    required_events = required_events,
    result = power_for(id, alpha_1s, required_events))
}

# ==== mtp/gsd: reachable states and unique calculations ======================

# Build a complete state catalogue and cache one calculation for every unique
# hypothesis-by-alpha combination. The state catalogue retains the applicable
# local-test context without duplicating identical univariate GSD tables.

local_test_for <- function(avail) {
  is_match <- vapply(mtp1$intersection_tests, function(x)
    setequal(x$hypotheses, avail), logical(1))
  if (!any(is_match)) return("singleton")
  test_id <- names(mtp1$intersection_tests)[which(is_match)[1]]
  sprintf("%s (%s)", mtp1$intersection_tests[[test_id]]$type, test_id)
}

state_records <- do.call(rbind, lapply(seq_along(reachable), function(i) {
  r <- reachable[[i]]
  a_1s_state <- as.numeric(local_alpha[r$row, r$avail])
  graph_alpha <- sides[r$avail] * a_1s_state
  data.frame(
    State = paste0("S", i),
    Available = paste(r$avail, collapse = ", "),
    Hypothesis = r$avail,
    GraphAlpha = graph_alpha,
    EffectiveAlpha = graph_alpha,
    Alpha1s = a_1s_state,
    LocalTest = local_test_for(r$avail),
    Status = ifelse(a_1s_state > alpha_tolerance,
                    "testable", "not currently testable"),
    stringsAsFactors = FALSE)
}))

state_records$Calculation <- "-"
configuration_keys <- character()
configurations <- list()
for (i in seq_len(nrow(state_records))) {
  if (state_records$Alpha1s[i] <= alpha_tolerance) next
  key <- paste(
    state_records$Hypothesis[i],
    format(state_records$Alpha1s[i], digits = 16, scientific = TRUE), sep = "|")
  config_index <- match(key, configuration_keys)
  if (is.na(config_index)) {
    configuration_keys <- c(configuration_keys, key)
    config_index <- length(configuration_keys)
    calc_id <- paste0("C", config_index)
    configurations[[calc_id]] <- list(
      id = state_records$Hypothesis[i],
      alpha_1s = state_records$Alpha1s[i],
      graph_alpha = state_records$GraphAlpha[i],
      effective_alpha = state_records$EffectiveAlpha[i])
  }
  state_records$Calculation[i] <- paste0("C", config_index)
}

for (calc_id in names(configurations)) {
  configurations[[calc_id]]$states <- unique(
    state_records$State[state_records$Calculation == calc_id])
  configurations[[calc_id]]$local_tests <- unique(
    state_records$LocalTest[state_records$Calculation == calc_id])
}
calculation_results <- lapply(configurations, function(x)
  solve_event_size(x$id, x$alpha_1s))

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

# Print one cached hypothesis-by-alpha calculation. State reuse is reported in
# the header instead of repeating identical tables under multiple states.
print_hypothesis <- function(calc_id, config, sizing) {
  id <- config$id
  alpha_1s <- config$alpha_1s
  x <- sizing$result
  nd <- mtp1$nodes[[id]]
  g <- gsd1[[id]]
  looks <- if (length(x$events) == 1) "FA" else
    c(paste0("IA", seq_len(length(x$events) - 1)), "FA")

  cat(sprintf("\ncalculation %s: %s\n", calc_id, labels[id]))
  cat(sprintf("  used by: %s | graph alpha %s | effective nominal alpha %s | target power %.2f%% | required final events %s | local test context(s): %s\n",
              paste(config$states, collapse = ", "),
              format_alpha(config$graph_alpha),
              format_alpha(config$effective_alpha),
              100 * (1 - g$beta), sizing$required_events,
              paste(config$local_tests, collapse = "; ")))
  cat(sprintf("  %s-sided, %s vs %s; looks by %s; %s approximation; futility: %s; calculation: sample size [events]\n",
              nd$test_sides,
              nd$intervention, nd$comparator, g$looks$trigger, g$approximation,
              if (has_futility(g$futility)) "configured" else "none"))

  # reference distribution, already looked up and gated by ref_curves above
  ref <- ref_curves[[id]]
  if (!is.null(ref))
    cat(sprintf("    reference (%s) %s: %s, median %s by cuts %s months\n",
                nd$comparator, nd$endpoint, ref$family,
                paste(ref$median, collapse = "/"), paste(ref$cut$at, collapse = "/")))

  boundary <- data.frame(
    Analysis       = looks,
    Boundary       = "efficacy",
    Z              = round(x$z, 3),
    CV_HR          = round(x$cv_hr, 3),
    `Nominal α`    = round(nd$test_sides * x$alpha_nom, 5),
    `Cumulative α` = round(nd$test_sides * x$spend_cum, 5),
    check.names = FALSE)
  cat("    -- type I error control --\n")
  print(boundary, row.names = FALSE)
  cat("    -- type II error / futility control --\n")
  print(beta_futility_summary(gsd1[[id]], looks, "sample-size"), row.names = FALSE)
  analyses <- data.frame(
    Analysis                = looks,
    `Information Fraction`  = round(x$events / max(x$events), 3),
    `N Events`               = round(x$events, 1),
    `Marginal Power`         = round(x$power_marginal, 4),
    `Cumulative Power`       = round(x$power, 4),
    check.names = FALSE)
  cat("    -- planned analyses --\n")
  print(analyses, row.names = FALSE)
}

cat(sprintf("\n================ %s ================\n", trial1$label))
ratio_text <- paste(sprintf("%s=%s:1", ids,
                            vapply(ids, allocation_ratio_for, numeric(1))),
                    collapse = ", ")
cat(sprintf("  allocation ratios (intervention:comparator) %s | stratification: %s | simulation: %s\n",
            ratio_text, if (isTRUE(trial1$strata$enabled)) "yes" else "none",
            if (isTRUE(trial1$sim$enabled)) "yes" else "none (closed form)"))

cat("\n-- reachable-state catalogue --\n")
state_table <- data.frame(
  State = state_records$State,
  Available = state_records$Available,
  Hypothesis = state_records$Hypothesis,
  `Graph alpha` = vapply(state_records$GraphAlpha, format_alpha, character(1)),
  `Effective alpha` = vapply(state_records$EffectiveAlpha, format_alpha, character(1)),
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
  print_hypothesis(calc_id, configurations[[calc_id]],
                   calculation_results[[calc_id]])
}
cat("\n")
