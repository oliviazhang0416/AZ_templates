# Standardized JSON Design-Input Data Dictionary

## Purpose and scope

This document defines the common JSON input schema for the `trial`, `mtp`, and `gsd` objects used by the reviewed R programs. After import, these are assigned to the R objects `trial1`, `mtp1`, and `gsd1`. The schema covers analytic and simulation-based power, sample-size and timeline calculations, time-to-event, continuous and binary endpoints, graphical and non-graphical multiplicity procedures, group-sequential efficacy monitoring, and formal or operational futility.

The schema is a union of supported fields. Unless a section explicitly defines
a field as conditional or omitted, every standardized field appears in its
corresponding JSON object. A field that is not applicable is set to JSON
`null`; an inactive feature with an explicit status field uses
`"enabled": false`.

## Conventions

- All calendar and follow-up times are in months. The schema intentionally has no `time_unit` field.
- In standardized filenames, `diff` identifies designs with difference-based effect measures (including mean difference and risk difference), and `fd` means a fixed design without group-sequential monitoring.
- Named IDs such as `arm1`, `PFS`, and `H1` are case-sensitive and must be unique within their JSON object.
- Canonical categorical values use lower-case, hyphenated strings unless an established method or R function name is case-sensitive, such as `"Schoenfeld"` or `"sfLDOF"`.
- `randomization_weight` replaces `allocation_weight`. It describes relative randomized allocation by arm. A comparison ratio is derived as intervention weight divided by comparator weight.
- JSON `null` means not applicable or not used by the selected calculation. It must not be interpreted as zero.
- Index arrays such as `active_at` use one-based analysis numbers.
- Spending-function names match the exported `gsDesign` names exactly. Spending parameters are stored inside the corresponding `family` object.
- Custom `family.nominal_values` are nominal levels on the hypothesis's `test_sides` scale. In efficacy they are alpha levels; in futility they are beta levels. A final array element of `null` means use the remaining available alpha or beta; `jsonlite` imports that array element as R `NA`.

## File layout and loading

The three standardized objects are properties of one JSON parameter file. Each
calculation script has an adjacent JSON file with exactly the same filename stem;
only the extension changes from `.R` to `.json`. For example:

```text
fixedseq-gsd-hr-gsdesign-3arm-PFS-OS-power-[AZ].json
```

Each parameter file represents one complete design configuration. A design with
operational futility and the corresponding design without futility use separate
same-stem R/JSON pairs; they are not combined as two power rules in one input.

Its top-level structure is:

```json
{
  "trial": {},
  "mtp": {},
  "gsd": {}
}
```

The R script and JSON file are stored together in the same folder. The reviewed
files currently use `/Users/oliviazhang/Desktop/AZ_templates` or
`/Users/oliviazhang/Desktop/Attachments`. A script reads its own same-stem JSON
file and retains the standardized R object names. For example:

```r
design_input <- jsonlite::fromJSON(
  "/Users/oliviazhang/Desktop/AZ_templates/fixedseq-gsd-hr-gsdesign-3arm-PFS-OS-power-[AZ].json",
  simplifyVector = TRUE,
  simplifyDataFrame = FALSE,
  simplifyMatrix = TRUE)
trial1 <- design_input$trial
mtp1   <- design_input$mtp
gsd1   <- design_input$gsd
```

The reader uses `simplifyVector = TRUE`, `simplifyDataFrame = FALSE`, and `simplifyMatrix = TRUE`. Consequently, JSON arrays become R vectors or matrices where appropriate, while object arrays such as `edges` remain lists. JSON contains no unevaluated R expressions. Store the original statistical assumption whenever the schema can derive a secondary quantity: for example, retain an exponential median and set its derived hazard rate to `null`, or retain landmark survival probabilities and derive interval hazards in R. Store an evaluated numeric result only when the source assumptions needed to reconstruct it are not represented by the schema.

## Power and sample-size counterpart convention

The calculation mode is identified by `power` or `samplesize` in the same-stem
R/JSON filename; no additional mode or objective field is required.

- A time-to-event power file fixes `gsd.<hypothesis_id>.max_events` and sets
  `beta` to `null`; achieved power is an output. For an event-count or
  information-fraction schedule, `max_events` is the planned final event
  target. For a calendar-time schedule, it is the planned expected cumulative
  events by the final calendar look; R scales the reference accrual shape to
  match that expectation and then calculates power. The only beta exception is
  beta-spending futility, where beta calibrates the prespecified lower boundary.
  An R engine that merely requires beta to instantiate an efficacy-only design
  must use an internal placeholder rather than store it in the JSON.
- A sample-size file requires a non-`null` hypothesis-specific `beta`; target
  power is `1 - beta`. For the time-to-event designs standardized here, the
  sizing estimand is the required final event count. `max_events` is `null` in
  the input and is solved in R.
- For event-driven sample-size calculations, `looks.trigger` is
  `"event-count"`, `looks.schedule.scale` is
  `"fraction-of-final-events"`, and the schedule values end at 1. The required
  final event count and the corresponding cumulative event counts are solved in
  R. These event fractions are not labelled information fractions.
- For an NPH sample-size calculation with calendar-time looks, the analysis
  times and reference accrual shape may remain inputs because they determine
  the AHR and information trajectory. The engine may derive an implied
  participant count internally, but the standardized sizing result is the
  required final event count; participant count is not reported as the solved
  design size.
## 1. `trial`: trial and data-generating assumptions

### Top-level fields

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `label` | string | Yes | Nonempty text | Human-readable design name. |
| `alpha` | object | Yes | See below | Overall type I error convention. |
| `arms` | object | Yes | At least two unique arms | Trial arms and randomized allocation. |
| `endpoints` | object | Yes | At least one endpoint | Endpoint definitions, assumptions, and tests. |
| `enrollment` | object | Yes | Fields may all be `null` | Enrollment model. |
| `dropout` | object | Yes | Fields may all be `null` | Independent dropout model. |
| `strata` | object | Yes | `"enabled": true` or `"enabled": false` | Stratification assumptions. |
| `sim` | object | Yes | `"enabled": true` or `"enabled": false` | Simulation controls. |

### `alpha`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `value` | number | Yes | `0 < value < 1` | Overall alpha on the declared sidedness scale. |
| `sidedness` | integer | Yes | `1` or `2` | One- or two-sided alpha convention. |

### Arm object: `trial.arms.<arm_id>`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `role` | string | Yes | `"control"`, `"experimental"`, `"component"` | Arm role; `"component"` identifies a component-treatment arm used in a component-of-combination comparison. |
| `component_of` | string or `null` | Conditional | Existing arm ID other than the current arm, or `null` | Combination arm containing this component. Required when `role = "component"`; otherwise `null`. |
| `randomization_weight` | number or `null` | Conditional | `> 0` or `null` | Relative randomized allocation weight. Required for randomized arms; may be `null` for a nonrandomized component arm. |

Example: weights `1, 1, 1` mean 1:1:1 randomization; weights `2, 1, 1` mean 2:1:1. A nonrandomized component arm uses `randomization_weight: null`; when an information ratio is needed for its assumed effect, R uses the randomized weight of the arm identified by `component_of`.

### Endpoint object: `trial.endpoints.<endpoint_id>`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `type` | string | Yes | `"time-to-event"`, `"continuous"`, `"binary"` | Endpoint data type. |
| `effect_measure` | string | Yes | Time-to-event: `"hazard-ratio"`; continuous: `"mean-difference"`; binary: `"risk-difference"`, `"risk-ratio"`, `"odds-ratio"` | Statistical effect scale. |
| `benefit_direction` | string | Yes | `"lower"`, `"higher"` | Direction favorable to the intervention. |
| `test` | object | Yes | See below | Primary analysis method. |
| `distributions` | object | Yes | One entry per arm; entries may be `null` | Absolute outcome assumptions by arm. |
| `effects` | object | Yes | One entry per arm; entries may be `null` | Relative treatment-effect assumptions by arm. |

The JSON preserves the parameterization used by the source R program. For any
intervention-versus-comparator comparison, at most two of the following are
specified: comparator distribution, intervention distribution, and
intervention effect. If the source supplies both arm distributions, the
intervention effect is `null`. If the source supplies the comparator
distribution and a treatment effect, the intervention distribution is `null`
and R derives it. A design may contain only an effect when absolute
distributions are unnecessary. A calculated third representation is not stored
as an additional input.

### Endpoint `test`

| Field | JSON type | Required | Allowed values | Meaning |
|---|---|---:|---|---|
| `method` | string | Yes | `"logrank"`, `"cox"`, `"t-test"`, `"chi-square"`, `"fisher-exact"` | Statistical test or model. |
| `effect_estimator` | string or `null` | Conditional | `"coxph"`, `null` | Estimator used for a reported observed treatment effect. For a time-to-event simulation, `"coxph"` requests a Cox-model hazard-ratio estimate while `method = "logrank"` may still provide the hypothesis-test p-value. |
| `alternative` | string | Yes | `"one-sided"`, `"two-sided"` | Test alternative. |
| `stratified` | boolean or `null` | Conditional | `true`, `false`, `null` | Whether the analysis uses the specified strata. |
| `equal_variance` | boolean or `null` | Conditional | `true`, `false`, `null` | Equal-variance assumption for a t-test. |
| `continuity_correction` | boolean or `null` | Conditional | `true`, `false`, `null` | Continuity correction for a chi-square test. |

### Distribution object: `trial.endpoints.<endpoint_id>.distributions.<arm_id>`

An entry is `null` when an absolute distribution is not specified directly. This can occur because absolute distributions are unnecessary, as in a fixed-event Schoenfeld power calculation, or because the arm distribution is derived by following that endpoint's effect chain to an arm with a non-`null` absolute distribution. When design assumptions differ by a declared stratum, the arm entry may instead be an object keyed directly by every ID in `trial.strata.definitions`; each value is a complete distribution object. A single ordinary distribution object applies across all strata.

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `family` | string | Yes when non-`null` | `"exponential"`, `"piecewise-exponential"`, `"landmark-survival"`, `"normal"`, `"bernoulli"` | Distribution family. |
| `cut` | object or `null` | Conditional | `"by": "follow-up-time"`; increasing `at` values `>= 0` | Piecewise interval definition. |
| `rate` | array of numbers or `null` | Conditional | Each value `>= 0` | Exponential hazards; length equals number of intervals. |
| `median` | number, array of numbers, or `null` | Conditional | Each value `> 0` | Median event time for an exponential model; an array supplies one median per piecewise interval. |
| `landmark_time` | array of numbers or `null` | Conditional | Strictly increasing and `>= 0` | Times for landmark survival assumptions. |
| `landmark_value` | array of numbers or `null` | Conditional | Each value in `[0, 1]`, nonincreasing | Survival probabilities corresponding to `landmark_time`. |
| `mean` | number or `null` | Conditional | Any finite value | Mean for a normal outcome. |
| `sd` | number or `null` | Conditional | `> 0` | Standard deviation for a normal outcome. |
| `probability` | number or `null` | Conditional | `[0, 1]` | Event/response probability for a Bernoulli outcome. |

Only parameters appropriate to the selected `family` are non-`null`. Store source-level assumptions rather than their calculated equivalents:

- `"exponential"`: exactly one of `median` and `rate` is non-`null`. Prefer `median` when the source specifies a median; R derives the rate as `log(2) / median`.
- `"piecewise-exponential"`: `cut` is non-`null` and exactly one of `rate` and `median` is non-`null`. The selected array has one value per interval. Prefer interval-specific medians when the source states medians; R derives each rate as `log(2) / median`.
- `"landmark-survival"`: `landmark_time` and `landmark_value` are non-`null`; `cut`, `rate`, and `median` are `null`. R may derive piecewise hazards when required by the calculation engine.
- `"normal"`: `mean` and `sd` are non-`null`; the survival and Bernoulli fields are `null`.
- `"bernoulli"`: `probability` is non-`null`; the survival and normal fields are `null`.

### Effect object: `trial.endpoints.<endpoint_id>.effects.<arm_id>`

As with distributions, an arm effect may be `null`, one common effect object,
or an object keyed directly by every ID in `trial.strata.definitions` when the
design effect differs by stratum.

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `vs` | string | Yes when non-`null` | Existing arm ID | Comparator arm. |
| `cut` | object or `null` | Conditional | `null`, or `"by": "follow-up-time"` or `"event-count"`; `"at": null` for a constant time-to-event effect, otherwise increasing values beginning at 0 | Piecewise-effect definition. Use `null` for a constant non-time-to-event effect. |
| `value` | number or array of numbers | Yes | Depends on `effect_measure` | Effect in each interval. A number when `cut.at` is `null`; otherwise array length equals the length of `cut.at`. |

Effect-value ranges:

- Hazard ratio, risk ratio, and odds ratio: `> 0`.
- Risk difference: `[-1, 1]`.
- Mean difference: any finite real number.

An effect chain is formed by following `effects.<arm_id>.vs`. For example, if arm C has an explicit distribution, arm B has an effect versus arm C, and arm A has an effect versus arm B, then the arm B and arm A distributions can be derived without storing calculated hazards. The calculation applies each link on the endpoint's `effect_measure` scale and uses the union of its piecewise cutpoints. Effect chains must not contain cycles.

Conversely, when both comparator and intervention distributions are supplied
directly, any effect between them is derived when needed and its JSON effect
entry remains `null`. This rule applies separately by stratum when the source
uses stratum-specific distributions or effects.

### `enrollment`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `model` | string or `null` | Yes | `null`, `"power"`, `"piecewise-rate"`, `"cumulative-count"` | Enrollment representation. |
| `n_accrual_max` | integer or `null` | Conditional | Positive integer | Maximum randomized enrollment. |
| `period` | number or `null` | Conditional | `> 0` | Total accrual duration in months. |
| `k` | number or `null` | Conditional | `> 0` | Shape parameter for the power accrual model. |
| `duration` | array of numbers with an optional final `null`, or `null` | Conditional | Numeric values `> 0`; only the final element may be `null` | Piecewise interval durations in months. A final `null` means the last rate continues until `n_accrual_max` is reached. |
| `rate` | array of numbers or `null` | Conditional | Each value `>= 0`; same length as `duration` | Enrollment rate in each interval. |
| `cumulative` | array of integers or `null` | Conditional | Nonnegative and nondecreasing | Cumulative randomized counts by month. |

Requirements by model:

- `"power"`: `n_accrual_max`, `period`, and `k` are required.
- `"piecewise-rate"`: `duration` and `rate` are required and have the same length; `n_accrual_max` may cap enrollment. When the final `duration` is `null`, `n_accrual_max` is required and the final rate continues until that cap is reached. With `simplifyVector = TRUE`, the final JSON `null` imports as R `NA`.
- `"cumulative-count"`: `cumulative` is required.

### `dropout`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `model` | string or `null` | Yes | `null`, `"exponential"` | Independent censoring model. |
| `probability` | number or `null` | Conditional | `[0, 1)` | Cumulative dropout probability by `window`. |
| `window` | number or `null` | Conditional | `> 0` | Dropout-probability window in months. |
| `applies_to` | string, array, object, or `null` | Conditional | `"all"` or valid arm and endpoint IDs | Scope of the dropout assumption. |

For `model = "exponential"`, both `probability` and `window` are required. The monthly censoring hazard is `-log(1 - probability) / window`.

### `strata`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `enabled` | boolean | Yes | `true`, `false` | Whether strata are part of generation or analysis. |
| `definitions` | object or `null` | Conditional | Named stratum definitions; see below | Meaning and population proportion of each stratum. |

When `enabled = true`, each property of `definitions` is a user-defined
stratum ID. Each definition contains:

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `values` | object | Yes | At least one named factor and nonempty level | Factor values defining the stratum. |
| `proportion` | number | Yes | `[0,1]`; proportions across definitions sum to 1 | Population proportion sampled from the stratum. |

For example, `"T2D-negative"` and `"T2D-positive"` may be used as the
definition IDs. A stratum can affect data generation even when the endpoint's
`test.stratified` field is `false`. The reviewed simulation engine implements
one convention, so it is not exposed as an input: arm-specific stratum counts
are fixed from the declared proportions. R rounds `n * proportion` for each
stratum other than the largest stratum, and the largest stratum receives the
remaining participants so that counts sum exactly to the arm size. When
`enabled = false`, `definitions` is `null`.

### `sim`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `enabled` | boolean | Yes | `true`, `false` | Whether simulation is used. |
| `n_sim` | integer or `null` | Conditional | Positive integer | Number of simulated trials. |
| `seed` | integer or `null` | Conditional | Valid JSON integer | Reproducibility seed. |
| `param_grid` | array of objects or `null` | Conditional | See below | One or more parameter grids at which operating characteristics are evaluated. |
| `outcome_models` | object or `null` | Conditional | Properties must match endpoint IDs in `trial.endpoints` | Source-parameterized distributions and effects used to generate simulated outcomes. |
| `empirical_quantities` | array of objects or `null` | Conditional | See below | Requested empirical summaries and operating characteristics. |

When `enabled = true`, `n_sim`, `seed`, `outcome_models`, and
`empirical_quantities` are required. When `enabled = false`, every other member
is `null`. The reviewed simulation programs generate participant outcomes from
the arm-specific outcome models, enrollment, and dropout inputs; a separate
simulation-unit field is therefore unnecessary.

#### Simulation parameter-grid object: `trial.sim.param_grid[i]`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `parameter` | string | Yes | `"calendar-time"`, `"event-count"`, `"enrollment-per-arm"` | Parameter varied over the grid. |
| `values` | array of numbers | Yes | Nonempty; positive and strictly increasing | Values at which the simulation result is evaluated. |

Each grid entry produces a separate series of evaluations. Multiple entries
are not interpreted as a Cartesian product. Calendar times are in months.
Event counts and enrollment-per-arm values must be integer-valued. A fixed
planned analysis schedule belongs in `gsd.<hypothesis_id>.looks`; `param_grid`
is used only when the program deliberately evaluates multiple candidate values
of a parameter.

#### Simulation outcome-model object: `trial.sim.outcome_models.<endpoint_id>`

| Field | JSON type | Required | Allowed values | Meaning |
|---|---|---:|---|---|
| `distributions` | object | Yes | Arm-specific object | Absolute distributions used for simulation; intervention entries may be `null` when derived from effects. |
| `effects` | object | Yes | Arm-specific object | Treatment effects used for simulation; entries are `null` when both arm distributions are supplied directly. |

The endpoint ID must exist in `trial.endpoints`. Each arm entry is one of:

1. `null`;
2. one ordinary distribution or effect object, which applies across all
   strata; or
3. an object keyed directly by IDs in `trial.strata.definitions`, with one
   complete distribution or effect object for every declared stratum.

There is no `by_stratum` wrapper. A common arm distribution may be stated once
instead of repeated for every stratum. When any assumption differs by stratum,
the arm entry must enumerate every defined stratum. Effect chains and parameter
ranges follow the rules for the corresponding design endpoint objects.

The design and simulation locations have distinct meanings:

- `trial.endpoints.<endpoint_id>.distributions/effects` contains the planning
  or design assumptions.
- `trial.sim.outcome_models.<endpoint_id>.distributions/effects` contains the
  complete assumptions used to generate simulated outcomes.

The design endpoint should contain only assumptions consumed by the design
calculation. For example, a fixed-event design may require a treatment effect
but no absolute event-time distribution, in which case all design distribution
entries are `null`. The effective simulation model must be sufficient to
generate the data. Its directly supplied blocks must provide either an absolute
distribution for every generated arm or a comparator distribution plus an
effect chain from which the remaining arm distributions can be derived. The
simulation model preserves whichever of these two parameterizations appears in
the source R program and does not store the derived third representation.
Derived endpoints are an exception when their treatment effect is induced by
the source observation and therefore is not an independent simulation input.
Directly generated outcomes are independent across participants, arms, and
endpoint models unless a dependency is established explicitly through
`derived_from`. Thus the reviewed body-weight and waist outcomes are generated
independently, while each derived responder outcome is calculated from the
same participant's body-weight observation. A future design requiring another
joint dependence structure would require an explicit schema extension.

#### Derived-Bernoulli simulation distribution

A binary observation deterministically constructed from another simulated
endpoint uses this simulation-only distribution form:

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `family` | string | Yes | `"derived-bernoulli"` | Binary value obtained by thresholding a simulated source endpoint. |
| `derived_from` | string | Yes | Existing endpoint ID | Source endpoint measured on the same participant and in the same arm. |
| `threshold_rule` | object | Yes | See below | Rule converting the source value to zero or one. |

The `threshold_rule` contains `operator`, one of `"<"`, `"<="`, `">"`, or
`">="`, and a finite numeric `value`. For example, body-weight percent change
below `-0.05` defines a five-percent responder. The source arm is implicitly
the same arm as the derived endpoint entry. The outcome must be calculated
from the already simulated source observation; it must not be generated as an
independent Bernoulli draw.

#### Empirical-quantity object: `trial.sim.empirical_quantities[i]`

| Field | JSON type | Required | Allowed values | Meaning |
|---|---|---:|---|---|
| `label` | string | Yes | Unique, nonempty text | User-defined label for the output column or row. |
| `grid` | string or `null` | Conditional | A `parameter` named in `trial.sim.param_grid`, or `null` | Parameter grid on which the quantity is evaluated. Required only when multiple grids exist and the quantity applies to a specific one. |
| `summary` | object | Yes | See below | Empirical summary function. |

The `summary` object contains:

| Field | JSON type | Required | Allowed values | Meaning |
|---|---|---:|---|---|
| `fn` | string | Yes | `"mean"`, `"median"`, `"quantile"`, `"probability"` | Function applied across simulation replicates. |
| `expression` | object | Yes | Quantity, comparison, Boolean, or internal-reference expression | Value or condition evaluated once within each simulated trial. |
| `probabilities` | array of numbers or `null` | Conditional | Strictly increasing values in `(0,1)`, or `null` | Quantile probabilities when `fn = "quantile"`; otherwise `null`. |

The calculation has one consistent order: evaluate `summary.expression`
within every simulated trial, then apply `summary.fn` across trials. There is
no separate `source`, `type`, `condition`, or rule field.

##### Quantity expression

A quantity expression identifies the hypothesis and replicate-level quantity:

```json
{
  "hypothesis_id": "H_PFS",
  "quantity": "event-count"
}
```

| Field | JSON type | Required | Allowed values | Meaning |
|---|---|---:|---|---|
| `hypothesis_id` | string | Yes | Existing `mtp.nodes` ID | Hypothesis providing the simulated value. |
| `quantity` | string | Yes | `"observed-effect-measure"`, `"p-value"`, `"test-statistic"`, `"event-count"`, `"calendar-time"`, `"participant-count"`, `"dropout-proportion"` | Replicate-level quantity evaluated. |

For `fn = "mean"` or `"median"`, the expression must return one numeric value
per trial. For `fn = "quantile"`, it returns a numeric value and
`probabilities` supplies the requested quantiles.

##### Fixed-threshold expression

For empirical probabilities at one or more fixed thresholds, the quantity
expression also contains `direction` and `thresholds`:

```json
{
  "hypothesis_id": "H_PFS",
  "quantity": "observed-effect-measure",
  "direction": ">",
  "thresholds": [0.1, 0.2]
}
```

`direction` is one of `"<"`, `"<="`, `">"`, or `">="`. `thresholds` is a
nonempty numeric array. With `fn = "probability"`, one empirical probability
is reported for each threshold.

##### Atomic decision expression

An atomic decision used within a marginal or joint power definition contains
a single `threshold`:

```json
{
  "hypothesis_id": "H_high_body_weight",
  "quantity": "p-value",
  "direction": "<=",
  "threshold": "current-local-alpha"
}
```

`threshold` is either a finite number or `"current-local-alpha"`. The latter
means the hypothesis's local alpha in the current reachable MTP state after
all applicable alpha transfers; it must not be replaced by the node's initial
alpha. For a fixed nominal test, use the numeric alpha directly.

##### Boolean expression

Nested `all_of` and `any_of` objects combine atomic decision expressions.
`all_of` requires every operand to be true in the same simulated trial;
`any_of` requires at least one. For example:

```json
{
  "all_of": [
    {
      "any_of": [
        {
          "hypothesis_id": "H_high_body_weight",
          "quantity": "p-value",
          "direction": "<=",
          "threshold": "current-local-alpha"
        },
        {
          "hypothesis_id": "H_high_5pct",
          "quantity": "p-value",
          "direction": "<=",
          "threshold": "current-local-alpha"
        }
      ]
    },
    {
      "hypothesis_id": "H_high_10pct",
      "quantity": "p-value",
      "direction": "<=",
      "threshold": "current-local-alpha"
    }
  ]
}
```

These expressions define reported marginal or joint operating characteristics;
they do not change testing, gating, rejection, or alpha transfer in `mtp`.

##### Internal-reference expression

An expression may be an internal `$ref` to a GSD rule:

```json
{
  "$ref": "#/gsd/H_PFS/futility"
}
```

A reference to `efficacy` evaluates whether its efficacy boundary is crossed;
a reference to `futility` evaluates whether its futility rule is met; and a
reference to the complete hypothesis-specific GSD object evaluates rejection
under both its efficacy and operational futility specifications. The referenced
object supplies the boundary construction and `active_at` looks, so the
empirical item does not repeat look or boundary values.

The empirical `grid` field is omitted when there is no parameter grid, when
exactly one grid exists, or when the quantity is evaluated on every grid. When
multiple grids exist and the quantity applies to only one, `grid` names that
grid's `parameter`.

## 2. `mtp`: hypotheses and multiplicity

### Top-level fields

| Field | JSON type | Required | Allowed values | Meaning |
|---|---|---:|---|---|
| `procedure` | object | Yes | See below | Multiplicity procedure. |
| `intersection_tests` | object or `null` | Yes | See below | Local tests indexed by intersection-hypothesis ID; `null` when there is no multiplicity. |
| `nodes` | object | Yes | At least one hypothesis | Hypothesis definitions. |
| `edges` | array of objects | Yes | May be empty | Alpha-transfer relationships. |

### `procedure`

| Field | JSON type | Required | Allowed values | Meaning |
|---|---|---:|---|---|
| `type` | string | Yes | `"none"`, `"single-step"`, `"graphical"`, `"holm"`, `"hochberg"`, `"gatekeeping"`, `"closed-testing"` | Multiplicity procedure class. |
| `strategy` | string or `null` | Conditional | `null`, `"fixed-sequence"`, `"parallel"`, `"parallel-hierarchies"`, `"alternating"`, `"hybrid"` | Procedure-specific organization. |
| `epsilon` | number or `null` | Conditional | `null` or `0 < epsilon < 0.5` | Small transfer weight used to encode an epsilon-edge AND gate. Required only when that construction is used. |

### Intersection-test object: `mtp.intersection_tests.<intersection_id>`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `hypotheses` | array of strings | Yes | At least two unique existing node IDs | Elementary hypotheses comprising the intersection null. |
| `type` | string | Yes | `"weighted-bonferroni"`, `"bonferroni"`, `"dunnett"` | Local test applied to this intersection hypothesis. |
| `correlation` | array of number arrays or `null` | Conditional | Square dimension equal to the number of `hypotheses`; symmetric, diagonal 1, entries in `[-1,1]`, positive semidefinite | Correlation matrix in the same order as `hypotheses`. Required for `"dunnett"`; otherwise `null`. |

The JSON object key is a readable intersection ID, such as `H123`, while `hypotheses` is the authoritative membership definition. This avoids having to parse node IDs from the intersection ID.

When `procedure.type = "graphical"`, the default local intersection test is `"weighted-bonferroni"`. Its weights are derived from the node alpha allocations and graph edges rather than repeated inside each intersection-test object.

For a closed procedure with three elementary hypotheses, the non-singleton intersections are `H12`, `H13`, `H23`, and `H123`. Each must have an entry unless the selected procedure supplies an equivalent shortcut that derives its local tests automatically.

Example for a three-hypothesis Dunnett intersection test:

```json
{
  "intersection_tests": {
    "H123": {
      "hypotheses": ["H1", "H2", "H3"],
      "type": "dunnett",
      "correlation": [
        [1.0, 0.5, 0.5],
        [0.5, 1.0, 0.5],
        [0.5, 0.5, 1.0]
      ]
    }
  }
}
```

### Node object: `mtp.nodes.<hypothesis_id>`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `label` | string | Yes | Nonempty text | Human-readable hypothesis label. |
| `intervention` | string | Yes | Existing arm ID | Intervention arm. |
| `comparator` | string | Yes | Existing arm ID and different from intervention | Comparator arm. |
| `endpoint` | string | Yes | Existing endpoint ID | Endpoint tested. |
| `initial_alpha` | number | Yes | `[0, trial.alpha.value]` | Alpha allocated before any rejection/recycling, on the declared alpha sidedness scale. |
| `test_sides` | integer | Yes | `1` or `2` | Sidedness for this hypothesis. |

For a confirmatory graph, the sum of initial node alpha must equal the available family alpha. A node may start at zero and receive alpha later. For a single-step Dunnett procedure, `initial_alpha` records the pre-adjustment allocation or weight on the declared family-alpha scale; it is not the Dunnett-adjusted marginal nominal level. The latter is derived from the family alpha and the correlation matrix in the applicable intersection test.

### Edge object: `mtp.edges[i]`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `from` | string | Yes | Existing node ID | Source hypothesis. |
| `to` | string | Yes | Existing node ID, different from `from` | Destination hypothesis. |
| `weight` | number | Yes | `[0,1]` | Fraction of released alpha transferred. |
| `on` | string | Yes | `"reject"` | Transfer event. |
| `condition` | string, object, or `null` | Conditional | Prespecified rejection condition or `null` | Additional condition for complex co-primary logic. |

Outgoing weights from a node must sum to no more than 1.

When an edge needs an additional condition, `condition` uses nested
`all_of`/`any_of` logic. In this location the condition changes alpha transfer
and therefore is part of the MTP, not merely an empirical output expression.
Hypothesis IDs must match `mtp.nodes` exactly.

For an epsilon-edge AND gate with two upstream nodes followed by an ordered pair
of downstream nodes, each upstream node transfers `epsilon` to the first
downstream node and `1 - epsilon` to the other upstream node. A unit edge from
the first downstream node to the second completes the sequence. Thus only a
negligible amount is available downstream after one upstream rejection;
rejecting both upstream hypotheses transfers the full available alpha to the
first downstream node, and its rejection transfers that alpha to the second.

## 3. `gsd`: hypothesis-specific sequential design

Each property is keyed by an ID in `mtp.nodes`.
When a design has no group-sequential monitoring, `gsd` is an empty object;
fixed final-analysis testing is then defined by the endpoint tests, MTP nodes,
and MTP procedure.

### GSD object: `gsd.<hypothesis_id>`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `looks` | object | Yes | See below | Planned analysis schedule. |
| `max_events` | number or `null` | Conditional | `> 0` | Final event target. For event-count or information-fraction power calculations it is the planned final event count. For calendar-time power calculations it is the planned expected cumulative events by the final look and is used to scale the reference accrual shape. It is `null` in sample-size files and solved in R. |
| `efficacy` | object | Yes | See below | Upper/efficacy boundary. |
| `beta` | number or `null` | Yes | `null` or `0 < beta < 1` | Type II error target for a sample-size solution; `null` for a power calculation except when required to calibrate beta-spending futility. Target power is `1 - beta`. |
| `futility` | object | Yes | See below | Formal lower boundary or operational stopping rule. |
| `approximation` | string or `null` | Yes | `null`, `"Schoenfeld"`, `"exact"`, `"simulation"` | Information/power approximation. |
| `engine` | string or `null` | Yes | `null`, `"gsDesign"`, `"gsDesign2"`, `"rpact"`, `"simulation"` | Calculation engine. |

### `looks`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `trigger` | string | Yes | `"calendar-time"`, `"event-count"`, `"information-fraction"` | Quantity that operationally triggers an analysis. |
| `schedule` | object | Yes | See below | One tagged representation of the planned look values. |

### `looks.schedule`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `scale` | string | Yes | `"calendar-time"`, `"count"`, `"fraction-of-final-events"`, `"information-fraction"` | Interpretation of `values`. |
| `values` | array of numbers | Yes | Strictly increasing; scale-specific restrictions below | Planned look values on the declared scale. |

Conditional rules:

- `trigger = "calendar-time"` requires `schedule.scale = "calendar-time"`.
  Values are positive months. Enrollment shape and, when applicable,
  dropout/distribution assumptions are required to obtain events and
  information. In a patient-level simulation, R performs each analysis at its
  stated calendar time; event counts and information fractions are simulation
  results. In an analytic power file, non-`null` `max_events` specifies the
  expected final event target and `n_accrual_max` is `null`; R derives the
  implied enrollment scale. In an analytic sample-size file, both are `null`,
  and beta determines the scale and required expected final events.
- `trigger = "event-count"` requires `schedule.scale = "count"` or
  `"fraction-of-final-events"`. Exact counts are positive integers; in a power
  file their final value equals `max_events`. Fractions lie in `(0,1]`, end at
  1, and are used when the final event count is being solved. After solving,
  R multiplies them by the required final event count to obtain operational
  cumulative event triggers.
- `trigger = "information-fraction"` requires
  `schedule.scale = "information-fraction"`. Values lie in `(0,1]` and end at
  1. This representation is used only when statistical information fractions
  are genuinely prespecified inputs.
- Information fraction is not a separate parallel field. For event-count and
  calendar-time schedules, the engine derives and reports information at each
  look. Under a proportional-hazards Schoenfeld approximation, event fractions
  are approximately information fractions. Under an AHR/NPH model they need
  not be equal and must not be equated automatically.

### `efficacy`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `enabled` | boolean | Yes | `true`, `false` | Whether efficacy testing is active. |
| `type` | string | Yes | `"alpha-spending"`, `"fixed"`, `"customized"` | Boundary construction. |
| `family` | object | Yes | See boundary-family definition below | Spending function, parameter, and customized nominal values. |
| `active_at` | array of integers | Yes | Unique values in `1:K` | Looks at which efficacy testing is active. |

For `type = "alpha-spending"`, `family.fn` is required and `family.nominal_values` is `null`. For `type = "customized"`, `family.fn = "customized"` and `family.nominal_values` is required. Custom values are always nominal and use the node's `test_sides`; there is no separate value-type, remainder, or alpha-scale argument.

### Boundary `family`

The same family structure is used for efficacy and formal futility spending. The parent object determines whether `nominal_values` represent alpha or beta.

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `fn` | string or `null` | Conditional | `"sfLDOF"`, `"sfLDPocock"`, `"sfHSD"`, `"sfPower"`, `"sfExponential"`, `"customized"`, `null` | Exact `gsDesign` spending-function name, or `"customized"` for supplied nominal values. |
| `param` | number or `null` | Conditional | Any finite value accepted by the selected function | Spending-function parameter, such as HSD gamma or the power-family exponent. |
| `nominal_values` | array of numbers or `null` | Conditional | Values in `[0,1)`; at most one final `null` | Customized nominal alpha or beta at each active look; `null` uses the remaining available alpha or beta. |

The total efficacy alpha is the alpha available to the node. It is not repeated inside `efficacy`. For beta spending, the total beta is `gsd.<hypothesis_id>.beta`.

### `futility`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `enabled` | boolean | Yes | `true`, `false` | Whether futility is used. |
| `binding` | boolean or `null` | Conditional | `true`, `false`, `null` | Whether crossing affects type I error accounting. Required when enabled. |
| `type` | string or `null` | Conditional | `"beta-spending"`, `"customized"`, `"threshold-rule"`, `null` | Futility representation. |
| `family` | object | Yes | See boundary-family definition above | Beta-spending function, parameter, and customized nominal values; members are `null` for a threshold rule or disabled futility. |
| `threshold_rule` | object | Yes | See below | Threshold-rule structure; its members are `null` unless `type = "threshold-rule"`. |
| `active_at` | array of integers or `null` | Conditional | Unique values in `1:K` | Looks at which the futility rule is evaluated. |

Formal beta-spending example:

```json
{
  "futility": {
    "enabled": true,
    "binding": false,
    "type": "beta-spending",
    "family": {
      "fn": "sfHSD",
      "param": 3,
      "nominal_values": null
    },
    "threshold_rule": {
      "metric": null,
      "threshold": null,
      "direction": null
    },
    "active_at": [1]
  }
}
```

In this example, the enclosing hypothesis uses `"beta": 0.10`; it is not repeated in the futility object.

Threshold-rule example:

```json
{
  "futility": {
    "enabled": true,
    "binding": false,
    "type": "threshold-rule",
    "family": {
      "fn": null,
      "param": null,
      "nominal_values": null
    },
    "threshold_rule": {
      "metric": "conditional-power",
      "threshold": 0.10,
      "direction": "<"
    },
    "active_at": [1]
  }
}
```

#### Futility `threshold_rule`

| Field | JSON type | Required | Allowed values or numeric range | Meaning |
|---|---|---:|---|---|
| `metric` | string or `null` | Conditional | `"conditional-power"`, `"predictive-probability"`, `"observed-effect-measure"`, `"p-value"`, `null` | Quantity compared with the threshold. |
| `threshold` | number or `null` | Conditional | Conditional/predictive power and p-value: `[0,1]`; effect estimate: range determined by the endpoint's `effect_measure`; otherwise `null` | Futility threshold. |
| `direction` | string or `null` | Conditional | `"<"`, `"<="`, `">"`, `">="`, `null` | Operator that causes a futility stop. |

For `"observed-effect-measure"`, the scale comes from the hypothesis endpoint's `effect_measure`. For example, it represents an observed hazard-ratio estimate for a hazard-ratio endpoint and an observed mean-difference estimate for a mean-difference endpoint.

## Cross-object validation

1. Every node's `intervention` and `comparator` must exist in `trial.arms`.
2. Every node's `endpoint` must exist in `trial.endpoints`.
3. Every non-`null` endpoint effect must reference an existing comparator through `vs`.
4. Each node must have a matching `gsd1` entry when group-sequential monitoring applies.
5. Initial confirmatory alpha must not exceed the overall alpha; a conventional graph allocates all available alpha initially, including zero allocations.
6. Edge node IDs must resolve, and graph weights must satisfy their sum constraints.
7. Every `intersection_tests` property must reference existing node IDs. A Dunnett intersection test requires a correlation matrix ordered as `hypotheses`; Bonferroni and weighted-Bonferroni intersection tests have `"correlation": null`.
8. Non-`null` `randomization_weight` values must be positive. A `null` value is allowed only when `role = "component"`; its `component_of` arm must have a positive randomized weight. Comparison ratios are derived and should not be duplicated in nodes.
9. `looks.schedule.scale` must agree with `looks.trigger`. Information-fraction and fraction-of-final-events schedules must end at 1; count schedules in power files must end at `max_events`; calendar-time schedules and all enrollment/dropout windows use months.
10. Piecewise cutpoints must be increasing, begin at zero, and have the same number of effect/rate values as intervals.
11. When `sim.enabled` is `true`, `n_sim`, `seed`, `outcome_models`, and `empirical_quantities` are required. Every outcome-model endpoint must exist in `trial.endpoints`, and its arm entries must match `trial.arms`.
12. When futility is disabled, `binding`, `type`, and `active_at` are `null`; all members of `family` and `threshold_rule` are also `null`.
13. Loose-alpha CoC/CoP evaluations that are not FWER-controlled should be represented as separate endpoint, hypothesis-node, and GSD calculation entries, not connected to the confirmatory alpha graph.
14. Every non-`null` `component_of` value must reference another arm; the referenced arm should have `role = "experimental"`.
15. For a calendar-time calculation, every relevant arm with a `null` distribution must have an acyclic `effects.<arm_id>.vs` path ending at an arm with a non-`null` absolute distribution.
16. For an epsilon-edge AND gate of this form, `procedure.epsilon` must agree with each small upstream-to-downstream edge, each reciprocal upstream-to-upstream edge must equal `1 - epsilon`, and every upstream row must sum to 1. The first downstream node transfers with weight 1 to the next downstream node in the prespecified order.
17. In a time-to-event sample-size file, every sized hypothesis must have non-`null` `gsd.<hypothesis_id>.beta` and `max_events = null`; R solves and reports the required final event count. Event-driven sizing uses `trigger = "event-count"` with `schedule.scale = "fraction-of-final-events"`, while calendar-time looks may be retained when needed to determine an NPH AHR and information trajectory. A genuine information-driven design uses `schedule.scale = "information-fraction"` instead.
18. In an analytic time-to-event power file, `max_events` must be non-`null` and `beta` must be `null`. For analytic calendar-time power, `n_accrual_max` is `null` and R scales the reference accrual shape so the expected cumulative events at the final look equal `max_events`. In a patient-level simulation with fixed enrollment and calendar-time looks, `n_accrual_max` is the simulated enrollment and `max_events` is `null`; event counts are random empirical results. Beta may remain non-`null` only when it calibrates an enabled beta-spending futility boundary.
19. Every stratum key used inside an outcome-model distribution or effect must exist in `trial.strata.definitions`. If an arm uses stratum-specific assumptions, every defined stratum must be present. Simulation uses fixed arm-specific stratum counts derived from the declared proportions.
20. A `derived-bernoulli` endpoint must reference an existing simulated endpoint through `derived_from`, use the same arm and participant record, and must not introduce an independent Bernoulli draw.
21. Every `hypothesis_id` in an empirical expression must exist in `mtp.nodes`. `all_of` and `any_of` expressions define reported operating characteristics and do not alter the MTP.
22. Every empirical-expression `$ref` must be an internal JSON Pointer that resolves to the applicable hypothesis-specific GSD or boundary object. A reference object cannot contain sibling properties.
23. An empirical `grid`, when present, must match a `parameter` in `trial.sim.param_grid`. It is omitted when no grid disambiguation is needed. An expression that references GSD futility or efficacy obtains its active looks, metric, direction, threshold, and boundary construction from that referenced object.
24. Within a design endpoint and within a simulation outcome model, a non-`null` intervention distribution and a non-`null` effect for that intervention are mutually exclusive. The JSON preserves the source R program's representation: both arm distributions with a `null` effect, or a comparator distribution plus an effect with a `null` intervention distribution. Planning and simulation blocks must not add a calculated third representation solely for convenience.
25. `effect_estimator = "coxph"` is valid only for a time-to-event endpoint and produces the observed hazard-ratio estimate; it does not replace the hypothesis-test method.
26. Directly generated simulation outcomes are independent unless their dependence is represented by `derived_from`. R must not silently introduce an unspecified correlation.

## Coverage of the reviewed R programs

| Design feature | Schema location |
|---|---|
| Piecewise control hazards and delayed HR | Endpoint `distributions` and `effects` |
| Dunnett adjustment and shared-control correlation | `mtp.intersection_tests.<intersection_id>` |
| Three-arm NPH and CoC comparisons | Arm weights, endpoint effects, MTP nodes |
| CoC/CoP GSD and nonbinding beta spending | `gsd.<hypothesis_id>.futility` and `gsd.<hypothesis_id>.efficacy.active_at` |
| Continuous and binary endpoints | Endpoint `type`, `test`, `distributions`, and `effects` |
| Parallel, fixed, alternating, and hybrid MTP logic | `mtp.procedure` and `mtp.edges` |
| Patient-level futility simulation | `trial.sim`, enrollment, dropout, endpoint `effect_estimator`, and GSD futility rule |
| Event- or calendar-driven simulation characterization | `trial.sim.param_grid` and `trial.sim.empirical_quantities` |
| Stratum-specific generating assumptions | `trial.strata.definitions` and stratum-keyed arm entries in `trial.sim.outcome_models` |
| Continuous outcomes converted to correlated binary responders | Simulation-only `derived-bernoulli` distributions |
| Marginal and joint empirical power | `trial.sim.empirical_quantities`, including nested `all_of`/`any_of` expressions |
| Sample-size search and QC | Hypothesis beta, GSD event targets, and enrollment |
