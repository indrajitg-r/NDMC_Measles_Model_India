# NDMC India measles model: TRC methods implementation

| Task | File in `simple/` |
| --- | --- |
| Run 1,000 bootstrap fits | [`run_bootstrap_1000.m`](simple/run_bootstrap_1000.m) |
| Read the main model equations | [`measles_model.m`](simple/measles_model.m) |
| Generate manuscript scenarios 0–11 | [`run_scenarios.m`](simple/run_scenarios.m) |

```matlab
% MATLAB Current Folder: INDIA_CODES_NDMC_MODEL/simple
fit = run_bootstrap_1000;
% Inspect outputs/india_1000/REPORT.txt first.
scenarios = run_scenarios;
```

The scenario command retains the identifiability and fit-adequacy check. 

## Original directory: requirements and quick start

Use MATLAB with Optimization Toolbox (`lsqcurvefit`). The code supports MATLAB R2016a syntax. The main workflow uses `randn`, so Statistics and Machine Learning Toolbox is not required. An optional Windows .NET Framework accelerator is compiled automatically and checked against the MATLAB reference solver. On other systems, or with backend `'matlab'`, the reference MATLAB solver is used. No Parallel Computing Toolbox or external spreadsheet is required.

Open MATLAB in **this folder**, `INDIA_CODES_NDMC_MODEL`, then run:

```matlab
addpath(fullfile(pwd, 'corrected', 'manuscript'))
test_india10_manuscript

% Diagnostic bootstrap first; inspect its REPORT.txt before a large run.
fit = run_india_manuscript(100, 'india_lsq_diagnostic');

% Final uncertainty run after reviewing the diagnostic analysis.
fit = run_india_manuscript(1000, 'india_lsq_1000');

% Separate command; refuses to run if calibration diagnostics have not passed.
scenarios = run_india_manuscript_scenarios('india_lsq_1000');
```

Choose a new run name for each calculation. A run name can contain letters, numbers, underscores and hyphens. Existing results are never overwritten. The calibration writes a checkpoint after the base fit and every bootstrap replicate; automatic checkpoint resumption is not implemented for this new workflow. Interrupted calculations must use a new run name, while the earlier checkpoint remains available for inspection.

For a short execution check only:

```matlab
check = run_india_manuscript(2, 'india_lsq_check', 'check');
```

```matlab
fit = run_india_manuscript(100, 'india_lsq_matlab', 'full', 'matlab');
```

Do not use `addpath(genpath(...))`: the archived scenario folders repeat function names. The entry points add only their shared dependencies. Restart MATLAB after editing the C# accelerator because loaded .NET assemblies cannot be unloaded.

## Methods and parameter order

| Item | Implementation |
| --- | --- |
| Age groups | 0–3, 3–6, 6–9, 9–12 months; 1–2, 2–3, 3–4, 4–5, 5–15 and 15–80 years |
| Compartments | M, S, F, I, R, V1, V2, VSIA; V0 added for scenarios 8–11 |
| Time | Biweeks; 26 complete intervals/year; historical simulation 1970–2023 |
| Observations | 24 annual reported case totals, 2000–2023 |
| Forecasts | 2024–2050; start from each fitted end-2023 state |
| Parameter vector | `[alpha_1_1, alpha_1_2, alpha_1_3, alpha_1_4, rho, kappa]` |
| Transmission | Four rates for under 1, 1–5, 5–15 and 15+ years; no additional 0.3686 multiplier |
| Force of infection | `alpha_age * (1 - 0.05*cos(2*pi*t/26)) * sum(I_j/N_j)` |
| Maternal infections | `kappa * [0.2, 0.6, 1, 0, ..., 0] .* M .* force_of_infection` |
| Parameter limits | Four alpha values up to 5; rho and kappa up to 1. A positive numerical floor of 1e-8 permits log-coordinate optimization. |
| Fitting | Unweighted ordinary nonlinear least squares using `lsqcurvefit`; all six parameters estimated |
| Bootstrap | Independent Normal draws with fitted mean and **SD** 20% of that mean; same least-squares estimator for each refit |
| Integration | RK4, four substeps/biweek; independent refinement checks at eight substeps |

The objective is `sum((predicted_annual_cases - observed_annual_cases).^2)`. A single common scaling constant is applied to both prediction and target to improve numerical conditioning; it does not change the minimizer. Predictions, rather than squared residuals, are passed to `lsqcurvefit`. The Normal error assumption defines the parametric bootstrap; fitting is not the earlier parameter-dependent maximum-likelihood objective.

Bootstrap samples are generated from the newly fitted mean, not an initial parameter guess or legacy saved fit. Normal observations remain continuous and may be negative; they are not rounded, truncated or replaced. Each dataset has its own deterministic seed. Eight starts are used for the observed data, and three for each bootstrap dataset, with a larger retry if convergence or agreement is inadequate. Failed or disagreeing fits remain visible and block scenarios.

## Population flows and manuscript clarifications

The manuscript is the methods reference; the following implementation details are explicit because some equations or table entries are incomplete:

- **Maternal waning:** the equations specify `M_i -> S_i` in the same age group, separately from ageing. The parameter table gives a single rate of `1/6` per biweek but no age-specific vector. The implementation applies that rate in every age group. The previous supplied convention `M3 -> S4` is not used by this entry point.
- **Ageing:** all compartments, including infected and maternally immune people, age simultaneously. The supplied rates are retained: 1/6 per biweek for the four infant bands, 1/26 for the four one-year bands, 1/260 for ages 5–15, and zero outgoing ageing from the final band. Six biweeks is the supplied approximation to three months, also used by the manuscript's waning-rate table.
- **Initialization:** preserve the specified total susceptible population, `0.20 * sum(population(4:9))`, and the supplied age weights where feasible. Cap allocations at age-group population capacity and redistribute the excess, as previously requested. Keep `F = 0.01*S`; account for infected seeds within maternal counts. This avoids negative recovered populations while retaining the original susceptible total.
- **MCV0 sequence:** successful MCV0 recipients remain immune, age and experience mortality, and can receive routine doses. The additional-dose equations show a V1-to-V2 loss but omit its receiving V2 term, while the prose allows routine doses after MCV0. The implementation conserves both `V0 -> V1` and `V1 -> V2` transfers. The latter is enabled only for the additional-dose scenarios, matching the distinction from the primary equations. These receiving terms should be made explicit in the manuscript equations before submission.
- **Coverage inputs:** routine coverage is converted to rates using the supplied 6- and 26-biweek divisors. Campaign coverage acts as a rate over one biweek, not an instantaneous vaccinated fraction. This inherited convention is explicit; an exact pulse-fraction interpretation would require a separate model change and refit.

These corrections change the numerical model. The old manuscript estimates, confidence intervals and figures have **not** been regenerated by editing the code and must not be presented as outputs of this corrected implementation.

## Identifiability before scenarios

Read `corrected/manuscript/outputs/<run_name>/REPORT.txt` and `parameters.csv` after calibration. Checks include:

- local log-sensitivity rank and conditioning;
- optimization convergence and agreement across starts;
- bootstrap log-parameter correlations, boundary frequencies and interval widths;
- population conservation, backend parity and timestep refinement;
- the observed standardized-residual RMS compared with the bootstrap refit distribution.

Scenario generation requires at least 100 bootstrap fits, all fits converged and agreeing, a converged base fit with at least two agreeing starts, full local rank, acceptable numerical checks, no flagged boundaries/strong coupling/wide intervals, and no fit-adequacy warning. Thresholds are defined in `india10_lsq_checks.m` and are practical screening rules, not a structural-identifiability proof. Short test runs cannot establish identifiability. A failed screen is a reason to examine profiles, constraints or model assumptions, not to suppress the warning. Acceptance of the previous two-parameter model does not override this gate.

## Scenarios

Routine coverage remains at the 2023 levels: MCV1 93%, MCV2 90%. Standard MCV1 efficacy is 95%; early MCV1, SIA and MCV0 efficacy are 70%; MCV2 efficacy is 98%. SIA targets ages 9 months–15 years at 16.5% coverage in every India SIA scenario. MCV0 targets the maternal compartment at 6–9 months, with 28% coverage once every three years. Campaigns begin in the first biweek of 2024. Birth and death trends continue from the historical inputs.

| Scenarios | Routine MCV1 | SIA intervals in scenario order | Additional MCV0 |
| --- | --- | --- | --- |
| 0, 1, 2, 3 | 9–12 months | None, 3 years, 2 years, 1 year | None |
| 4, 5, 6, 7 | 6–9 months | None, 3 years, 2 years, 1 year | None |
| 8, 9, 10, 11 | 9–12 months | None, 3 years, 2 years, 1 year | Every 3 years |

The same bootstrap parameter vector is paired across all twelve scenarios. Cases averted are computed against scenario 0 within each replicate; negative values are retained. Interval endpoints are linearly interpolated 2.5th and 97.5th percentiles of fitted model means, not future observation prediction intervals. The code does not infer formal elimination status from a deterministic expected case count.

## Files and outputs

| Path | Purpose |
| --- | --- |
| `corrected/manuscript/` | Current six-parameter least-squares entry points, diagnostics and tests |
| `corrected/india10_data.m`, `india10_rhs.m`, `india10_simulate.m` | Shared historical inputs, differential equations and integration |
| `corrected/India10Kernel.cs` | Optional historical accelerator; manuscript and earlier conventions have separate entry methods |
| `corrected/scenarios/india10_scenario_*.m` | Shared scenario definitions and reference solver |
| `corrected/two_parameter/` | Earlier shared-alpha likelihood analysis; not the current manuscript workflow |
| `corrected/run_india_bootstrap.m` | Earlier six-parameter likelihood comparison; not least squares |
| `INDIA_CODES/` | Unmodified supplied research archive, with known defects; not the corrected entry point |
| `CHECKS.md` | Dated check record and earlier analyses |

Calibration saves `calibration.mat`, `checkpoint.mat`, `REPORT.txt`, labelled parameter/fit/initial-condition CSV files, every synthetic observation and refitted mean, optimizer records, seeds, software version and a source snapshot. `calibration.mat` also contains full sensitivity and correlation diagnostics. Scenario generation verifies the source snapshot and saves `scenarios/scenarios.mat`, `annual_scenarios.csv` and `percent_cases_averted.csv` beneath the calibration run. The portable reference solver is used for six-parameter forecasts.

Historical cases and vaccination inputs are embedded in `india10_data.m`; consult the manuscript's Data description for their WHO, Auzenbergs et al. and population-pyramid sources. Input dates, parameter conventions and initial assumptions are fixed, not automatically refreshed from the web.

## Validation status

MATLAB R2016a equation, regression and all twelve scenario-configuration tests passed. A base calibration and two bootstrap refits verified fitting, deterministic 20%-SD generation, CSV export and the scenario guard. This small run found a deficient local sensitivity rank (4/6) and boundary estimates; it did not establish identifiability. No new production bootstrap/scenario ensemble has been run. See [CHECKS.md](CHECKS.md) for measured numerical tolerances and remaining manuscript clarifications.

## GitHub upload

Upload this folder with its subdirectories and the included `.gitignore`. Working `outputs/`, compiled libraries, temporary files and editor backups are ignored. Previously saved two-parameter result packages elsewhere in the parent repository are historical analyses and should not be labelled as results of this six-parameter method. If publishing selected new results, copy them into a separately named results folder and retain their exact source snapshot and run report.

The supplied archive is retained for provenance; its old `Fit_measles_model_bootstrap.m`, `plot_results.m`, `Phatss.mat` and aggregate-path scripts are not corrected analysis entry points. Do not use the archived saved parameter draws with the new equations. No license is invented for the research code; the authors should select a license before redistribution. This preparation does not upload or publish to GitHub.
