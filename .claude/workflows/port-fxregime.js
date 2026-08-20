export const meta = {
  name: 'port-fxregime',
  description: 'Port the R package fxregime to Julia, module by module, verified against R',
  whenToUse: 'Rebuilding or extending FXRegime.jl against the R reference implementation',
  phases: [
    { title: 'Implement', detail: 'one agent per module group, each verified against Rscript' },
    { title: 'Integrate', detail: 'assemble the module, make it load and run end to end' },
    { title: 'Fixtures', detail: 'generate R reference outputs as JSON/CSV test fixtures' },
    { title: 'Tests', detail: 'Julia test suite asserting equality with the R fixtures' },
    { title: 'Docs', detail: 'README explaining the method and the package' },
  ],
}

const REPO = '/home/ayush/REPOSITORIES/FXRegime.jl'
const SCRATCH = '/tmp/claude-1000/-home-ayush-REPOSITORIES-FXRegime-jl/567d1d45-361a-40a5-aa9a-6d66cace1f85/scratchpad'
const CONTRACT = `${SCRATCH}/CONTRACT.md`
const RSRC = `${SCRATCH}/src`

const PREAMBLE = `You are implementing part of FXRegime.jl, a Julia port of the R package fxregime.

READ FIRST, in this order:
1. ${CONTRACT} — the binding API contract. Follow it exactly; do not invent different names or signatures.
2. The R sources it points at under ${RSRC}/fxregime/R/ and ${RSRC}/strucchange/R/. These are the spec.
3. ${REPO}/src/series.jl and ${REPO}/src/data.jl — already written. Use FXSeries; do not redefine it.

R is installed. \`Rscript -e '...'\` works and \`library(fxregime)\`, \`library(strucchange)\`,
\`library(zoo)\`, \`library(sandwich)\` all load. (\`car\` may be missing — work around it.)
NEVER guess what an R function does: run it and look.

Julia 1.10 is at \`julia\`. Use \`julia --project=${REPO}\` to get the project's dependencies
(Distributions, and the stdlibs). Available modules inside the package: Statistics, LinearAlgebra,
Dates, Printf, DelimitedFiles, Distributions. Do not add dependencies.

Your files are part of a module, so write bare top-level definitions — no \`module\`, no \`using\`
statements (the module file supplies those), no \`include\` of sibling files.

HOW TO VERIFY (required — do not skip):
Write a scratch driver script under ${SCRATCH}/verify/ that does
\`include\`s of the sibling src files you depend on plus your own, in dependency order, wrapped in
\`using Statistics, LinearAlgebra, Dates, Printf, DelimitedFiles, Distributions\`, then compares
your output against numbers you obtained from Rscript on the same input. Iterate until they agree
to at least 1e-8 relative. Report the achieved agreement in your final answer.

Use FXRatesCHF() as the test input — it is the same data R has, so R and Julia see identical bytes.

Style: docstrings on every exported function naming the R equivalent, explicit Float64,
no unnecessary allocation in the O(n^2) loops.

Your final message must be a short report: files written, what you verified against R, the numeric
agreement achieved, and anything you could not verify. Do not paste the code.`

phase('Implement')

const MODULES = [
  {
    key: 'breakpoints',
    label: 'recresid+breakpoints',
    task: `Implement ${REPO}/src/recresid.jl and ${REPO}/src/breakpoints.jl.

This is the numerical core and the part that most needs to be fast: it replaces an O(n^3 k^2)
brute force with strucchange's O(n^2 k^2) recursive-residual algorithm.

recresid.jl: port strucchange's \`recresid_r\` (${RSRC}/strucchange/R/recresid.R) including the
periodic full-refit check against the rank-1 recursion and the rank-deficiency handling.

breakpoints.jl: port \`RSSi\`/\`RSS.triang\`/\`extend.RSS.table\`/\`extract.breaks\` from
${RSRC}/strucchange/R/breakpoints.R and \`RSS2obj.triang\` + the whole dynamic program and the
information criteria from ${RSRC}/fxregime/R/gbreakpoints.R. Follow the contract's index sets
exactly — an off-by-one here silently changes every breakdate.

Include the intercept-only fast path R has. Thread the RSS triangle rows with Threads.@threads
(rows have very unequal cost, so use :dynamic scheduling or reverse the iteration order).

VERIFY against R:
- \`recresid\` vs \`strucchange::recresid(X, y)\` on random data and on real FX returns.
- The RSS triangle vs \`strucchange:::breakpoints.formula(...)$RSS.triang\` (access it via
  \`bp <- breakpoints(y ~ ., data=d, h=...); bp$RSS.triang\` — check the actual field name in R).
- The full dating vs \`fxregime::fxregimes(CNY ~ USD + JPY + EUR + GBP, data=cny, h=20, breaks=10)\`
  on the daily CNY data from the package vignette (start 2005-07-25, end 2009-07-31). R gives
  breakdates 2006-03-14, 2008-08-22, 2008-12-31 under LWZ. Your breakpoints must be the same
  observation numbers, and the NLL/BIC/LWZ values must match to 1e-8.
  (That R run takes ~3 minutes; yours should be far faster — report the timing of both.)`,
  },
  {
    key: 'model',
    label: 'fxreturns+fxlm',
    task: `Implement ${REPO}/src/returns.jl and ${REPO}/src/fxlm.jl.

returns.jl: port \`fxreturns\` (${RSRC}/fxregime/R/fxreturns.R). The Friday anchoring, the
"last non-missing in the week, per column, independently" aggregation, zoo::na.locf semantics
(including what it does with leading NAs) and the trimming are all easy to get subtly wrong —
check every one of them against R on real data with real gaps.

fxlm.jl: port \`fxlm\`, \`coef.fxlm\`, \`estfun.fxlm\`, \`bread.fxlm\` (${RSRC}/fxregime/R/fxlm.R)
and \`fxpegtest\` (${RSRC}/fxregime/R/fxtools.R). \`bread.fxlm\` builds a (k+1)x(k+1) matrix from
sandwich's \`bread.lm\` plus the beta/variance cross terms — inspect \`sandwich:::bread.lm\` in R
and reproduce it numerically rather than reasoning about it on paper. Also provide the summary
table (OLS standard errors, t, p) and an HC3 sandwich option for vcov.

If \`car\` is unavailable in R, implement fxpegtest as the standard Wald F test of
"peg slope = 1, other slopes = 0" (intercept and variance unrestricted) using the df-adjusted OLS
covariance, and verify it by hand-computing the same F statistic in R with plain matrix algebra
(R = restriction matrix, F = (Rb-q)'(R V R')^{-1}(Rb-q)/q). Say clearly in your report that it was
verified that way rather than against car.

VERIFY against R:
- \`fxreturns("INR", frequency="weekly", start=as.Date("1993-04-01"), end=as.Date("2008-01-04"),
   other=c("USD","JPY","DUR","GBP"), data=FXRatesCHF)\` — same dimensions, same index, same values.
- The same for the daily CNY series of the CNY vignette.
- \`coef(fxlm(...))\` including "(Variance)", \`estfun\`, \`bread\` — all to 1e-10.`,
  },
  {
    key: 'inference',
    label: 'confint+monitor+efp',
    task: `Implement ${REPO}/src/confint.jl, ${REPO}/src/monitor.jl and ${REPO}/src/efp.jl.

confint.jl: port strucchange's unexported \`pargmaxV\` (get the source with
\`Rscript -e 'strucchange:::pargmaxV'\`) and \`confint.fxregimes\`
(${RSRC}/fxregime/R/fxregimes.R). Note the trap the earlier prototype fell into: R forms
\`delta\`, \`Q\` and \`Omega\` over the FULL coefficient vector INCLUDING the variance parameter,
with \`Q = solve(bread(fit))\`, not over the regressors with X'X/n. Reproduce R's expanding
\`ub += 1000\` / \`lb -= 1000\` bracket search and its uniroot, and cap the bracket search so a
non-convergent case cannot loop forever (R will hang; you should raise instead).
Depends on FXRegimes/refit from src/regimes.jl which another agent is writing in parallel —
code against the contract's signature for it and note the dependency; do not write regimes.jl.

monitor.jl: port \`fxmonitor\` (${RSRC}/fxregime/R/fxmonitor.R) including the Table III critical
value matrix verbatim, the two-way \`approx\` interpolation, the Bonferroni correction, and
\`breakpoints\`/\`breakdates\` for the first boundary crossing (read the rest of fxmonitor.R for
those methods).

efp.jl: \`gefp(m::FXLM)\` matching \`strucchange::gefp(fxlm_obj, fit = NULL)\` — note it uses
strucchange's \`root.matrix\` (symmetric square root via eigen decomposition, NOT Cholesky).
\`sctest\` with the double-maximum functional \`maxBB\`: statistic max|process|, p-value
1 - (1 + 2*sum_{i=1..100} (-1)^i exp(-2 i^2 x^2))^k (see pvalue.efp in
${RSRC}/strucchange/R/efp.R, "Brownian bridge"/"max" branch).

VERIFY against R:
- pargmaxV over a grid of (x, xi, phi1, phi2) including negative x — 1e-10.
- \`gefp\` process values and \`sctest(cny_efp)\` statistic and p-value for the CNY vignette model.
- \`fxmonitor(CNY ~ USD+JPY+EUR+GBP, data=window(cny, end=as.Date("2006-05-31")),
   start=as.Date("2005-11-01"), end=4)\` — process, critval, and the reported break date.
- For confint, use the CNY dating result: R gives, at level 0.9, breakpoint 158 with interval
  [143, 159], 778 with [762, 779], 866 with [865, 880]. You will need a segmentation to test
  against; if src/regimes.jl is not ready, hard-code those breakpoints in your scratch driver.`,
  },
]

const impl = await parallel(MODULES.map(mod => () =>
  agent(`${PREAMBLE}\n\n=== YOUR TASK ===\n${mod.task}`, {
    label: mod.label,
    phase: 'Implement',
  })))

phase('Integrate')

const integration = await agent(
  `${PREAMBLE}

=== YOUR TASK ===
Three agents have just written, in parallel:
  src/recresid.jl, src/breakpoints.jl   -> ${impl[0]}
  src/returns.jl, src/fxlm.jl           -> ${impl[1]}
  src/confint.jl, src/monitor.jl, src/efp.jl -> ${impl[2]}

Now write ${REPO}/src/regimes.jl (the fxregimes / refit / coef / fitted / residuals layer from
${RSRC}/fxregime/R/fxregimes.R, per the contract) and ${REPO}/src/FXRegime.jl (the module file:
usings, includes in dependency order, exports), plus port the Shapley R^2 decomposition
(\`calculate_r2_decomposition\`) from ~/Downloads/fxreg/err_functions.jl as \`r2_decomposition\`
— that one is an addition of ours, not in R.

Then make the whole thing actually work:
- \`julia --project=${REPO} -e 'using FXRegime'\` must load cleanly with no warnings.
- Reconcile any disagreements between the three agents' files (duplicate definitions, mismatched
  signatures, names that drifted from the contract). The contract wins.
- Run the full CNY vignette analysis end to end in Julia and check it against R:
  fxreturns -> fxlm -> gefp/sctest -> fxregimes -> breakdates -> coef -> confint -> refit,
  and the INR weekly vignette analysis (h=20, breaks=10) as well.
  R's INR answer: run it yourself to get it, it takes a few minutes.
- Report the end-to-end agreement and the Julia vs R timings.

Fix whatever is broken rather than reporting it as someone else's problem.`,
  { label: 'integrate', phase: 'Integrate' })

phase('Fixtures')

const fixtures = await agent(
  `${PREAMBLE}

=== YOUR TASK ===
The Julia package now exists and loads. Build the R-reference test fixtures.

Write ${REPO}/test/generate_fixtures.R — a self-contained R script that regenerates every fixture
from the FXRatesCHF data shipped with the R package, and writes them to ${REPO}/test/fixtures/ as
plain CSV (one file per fixture, with a header; dates as ISO yyyy-mm-dd). Committing the CSVs
means the Julia tests run without R installed, while the script documents exactly how each number
was produced and lets anyone regenerate them.

Cover at least:
- fxreturns: weekly INR series and daily CNY series (index + all columns).
- fxlm: coefficients incl. (Variance), residuals, fitted, estfun, bread, and the summary table
  (se/t/p) for the CNY history model and the full-sample INR model.
- recresid on both series.
- The RSS triangle for a small subset (say the first 120 CNY observations, h=20) — full triangle.
- fxregimes on CNY (h=20, breaks=10) and INR (h=20, breaks=10): breakpoints for every m from 1 to
  10, and the NLL/BIC/LWZ for every m, plus the selected m under both BIC and LWZ.
- confint at levels 0.9 and 0.95 for both.
- refit/coef per segment for both.
- gefp process + sctest statistic and p-value for both.
- fxmonitor for the CNY vignette setup: process, critval, break date.
- A couple of small synthetic regressions with known structural breaks, so a failure can be
  localised without a 1000-row FX series.

Each fixture CSV must be accompanied by a row in ${REPO}/test/fixtures/MANIFEST.csv giving
file, what it is, the R call that produced it, and the R package versions used.

Run the script. Verify each CSV is non-empty and parses. Report the list of fixtures and the total
size on disk. Keep the whole fixtures directory under 5 MB.`,
  { label: 'fixtures', phase: 'Fixtures' })

phase('Tests')

const tests = await agent(
  `${PREAMBLE}

=== YOUR TASK ===
Write the Julia test suite in ${REPO}/test/, driven by the R fixtures in ${REPO}/test/fixtures/
(see ${REPO}/test/fixtures/MANIFEST.csv and ${REPO}/test/generate_fixtures.R — fixture report
follows below).

${fixtures}

Structure:
  test/runtests.jl        — includes the rest
  test/test_series.jl     — FXSeries, window, data loading
  test/test_returns.jl    — fxreturns vs R
  test/test_fxlm.jl       — fxlm, coef, estfun, bread, vcov, summary, fxpegtest vs R
  test/test_recresid.jl   — recursive residuals vs R
  test/test_breakpoints.jl— RSS triangle, dynamic program, ICs, breakpoints for every m vs R
  test/test_regimes.jl    — fxregimes end to end, refit, coef, breakdates vs R
  test/test_confint.jl    — confidence intervals vs R at both levels
  test/test_monitor.jl    — fxmonitor vs R
  test/test_efp.jl        — gefp + sctest vs R
  test/test_edge.jl       — edge cases: h at its limits, zero breaks, rank-deficient segments,
                            all-missing columns, single-column input, n < 2h, NaN handling.

Rules:
- Every numeric comparison states its tolerance explicitly and justifies anything looser than
  1e-8 with a comment (e.g. recursive vs direct RSS accumulates float error).
- Test the exact index conventions: that breakpoint i is the LAST observation of segment i, and
  that breakdates line up with R's to the day.
- Use @testset names that identify the R function being matched.
- The suite must not need R or the network to run.
- Keep the whole suite under 60 seconds. The full INR/CNY datings are the expensive part — if one
  is slow, keep the CNY one (it is the smaller) in the fast path and mark the INR one to run only
  when ENV["FXREGIME_SLOW_TESTS"] == "true", and make sure the fast path still covers the dating
  logic properly.
- Add ${REPO}/test/Project.toml with Test and whatever else the suite needs.

Then RUN it: \`julia --project=${REPO} -e 'using Pkg; Pkg.test()'\`. Every test must pass. If a test
fails because the Julia implementation is wrong, FIX THE JULIA IMPLEMENTATION in src/ — do not
loosen a tolerance or delete a test to make it green. If you genuinely cannot make one match R,
mark it @test_broken with a comment explaining precisely what differs and why.

Report: number of tests, runtime, anything left @test_broken, and any src/ bug the tests caught.`,
  { label: 'tests', phase: 'Tests' })

phase('Docs')

const docs = await agent(
  `${PREAMBLE}

=== YOUR TASK ===
Write ${REPO}/README.md. Audience: an economist or quantitative researcher who has heard of
"de facto exchange rate regime classification" but does not know this package.

It must explain, in prose that stands on its own:
- What the Frankel-Wei regression is and why cross-currency returns against an independent
  numeraire identify the implicit currency basket weights. Give the model equation.
- Why the variance matters as much as the coefficients: sigma is the flexibility of the regime,
  the basket weights are its anchor, and the package treats the variance as a parameter of the
  model rather than a nuisance — which is exactly why it uses a normal log-likelihood rather than
  plain least squares to date the breaks.
- What the dating procedure does: Bai-Perron style dynamic programming over segment
  log-likelihoods, minimal segment size h, model selection by LWZ or BIC, and confidence
  intervals for the break dates from the argmax of a two-sided Brownian motion with drift.
- What the monitoring procedure does and how it differs from dating (prospective vs retrospective).
- A worked example, copy-pasteable, that reproduces the CNY analysis: the code and its actual
  output (run it and paste the real output, do not invent it).
- A short "relation to R" section: this is a port of the R package fxregime by Achim Zeileis,
  Ajay Shah and Ila Patnaik; the test suite asserts numerical agreement with it; a table mapping
  R function -> Julia function; and what is NOT ported (say so honestly — check the R NAMESPACE
  against our exports and list the gaps, e.g. plotting methods).
- Installation, and how to run the tests, and how to regenerate the R fixtures.
- Citation: the CSDA 2010 paper (Zeileis A, Shah A, Patnaik I (2010), "Testing, Monitoring, and
  Dating Structural Changes in Exchange Rate Regimes", Computational Statistics & Data Analysis
  54(6), 1696-1706, doi:10.1016/j.csda.2009.12.005) and the R package.
- The FXRatesCHF data provenance (US Federal Reserve H.10, public domain) and the note that DUR
  is DEM adjusted to EUR before 1999.

Also write ${REPO}/docs/src/index.md as a thin pointer to the README, and check that the
Documenter setup PkgTemplates left in ${REPO}/docs/ still builds or is at least not broken by
what we added.

Run every code block you put in the README and paste real output. A README that does not run is
worse than no README.`,
  { label: 'readme', phase: 'Docs' })

return {
  implement: impl,
  integrate: integration,
  fixtures,
  tests,
  docs,
}
