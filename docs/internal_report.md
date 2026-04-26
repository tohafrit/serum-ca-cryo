# Internal Project Report
## InnoCentive: Post-thaw Ca deficit in serum QC standards

Anton Pakhunov | github.com/tohafrit/serum-ca-cryo | April 2026

This is a working document. It summarizes what I found, what the model predicts,
where the uncertainty is, and how to turn it into a submission. Not for publication.

---

## Section 1 — What the Seeker reported

The Seeker makes liquid human-serum QC standards for clinical analyzers. After
storing vials at -20°C for more than 6 months, they see a calcium drop of 4% or
more when the vials are thawed. The effect is reversible: it goes away after
mixing or after 24-48 hours at 2-8°C. It only appears in some vials, not all.
Vials stored less than 6 months are fine.

The Seeker wants to know: why does this happen, and how to stop it?

**Constraints they set:**
- Cannot change the formulation chemistry (no new additives)
- Can only change freezing, degassing, or thawing protocols
- Must stay ISO 13485 compliant
- No REACH/PFAS restricted substances

**Prize:** $40,000 total, $15,000 guaranteed minimum. Up to 3 solutions can be
submitted (only last 3 count). Deadline: June 1, 2026.

---

## Section 2 — What the model found: mechanism

The mechanism has four steps. All of them are standard physical chemistry, no
exotic assumptions.

### Step 1: Cryoconcentration during freezing (Module 2, fig01)

When serum freezes, most water turns to ice. The ions stay in a shrinking pool
of liquid. This is called cryoconcentration. For glycerol 15% w/w (typical
cryoprotectant), the concentration factor k = 5.58 at -20°C. Derived from
Raoult's law, consistent with Fennema 1973 experimental data.

What this means in numbers:
- Ca in the pool: 2.75 mM -> 15.4 mM
- Inorganic phosphate: 1.8 mM -> 10.1 mM
- Ion product Ca x Pi: increases 31-fold vs physiological
- Ionic strength: 0.14 -> 0.80 mol/kg

Also: CO2 escapes during freeze. Local pH rises from 7.4 to 7.8-8.8 depending
on how much CO2 is lost.

Different cryoprotectants give different k: DMSO 10% gives k=6.8, no
cryoprotectant gives k=35.5. All three scenarios give the same qualitative
result (serum is supersaturated with respect to calcium phosphate).

### Step 2: Supersaturation and nucleation (Module 3, fig02)

At the cryoconcentrated state, the saturation index for hydroxyapatite is:

**SI(HAp) = +7.55 +/- 0.9**

This was calculated with the Davies activity model (standard for I < 0.5 mol/L,
approximate for I = 0.8 mol/L). I cross-checked it against WATEQ extended
Debye-Huckel, which is the exact activity model used in PHREEQC (USGS
standard). Agreement: delta-SI < 0.35 for all four phases at the cryoconc
state (Module 8, fig11). The uncertainty from activity model choice is one
order of magnitude smaller than the uncertainty from Ksp choice, so this is
not a problem.

All four phases are supersaturated at the cryoconcentrated state:

| Phase | SI (Davies) | SI (WATEQ) | Delta |
|-------|-------------|------------|-------|
| Hydroxyapatite | +15.75 | +15.49 | 0.26 |
| Calcite | +2.58 | +2.27 | 0.32 |
| Brushite | +0.63 | +0.76 | 0.13 |
| Monetite | +0.94 | +1.07 | 0.13 |

Note: the SI values in Module 8 are higher than in Module 3 (+7.55) because
Module 8 does not include albumin binding (PHREEQC does not model serum
proteins). Module 3 is the realistic estimate for serum. Both are far above 0.

Nucleation happens on glass surface heterogeneities. The delay is stochastic:
modeled as Gamma(shape=3, scale=30 days) per vial. This is the explanation
for vial-to-vial variability.

### Step 3: Ostwald ripening at -20°C (Module 5, fig04)

ACP converts to OCP, then OCP to HAp. Rate constants from Boskey & Posner 1973,
Heughebaert & Nancollas 1984, Christoffersen et al. 1989. Extrapolated to -20°C
via Arrhenius (Ea ~ 60-80 kJ/mol) with viscosity correction for glycerol
(eta x10, Stokes-Einstein).

Phase distribution over time:

| Time | ACP | OCP | HAp |
|------|-----|-----|-----|
| 1 month | 90% | 9% | 1% |
| 3 months | 72% | 25% | 3% |
| 6 months | 52% | 41% | 7% |
| 12 months | 22% | 46% | 32% |

The 6-month point is key: 7% of precipitated calcium is in HAp form. That is
enough to produce a measurable deficit in a 60-min thaw.

### Step 4: Thawing — not enough time for HAp to dissolve (Module 5, fig05)

Dissolution follows Noyes-Whitney with Ostwald-Freundlich size correction.
At 22°C, quiescent (no mixing):

- ACP (radius ~50 nm): dissolves in < 10 min. No deficit if storage < 3 months.
- OCP (radius ~100 nm): dissolves in 20-40 min. Partial deficit.
- HAp (radius ~200-400 nm at 6 months): half-life ~3.2 hours. Only ~35%
  dissolves in 60 min.

Predicted calcium deficit at 60 min (quiescent 22°C):

| Storage | Deficit |
|---------|---------|
| 1 month | 1.1% |
| 3 months | 3.1% |
| **6 months** | **5.4%** |
| 12 months | 8.5% |

**The 5.4% at 6 months matches the Seeker's reported >= 4% threshold.**
No parameter fitting on Seeker data. The number comes from literature kinetics.

### Why the three observations are explained

**6-month threshold.** At 3 months, 72% of precipitate is ACP — dissolves in
minutes, no deficit. At 6 months, 7% is HAp with hours-scale dissolution time.
That threshold is Ostwald ripening kinetics, not a chosen parameter.

**Reversibility.** HAp is thermodynamically metastable at physiological ionic
strength (SI = +5.1 means dissolution is thermodynamically favorable). It is
only kinetically stuck. Give it 48 hours or vigorous mixing and it dissolves.
No irreversible chemistry happens.

**Vial-to-vial variability.** Nucleation delay is stochastic (glass surface
heterogeneities). Sobol total-effect index for nucleation delay is ST = 0.71,
the single dominant parameter. Vials that nucleate early have more HAp at
6 months. Vials that nucleate late are still mostly ACP. This is why some vials
show the effect and some do not.

---

## Section 3 — Interventions

Three protocol changes, no formulation chemistry modification. All tested in
Monte Carlo simulation of 10,000 vials (Module 7).

### Intervention A: Pre-freeze vacuum degassing (10% residual CO2)

Removes CO2 before freezing. Less CO2 means less pH rise during freeze, which
means lower supersaturation, which means fewer and later nuclei.

In the model: k_sig (nucleation driving force parameter) decreases from 0.15
to 0.10.

### Intervention B: Controlled-rate freezing (2 degrees C per min)

Slower freeze creates more uniform cryoconcentration. Local peaks of [Ca] and
[Pi] are smaller. Fewer nucleation sites are activated.

In the model: k_sig = 0.05, nucleation multiplier = 2.5 (later nucleation).
Requires a controlled-rate freezer (~$15,000-30,000 capital cost).

### Intervention C: Double-pulse vortex thaw (30 s at 5 min + 60 s at 25 min)

The main mechanism is Noyes-Whitney: vortexing reduces the stagnant diffusion
layer around each crystal from ~10 micrometers (quiescent) to ~1 micrometer.
Mass transfer rate increases about 10 times. Using 90-min window instead of
60-min gives additional time.

In the model: thaw_h = 1 micrometer, thaw_min = 90. No change to any
chemistry parameter. Test confirms formulation neutrality:
`test_formulation_chemistry_unchanged_by_vortex` passes.

I think the vortex protocol is the most practical of the three. It requires no
capital equipment and can be validated in a few weeks.

### Quantitative results (10,000 vials, seed=42)

| Scenario | Mean deficit 6mo | >4% at 6mo | Mean deficit 12mo | >4% at 12mo | Complexity | Thaw time |
|----------|-----------------|------------|------------------|-------------|------------|-----------|
| Baseline (pH 8.0) | 3.97% | 57.9% | 7.6% | 80.6% | 1 | 60 min |
| +Degas 10% CO2 | 2.88% | 42.2% | 6.1% | 77.6% | 2 | 60 min |
| +CRF 2°C/min | 2.31% | 32.8% | 5.5% | 60.8% | 5 | 60 min |
| +Vortex 30s | 2.49% | 27.6% | 5.2% | 77.2% | 2 | 60 min |
| +Vortex 60s | 2.05% | 6.3% | 4.5% | 74.3% | 2 | 60 min |
| +Combined | 1.04% | 0.4% | 2.95% | 46.3% | 5 | 60 min |
| **+Combined+** | **0.32%** | **0.0%** | **1.15%** | **0.0%** | **5** | **90 min** |
| Seeker workaround (48h, 4°C) | 0.16% | 0.0% | 0.57% | 0.0% | 1 | 2880 min |

**Combined+** = degassing + CRF 2°C/min + double-pulse vortex + 90-min window.
Result: 0.0% vials above 4% at 12 months. Same result as Seeker's existing
48-hour workaround, but 32 times faster.

The model predicts the 48-hour workaround result independently, without fitting.
This is the main validation of the mechanism — if the physics is right, the
workaround should work, and it does.

---

## Section 4 — Risk assessment

| Intervention | pH change | Osmolality | Metal balance | ISO 13485 | Shear risk | Capex |
|---|---|---|---|---|---|---|
| Degassing | -0.3 units | none | none | OK | none | low (vacuum line) |
| CRF 2°C/min | none | none | none | OK | none | high (~$20K) |
| Vortex thaw | none | none | none | OK | **YES** | none |

**The shear risk is the main concern for vortex.** At 1500 rpm for 30-90 seconds,
literature reports 5-15% decrease in LDH and CK activity in plasma. Serum is
less sensitive than plasma, but still needs validation.

Required before deployment: test panel of 12-15 analytes before and after
vortex protocol. At minimum: LDH, CK, ALP, ALT, albumin, total protein.

Fallback: 500 rpm for 120 seconds gives boundary layer ~3 micrometers (3x
better than quiescent). Less efficacious but probably safe for enzymes.

Degassing and CRF have no known risks related to analyte stability. Their
main validation need is process control documentation (residual CO2 measurement
for degassing, freeze profile recording for CRF).

---

## Section 5 — Epistemic honesty: what is solid and what is not

**Solid:**
- Cryoconcentration factor k=5.58 for glycerol 15% (derived from Raoult's law,
  matches Fennema 1973). Tested for 3 cryoprotectant scenarios.
- SI(HAp) > 0 at cryoconcentrated state is robust across all parameter
  combinations. Serum is definitely supersaturated.
- Qualitative mechanism (ACP -> OCP -> HAp, slow HAp dissolution) is consistent
  with biomineralization literature for 50+ years.
- Nucleation delay as dominant variability source: Sobol ST=0.71, confirmed
  without any parameter fitting to Seeker data.
- Activity model: Davies cross-checked with WATEQ, delta-SI < 0.35 (Module 8).

**Uncertain:**
- Ksp(ACP): Boskey 1973 vs Christoffersen 1990 gives delta-SI = 2.8. I report
  as a band, not a point. Mitigation: does not affect qualitative conclusions.
- Arrhenius extrapolation to -20°C: factor-3 uncertainty on rate constants.
  Mitigation: 6-month threshold is robust across the full range.
- Cryoprotectant identity unknown: glycerol 15% is assumed. If it is DMSO 10%,
  k changes to 6.8 and SI changes by ~0.5 log-unit. Still clearly supersaturated.
- Albumin-Ca binding constant: Fogh-Andersen 1995 used, +/-15% range across
  literature. Small effect on final result.
- Glass surface parameters: Sobol ST=0.46, but I used a best estimate from
  literature, not measured for the specific vials.

What would improve the model most: (1) ICP-MS measurement of Ca in solution
after centrifuging 6-month vials -> gives the real precipitation fraction.
(2) Cryoprotectant specification from the Seeker.

---

## Section 6 — Three experiments to offer the Seeker

These are designed so the Seeker can run them in their own lab with standard
equipment. Each one directly confirms or rules out the mechanism.

**Experiment 1: DLS or NTA on freshly thawed vs equilibrated vials**

Take 6-month vials. Split each into two aliquots. One measured right after
60-min thaw at 22°C. One measured after 48-h at 4°C. Look for particles.

Prediction: 50-500 nm particles in the fresh aliquot, none in the equilibrated
one. Control: same test on 1-month vials, no particles in either case.

If no particles are found in 6-month fresh aliquots, the mechanism is wrong.

Cost: Malvern Zetasizer at university core facility, ~$500/day.

**Experiment 2: Ionic Ca (ISE) and total Ca (ICP-MS) measured simultaneously**

Same vials, same time points. Ionic Ca by calcium-selective electrode. Total Ca
by ICP-MS on the same aliquot.

Prediction: both drop proportionally -> precipitation (our mechanism). If only
ionic Ca drops and total stays the same -> it is binding, not precipitation,
and we are wrong.

After 48-h equilibration: both return to 100%.

Cost: ISE setup ~$200, ICP-MS at core facility ~$100/sample.

**Experiment 3: Cryo-SEM at 1 month vs 6 months**

Freeze vials at -20°C. After 1 month and 6 months, vitrify and image by
cryo-SEM.

Prediction: 1 month - amorphous aggregates < 50 nm (ACP). 6 months -
crystalline particles > 200 nm with defined faces (HAp). This directly shows
Ostwald ripening happening in the vials.

Cost: cryo-SEM session at university facility, ~$300/session.

---

## Section 7 — Submission strategy

Three variants, each with different positioning:

**Variant A (primary):** Full mechanism + Combined+ protocol.
- Best for evaluators who want to understand the whole problem.
- Headline: 0.0% vials above 4% at 12 months, 32x faster than the 48-hour
  workaround you already use.
- Submit first (~May 27).

**Variant B:** Mechanism + vortex-only intervention.
- Best for evaluators who want a cheap, fast solution.
- Headline: 30-second vortex reduces problem vials from 58% to 6% at 6 months.
  No new equipment. Validate in 4 weeks.
- Submit ~May 29.

**Variant C:** Mechanism + glass supplier specification.
- Long-term root cause solution.
- Sobol ST(glass surface) = 0.46. Low-silanol glass or siliconized vials would
  shift the nucleation delay distribution later and reduce the number of
  affected vials upstream.
- Submit by May 31. Do not submit on the deadline day (server load risk).

---

## Section 8 — What goes in each form field

**Field 1 — Problem & Opportunity (500 words)**
Start with: clinical labs depend on QC standards, calcium errors corrupt
calibration. Then explain why this is puzzling: 6-month threshold, reversibility,
batch variability — three observations that do not have an obvious explanation.
End with: my model explains all three from first principles, and the model
independently predicts the 48-hour workaround. Include the GitHub link and
that it reproduces with `make all` in under 2 minutes.

**Field 2 — Solution Overview (500 words)**
Part 1 (~250 words): the mechanism, four steps, key numbers (k=5.58, IAP x31,
SI(HAp)=+7.5, 5.4% at 6 months, 42% vials affected). Part 2 (~200 words):
three interventions, Combined+ result (0%, 90 min, 32x faster). End with
three falsifiable experiments.

**Field 3 — Solution Feasibility (500 words)**
GitHub link, `make all`. WATEQ cross-check (delta-SI < 0.35). Literature
base: Boskey 1973, Christoffersen 1989, Carpenter/Pikal for cryoconcentration,
Combes & Rey 2010 review. Same calcium phosphate chemistry as bone and teeth,
just in a different liquid. The novelty is applying this to clinical QC standards.

**Field 4 — Experience (500 words)**
Honest: independent researcher, physical chemistry degree (Mendeleev/MIPT),
DevOps engineer, Israel. Not a clinical formulator. What compensates for that:
139 tests, 11 figures, full reproducibility, WATEQ cross-check. The Seeker can
clone the repo and verify every number.

**Field 5 — Solution Risks (500 words)**
Per claim: Ksp(ACP) reported as band, not point. Activity model tested with
WATEQ. Cryoprotectant swept. Vortex shear risk: validation panel required
before deployment, specific analytes listed. Glass lot variability: recommend
silanol density spec to glass supplier.

**Field 6 — Timeline, Capability and Costs (500 words)**
Phase 1 (4 weeks, $8-12K): DLS/NTA + ISE/ICP-MS + cryo-SEM on Seeker samples.
Phase 2 (8 weeks, $15-25K): vortex dose-response + analyte panel + CRF pilot.
Phase 3 (4 weeks, $5-8K): SOP + ISO 13485 docs. Total: 16 weeks, $30-45K.
Contractor model from Israel, can travel if needed.
