# Internal Project Report
## InnoCentive: Post-thaw Ca deficit in serum QC standards

Working document, not for submission. Summarizes what the model shows, what is
robust, what is uncertain, and how it maps to the submission. All numbers are
reproduced by `make all` (135 unit tests).

---

## Section 1 — What the Seeker reported (public description)

A liquid, human-serum-based QC standard shows a **reversible** post-thaw drop in
measured calcium: a freshly-thawed sample reads lower than the added amount; the
effect is **reversible with additional mixing**; it is **batch-to-batch and
vial-to-vial dependent** ("in some samples"). Goal: root cause + mechanism +
an intervention to prevent it. Constraints (typical for this challenge):
formulation-neutral, ISO 13485 / REACH / PFAS compliant, process changes only
(freezing / degassing / thawing).

Note: the public page does not state a numeric threshold, a storage temperature,
or a 6-month onset. We therefore avoid claiming those specifics and frame the
mechanism around what is stated. (If the full brief gives them, they can be
re-added.)

---

## Section 2 — Mechanism (what the model shows)

Four steps; all standard physical chemistry.

1. **Cryo-concentration.** On freezing, ions concentrate into the unfrozen pool.
   For 15% glycerol at −20°C, k ≈ 5.58: Ca 2.75→15.4 mM, Pi 1.8→10.1 mM, ionic
   strength 0.14→0.80 mol/kg. (Module 2.)

2. **Supersaturation.** At the pool composition SI(HAp) ≈ +7.5. Computed with the
   Davies model and cross-checked against WATEQ extended Debye-Hückel (PHREEQC):
   max ΔSI = 0.32 at cryo conditions, 0.03 at physiological — so the activity
   model is not the dominant uncertainty. Serum is forced past the inhibitors
   that keep it metastable, and amorphous calcium phosphate (ACP) nucleates on
   glass-surface heterogeneities. (Modules 3, 4, 8.)

3. **The precipitate stays amorphous (key correction).** We tested whether ACP
   ripens to crystalline hydroxyapatite during frozen storage. It does not: the
   cryo-pool is ~84% glycerol with viscosity ~4100 mPa·s at −20°C (Cheng 2008).
   Solution-mediated transformation is diffusion-limited (Stokes-Einstein), so it
   is ~450× slower than a naïve dilute-serum estimate. Reaching even a few percent
   HAp would take centuries; over a year the precipitate stays >95% ACP. This
   matches Combes & Rey (2010): ACP is kinetically stable for months below 0°C.
   **An earlier version of this model used the dilute-serum viscosity (~9 mPa·s)
   and wrongly produced fast ripening to HAp; that was a 457× error and is now
   corrected.** (Module 5.)

4. **Under-sampling at thaw, reversible by mixing.** The amorphous deposit sits on
   the glass surface (it nucleated there) and aggregates into ~µm-scale clusters.
   A quiescent freshly-thawed draw under-samples it → low reading. Mixing thins
   the diffusion boundary layer (~10 µm → 1–2 µm), re-dispersing/re-dissolving the
   deposit → reading returns. Extended standing does the same by diffusion alone.
   (Modules 5–7.)

**Why the three observations follow:**
- *Reversible with mixing* — an amorphous, surface-bound deposit re-disperses;
  a crystalline phase would not. Reversibility ⇒ amorphous.
- *In some samples / batch dependence* — nucleation is stochastic and
  surface-catalysed (Sobol ST(nucleation)=0.71, ST(glass)=0.46); only some vials
  carry an appreciable deposit, and the fraction grows with storage.
- *Develops over storage* — the affected *fraction* grows as more vials cross
  their nucleation induction time (not because the per-vial amount ripens).

---

## Section 3 — Numbers and honesty about magnitude

Per-vial deficit at a quiescent 60-min thaw ≈ **1.6%** for representative
parameters. This is **not** a first-principles prediction: it scales with the
precipitated fraction F_PRECIP (mass-balance band 0.07–0.97; albumin buffering
dominant) and the deposit size (~1–50 µm). The honest band is **~0.5–15%**, which
brackets a plausible observed deficit. The model robustly predicts the
*direction* and *drivers*; the magnitude is set by quantities the experiments
measure.

Affected-vial fraction (Module 6, 10,000 vials): ~25% early, rising to ~70% and
then ~90% over extended storage — the nucleation-driven "in some samples".

Interventions (Module 7, fraction with a measurable deficit at the mid-storage
point): baseline 0.71; loose seal (higher pH) 0.87; **deep-freeze ≤ −80°C → ~0
(prevents)**; **mixing step → ~0 (neutralizes)**; CRF 0.48 (fewer vials
nucleate); extended standing → ~0.

---

## Section 4 — Risk assessment

- **Mixing & shear-sensitive analytes** (main deployment risk): vortexing can
  reduce LDH/CK activity by single-digit % in plasma; serum less so. Validate a
  full analyte panel; use the gentlest mixing that recovers calcium (the deposit
  is amorphous and loosely bound, so gentle inversion may suffice).
- **Magnitude uncertainty**: explicit band 0.5–15%; Experiment 2 measures it.
- **Cryoprotectant assumption**: a low-viscosity pool would not suppress ripening;
  the qualitative mechanism holds across cryoprotectants, the kinetics do not.
- **Compliance**: all interventions are process-only; no substance added/removed
  → formulation unchanged (ISO 13485), REACH/PFAS status unchanged.
- **Seeker-specified factors**: pH/ionic strength revert on thaw; osmolality
  unaffected; total calcium conserved (only the *sampled* fraction changes).

---

## Section 5 — What is solid vs uncertain

**Solid:** cryo-concentration and supersaturation (thermodynamics, WATEQ-checked);
ripening suppressed → amorphous deposit (viscosity + literature); nucleation-driven
vial-to-vial variability; reversibility by mixing (mass transfer).

**Uncertain:** absolute deficit magnitude (F_PRECIP × deposit size); cryoprotectant
identity; exact deposit morphology/location. All are directly measurable.

What would improve the model most: (1) ISE+ICP-MS on fresh vs mixed aliquots
(gives F_PRECIP and confirms precipitation vs binding); (2) DLS/NTA + imaging
(gives deposit size/location); (3) the actual cryoprotectant.

---

## Section 6 — Three experiments to offer the Seeker

1. DLS/NTA on freshly-thawed vs mixed aliquot (particles present then cleared;
   none in young vials). Sizes the deposit.
2. Simultaneous ISE (ionic Ca) + ICP-MS (total Ca) before/after mixing. Both drop
   and recover → precipitation; only ionic drops → binding (falsifies us).
   Measures the precipitated fraction.
3. Deposit location: bottom/wall rinse vs bulk, or surface imaging.

---

## Section 7 — Submission strategy

Single strong submission, **two-tier solution** (`docs/submission_final.md`):
**prevent** (deep-frozen / vitrified storage ≤ −80°C → nucleation arrested) or
**neutralize** (defined mixing step, no new cold chain), with the mechanism as
rationale, an honest magnitude band, decisive experiments, and a reproducible
repo. Differentiators vs a crowded field: the deep-freeze prevention falls out
of our rigorous viscosity analysis (the same physics that *ruled out* crystalline
ripening), we are explicit about the magnitude band, and every number is
reproducible with one command.
