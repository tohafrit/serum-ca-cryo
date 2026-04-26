# Submission Text — Variant B
## Mechanism + vortex-only intervention (cheapest deployment)
## InnoCentive: Post-thaw Ca deficit in serum QC standards

Positioning: practical quick-win. Same mechanism as Variant A, intervention
section focuses on vortex only. Lower cost, no capital equipment, faster to
validate. 12-month results are not as good as Combined+, but 6-month result
is near-elimination. Good for a Seeker who wants to move fast.

Figures to attach: fig04, fig05, fig06, fig08, fig09.

---

## FIELD 1: PROBLEM AND OPPORTUNITY

Clinical laboratories depend on QC standards to know their calcium results are
correct. A post-thaw calcium deficit in the standard looks like an analyzer
problem. Finding the root cause, and correcting it with a protocol change that
takes less than 5 minutes and requires no new equipment, is the opportunity here.

The deficit appears after more than 6 months of storage at -20°C. Before 6
months, vials are fine. After 6 months, some vials show 4% or more calcium drop.
Not all vials, not all batches. Mixing or 48 hours of cold equilibration reverses
the effect completely.

Three observations that need one explanation: a threshold at 6 months, a
reversal with simple physical steps, and variability between vials from the same
batch. Standard explanations like protein binding or freezing artifacts do not
account for all three at once.

The proposed cause is cryoconcentration during freezing, which drives calcium
and phosphate to supersaturation levels that cause a small amount of calcium
phosphate to precipitate. Over months, this precipitate transforms into
hydroxyapatite crystals. A standard 60-minute thaw is not long enough for those
crystals to dissolve. Mixing or extended equilibration gives them more time.

The model, built from published physical chemistry without parameter fitting to
the Seeker's data, predicts a 5.4% calcium deficit at 6 months and 42% of vials
above the 4% threshold. It also reproduces the 48-hour workaround independently.
The full code is at github.com/tohafrit/serum-ca-cryo.

---

## FIELD 2: SOLUTION OVERVIEW

### Part 1 — Mechanism

Freezing concentrates ions in the remaining liquid. For glycerol 15% w/w at
-20°C, calcium goes from 2.75 mM to 15.4 mM and phosphate goes from 1.8 mM to
10.1 mM. The ion product increases 31 times. Hydroxyapatite saturation index
rises to +7.5, well above the precipitation threshold.

Amorphous calcium phosphate (ACP) nucleates on glass surface heterogeneities.
The nucleation delay varies between vials: this is why some vials show the
deficit and others do not. Sobol sensitivity analysis shows nucleation delay
explains 71% of the variability in outcome.

During storage, ACP transforms to octacalcium phosphate and then to
hydroxyapatite (Ostwald ripening). At 6 months, 7% of the precipitated calcium
is hydroxyapatite. At 12 months, it is 32%. HAp crystals at 6 months are about
200-400 nm in radius. Their dissolution half-life at 22°C quiescent is
about 3 hours.

In a standard 60-minute quiescent thaw, only about 35% of HAp dissolves. The
rest stays as crystals. The analyzer measures dissolved calcium, so it reads low.

### Part 2 — Recommended intervention: vortex thaw

The Noyes-Whitney equation for dissolution shows that mass transfer rate is
proportional to 1/h, where h is the thickness of the stagnant diffusion layer
around each crystal. Quiescent conditions give h about 10 micrometers. Vortexing
at 1000-1500 rpm reduces h to about 1-3 micrometers, increasing dissolution
rate by 3 to 10 times.

A 30-second vortex pulse applied at 30 minutes post-thaw reduces the fraction
of vials above 4% from 57.9% to 27.6% at 6 months. Extending to a 60-second
pulse brings it to 6.3%. These numbers come from Monte Carlo simulation of
10,000 vials (Module 7 in the repository).

I think the vortex protocol is the most practical starting point. It requires
no capital equipment. A standard lab vortex is enough. The SOP change is:
vortex 30-60 seconds at 30 minutes after starting the thaw. The only validation
required is a panel of shear-sensitive analytes before and after.

The vortex alone does not solve the problem at 12 months (77% still above
threshold at 12 months with 30-second pulse). For 12-month storage, the
combined protocol from Variant A is needed. But for most products with 6-month
shelf life, vortex alone is a fast and cheap first step.

Three experiments to confirm the mechanism before deploying: DLS/NTA particle
counting on freshly thawed vs equilibrated vials, simultaneous ionic and total
Ca by ISE and ICP-MS, cryo-SEM of frozen samples at 1 month and 6 months.

---

## FIELD 3: SOLUTION FEASIBILITY

*(Same as Variant A — see submission_variant_a.md Field 3)*

The chemistry is the same. The only difference is the intervention scope.
Vortex thaw has direct precedent in pharmaceutical particle dissolution
(Noyes-Whitney kinetics are standard in USP dissolution testing) and in
blood banking (red blood cells are routinely resuspended by vortex).

---

## FIELD 4: EXPERIENCE

*(Same as Variant A — see submission_variant_a.md Field 4)*

---

## FIELD 5: SOLUTION RISKS

The main risk in this variant is the vortex protocol and its effect on
shear-sensitive analytes. This is addressed in full in submission_variant_a.md
Field 5.

Additional risk specific to this variant: the vortex-only intervention does
not solve the problem at 12 months. If the Seeker's product has 12-month shelf
life, this solution alone is not enough. The recommendation is to start with
vortex validation, and if the 6-month result is confirmed experimentally, add
degassing and controlled-rate freezing later.

---

## FIELD 6: TIMELINE, CAPABILITY AND COSTS

**Phase 1 — Mechanism validation (4 weeks, $8,000-12,000)**

Same as Variant A: DLS/NTA, ISE+ICP-MS, cryo-SEM on Seeker-provided vials.

**Phase 2 — Vortex validation only (4 weeks, $8,000-12,000)**

Vortex dose-response: three speeds (500, 1000, 1500 rpm), three durations
(15, 30, 60 seconds), three storage durations (3, 6, 12 months). Endpoint:
calcium recovery by ISE and ICP-MS.

Full analyte panel before and after vortex: LDH, CK, ALP, ALT, albumin, total
protein, and product-specific critical analytes. This is the shear validation.

**Phase 3 — SOP and documentation (2 weeks, $3,000-5,000)**

Shorter than Variant A because only one intervention needs documentation.

Total: 10 weeks, $19,000-29,000. No capital equipment required.

---

*Editing notes:*
- *Shorten Field 2 if over limit — the vortex mechanism explanation can be simplified*
- *Field 6: if the Seeker's product is 6-month shelf life only, remove the 12-month caveat*
