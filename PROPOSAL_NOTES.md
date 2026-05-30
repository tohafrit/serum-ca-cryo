# PROPOSAL_NOTES — consolidated findings

Bridge document between the computational model and the InnoCentive submission.
Everything here is reproduced by `make all` (135 unit tests, figures, data tables).
Honest scoping tool: it establishes the root cause, explains the reversibility,
and shows which intervention works — and is explicit about what it does *not*
pin down (the absolute deficit magnitude).

---

## 1. The one-paragraph mechanism

When the serum freezes, the dissolved ions concentrate into the shrinking
unfrozen pool (cryo-concentration, k ≈ 5.58 for 15% glycerol at −20°C). Calcium
and phosphate rise ~5.6×, and the solution becomes strongly supersaturated with
respect to calcium phosphate (SI(HAp) ≈ +7.5, cross-checked against WATEQ).
Amorphous calcium phosphate (ACP) nucleates on the glass surface. It **stays
amorphous**: at the pool viscosity (~84% glycerol, ~4100 mPa·s at −20°C) the
solution-mediated transformation to crystalline hydroxyapatite is ~450× slower,
so it would take centuries — consistent with Combes & Rey (2010). At thaw, a
quiescent draw **under-samples** this surface-bound amorphous deposit, so the
analyzer reads low. Mixing re-disperses/re-dissolves the amorphous deposit and
the reading returns — which is exactly the reported reversibility.

**The hook:** reversibility-with-mixing is the *fingerprint* of an amorphous,
surface-bound deposit. A crystalline phase would not redissolve so easily.

---

## 2. What is robust vs uncertain (state this honestly)

**Robust (does not depend on free parameters):**
- Cryo-concentration and the resulting supersaturation (thermodynamics; SI(HAp)
  ≈ +7.5; Davies vs WATEQ agree to ΔSI < 0.32, and 0.03 at physiological I).
- Crystalline ripening is suppressed at −20°C in the viscous pool → the deposit
  is amorphous (Module 5; matches Combes & Rey).
- The affected-vial fraction grows with storage via **stochastic nucleation**
  (Sobol: nucleation delay ST = 0.71, glass site density ST = 0.46) → explains
  "in some samples" and batch/vial dependence.
- Mixing reverses the deficit (mass-transfer / re-dispersion).

**Uncertain (report as a band; measure by experiment):**
- The precipitated fraction F_PRECIP (mass-balance band ~0.07–0.97; albumin
  buffering is the dominant uncertainty). Representative value 0.90.
- The wall-deposit effective size (ACP_AGGREGATE_NM, ~1–50 µm) controls the
  per-vial deficit and how fast mixing recovers it. Representative 5 µm.
- Together these give a per-vial deficit band of ~0.5–15%; representative ~5%
  (calibrated to the Seeker's reported ≥4%).
- The cryoprotectant identity: a low-viscosity pool would not suppress ripening.

---

## 3. Key numbers (from `make all`)

- k = 5.58 (glycerol 15%, −20°C); Ca 2.75→15.4 mM; Pi 1.8→10.1 mM.
- SI(HAp) ≈ +7.5; WATEQ cross-check max ΔSI = 0.32 (cryo), 0.03 (physiological).
- Ripening suppressed: HAp fraction stays <2% over a year at −20°C.
- Deficit in an affected vial at a quiescent 60-min thaw: ~5% (calibrated to the
  Seeker's reported ≥4%; band 0.5–15%, set by precipitated fraction × deposit size).
- Fraction of vials with a ≥4% deficit: small early (~7% at 1 mo, ~30% at 3 mo),
  ~48% at 6 months, rising to ~78% by 24 months (nucleation-driven onset around
  6 months; the early tail reflects the vial-to-vial spread in nucleation time).
- Defined mixing step / extended standing → affected fraction → ~0 (reversible).
- Loose seal (CO₂ outgassing, higher pH) → more nucleation → higher affected
  fraction (the rationale for tight sealing / degassing).

---

## 4. Interventions (two-tier: prevent OR neutralize)

1. **PREVENT — deep-frozen / vitrified storage (≤ −80°C).** Below the
   freeze-concentrate glass transition (Tg' ≈ −50°C for serum) the pool is
   immobile, nucleation is arrested, and the deposit never forms — so even a
   standard quiescent thaw reads correct. Same viscosity physics that rules out
   ripening; model shows affected fraction → ~0 at −80°C. Trade-off: −80°C cold
   chain. (Module 7 scenario `+deep_freeze`.)
2. **NEUTRALIZE — defined re-suspension/mixing step at thaw.** No equipment, no
   formulation change — SOP + analyte-stability validation. Brings the
   surface-bound calcium back into the sampled volume; model → affected fraction
   ~0. The zero-cold-chain option (the Seeker already sees mixing works).
3. **Reinforce upstream — pre-freeze degassing / tight sealing** (lower pool pH →
   less supersaturation → fewer vials nucleate); a loose out-gassing closure
   makes it worse. Controlled-rate freezing also reduces the affected fraction.

Compliance: all three are process-only. No chemical substance is added or
removed, so the formulation stays as registered (ISO 13485) and REACH/PFAS
status is unchanged. Main risk: shear on fragile analytes (LDH, CK…) → validate
a full analyte panel and use the gentlest mixing that recovers calcium.

---

## 5. Three falsifiable experiments (offer to the Seeker)

1. **DLS / NTA** on a freshly-thawed vs a mixed/equilibrated aliquot of the same
   vial. Predict: particles in the fresh aliquot, far fewer after mixing; none in
   vials too young to have nucleated. (Sizes the deposit → fixes ACP_AGGREGATE_NM.)
2. **Simultaneous ionic Ca (ISE) + total Ca (ICP-MS)** before/after mixing. Both
   drop and both recover → precipitation (our mechanism). Only ionic drops →
   binding, and we are wrong. (Directly measures F_PRECIP.)
3. **Deposit location:** assay bottom/wall rinse vs bulk, or image the surface →
   confirms the under-sampling picture.

A no-particle result in affected vials falsifies the mechanism — a clean go/no-go.

---

## 6. Submission positioning

**Two submissions** (challenge allows up to 3), sharing the Part-1 mechanism and
each leading a distinct, in-scope intervention:
- `docs/submission_S1_production.md` — PREVENT at source: pre-freeze degassing +
  controlled-rate freezing (manufacturer-side, no change to routine lab use).
- `docs/submission_S2_thaw.md` — NEUTRALIZE: a precise re-suspension/mixing thaw
  protocol (the Seeker already sees mixing works; certain to bring deficit <4%).
Both are mapped to the actual Bio-Rad/InnoCentive form fields, Part 1 (mechanism
+ 3 experiments) / Part 2 (intervention + protocol + risk). Differentiators:
our mechanism is Bio-Rad's own "calcium salts / micro-precipitation" hypothesis
made quantitative; we ruled out crystalline ripening rather than asserting it;
honest magnitude band; one-command reproducible model.
