# SUBMISSION S2 — Neutralise at thaw (precise re-suspension / mixing protocol)

Paste each section into the matching field of the InnoCentive form. Character budget per field: 3000. Title: ≤100.

---

## Title (≤100 chars)

Neutralise the post-thaw Ca deficit with a defined re-suspension (mixing) thaw protocol

## Field 1 — Participation Type

Solver (Individual)

## Field 2 — Solution Level (TRL)

TRL 4. The core effect (mixing restores the reading) is already observed by the Seeker; this submission turns it into a defined, validated protocol.

## Field 3 — Problem & Opportunity (≤3000)

You have already seen the key fact: the deficit reverses with additional mixing, or with 24–48 h of equilibration. So the calcium is still in the vial. It is held as a solid that a quiescent draw under-samples, and a precise thawing protocol can recover it every time. The opportunity is to turn that incidental observation into a specified, validated re-suspension step that brings recovery below the 4% threshold reliably, with no new equipment, using only the permitted thawing protocol.

The contribution is the mechanism behind it, worked out quantitatively in a public, reproducible model. Your preliminary studies already point to "calcium salts" and micro-precipitation. We show the precipitate is amorphous calcium phosphate: crystalline hydroxyapatite cannot grow in the viscous cryo-pool, where it would take years to form. An amorphous solid re-disperses easily, which is why gentle mixing recovers the calcium and a crystalline deposit would not. Knowing this also tells us how gentle the mixing can be, which protects the shear-sensitive analytes.

The full model is at github.com/tohafrit/serum-ca-cryo (`make all`, 135 tests).

## Field 4 — Solution Overview (≤3000)

PART 1 — MECHANISM (Solution Requirement 1)

1. Cryo-concentration. On freezing, solutes concentrate into a shrinking unfrozen pool (~5–6× at −20 °C with glycerol); calcium and phosphate become strongly supersaturated for calcium phosphate (SI(HAp) ≈ +7.5, WATEQ-cross-checked). Higher QC calcium/phosphate levels are even more supersaturated.
2. Surface nucleation. Amorphous calcium phosphate nucleates on amber-glass heterogeneities. Nucleation is stochastic and surface-dependent, so only some vials are affected and the fraction grows with storage (the batch/vial variability and the ≥6-month onset).
3. It stays amorphous. The viscous pool suppresses crystalline ripening (~450× slower; Combes & Rey 2010). The deposit remains a loosely-bound amorphous solid.
4. Under-sampling at thaw. A quiescent freshly-thawed draw under-samples the wall-bound deposit, so the analyzer reads low. Mixing or 24–48 h equilibration re-disperses it and the value returns (your observation).

PART 2 — INTERVENTION (Solution Requirement 2): a defined re-suspension thaw protocol. It is formulation-neutral and the only change is to the permitted thaw step. It also has one advantage prevention lacks: it works on already-manufactured stock, so it fixes vials in the field today as well as future lots.

Protocol: thaw at ambient ≤60 min as today, then apply a defined, gentle re-suspension before sampling (for example a fixed number of slow inversions, or a brief low-speed vortex, then a short hold). Choose the gentlest action that restores calcium within acceptance. In the model, moving from a quiescent draw to a thin mixed boundary layer takes the affected fraction to ~0 and the deficit from ≥4% to about 1%, well under the 4% target, and it reproduces your reported reversibility with no parameter fitted to it. Since the deposit is amorphous and loosely bound, gentle inversion may be enough, which avoids high shear. The step uses only the permitted thaw protocol and needs no equipment the lab does not already have. It complements the manufacturing-side prevention in our companion submission.

## Field 5 — Solution Feasibility (≤3000)

This is the most certain route, because the Seeker has already observed that mixing or equilibration reverses the deficit. We are specifying and bounding something that is already known to work.

The supporting physics is standard and cross-validated. Saturation indices (Davies) were checked against the WATEQ extended Debye-Hückel model in USGS PHREEQC: ΔSI ≤ 0.32 at cryo conditions, 0.03 at physiological strength. Re-dispersion is modelled as mass-transfer-limited (Noyes-Whitney): a thinner diffusion boundary layer under mixing speeds recovery. It is the same physics that makes 24–48 h quiescent equilibration work, just faster. Everything is reproducible: `make all` runs 135 unit tests and regenerates all figures and tables.

PART 1 EXPERIMENTS (to prove/disprove the mechanism):

1. DLS or NTA on a freshly-thawed vs a mixed/equilibrated aliquot from the same ≥6-month vial. TRUE → particles present in the fresh aliquot, far fewer after mixing; none in <6-month vials. FALSE → no particles.
2. Simultaneous ionic Ca (ISE) + total Ca (ICP-MS) before and after mixing. TRUE (precipitation) → both drop and both recover. FALSE (binding) → only ionic Ca drops. Also measures the precipitated fraction directly.
3. Deposit location: bottom/wall rinse vs bulk, or surface imaging. TRUE → calcium concentrated at the wall/bottom of fresh, unmixed vials.

A mixing dose-response (speed × duration vs calcium recovery and analyte stability) then fixes the gentlest validated setting.

## Field 6 — Experience (≤3000)

My background is physical chemistry. I hold a specialist degree from the Mendeleev University of Chemical Technology, Moscow (Institute of Physical Chemistry, 2008), specialising in oxide single-crystal growth: nucleation from supersaturated solution, control of growth rate, and characterisation by UV-Vis, IR, ICP and AAS. The physics here is the physics I worked with in the laboratory: nucleation from a supersaturated ionic solution, and an amorphous solid that re-dissolves, with different ions.

For the last several years I have worked as a senior DevOps engineer at a large enterprise (cloud, infrastructure-as-code, CI). That is why this submission ships as a reproducible repository with 135 unit tests and a one-command build, and why I think in terms of a validated, staged protocol with a clear acceptance criterion and a rollback, which is exactly what an SOP change to a regulated product needs.

I am based in Israel with access to Technion core facilities (DLS, ICP-MS, electron microscopy). I am not a clinical IVD formulator; I rely on reproducibility and first-principles physical chemistry, and on letting you verify every number yourself.

## Field 7 — Partnering

Yes.

## Field 8 — Solution Risks (≤3000)

Main deployment risk: shear on fragile analytes. Vigorous vortexing can lower some enzyme activities (e.g. LDH, CK) by a few percent. Mitigation: because the deposit is amorphous and loosely bound, use the gentlest action that recovers calcium (slow inversions may suffice), and validate a full analyte panel (at minimum LDH, CK, ALP, ALT, albumin, total protein, plus product-critical analytes) across mixing intensities to set the safe, effective setting. This directly addresses "avoid changes to routine laboratory use" beyond the thaw step, and protects the other 40+ analytes in the product.

Seeker-required risk factors:

- pH / ionic strength: a mechanical mixing step changes neither; the thawed, re-diluted product is chemically identical, just homogeneous.
- Metal balance: total calcium (and other metals) is conserved; mixing only returns the wall-bound calcium to the sampled volume.
- Osmolality: unaffected (no solute added or removed).

Compliance: the step is a "precise protocol of sample thawing", which is explicitly permitted. No excipient or substance is added or removed, so it is ISO 13485 formulation-neutral and REACH/PFAS status is unchanged.

Trade-offs: a thaw step adds a small action at the bench, and we minimise it to the gentlest validated move. It is fully effective on already-manufactured stock, where prevention only helps future lots, so it is the immediate remedy and complements the manufacturing fix. The magnitude and morphology uncertainty (band ~0.5–15%) is resolved by Experiment 2 and the dose-response.

## Field 9 — Timeline, capability and costs (≤3000)

Phase 1 — confirm the mechanism (≈4 weeks). On your vials of different storage ages: DLS/NTA on fresh vs mixed aliquots; simultaneous ISE + ICP-MS; wall-vs-bulk calcium. Deliverable: mechanism confirmed and a measured precipitated fraction. Go/no-go before further spend.

Phase 2 — validate the thaw protocol (≈6 weeks). Mixing dose-response: speed × duration × storage age vs calcium recovery, with a full analyte-stability panel at each setting, to fix the gentlest action that brings recovery <4% within your acceptance criteria. Confirm robustness across batches and fill volumes.

Phase 3 — documentation (≈4 weeks). Updated thaw SOP with the defined step and endpoint, ISO 13485 change-control, and a verification study design.

Indicative total ≈ 14 weeks of instrument time and materials. No capital equipment (mixing devices are already present in labs). I work as a contractor from Israel, can run Phase 1 at Technion core facilities, and am open to remote work for analysis and on-site visits for the validation. Genuinely interested in seeing this through.

## Field 10 — Online References (≤3000)

Full model and code: github.com/tohafrit/serum-ca-cryo (`make all` reproduces all figures and 135 unit tests; no proprietary software).

1. Combes C, Rey C (2010). Amorphous calcium phosphates. Acta Biomater 6:3362.
2. Boskey AL, Posner AS (1973). ACP→hydroxyapatite conversion. J Phys Chem 77:2313.
3. Christoffersen J et al. (1989). Dissolution kinetics of calcium hydroxyapatite. J Crystal Growth 94:767.
4. Cheng N-S (2008). Viscosity of glycerol-water mixtures. Ind Eng Chem Res 47:3285.
5. Carpenter JF, Crowe JH (1988). Mechanism of cryoprotection. Cryobiology 25:244.
6. Parkhurst DL, Appelo CAJ (1999). PHREEQC v2 (WATEQ activity model). USGS 99-4259.
7. Fogh-Andersen N et al. (1995). Calcium binding to serum albumin. Clin Chem 41:1522.
8. Fennema O (1973). Solid-liquid equilibria (freezing-point depression). M. Dekker.

## Solution Summary (optional field)

The calcium is still in the vial. It micro-precipitates as amorphous calcium phosphate on the glass and is under-sampled until mixed. A defined, gentle re-suspension thaw step returns it to solution and brings recovery below 4% reliably, and it works on existing stock as well as future lots. Mechanism, numbers and the three confirming experiments are fully reproducible.
