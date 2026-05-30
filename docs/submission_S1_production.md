# SUBMISSION S1 — Prevent at source (manufacturing: degassing + controlled freezing)

Paste each section into the matching field of the InnoCentive form.
Character budget per field: 3000. Title: ≤100.

---

## Title (≤100 chars)
Prevent Ca micro-precipitation at source: pre-freeze degassing + controlled-rate freezing

## Field 1 — Participation Type
Solver (Individual)

## Field 2 — Solution Level (TRL)
TRL 3. Analytical and experimental proof of concept: the mechanism is quantified
in a reproducible model, and the confirming experiments are specified but not yet
run.

## Field 3 — Problem & Opportunity (≤3000)

The freshly-thawed vial reads low because a small amount of calcium has come out
of solution as a solid. The calcium is still in the vial. Your own preliminary
studies already point to "calcium salts" and list micro-precipitation/phase
separation as a candidate. This submission takes that hypothesis, makes it
quantitative, and shows how to keep the solid from forming in the first place.

Two things set the approach apart. The mechanism is worked through in a fully
reproducible thermodynamic and kinetic model (public repository, one-command
rebuild), so every number can be checked. And it explains the part most accounts
skip: why the loss reverses. The precipitate is amorphous calcium phosphate, not
a crystalline phase. In the cryo-concentrated, glycerol-rich pool the viscosity
reaches about 4000 mPa·s at −20 °C, so ripening to crystalline hydroxyapatite
would take years. An amorphous deposit re-disperses easily when the vial is
mixed. That is why the reading comes back, and it is the clearest sign the
deposit is amorphous.

The deposit forms only while the pool is supersaturated during freezing and
storage. The cheapest place to act is therefore at manufacture, with two levers
you already control: pre-freeze degassing and a controlled freezing protocol.
Both are formulation-neutral and ISO 13485-compatible, and both keep the calcium
in solution so the analyzer reads the labelled value. Because they act at
manufacture, they leave routine laboratory use completely unchanged.

The full model is at github.com/tohafrit/serum-ca-cryo (`make all`, 135 tests).

## Field 4 — Solution Overview (≤3000)

PART 1 — MECHANISM (Solution Requirement 1)
1. Cryo-concentration. As the vial freezes, water turns to ice and the solutes
   concentrate into a shrinking unfrozen pool. At −20 °C with a glycerol
   cryoprotectant the pool is ~5–6× concentrated; calcium and phosphate rise into
   the mM-tens range and the solution becomes strongly supersaturated for calcium
   phosphate (saturation index SI(HAp) ≈ +7.5, cross-checked against the USGS
   WATEQ activity model). At your higher QC calcium/phosphate levels the product
   is even more supersaturated.
2. Surface nucleation. Amorphous calcium phosphate nucleates on amber-glass
   surface heterogeneities. Nucleation is stochastic and surface-dependent, so
   only some vials carry an appreciable deposit at a given time, and that
   fraction grows with storage. This gives the batch-to-batch and vial-to-vial
   variability, and the ≥6-month onset set by the induction time at −20 °C in the
   viscous pool.
3. It stays amorphous. The viscous pool suppresses crystalline ripening (~450×
   slower; matches Combes & Rey 2010, where ACP is stable for months below 0 °C).
4. Under-sampling at thaw. The amorphous deposit sits on the glass. A quiescent
   freshly-thawed draw under-samples it → the analyzer reads low. Mixing or
   24–48 h equilibration re-disperses it → the value returns (your observation).

The model reproduces ≥4% deficits in affected vials and a fraction that grows
with storage; magnitude depends on how much precipitates (an explicit, bounded
parameter the experiments measure).

PART 2 — INTERVENTION (Solution Requirement 2; manufacturing-side, so it makes
NO change to routine laboratory use)
A. Pre-freeze degassing. Removing dissolved CO2 before freezing prevents the pH
   rise that otherwise raises supersaturation in the pool. Lower supersaturation
   → fewer and later nuclei → less precipitate. (Your model run shows a loosely
   sealed, out-gassing vial is markedly worse, so controlling gas and seal helps.)
B. Controlled-rate freezing. A defined, faster, uniform freeze reduces local
   cryo-concentration peaks and shortens time in the precipitation-prone window,
   lowering the affected fraction.
Protocol: (i) inline vacuum degassing of the bulk to a defined residual CO2
before fill; (ii) tight, validated closure; (iii) controlled-rate freezing on a
defined ramp (e.g. ~1–2 °C/min) with recorded thermal profile. All steps are at
the manufacturer; the lab thaws exactly as today. Optional stronger form:
deep-frozen (≤ −80 °C) storage vitrifies the pool and arrests nucleation
entirely. That is the most complete prevention, with a cold-chain trade-off.

## Field 5 — Solution Feasibility (≤3000)

The chemistry is standard and well-supported. Cryo-concentration in frozen
biologicals is established (Carpenter & Crowe 1988; Pikal 1990). Calcium
phosphate metastability and the stability of the amorphous phase at low
temperature are textbook (Combes & Rey 2010). Degassing and controlled-rate
freezing are routine, validated unit operations in regulated manufacturing.

Thermodynamics is the solid core and is cross-validated: saturation indices use
the Davies model and were checked against the WATEQ extended Debye-Hückel model
(USGS PHREEQC). Agreement is ΔSI ≤ 0.32 at the cryo-concentrated condition and
0.03 at physiological strength, so the supersaturation conclusion does not rest
on one activity model. Everything is reproducible: from a clean clone, `make all`
runs 135 unit tests and regenerates every figure and table in minutes.

The absolute deficit is the uncertain part. It scales with the precipitated
fraction and the deposit morphology, both of which the experiments below measure.
We therefore report the magnitude as a band and let the data fix the value.

PART 1 EXPERIMENTS (to prove/disprove the mechanism):
1. DLS or NTA on a freshly-thawed vs a mixed/equilibrated aliquot from the same
   ≥6-month vial. TRUE → 50–500 nm+ particles in the fresh aliquot, far fewer
   after mixing, and none in <6-month vials. FALSE → no particles.
2. Simultaneous ionic Ca (ISE) + total Ca (ICP-MS), before and after mixing.
   TRUE (precipitation) → both drop and both recover. FALSE (binding) → only
   ionic Ca drops. This also measures the precipitated fraction directly.
3. Deposit location: assay a bottom/wall rinse vs the bulk, or image the glass
   surface. TRUE → calcium concentrated at the wall/bottom of fresh vials.

A no-particle result falsifies the mechanism. That gives a clean go/no-go before
any further spend.

## Field 6 — Experience (≤3000)

My background is physical chemistry. I hold a specialist degree from the
Mendeleev University of Chemical Technology, Moscow (Institute of Physical
Chemistry, 2008), specialising in oxide single-crystal growth: nucleation from
supersaturated solution, control of growth rate, and phase verification by
UV-Vis, IR, ICP and AAS. The physics in these vials is the physics I worked with
in the laboratory: nucleation from a supersaturated ionic solution and an
amorphous-to-crystalline pathway, with different ions.

For the last several years I have worked as a senior DevOps engineer at a large
enterprise (cloud, infrastructure-as-code, CI). That is why this submission ships
as a reproducible repository with 135 unit tests and a one-command build rather
than as assertions: every claim is checked and version-controlled, the same
discipline a validated manufacturing change needs.

I am based in Israel with access to Technion core facilities (DLS, ICP-MS,
electron microscopy) for the validation proposed here. I am not a clinical IVD
formulator; I lean on reproducibility and first-principles physical chemistry,
and on letting you verify every number yourself.

## Field 7 — Partnering
Yes.

## Field 8 — Solution Risks (≤3000)

Seeker-required risk factors:
- pH / ionic strength: degassing slightly lowers the CO2-driven pH rise during
  freezing; both pH and ionic strength revert fully on thaw. No change to the
  thawed product the lab measures.
- Metal balance: total calcium (and Mg, Na, K…) is conserved in the vial. The
  intervention only keeps calcium in solution rather than on the glass; nothing
  is added or removed.
- Osmolality: degassing removes a dissolved gas, not solutes; the effect on
  osmolality is negligible (<5 mOsm/kg).

Compliance: both levers are process-only and act exactly within the permitted
scope (initial freezing and degassing). No excipient is added or removed, so the
formulation stays as registered (ISO 13485); no new chemical substance is
introduced, so REACH and PFAS status is unchanged.

Efficacy and trade-offs:
- Degassing and controlled freezing lower the supersaturation and the number of
  affected vials, and they bring the mean deficit below 4% in the model. A
  residual fraction of vials can still precipitate. To guarantee <4% in every
  vial, combine them with the defined thaw step in our companion submission.
- The optional deep-freeze (≤ −80 °C) form prevents the deposit entirely, at the
  cost of a colder cold chain. It suits products that can carry that logistics
  burden; otherwise the two levers above are enough for most lots.
- Magnitude uncertainty: the absolute deficit depends on the precipitated
  fraction/morphology (band ~0.5–15%); Experiment 2 measures it directly, so the
  pilot is self-correcting.
- Cryoprotectant dependence: a much less viscous pool would change the kinetics;
  confirming your cryoprotectant sharpens the model immediately.

## Field 9 — Timeline, capability and costs (≤3000)

Phase 1 — confirm the mechanism (≈4 weeks). On your vials of different storage
ages: DLS/NTA on fresh vs mixed aliquots; simultaneous ISE + ICP-MS; a
wall-vs-bulk calcium check. Deliverable: yes/no on the mechanism and a measured
precipitated fraction. Clear go/no-go before further spend.

Phase 2 — manufacturing pilot (≈6–8 weeks). Degassing pilot (residual CO2 by
headspace GC; pH and Ca recovery vs control). Controlled-rate freezing pilot
(thermal profile; deficit distribution at ≥6 months vs current process).
Optional deep-freeze arm (−20 vs ≤ −80 °C storage) for comparison.

Phase 3 — documentation (≈4 weeks). SOP for the degassing + freezing steps,
ISO 13485 change-control package, and a stability/verification study design.

Indicative total ≈ 14–16 weeks of instrument time and materials. Degassing is
low capital (inline degasser); controlled-rate freezing may need a unit if not
already present. I work as a contractor from Israel, can run Phase 1 at Technion
core facilities, and am open to remote work for analysis and on-site visits for
the manufacturing pilots. I am genuinely interested in seeing this through.

## Field 10 — Online References (≤3000)

Full model and code: github.com/tohafrit/serum-ca-cryo (`make all` reproduces all
figures and 135 unit tests; no proprietary software).

1. Combes C, Rey C (2010). Amorphous calcium phosphates. Acta Biomater 6:3362.
2. Carpenter JF, Crowe JH (1988). Mechanism of cryoprotection. Cryobiology 25:244.
3. Pikal MJ (1990). Freeze-drying of proteins, Part I. Pharm Biotechnol 2:120.
4. Boskey AL, Posner AS (1973). ACP→hydroxyapatite conversion. J Phys Chem 77:2313.
5. Cheng N-S (2008). Viscosity of glycerol-water mixtures. Ind Eng Chem Res 47:3285.
6. Fennema O (1973). Solid-liquid equilibria (freezing-point depression). M. Dekker.
7. Parkhurst DL, Appelo CAJ (1999). PHREEQC v2 (WATEQ activity model). USGS 99-4259.
8. Fogh-Andersen N et al. (1995). Calcium binding to serum albumin. Clin Chem 41:1522.

## Solution Summary (optional field)
The calcium is still in the vial. It micro-precipitates as amorphous calcium
phosphate on the glass and is under-sampled until mixed. Prevent it at source with pre-freeze
degassing and controlled-rate freezing (no change to routine lab use); a defined
thaw step (companion submission) guarantees <4%. Mechanism, numbers and the three
confirming experiments are fully reproducible.
