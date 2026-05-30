# Title
The reversible deficit is amorphous calcium phosphate on the vial surface,
under-sampled at a quiescent thaw — prevent it by deep-frozen (vitrified)
storage, or neutralize it with a defined mixing step

# Solution Summary
The deficit is reversible because the calcium is not chemically lost — during
freezing the serum cryo-concentrates, becomes strongly supersaturated, and
precipitates amorphous calcium phosphate on the glass surface. A quiescent
freshly-thawed draw under-samples this surface-bound fraction, so the analyzer
reads low. Mixing re-disperses and re-dissolves it, and the reading returns.
The precipitate stays amorphous (we show crystalline hydroxyapatite cannot form
at −20 °C in a viscous cryoprotected pool) — and that is exactly why simple
mixing reverses it.

Two-tier solution, both process-only (no formulation change):
- **Prevent it at the root — deep-frozen / vitrified storage (≤ −80 °C).** The
  precipitate forms only while the cryo-concentrated pool is still mobile. Below
  the glass transition of the freeze-concentrate (Tg' ≈ −50 °C for serum) the
  pool vitrifies and nucleation is arrested — the deposit never forms, so even a
  standard thaw reads correct. This follows directly from the same viscosity
  physics that rules out crystalline ripening, and our model reproduces it (the
  affected fraction falls to ~0 at −80 °C).
- **Or neutralize it with no new cold chain — a defined mixing/re-suspension
  step at thaw** (invert/vortex or a short equilibration), which brings the
  surface-bound calcium back into the sampled volume.

Three lab experiments confirm the mechanism and size the effect.

Reproducible model: github.com/tohafrit/serum-ca-cryo (`make all`, 135 tests).

## FIELD 1: PROBLEM AND OPPORTUNITY

Clinical laboratories trust a liquid QC standard to read its labelled calcium.
When a freshly-thawed vial reads a few percent low, it looks like analyzer
drift and can trigger an unnecessary recalibration or a service call. The
useful insight here is that the calcium is almost certainly still in the vial —
just not where the probe samples it. That reframes an alarming "lost calcium"
into a solvable sampling-and-re-dissolution problem.

The Seeker's own observations point the way. The deficit is reversible with
additional mixing. It appears in some vials and not others, and varies batch to
batch. A process that reverses on mixing is not chemical degradation (that would
be irreversible) and is not a uniform binding shift (that would hit every vial
equally). It is something that forms a physical, re-dispersible solid in some
vials.

That is calcium phosphate. When serum freezes, almost all the water turns to
ice and the dissolved ions concentrate into the shrinking unfrozen pool. In our
model, with 15 % glycerol at −20 °C the pool concentrates about 5.6-fold:
calcium rises from ~2.75 to ~15 mM, inorganic phosphate from ~1.8 to ~10 mM, and
the solution becomes strongly supersaturated with respect to calcium phosphate
(saturation index for hydroxyapatite about +7.5, cross-checked against the USGS
WATEQ activity model). The serum is forced past the inhibitors that normally
keep it metastable, and amorphous calcium phosphate nucleates — preferentially
on glass-surface heterogeneities.

One point matters for why it is reversible. We checked whether this precipitate
ripens into crystalline hydroxyapatite during storage, and found it cannot: the
cryo-concentrated pool is roughly 84 % glycerol, with a viscosity near 4000 mPa·s
at −20 °C, which slows the solution-mediated transformation by more than two
orders of magnitude. Reaching even a few percent crystalline apatite would take
many years, not months. This matches the literature (Combes & Rey 2010: amorphous
calcium phosphate is kinetically stable for months below 0 °C). So the solid
stays amorphous — and amorphous, surface-bound calcium phosphate re-disperses and
re-dissolves easily when you mix the vial. A crystalline phase would not. The
reversibility you see is the fingerprint of an amorphous deposit.

Because nucleation is stochastic and surface-dependent, only some vials carry an
appreciable deposit at any given time, and the affected fraction grows with
storage. That accounts for the vial-to-vial and batch-to-batch variability.

The full analysis is at github.com/tohafrit/serum-ca-cryo; `make all` reproduces
every figure and 135 unit tests in a few minutes.

## FIELD 2: SOLUTION OVERVIEW

There are two ways to act, and they follow directly from the mechanism. Choose by
whether a deep-frozen cold chain is acceptable for the product.

**Option A — PREVENT it (root cause): deep-frozen / vitrified storage (≤ −80 °C).**
The deposit only forms while the cryo-concentrated pool is still a mobile liquid.
That pool is ~84 % glycerol; its viscosity climbs steeply as it cools, and below
the glass transition of the freeze-concentrate (Tg' ≈ −50 °C for serum) it
vitrifies — molecular diffusion stops and nucleation is arrested. Store below
that point (≤ −80 °C is standard practice for sensitive serum analytes) and the
precipitate never forms, so even an ordinary quiescent thaw reads correct. This
is the same viscosity physics that rules out crystalline ripening, used the other
way round: our model's affected fraction falls from the −20 °C value to ~0 at
−80 °C. Fast freezing through the −5…−40 °C window helps by minimising time in
the mobile-concentrated state. Trade-off: a −80 °C cold chain (freezers,
logistics) — which is why Option B exists.

**Option B — NEUTRALIZE it (no new cold chain): a defined re-suspension step at
thaw.** The deposit sits on the glass and re-disperses with mixing — which you
already see works. A specified step (thaw, then invert or briefly vortex, then a
short hold before sampling) brings the surface-bound calcium back into the
sampled volume. No equipment, no formulation change — only an SOP and an
analyte-stability validation. In our re-dispersion model a quiescent draw leaves
a measurable deficit while a defined mixing step brings recovery essentially to
completion, reproducing your reported reversibility with no parameter fitted to it.

Either can be reinforced upstream by **pre-freeze degassing / tight sealing**
(lower pool pH → less supersaturation → fewer vials nucleate). Our model shows a
loose, out-gassing closure makes the problem worse — so controlling that is a
cheap risk reduction.

On magnitude, we are deliberately honest: the per-vial deficit depends on how
much calcium precipitates and on the deposit's size/morphology, both of which we
cannot fix from theory alone. For representative values our model gives a
few-percent deficit, and the plausible band spans roughly 0.5 % to 15 %. The
model robustly predicts the *direction* and *drivers* — that precipitation
occurs, that the affected fraction grows with storage through nucleation, and
that mixing reverses it — and brackets the magnitude. The experiments below pin
it down.

Three experiments the Seeker can run to confirm and quantify:
1. Particle counting (DLS or NTA) on a freshly-thawed vs a mixed/equilibrated
   aliquot from the same vial — expect particles in the fresh aliquot, far fewer
   after mixing; none in vials too young to have nucleated.
2. Simultaneous ionic calcium (ISE) and total calcium (ICP-MS) before and after
   mixing — if both drop and both recover, it is precipitation (our mechanism);
   if only ionic drops, it is binding and we are wrong. This experiment directly
   measures the precipitated fraction.
3. Look for the deposit's location — assay the bottom/wall rinse vs the bulk, or
   image the surface — to confirm the under-sampling picture.

## FIELD 3: SOLUTION FEASIBILITY

The chemistry is standard and old. Cryo-concentration in frozen biological
systems is well established (Carpenter & Crowe 1988; Pikal 1990). Calcium
phosphate precipitation, the metastability of serum, and the stability of the
amorphous phase at low temperature are textbook biomineralization (Combes & Rey
2010). We did not invent a mechanism; we applied a known framework to a QC vial
and removed the parts that do not hold up.

The thermodynamics is the solid core and is cross-validated. Saturation indices
use the Davies activity model; because the cryo-pool ionic strength (~0.8 mol/kg)
is at the edge of Davies' range, we cross-checked every phase against the WATEQ
extended Debye-Hückel model used in USGS PHREEQC. Agreement is within ΔSI = 0.32
at cryo conditions and 0.03 at physiological conditions — so the supersaturation
conclusion does not rest on one activity model.

We are also explicit about what the model does *not* do. It does not predict the
absolute deficit from first principles — that depends on the precipitated
fraction and deposit morphology, which is why we propose measuring them. And the
re-dissolution step is modelled as mass-transfer limited (Noyes-Whitney), which
is appropriate for an amorphous, undersaturated-on-redilution solid and is why
mixing helps; it is not a claim about a fixed equilibrium endpoint.

The whole model is at github.com/tohafrit/serum-ca-cryo. From a clean clone,
`make all` installs dependencies, runs 135 unit tests, and regenerates every
figure and data table in a few minutes. No proprietary software is needed. Each
physical step — cryo-concentration, saturation index, nucleation statistics,
suppressed ripening, re-dispersion kinetics — is a separately testable module.
You can change an assumption (cryoprotectant, pH, deposit size) and watch every
downstream number move.

One honest dependency: the cryoprotectant matters. A high-glycerol pool suppresses
ripening (our case); a low-viscosity pool (less glycerol, or a salt-only freeze)
would not. Telling us the actual cryoprotectant would sharpen the model
immediately.

## FIELD 4: EXPERIENCE

My background is physical chemistry. I hold a specialist degree from the
Mendeleev University of Chemical Technology in Moscow (Institute of Physical
Chemistry, 2008), specializing in oxide single-crystal growth. The diploma work
was on preparing laser-grade oxide crystals up to the pre-pressing stage:
nucleation from supersaturated solution, control of growth rate, and phase
verification by UV-Vis, IR, ICP and AAS. The physics in these vials —
nucleation from a supersaturated ionic solution, an amorphous-to-crystalline
pathway, dissolution kinetics — is the physics I worked with in the lab, with
different ions. Nucleation control and phase verification in a supersaturated
melt or solution are not adjacent to this problem — they are the same problem.

For the last several years I have worked as a senior DevOps engineer at a large
enterprise (AWS/GCP/Azure, Kubernetes, Terraform, CI). That is why this
submission ships as a reproducible repository with 135 unit tests and a one-command
build, rather than as assertions: I treat a model the way I treat production code
— every claim is checked and version-controlled. The same discipline applies to
a process change like a thaw-protocol SOP: validate it, stage it, and keep a
clear rollback criterion.

I am based in Israel and have access to Technion core facilities (DLS, ICP-MS,
electron microscopy) for the validation work proposed here. I am not a clinical
IVD formulator and have no GMP track record — which is partly why I leaned on
reproducibility instead of authority, and why I built the analysis so you can
check every number yourself.

## FIELD 5: SOLUTION RISKS

Mixing and shear-sensitive analytes. The main risk of a mixing step is shear on
fragile analytes. Vortexing can reduce some enzyme activities (e.g. LDH, CK) by
single-digit percentages in plasma; serum is generally less sensitive, but this
must be verified on your matrix and acceptance criteria. Mitigation: validate a
full analyte panel (at minimum LDH, CK, ALP, ALT, albumin, total protein, plus
any product-critical analytes) across mixing intensities, and choose the gentlest
mixing that recovers calcium — gentle inversion or a brief low-speed step may
suffice, since the deposit is amorphous and loosely bound.

Magnitude uncertainty. We are transparent that the absolute deficit is uncertain
(band ~0.5–15 %), because it scales with the precipitated fraction and deposit
size. This is not hidden in the model — it is an explicit parameter, and
Experiment 2 (ISE + ICP-MS) measures it directly. The decision-useful predictions
(reversible-by-mixing, grows with storage, vial-dependent) are robust to it.

Cryoprotectant assumption. We modelled 15 % glycerol. A different cryoprotectant
changes the pool viscosity and therefore whether any crystallization could occur
and how fast precipitation proceeds. The qualitative mechanism (supersaturation →
amorphous precipitate → reversible by mixing) holds across cryoprotectants; the
quantitative kinetics do not. Confirming the cryoprotectant removes this risk.

Glass/surface dependence. Nucleation is surface-catalysed, so vial-lot and
silanol-density variation drive the vial-to-vial spread (the dominant variance
term in our sensitivity analysis). This is consistent with the batch dependence
you see and is testable by comparing glass lots.

Required Seeker-specified factors. pH/ionic strength: freezing concentrates ions
and can raise pool pH; both revert fully on thaw. Osmolality: precipitation at
this scale does not measurably change it. Metal balance: total calcium is
conserved in the vial — only the *sampled* fraction changes, which is exactly why
mixing recovers it and why ISE+ICP-MS is the right check.

## FIELD 6: TIMELINE, CAPABILITY AND COSTS

Phase 1 — confirm and size the mechanism (about 4 weeks). On Seeker-provided
vials of different storage ages: DLS/NTA on fresh vs mixed aliquots; simultaneous
ISE + ICP-MS before/after mixing; a bottom/wall-vs-bulk calcium check. Deliverable:
yes/no on the mechanism and a measured precipitated fraction and deposit size.
A clear go/no-go before any further spend — if no particles appear in affected
vials, the mechanism is wrong and we stop.

Phase 2 — validate the chosen tier (about 6–8 weeks).
- If pursuing PREVENTION: a deep-freeze pilot — store matched lots at −20 °C vs
  ≤ −80 °C and compare particle counts and Ca recovery over time. Expect the
  −80 °C arm to show no deposit and no deficit. This directly tests the
  vitrification prevention.
- If pursuing NEUTRALIZATION: dose-response of mixing intensity vs calcium
  recovery, with a full analyte-stability panel at each setting, to find the
  gentlest step that recovers calcium within your acceptance criteria.
- Optional: degassing/sealing and controlled-rate freezing pilots to reduce
  precipitation at source.

Phase 3 — documentation (about 4 weeks). SOP, ISO 13485 change-control package,
and a short stability/verification study design.

Indicative total: roughly 14–16 weeks of instrument time and materials. The
primary fix (a mixing step) needs no capital equipment; upstream options do.
I work as a contractor from Israel, can run Phase 1 at Technion core facilities,
and am open to remote collaboration for the analysis phases and on-site visits
for the lab phases.

I am genuinely interested in seeing this through validation, not in a one-off
write-up.

## FIELD 10: ONLINE REFERENCES

Full model and code: github.com/tohafrit/serum-ca-cryo
`make all` reproduces all figures and 135 unit tests in a few minutes on a normal
computer; no proprietary software required.

1. Combes C, Rey C (2010). Amorphous calcium phosphates: synthesis, properties
   and uses in biomaterials. Acta Biomaterialia 6(9):3362–3378.
2. Carpenter JF, Crowe JH (1988). The mechanism of cryoprotection of proteins by
   solutes. Cryobiology 25(3):244–255.
3. Pikal MJ (1990). Freeze-drying of proteins, Part I: process design.
   Pharm Biotechnol 2:120–160.
4. Boskey AL, Posner AS (1973). Conversion of amorphous calcium phosphate to
   microcrystalline hydroxyapatite. J Phys Chem 77(19):2313–2317.
5. Fennema O (1973). Solid-liquid equilibria. In: Low Temperature Preservation of
   Foods and Living Matter. Marcel Dekker.
6. Cheng N-S (2008). Formula for the viscosity of a glycerol-water mixture.
   Ind Eng Chem Res 47(9):3285–3288.
7. Fogh-Andersen N et al. (1995). Ionic binding, net charge and Donnan effect of
   human serum albumin. Clin Chem 41(12):1522–1525.
8. Parkhurst DL, Appelo CAJ (1999). User's guide to PHREEQC (v2). USGS WRIR
   99-4259. (WATEQ activity model used for the cross-check.)
9. Davies CW (1962). Ion Association. Butterworths, London.
10. Sahai N (2005). pH-dependent mineral-water interfacial energetics: calcium
    phosphate heterogeneous nucleation on silica. Am J Sci 305:661–672.
