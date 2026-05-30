# serum-ca-cryo

Computational backbone for an InnoCentive submission:
**"Preventing post-thaw decrease in calcium concentration in serum-based QC standards"**

The model is a *scoping* tool: it establishes the thermodynamic root cause,
explains why the deficit is reversible, and shows which intervention works and
why. It is fully reproducible. The absolute deficit magnitude carries an explicit
uncertainty band (see "Precipitated fraction" below) — the proposed experiments
pin it down.

**The fix — two tiers, both process-only.** The calcium is not lost: it is
amorphous calcium phosphate on the vial surface that a quiescent draw
under-samples. **Prevent** it by deep-frozen / vitrified storage (≤ −80°C →
nucleation arrested, deposit never forms), or **neutralize** it with no new cold
chain by a defined mixing/re-suspension step at thaw.

---

## Headline results

| Result | Value |
|--------|-------|
| Cryo-concentration → supersaturation (glycerol 15%, −20°C) | k = 5.58; SI(HAp) ≈ **+7.5** (WATEQ-checked, ΔSI < 0.32) |
| Crystalline ripening at −20°C in the viscous pool | **suppressed** (~450× slower) → precipitate stays **amorphous** (fig04) |
| Deficit in an affected vial (quiescent thaw) | **~5%** (calibrated to the Seeker's reported ≥4%); band **~0.5–15%** (fig06) |
| Fraction of vials with a ≥4% deficit | small early (~7% at 1 mo) → **~48% at 6 mo → ~78% at 24 mo**, via stochastic nucleation onset (fig06/07) |
| Defined mixing step, or standing | deficit → **~0** — reproduces the Seeker's "reversible with mixing" (fig09) |
| **Deep-frozen storage (≤ −80°C)** | nucleation arrested → affected fraction → **~0** — root-cause **prevention** (fig09) |

The key insight: the precipitate is **amorphous**, not crystalline — and that is
*why* simple mixing reverses it. We show crystalline hydroxyapatite cannot form
at −20°C in the viscous cryoprotected pool, so an easily re-dispersible amorphous
deposit is what remains. The reversibility is the fingerprint of an amorphous
phase.

---

## Mechanism (textual schematic)

```
FREEZING
  → cryoconcentration (k=5.58 for glycerol 15%): [Ca] 2.75→15.4 mM, [Pi] 1.8→10.1 mM
  → SI(HAp) ≈ +7.5 (WATEQ-checked) — serum forced past its metastable inhibitors

NUCLEATION (stochastic; glass surface-catalysed)
  → amorphous calcium phosphate (ACP) precipitates on the glass surface
  → nucleation delay ~ Lognormal(median ≈ 90 days) per vial, pH/supersaturation-
    dependent → only SOME vials carry a deposit at a given time (vial/batch spread)
  → Sobol: nucleation delay ST = 0.71 (dominant), glass site density ST = 0.46

FROZEN STORAGE (months, −20°C)
  → pool is ~84% glycerol, η ≈ 4100 mPa·s → solution-mediated ripening is
    ~450× slower → ACP → HAp would take centuries → precipitate STAYS AMORPHOUS
  → individual particles stay ~tens of nm; they aggregate into ~µm wall deposits
  → consistent with Combes & Rey (2010): ACP is kinetically stable below 0°C

THAW, QUIESCENT DRAW (22°C, 60 min, no mixing)
  → the surface-bound amorphous deposit is under-sampled by a bulk draw
  → analyzer reads low: deficit = F_PRECIP × (fraction not re-dispersed) ~ a few %

THAW WITH MIXING (or extended standing)
  → mixing thins the diffusion boundary layer (~10 µm → 1–2 µm) → fast re-dispersion
  → the amorphous deposit re-disperses/re-dissolves → deficit → 0 (reversible)
```

Three observations, one mechanism:

- **Reversibility with mixing** — the calcium is not chemically lost; an amorphous,
  surface-bound deposit re-disperses on mixing. A crystalline phase would not —
  so reversibility itself tells us the deposit is amorphous.
- **In some samples / batch-to-batch** — nucleation is stochastic and surface-
  catalysed, so only some vials carry an appreciable deposit; the affected
  fraction grows with storage (Sobol ST(nucleation) = 0.71).
- **Develops over storage** — the *affected fraction* (not the per-vial amount)
  grows as more vials cross their nucleation induction time.

---

## Precipitated fraction (the honest uncertainty)

The deficit (as % of total Ca) = `F_PRECIP × (1 − recovery)`, where `F_PRECIP`
is the fraction of total calcium that precipitates at peak cryoconcentration.
This is bounded by **mass balance**, not by albumin binding:

- phosphate stoichiometry and the extreme supersaturation (SI(HAp) ≈ +7.5)
  permit up to ~0.97 of the calcium to precipitate;
- albumin buffering of free Ca²⁺ lowers the actual figure, with a wide plausible
  band of **~0.07 … 0.97**.

The model adopts a representative `F_PRECIP = 0.90` and reports the deficit
**magnitude** as uncertain to this factor. The predictions that are *robust* to
`F_PRECIP` are the onset timescale and the relative effect of mixing. `F_PRECIP`
is exactly what the proposed simultaneous ISE + ICP-MS experiment measures.

The re-dispersion step is modelled as sink-limited Noyes-Whitney (the
back-reaction term is dropped); this is appropriate for the amorphous deposit,
which is undersaturated once the vial is re-diluted and so fully re-disperses
given mixing or time. See `src/ripening_kinetics.py` for the full caveat.

---

## Repository structure

```
src/
  freezing_trajectory.py   Module 2: cryoconcentration & pH trajectory
  saturation_indices.py    Module 3: SI for 8 mineral phases with albumin binding
  supersaturation_map.py   Module 4: SI heatmap over k × pH grid
  ripening_kinetics.py     Module 5: ripening (shown suppressed) + re-dispersion kinetics
  vial_simulation.py       Module 6: Monte Carlo (10,000 vials, 50 batches)
  interventions.py         Module 7: intervention scenarios (sealed baseline)
  phreeqc_runner.py        Module 8: Davies vs WATEQ thermodynamic cross-check
  plot_fig*.py             Figure generators

tests/                     135 tests, 7 modules
data/                      7 CSV tables (auto-generated by `make figures`)
figures/                   11 figures + 4 supplementary (auto-generated)
phreeqc/                   WATEQ parameter database (documentation)
```

## Intervention scenarios (Module 7)

Metric: **fraction of vials at or above the Seeker's 4% threshold at 6 months**.
Interventions act either on the *affected fraction* (via nucleation) or on
*recovery* (via mixing). Per affected vial the deficit is ~5% (band 0.5–15%).

| Scenario | What it changes | Frac ≥4% @6mo | Complexity |
|----------|-----------------|----------------|------------|
| baseline (sealed) | — | 0.48 | — |
| loose_seal (risk) | CO₂ vents, pH→8.0 → more nucleation | 0.72 | — |
| **+crf_2C** | **degassed + controlled freezing (PREVENT, S1)** | 0.22 | manufacturer |
| **+vortex_30s** | **thaw mixing only (NEUTRALIZE, S2)** | **~0.00** | **SOP only** |
| +vortex_60s | longer mixing | ~0.00 | SOP only |
| +combined | CRF + mixing | ~0.00 | high |
| +combined_plus | double-pulse mixing + 90 min | ~0.00 | high |
| +extended_mixing | long standing / cold soak | ~0.00 | none (slow) |
| +deep_freeze | store ≤ −80°C (full prevent; cold-chain trade-off) | ~0.00 | cold chain |

Two tiers: **deep-freeze prevents** the deposit from forming (root cause);
**a mixing step neutralizes** it with no new cold chain. Degassing/tight sealing
and controlled freezing reduce the affected fraction upstream.

---

## Reproduce everything

```bash
git clone https://github.com/tohafrit/serum-ca-cryo
cd serum-ca-cryo
make all          # installs deps, runs 135 tests, regenerates all figures + tables
                  # a few minutes on a modern machine (after deps cached)
```

`make all` regenerates **every** data table and figure from source, including
the Module 5 and Module 7 CSVs that the headline numbers come from. Individual
targets: `make setup`, `make test`, `make figures`, `make clean`.

---

## How this maps to the InnoCentive challenge

| Seeker requirement | Model response | Reference |
|--------------------|---------------|-----------|
| Mechanistic explanation | Cryoconcentration → CaP supersaturation → amorphous CaP on glass → under-sampled at quiescent thaw | fig01–fig04 |
| Reversible with mixing | Amorphous (not crystalline) surface deposit re-disperses on mixing → deficit → 0 | fig09 (mixing scenarios) |
| Vial-to-vial / batch variability | Stochastic surface nucleation; Sobol ST(delay)=0.71, ST(glass)=0.46 | fig06, fig07, figS4 |
| Intervention to prevent / neutralize | Prevent: deep-freeze ≤ −80°C (nucleation arrested). Neutralize: defined mixing step at thaw. Both take affected fraction 0.71 → ~0 | fig09, fig10 |
| Formulation neutrality | Process changes only; no excipient added/removed | test_formulation_chemistry_* |
| ISO 13485 / REACH / PFAS | No new chemical substance introduced | PROPOSAL_NOTES.md §Part 2 |
| Independent thermodynamic check | Davies vs WATEQ: max ΔSI = 0.32; physiological ΔSI = 0.03 | fig11 |

---

## Dependencies

- Python ≥ 3.11; numpy, scipy, matplotlib, pandas, SALib
- No PHREEQC binary required — Module 8 implements the WATEQ extended
  Debye-Hückel activity model in pure Python.

---

## Bibliography

1. Boskey AL & Posner AS (1973) "Conversion of amorphous calcium phosphate to microcrystalline hydroxyapatite." *J Phys Chem* 77:2313–2317.
2. Eanes ED & Posner AS (1965) "Kinetics and mechanism of conversion of noncrystalline calcium phosphate to crystalline hydroxyapatite." *Trans NY Acad Sci* 28:233–241.
3. Christoffersen J et al. (1989) "Kinetics of dissolution of calcium hydroxyapatite." *J Crystal Growth* 94:767–777.
4. Heughebaert JC & Nancollas GH (1984) "Kinetics of crystallization of octacalcium phosphate." *J Phys Chem* 88:2478–2481.
5. Combes C & Rey C (2010) "Amorphous calcium phosphates: synthesis, properties and uses in biomaterials." *Acta Biomater* 6:3362–3378.
6. Fogh-Andersen N et al. (1995) "Ionic binding, net charge, and Donnan effect of human serum albumin." *Clin Chem* 41:1522–1525.
7. Pedersen KO (1972) "Binding of calcium to serum albumin: optimized normal reference values." *Scand J Clin Lab Invest* 30:321–327.
8. Davies CW (1962) *Ion Association*. Butterworths, London.
9. Plummer LN & Busenberg E (1982) "Solubilities of calcite, aragonite and vaterite in CO₂-H₂O solutions." *Geochim Cosmochim Acta* 46:1011–1040.
10. Carpenter JF & Crowe JH (1988) "The mechanism of cryoprotection of proteins by solutes." *Cryobiology* 25:244–255.
11. Pikal MJ (1990) "Freeze-drying of proteins: process, formulation, and stability." ACS Symposium Series 567.
12. Fennema OR (1973) *Low Temperature Preservation of Foods and Living Matter*. Marcel Dekker.
