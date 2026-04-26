# serum-ca-cryo

Computational backbone for an InnoCentive submission:  
**"Preventing post-thaw decrease in calcium concentration in serum-based QC standards"**

## The phenomenon

Human-serum quality-control standards stored at −20 °C for >6 months show a reversible post-thaw
calcium deficit of ≥4%. The deficit disappears after mixing or 24–48 h equilibration at 2–8 °C,
and is absent in samples stored <6 months.

## Proposed mechanism (tested here)

```
Freezing
  → cryoconcentration (unfrozen interstitial pool reaches k×10–50 normal concentration)
  → local supersaturation in Ca–phosphate and Ca–carbonate systems
  → nucleation of amorphous calcium phosphate (ACP) and possibly CaCO₃ phases

Storage at −20 °C (6+ months)
  → Ostwald ripening: ACP → octacalcium phosphate (OCP) → hydroxyapatite (HAp)
  → particle growth into microcrystals (slow dissolution kinetics)

Thaw
  → ACP (formed early) dissolves in minutes → no deficit if <6 months
  → HAp microcrystals dissolve in days → deficit persists unless mixed/equilibrated
```

Secondary: CO₂ outgassing during freeze raises local pH, promoting CaCO₃ co-precipitation.

## Modules

| Module | File | Figure |
|--------|------|--------|
| 1 | Project scaffold | — |
| 2 | Cryoconcentration trajectory | fig01 |
| 3 | Saturation indices | fig02 |
| 4 | Supersaturation map | fig03 |
| 5 | Ostwald ripening kinetics | fig04, fig05 |
| 6 | Monte Carlo vial variability | fig06 |
| 7 | Intervention modeling | fig07 |
| 8 | PHREEQC validation (optional) | fig03 overlay |
| 9 | Final assembly | all |

## Reproduce all figures

```bash
# 1. Clone and set up
git clone <repo>
cd serum-ca-cryo
make setup          # creates .venv and installs deps

# 2. Run tests
make test

# 3. Generate all figures (~5–8 min on a modern Mac)
make figures

# 4. Open notebooks for interactive exploration
make notebooks
```

## Dependencies

- Python ≥ 3.11
- numpy, scipy, matplotlib, pandas, SALib
- Optional: phreeqpy + PHREEQC binary (`brew install phreeqc`) for Module 8 validation

## Bibliography

Full citations are inline in each source file. Key references:

- Ksp for HAp: Dorozhkin SV (2002) *J Mater Sci* 37:4871–4880
- ACP solubility: Boskey AL & Posner AS (1973) *J Phys Chem* 77:2313–2317
- Ostwald ripening rates: Nancollas GH & Mohan MS (1970) *Arch Oral Biol* 15:731–745
- Cryoconcentration: Fennema OR (1973) *Low Temperature Preservation of Foods and Living Matter*
- Albumin–Ca binding: Fogh-Andersen N et al. (1995) *Clin Chem* 41:1522–1525
- Davies activity correction: Davies CW (1962) *Ion Association*, Butterworths
