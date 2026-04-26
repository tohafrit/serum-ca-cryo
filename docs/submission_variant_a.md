# Submission Text — Variant A
## Full mechanism + Combined+ protocol
## InnoCentive: Post-thaw Ca deficit in serum QC standards

This is the primary submission. Lead with the complete mechanistic explanation
and the combined intervention.

Figures to attach: fig04, fig05, fig06, fig09, fig11 (five files, one caption each).

---

## FIELD 1: PROBLEM AND OPPORTUNITY

Clinical laboratories use serum-based QC standards to verify that their
analyzers give correct results. Calcium is one of the most critical analytes,
because ionized calcium is used in diagnosis of parathyroid disorders, renal
failure, and cardiac conditions. When the QC standard shows a calcium reading
that is too low, it is hard to know if the analyzer is wrong or the standard
is wrong.

The deficit appears only after more than 6 months of storage at -20°C. Vials
stored for shorter times are fine. The effect is different from vial to vial,
and sometimes from batch to batch. And if you leave the vial at 2-8°C for
48 hours, the calcium reading comes back to normal. These three facts together
are puzzling. A simple chemical explanation like binding or protein degradation
does not account for the threshold at 6 months or for the variability between
vials from the same batch.

The proposed mechanism is cryoconcentration-driven calcium phosphate
precipitation, followed by slow crystallization during storage, followed by
incomplete dissolution at thaw. The key insight is that freezing concentrates
all solutes in the remaining liquid, increasing the calcium and phosphate ion
product by about 31 times. This causes amorphous calcium phosphate to form on
glass surfaces. Over 6 or more months, it converts to hydroxyapatite, which
is much harder to dissolve. A standard 60-minute thaw at room temperature is
not enough time for hydroxyapatite crystals to go back into solution.

The model was built entirely from published physical chemistry data, without
fitting any parameters to the Seeker's observations. It predicts a 5.4%
calcium deficit at 6 months, matching the reported threshold of 4% or more.
It also predicts that 42% of vials will show the effect, which is consistent
with "in some samples." The 48-hour cold equilibration workaround that the
Seeker already uses is independently reproduced by the model.

The full analysis is available at github.com/tohafrit/serum-ca-cryo. Running
`make all` reproduces all figures and all 139 tests in under 2 minutes on any
modern computer. The thermodynamic calculations were cross-checked against the
WATEQ database, which is the standard used by USGS in PHREEQC.

---

## FIELD 2: SOLUTION OVERVIEW

### Part 1 — Mechanism

The deficit is caused by a four-step process.

First, freezing concentrates ions in the unfrozen liquid. For glycerol 15% w/w,
the concentration factor at -20°C is k = 5.58. Calcium goes from 2.75 mM to
15.4 mM. Phosphate goes from 1.8 mM to 10.1 mM. The ion activity product for
hydroxyapatite increases 31 times. The saturation index reaches +7.5 log units,
compared to +5.1 in physiological serum. CO2 also escapes during freezing,
which raises local pH from 7.4 to about 7.8-8.8.

Second, amorphous calcium phosphate (ACP) nucleates on glass surface
heterogeneities. The nucleation delay varies from vial to vial, modeled from
glass surface chemistry data. This variability is the main reason some vials
show the deficit and others do not. Sobol sensitivity analysis shows that
nucleation delay explains 71% of the variance in outcome.

Third, during months of storage at -20°C, ACP transforms to octacalcium
phosphate (OCP) and then to hydroxyapatite (HAp). This is Ostwald ripening,
the same process that forms bone and dental enamel. At 6 months, the
precipitated calcium is roughly 52% ACP, 41% OCP, and 7% HAp.

Fourth, at thaw, ACP and OCP dissolve quickly. HAp crystals, which have grown
to about 200-400 nm in radius, dissolve much more slowly. In a quiescent 60-min
thaw at 22°C, only about 35% of HAp dissolves. The rest stays in solution as
tiny crystals, not registered as calcium by the analyzer. This gives 5.4%
deficit at 6 months.

### Part 2 — Interventions

Three process changes, no formulation modification.

Pre-freeze degassing to 10% residual CO2 reduces the pH rise during freezing,
which lowers supersaturation and delays nucleation. Controlled-rate freezing
at 2°C/min creates more uniform cryoconcentration and fewer nucleation sites.
A double-pulse vortex thaw (30 seconds at 5 minutes, 60 seconds at 25 minutes,
90-minute total window) reduces the diffusion layer around each crystal from
about 10 micrometers to about 1 micrometer, increasing dissolution rate by
about 10 times.

The combined protocol of all three gives 0.0% of vials above the 4% threshold
at 12 months, in 90 minutes of total thaw time. The Seeker's existing 48-hour
cold equilibration achieves the same result but takes 32 times longer.

To validate the mechanism before deploying the intervention, three experiments
are possible: (1) DLS or NTA particle counting on freshly thawed vs
equilibrated vials (predict 50-500 nm particles in fresh, none in equilibrated),
(2) simultaneous ionic Ca by ISE and total Ca by ICP-MS (proportional drop
confirms precipitation, not binding), (3) cryo-SEM of frozen samples at 1
month vs 6 months (predict visible crystal growth).

---

## FIELD 3: SOLUTION FEASIBILITY

The model builds on published physical chemistry that has been validated for
decades in the context of bone and dental mineralization. The same calcium
phosphate phases (ACP, OCP, hydroxyapatite) form in bone via the same
transformation pathway. The same Noyes-Whitney dissolution kinetics apply.
This is not a novel theory. It is an established framework applied to a new
context.

The key references for the kinetic parameters are Boskey and Posner (1973) for
ACP-to-HAp transformation rates, Heughebaert and Nancollas (1984) for OCP
crystallization, and Christoffersen et al. (1989) for hydroxyapatite dissolution
kinetics. For cryoconcentration in pharmaceutical systems, Carpenter and Crowe
(1988) and Pikal (1990) established that the same concentration mechanism
operates in freeze-drying. The review by Combes and Rey (2010) in Acta
Biomaterialia covers the full ACP-to-HAp transformation literature.

The thermodynamic model was cross-checked against WATEQ extended Debye-Huckel,
which is the activity model used in PHREEQC, the standard geochemical
speciation code from USGS. At physiological ionic strength, the two models
agree to within delta-SI = 0.03 (three decimal places). At the cryoconcentrated
state (ionic strength 0.8 mol/kg, at the edge of validity for Davies equation),
the difference is delta-SI = 0.32. This is smaller by one order of magnitude
than the uncertainty from the choice of ACP solubility constant (Boskey 1973
vs Christoffersen 1990 endpoint: delta-SI = 2.8). So the activity model is
not the dominant source of uncertainty.

The complete computational model is open-source at
github.com/tohafrit/serum-ca-cryo. It includes 139 passing unit tests,
11 figures, and 6 data tables. Running `make all` from a fresh clone reproduces
every result in under 2 minutes. No proprietary software is required.

---

## FIELD 4: EXPERIENCE

I have a degree in physical chemistry from Mendeleev University of Chemical
Technology, Moscow (Institute of Physical Chemistry). My focus was solution
thermodynamics and kinetics, which is exactly the relevant background here.

Currently I work as a DevOps engineer, based in Israel. I have been working
in software infrastructure for several years, and physical chemistry is not my
daily job. I know this is an unusual combination. I will explain why I think
it is an asset in this case.

Formulation scientists in the IVD industry often work with empirical approaches
because they are fast and reliable for product development. This problem seems
to have been treated empirically too: the 48-hour workaround was found by
observation, not by understanding the chemistry. A fresh look from first
principles, without the assumption that the cause must be something specific
to serum proteins or the analyzer, is what allowed me to identify cryoconcentration
as the likely cause.

The model has 139 unit tests that verify each physical assumption independently.
The Sobol sensitivity analysis identifies which parameters actually matter
(nucleation delay and glass surface area dominate; most other parameters have
low sensitivity). The thermodynamic cross-check with WATEQ is there to make
sure the activity model choice does not change the conclusions.

The Seeker can clone the repository and run `make all`. Every number in this
submission can be reproduced that way. I think transparent, reproducible
computation is a reasonable substitute for a wet-lab track record in this case.
The model predicts the Seeker's own workaround. That is the validation.

---

## FIELD 5: SOLUTION RISKS

The main model risks and how I handled them:

**ACP solubility.** The solubility of amorphous calcium phosphate has two
commonly used values in the literature: Boskey and Posner 1973, and
Christoffersen 1990. They differ by about 2.8 log-units in SI. I do not
claim to know which one is correct for the specific ACP that forms in these
vials. The saturation index for HAp at the cryoconcentrated state is reported
as a range (+6.1 to +9.0), not a point estimate. The qualitative conclusion
(HAp is supersaturated and will form) is the same for any value in that range.

**Activity model at high ionic strength.** Davies equation is standard for
ionic strength below 0.5 mol/L. The cryoconcentrated pool has I = 0.8 mol/kg.
I cross-checked the Davies results with WATEQ extended Debye-Huckel across
all four phases at both conditions. Maximum difference was 0.32 log-units. This
is inside the uncertainty from Ksp choice, so the activity model is not the
limiting factor.

**Unknown cryoprotectant.** I used glycerol 15% w/w as the assumed
cryoprotectant, because it is the most common for serum standards. If it is
DMSO 10%, the concentration factor changes to k=6.8 (higher, so more
supersaturated). If there is no cryoprotectant, k=35.5 (much higher). All
three scenarios predict the same problem and the same interventions work.

**Vortex and shear-sensitive analytes.** Vortexing at 1500 rpm for 30-90
seconds can reduce LDH and CK activity by 5-15% based on published data for
plasma. Serum is less sensitive, but this needs to be confirmed. Before
deploying the vortex protocol, a validation panel of at least 12 analytes
(LDH, CK, ALP, ALT, albumin, total protein, and others critical for the
specific product) must be measured before and after vortex. If the reduction
is unacceptable, the protocol can be modified to 500 rpm for 120 seconds,
which is less effective but likely safer.

**Glass lot variability.** The Sobol analysis shows that glass surface site
density explains 46% of vial-to-vial variance (second after nucleation delay).
Different glass suppliers and different manufacturing lots may have different
surface silanol densities, which would change the nucleation rate. I did not
have data on the specific glass used in the Seeker's vials.

---

## FIELD 6: TIMELINE, CAPABILITY AND COSTS

The proposed work has three phases. Phase 1 validates the mechanism. Phase 2
validates the intervention. Phase 3 produces documentation.

**Phase 1 — Mechanism validation (4 weeks, $8,000-12,000)**

The Seeker needs to provide vials at 1, 3, 6, and 12 months of storage, both
freshly thawed and 48-hour equilibrated. Three measurements:

DLS or NTA particle counting (Malvern Zetasizer or NanoSight at university
core facility, about $500/day). This will confirm or rule out particles in the
50-500 nm size range predicted for HAp crystals.

Simultaneous ionic Ca by ISE and total Ca by ICP-MS. This distinguishes
precipitation from binding. The ISE setup costs about $200. ICP-MS at a core
facility is about $100 per sample.

Cryo-SEM of vitrified samples at 1 month and 6 months. Predict particle growth
from under 50 nm (ACP) to over 200 nm (HAp). Core facility session about $300.

**Phase 2 — Intervention validation (8 weeks, $15,000-25,000)**

Vortex dose-response study: three rotation speeds (500, 1000, 1500 rpm), three
pulse durations, across three storage durations. Full analyte panel before and
after vortex. This closes the shear-sensitivity risk.

Controlled-rate freeze pilot: compare standard uncontrolled freezing with
1°C/min and 2°C/min profiles. Measure deficit distribution at 6 months
post-freeze for each condition.

Degassing pilot: measure residual CO2 by infrared spectroscopy before and after
vacuum degassing. Confirm target of 10% residual.

**Phase 3 — Documentation (4 weeks, $5,000-8,000)**

SOP writing for the combined protocol (degassing, controlled-rate freeze, vortex
thaw, 90-minute window). ISO 13485 traceability documentation. Stability study
design for the modified protocol.

Total: 16 weeks, estimated $30,000-45,000 for instrument time and materials.
This does not include capital equipment for a controlled-rate freezer, if the
Seeker does not have one. I work as a contractor from Israel and can travel
for on-site work if the Seeker's laboratory is in Europe or North America.

---

*Editing notes:*
- *Field 4: add one personal sentence about a specific project or experience*
- *Field 5 vortex section: if you have internal data on LDH stability, add it*
- *All fields: read aloud, replace any sentence that sounds like a presentation with something more direct*
- *Check all numbers against data/module7_intervention_outcomes.csv before submitting*
