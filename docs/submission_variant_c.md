# Submission Text — Variant C
## Mechanism + glass supplier specification (long-term root cause)
## InnoCentive: Post-thaw Ca deficit in serum QC standards

Positioning: long-term root cause solution. Complements Variant B. Addresses
the glass surface as the upstream cause of nucleation variability. No protocol
change needed if the glass specification is right.

Sobol ST(glass surface site density) = 0.46. This is the second most important
parameter after nucleation delay, and it is the one the manufacturer can control
by specifying the glass.

Figures to attach: figS4 (Sobol chart), fig06, fig07, fig04, fig09.

---

## FIELD 1: PROBLEM AND OPPORTUNITY

Serum QC standards for clinical analyzers show a reversible calcium deficit
after more than 6 months of storage at -20°C. The deficit appears in some vials,
not all. This variability between vials from the same batch has no obvious
explanation if you assume the vials are identical.

They are not identical. Glass vials differ in surface chemistry from piece to
piece, from lot to lot, and from supplier to supplier. The density of silanol
groups on the inner surface of a borosilicate vial determines how many sites
are available for heterogeneous nucleation. More nucleation sites means earlier
and more extensive calcium phosphate precipitation. Fewer sites means less
precipitation and a lower chance of seeing the 4% deficit.

The proposed mechanism is cryoconcentration driving calcium phosphate
supersaturation, followed by glass-surface nucleation of amorphous calcium
phosphate, followed by months of crystallization to hydroxyapatite, followed
by incomplete dissolution in a standard 60-minute thaw. The 6-month threshold
and the reversal with extended equilibration are both explained by the
kinetics of this pathway.

The opportunity here is upstream control. Instead of correcting for the problem
after it forms, specify the glass to prevent it from forming. Low-silanol glass
or siliconized borosilicate reduces the nucleation rate and shifts the
distribution of affected vials toward zero. This is a one-time specification
change with the vial manufacturer, not an ongoing protocol change in the lab.

---

## FIELD 2: SOLUTION OVERVIEW

### Part 1 — Mechanism

*(Same core mechanism as Variant A: k=5.58, IAP x31, SI(HAp)=+7.5, ACP->OCP->HAp,
5.4% deficit at 6 months. Write in same plain style.)*

### Part 2 — Root cause: glass surface chemistry

Sobol sensitivity analysis (Monte Carlo, 10,000 vials) shows that nucleation
delay is the single dominant source of vial-to-vial variability, explaining 71%
of variance in calcium deficit. Glass surface site density is the second
parameter, explaining 46%.

Nucleation delay depends on two things: the supersaturation driving force (set
by the freezing and formulation) and the density of heterogeneous nucleation
sites (set by the glass surface). You can reduce supersaturation by changing
the freezing protocol, or you can reduce nucleation site density by changing
the glass specification.

Silanol groups (Si-OH) on the inner surface of borosilicate glass are the
primary nucleation sites for calcium phosphate. High-silanol glass nucleates
more ACP, earlier, in more vials. Low-silanol glass or siliconized glass
(treated with dimethyldichlorosilane or similar) has far fewer active sites.
The literature on heterogeneous nucleation in biomineralization (Combes and
Rey 2010) shows that surface chemistry differences of less than one order of
magnitude in site density can change nucleation times by months.

Practical recommendation: specify silanol density in the vial purchase
specification. Request Type I borosilicate with controlled surface silanol
density below a defined threshold. Ask the supplier for surface hydroxyl density
data (measured by deuterium exchange or XPS). Alternatively, specify
siliconized glass. This is already done in some pharmaceutical packaging
applications for the same reason.

The Sobol model suggests that reducing glass surface site density by 5x would
shift the median nucleation delay from 70 days to over 200 days, which would
cut the fraction of vials showing deficit at 6 months roughly in half, with no
protocol change at all.

Three experiments to validate: DLS/NTA particle counting on vials from different
glass lots (predict correlation between silanol density and particle count),
ISE+ICP-MS to confirm precipitation rather than binding, cryo-SEM to observe
nucleation density differences between glass types.

---

## FIELD 3: SOLUTION FEASIBILITY

*(Same core references as Variant A.)*

The glass surface chemistry aspect has strong precedent. Heterogeneous
nucleation of calcium phosphate on silica surfaces is well documented in
biomineralization research (Sahai 2005, Journal of Colloid and Interface
Science; Combes and Rey 2010). The effect of surface silanol density on
nucleation kinetics in pharmaceutical systems is covered in the freeze-drying
formulation literature (Pikal 1990; Nail and Jiang 2002).

Siliconized glass is already a standard specification option for pharmaceutical
vials. USP Type I borosilicate glass with silicone coating is available from
multiple suppliers and is routinely used for products sensitive to glass
surface interactions.

The model is at github.com/tohafrit/serum-ca-cryo, fully reproducible with
`make all`. The Sobol sensitivity analysis (figS4 in the repository) shows
the glass surface site density parameter clearly.

---

## FIELD 4: EXPERIENCE

*(Same as Variant A — see submission_variant_a.md Field 4)*

---

## FIELD 5: SOLUTION RISKS

**Glass specification change.** Changing the vial supplier or the glass
specification requires validation under ISO 13485. At minimum: qualification
study comparing calcium recovery from old glass vs new glass after 6 and 12
months of storage. This adds 6-12 months to the implementation timeline.

**Model uncertainty for glass parameter.** The glass surface site density
parameter in the model is a literature estimate, not a measurement on the
specific vials the Seeker uses. The Sobol index of 0.46 was computed for a
range based on published silanol density data for borosilicate glass. The
actual sensitivity may be different if the vials are siliconized or if a
different glass type is used.

**Supplier qualification.** Not all glass suppliers provide silanol density
data. Getting this specification into a purchase order and having it measured
reproducibly adds complexity to the supplier relationship.

**This solution does not give immediate results.** It requires at least one
storage cycle (6-12 months) to confirm improvement. For immediate relief,
Variant B (vortex protocol) can be implemented in parallel.

The recommended path is: Variant B first (fast, no capital, 4-week validation),
Variant C second (glass spec, longer lead time, permanent fix). They are
complementary, not alternatives.

---

## FIELD 6: TIMELINE, CAPABILITY AND COSTS

**Phase 1 — Glass characterization (4 weeks, $5,000-8,000)**

Obtain silanol density measurements for current vial lot (XPS or deuterium
exchange, core facility). Obtain 2-3 alternative glass types from suppliers
(Type I borosilicate standard, Type I siliconized, one low-silanol alternative).
Measure nucleation delay distribution in each type using the DLS/NTA assay
from the mechanism validation.

**Phase 2 — Comparative storage study (6-12 months, $10,000-20,000)**

Store samples in each glass type at -20°C. Measure calcium deficit distribution
at 3, 6, and 12 months. This is the key evidence that the glass specification
matters.

This phase cannot be shortened. It requires actual storage time. This is why
implementing Variant B in parallel is recommended — it gives immediate
improvement while this study runs.

**Phase 3 — Supplier qualification and SOP (4 weeks, $5,000-8,000)**

Update vial specification with silanol density requirement. Qualify one or two
suppliers to the new spec. Document the change under ISO 13485.

Total timeline: 8-14 months (dominated by storage study). Total cost:
$20,000-36,000. This does not include the comparative storage study consumables,
which depend on the Seeker's product volume.

---

*Editing notes:*
- *Field 2 Part 1 mechanism: copy from Variant A and shorten to 200 words*
- *Field 2 Part 2: this is the differentiating section — keep it detailed*
- *Field 6: be clear that this is a long-term solution; do not oversell timeline*
