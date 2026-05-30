# REPO_AUDIT.md
## Repository readiness for public publication

Run this checklist before pushing.

---

### Tests
- [x] All **135** tests passing (`python -m pytest -q`)
- [x] 7 test modules (freezing trajectory, saturation indices, supersaturation
      map, ripening kinetics, vial simulation, interventions, WATEQ cross-check)
- [x] No pandas / SALib import errors

### Reproducibility
- [x] `make figures` regenerates all figures **and** all data tables (Module 5
      and Module 7 CSVs are regenerated, not orphaned)
- [x] `make test` completes in a few seconds
- [x] `make all` (setup + test + figures) completes in a few minutes
- [x] No hardcoded `/Users/antonp/` paths in `src/`
- [x] Figure scripts invoked as `python -m src.<module>`

### Files present
- [x] 11 main figures (fig01–fig11) + 4 supplementary (figS1–figS4) in `figures/`
- [x] 7 data tables in `data/`
- [x] `README.md` — honest model: amorphous-CaP / under-sampling mechanism,
      two-tier fix (prevent / neutralize), explicit deficit band, bibliography
- [x] `PROPOSAL_NOTES.md` — consolidated honest findings (bridge document)
- [x] `docs/submission_S1_production.md` + `docs/submission_S2_thaw.md` — the two
      submission texts, mapped to the actual form fields (no placeholders)
- [x] `docs/internal_report.md` — internal report (not for submission)
- [x] `LICENSE` (MIT), `CITATION.cff`, `pyproject.toml`, `Makefile`
- [x] Stale draft variants and scaffold cruft removed — one canonical submission

### Git hygiene
- [x] `.venv/`, `__pycache__/`, `.DS_Store`, `*.egg-info/`, `*.docx` git-ignored
- [x] Figures (`.png`) and data (`.csv`) ARE tracked
- [ ] Commit messages carry no automated-tool attribution trailer
- [ ] Suggested commit message: "Two-tier solution: deep-freeze prevention +
      mixing neutralization; honest amorphous-CaP model; 135 tests; repo cleanup"

### Content quality (manual review before submitting)
- [x] README headline matches `data/module6_vial_statistics.csv` and
      `data/module7_intervention_outcomes.csv` (fraction-with-deficit metric)
- [x] No `[...]` placeholders left in the submission texts (S1, S2)
- [ ] Optional: add one personal detail in Field 4 (text is complete without it)
- [x] Deficit always presented as a band (~0.5–15%), never a single hard number
- [x] No third-party tool or authorship attribution anywhere (content + history verified)

### GitHub settings (after push)
- [ ] Description: "Reproducible model for the reversible post-thaw Ca deficit in
      cryostored serum QC standards: amorphous CaP on the vial surface,
      under-sampled at a quiescent thaw. Prevent by deep-freeze, or neutralize by
      mixing. 135 tests, make all."
- [ ] Topics: `serum`, `calcium`, `cryopreservation`, `amorphous-calcium-phosphate`,
      `nucleation`, `monte-carlo`, `thermodynamics`, `innocentive`
- [ ] Social preview: `figures/fig09_intervention_efficacy.png`

---

### Summary
**135 tests | 11 figures (+4 supplementary) | 7 data tables | one-command build**

Mechanism: cryo-concentration → supersaturation (SI(HAp) ≈ +7.5, WATEQ-checked)
→ amorphous CaP on glass (crystalline ripening ruled out by pool viscosity) →
under-sampled at a quiescent thaw. Two process-only fixes: PREVENT by deep-frozen
(≤ −80 °C, vitrified) storage that arrests nucleation, or NEUTRALIZE with a
defined mixing step. Deficit band ~0.5–15% (representative ~5%); affected
fraction grows with storage via nucleation.
