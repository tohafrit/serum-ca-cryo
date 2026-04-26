# REPO_AUDIT.md
## Repository readiness for public publication

Run this checklist manually before `git remote add origin && git push`.

---

### Tests

- [x] All 139 tests passing (`python -m pytest tests/ -v`)
- [x] Module 1 (scaffold): 6 tests
- [x] Module 2 (freezing trajectory): 6 tests
- [x] Module 3 (saturation indices): 9 tests
- [x] Module 4 (supersaturation map): included in Module 3 test suite
- [x] Module 5 (ripening kinetics): 23 tests
- [x] Module 6 (vial simulation): 16 tests
- [x] Module 7 (interventions): 16 tests
- [x] Module 8 (WATEQ cross-check): 13 tests
- [x] No pandas/SALib import errors (both installed)

### Reproducibility

- [x] `make figures` completes cleanly in ~13 s on Apple Silicon
- [x] `make test` completes cleanly in ~3 s
- [x] `make all` (setup + test + figures) completes in < 2 min after pip cache warm
- [x] No hardcoded `/Users/antonp/` paths in any `src/` file
- [x] All figure scripts invoked as `python -m src.<module>` (not `python src/<module>.py`)
- [x] Makefile uses `python` (active env) not `.venv/bin/python` (venv-dependent path)

### Files present

- [x] 11 figures in `figures/`: fig01–fig11 + figS1–figS4
- [x] 6 data tables in `data/`: module3, module5, module6 (×2), module7 (×2), module8
- [x] `README.md` with executive summary, mechanism, module table, InnoCentive mapping, bibliography
- [x] `PROPOSAL_NOTES.md` with 7 parts (mechanism, interventions, feasibility, risks, experience, timeline, variants)
- [x] `LICENSE` (MIT)
- [x] `CITATION.cff` with author info and key references
- [x] `REPO_AUDIT.md` (this file)
- [x] `pyproject.toml` with all runtime deps declared (numpy, scipy, matplotlib, pandas, SALib)
- [x] `Makefile` with setup / test / figures / clean / all targets

### Git hygiene

- [ ] Verify `.venv/` is NOT staged (`git status` should show it excluded)
- [ ] Verify `__pycache__/` is NOT staged
- [ ] Verify `.DS_Store` is NOT staged
- [ ] Verify `serum_ca_cryo.egg-info/` is NOT staged (add to .gitignore if needed)
- [ ] All figures (`.png`) and data (`.csv`) ARE staged (removed from .gitignore)
- [ ] Commit message for initial push: "Initial public release: Modules 1–8, 139 tests, 11 figures"

### Content quality (manual review)

- [ ] README executive summary paragraph verified accurate against current model output
- [ ] Three headline numbers (5.4%, 42%, 0.0%) match `data/module6_vial_statistics.csv` and `data/module7_intervention_outcomes.csv`
- [ ] PROPOSAL_NOTES.md Part 1 experiment predictions are specific enough to be falsifiable
- [ ] PROPOSAL_NOTES.md Part 6 cost estimates are defensible (instrument time rates from your experience)
- [ ] No claim in PROPOSAL_NOTES.md that isn't backed by a figure or data file

### GitHub repository settings (after push)

- [ ] Repository description: "Computational model for post-thaw Ca deficit in cryostored serum QC standards. 139 tests, 11 figures, make all < 2 min."
- [ ] Topics: `serum`, `calcium`, `cryopreservation`, `hydroxyapatite`, `ostwald-ripening`, `monte-carlo`, `thermodynamics`, `innocentive`
- [ ] License: MIT (auto-detected from LICENSE file)
- [ ] Social preview: upload `figures/fig09_intervention_efficacy.png` (most visually striking)

---

### Summary

**139 tests | 11 figures | 6 data tables | make all < 2 min**

Ready for: `git remote add origin https://github.com/antonphilippov/serum-ca-cryo && git push -u origin main`
