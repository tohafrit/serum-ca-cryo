"""
Tests for Module 8: WATEQ vs Davies thermodynamic cross-check.
All assertions are physical sanity checks or model-agreement criteria.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
import numpy as np
from src.phreeqc_runner import (
    compute_si, run_validation, _pass_summary, _speciate,
    BASE_COMP_mM, SCENARIOS, PHASES, AGREE_TOL,
)


@pytest.fixture(scope="module")
def validation_rows():
    return run_validation()


# ── Physical sanity: HAp is supersaturated in both scenarios ──────────────────

def test_hap_supersaturated_physiological_davies():
    """HAp SI > 0 at k=1, pH 7.4 using Davies (Module 3 model)."""
    comp = {el: v for el, v in BASE_COMP_mM.items()}
    si = compute_si(comp, pH=7.40, model="davies")
    assert si["Hydroxyapatite"] > 0, \
        f"HAp SI = {si['Hydroxyapatite']:.2f} at physiological; expected > 0"

def test_hap_supersaturated_physiological_wateq():
    """HAp SI > 0 at k=1, pH 7.4 using WATEQ (PHREEQC model)."""
    comp = {el: v for el, v in BASE_COMP_mM.items()}
    si = compute_si(comp, pH=7.40, model="wateq")
    assert si["Hydroxyapatite"] > 0


def test_hap_si_increases_with_cryoconcentration():
    """Cryoconcentration (k=5.58) raises SI(HAp) vs physiological (both models)."""
    for model in ("davies", "wateq"):
        comp_1  = {el: v       for el, v in BASE_COMP_mM.items()}
        comp_k  = {el: v * 5.58 for el, v in BASE_COMP_mM.items()}
        si_1  = compute_si(comp_1,  pH=7.40, model=model)
        si_k  = compute_si(comp_k,  pH=7.81, model=model)
        assert si_k["Hydroxyapatite"] > si_1["Hydroxyapatite"], \
            f"{model}: cryo SI(HAp) not > physiological SI(HAp)"


# ── Davies vs WATEQ agreement criteria ────────────────────────────────────────

def test_physiological_models_agree_on_all_phases(validation_rows):
    """At physiological I≈0.15M, Davies and WATEQ must agree within ±1.0 for all 4 phases."""
    phys_rows = [r for r in validation_rows if r["scenario"] == "physiological"]
    for r in phys_rows:
        assert r["abs_delta_SI"] <= AGREE_TOL, (
            f"Physiological {r['phase']}: |ΔSI| = {r['abs_delta_SI']:.3f} > {AGREE_TOL}"
        )

def test_physiological_models_agree_very_tightly(validation_rows):
    """At I≈0.15M (physiological), Davies and WATEQ should agree within ±0.15 log units."""
    phys_rows = [r for r in validation_rows if r["scenario"] == "physiological"]
    for r in phys_rows:
        assert r["abs_delta_SI"] < 0.15, (
            f"Physiological {r['phase']}: |ΔSI| = {r['abs_delta_SI']:.3f}; "
            "expected <0.15 at low I"
        )

def test_cryo_models_agree_on_at_least_3_phases(validation_rows):
    """
    At cryoconcentrated I≈0.80M, Davies and WATEQ agree within ±1.0 for ≥ 3/4 phases.
    This is the headline Module 8 pass criterion (same spec as original PHREEQC target).
    """
    passes = _pass_summary(validation_rows)
    assert passes["cryoconcentrated"], (
        "Cryoconcentrated scenario: <3 of 4 phases agree within |ΔSI|≤1.0 between Davies and WATEQ"
    )

def test_cryo_models_agree_within_half_unit(validation_rows):
    """At I≈0.80M (cryo), all phases should agree within ±0.5 log units."""
    cryo_rows = [r for r in validation_rows if r["scenario"] == "cryoconcentrated"]
    for r in cryo_rows:
        assert r["abs_delta_SI"] < 0.5, (
            f"Cryoconc {r['phase']}: |ΔSI| = {r['abs_delta_SI']:.3f}; expected <0.5"
        )


# ── Monotonicity: SI increases with ionic strength ────────────────────────────

def test_calcite_si_increases_with_k(validation_rows):
    """SI(Calcite) at k=5.58 > SI(Calcite) at k=1 — both models."""
    for model_col in ("SI_davies", "SI_wateq"):
        phys = next(r[model_col] for r in validation_rows
                    if r["scenario"] == "physiological" and r["phase"] == "Calcite")
        cryo = next(r[model_col] for r in validation_rows
                    if r["scenario"] == "cryoconcentrated" and r["phase"] == "Calcite")
        assert cryo > phys, f"{model_col}: Calcite SI not higher at cryo vs physiological"

def test_brushite_supersaturated_at_cryo():
    """At k=5.58, Brushite (CaHPO4·2H2O) should be supersaturated (SI > 0)."""
    comp = {el: v * 5.58 for el, v in BASE_COMP_mM.items()}
    for model in ("davies", "wateq"):
        si = compute_si(comp, pH=7.81, model=model)
        assert si["Brushite"] > 0, \
            f"{model}: Brushite SI = {si['Brushite']:.2f} at cryo; expected > 0"


# ── Speciation sanity ─────────────────────────────────────────────────────────

def test_ca_free_less_than_total():
    """Free Ca²⁺ must be < total Ca due to complex formation."""
    comp = {el: v for el, v in BASE_COMP_mM.items()}
    sp = _speciate(comp, pH=7.40, model="davies")
    total_Ca = comp["Ca"] * 1e-3
    assert sp["Ca_free"] < total_Ca, \
        f"Ca_free ({sp['Ca_free']*1e3:.3f} mM) >= total Ca ({comp['Ca']} mM)"

def test_ionic_strength_physiological_range():
    """Physiological I should be in range 0.10–0.20 mol/L."""
    comp = {el: v for el, v in BASE_COMP_mM.items()}
    sp = _speciate(comp, pH=7.40)
    assert 0.10 <= sp["I"] <= 0.20, \
        f"Physiological I = {sp['I']:.3f} mol/L; expected 0.10–0.20"

def test_ionic_strength_cryo_range():
    """Cryoconcentrated I should be in range 0.50–1.20 mol/L."""
    comp = {el: v * 5.58 for el, v in BASE_COMP_mM.items()}
    sp = _speciate(comp, pH=7.81)
    assert 0.50 <= sp["I"] <= 1.20, \
        f"Cryo I = {sp['I']:.3f} mol/L; expected 0.50–1.20"

def test_csv_exists():
    """Module 8 CSV must exist after run_validation()."""
    csv_path = Path(__file__).parent.parent / "data" / "module8_validation.csv"
    assert csv_path.exists(), "data/module8_validation.csv not found; run main()"
