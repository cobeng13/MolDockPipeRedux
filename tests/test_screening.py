import pytest

from moldockpipe.models import ScreeningPolicy
from moldockpipe.services.screening import screen_smiles


def test_annotate_only_keeps_rule_failures_eligible() -> None:
    pytest.importorskip("rdkit")
    result = screen_smiles("CCCCCCCCCCCCCCCCCCCC", {"lipinski": True, "veber": True, "egan": False, "ghose": False}, ScreeningPolicy.ANNOTATE_ONLY)
    assert result.decision == "warning"
    assert result.rules["lipinski"] is False


def test_invalid_smiles_is_failure() -> None:
    pytest.importorskip("rdkit")
    result = screen_smiles("not smiles", {"lipinski": True}, ScreeningPolicy.ANNOTATE_ONLY)
    assert result.decision == "fail"


def test_boiled_egg_uses_reference_yolk_thresholds() -> None:
    pytest.importorskip("rdkit")
    result = screen_smiles("c1ccccc1", {"boiled_egg": True}, ScreeningPolicy.EXCLUDE_FAILING_ANY)
    assert result.descriptors["boiled_egg"] == "YOLK"
    assert result.rules["boiled_egg"] is True
    assert result.decision == "pass"


def test_boiled_egg_rejects_non_yolk_compounds_when_selected() -> None:
    pytest.importorskip("rdkit")
    result = screen_smiles("OCCO", {"boiled_egg": True}, ScreeningPolicy.EXCLUDE_FAILING_ANY)
    assert result.descriptors["boiled_egg"] != "YOLK"
    assert result.rules["boiled_egg"] is False
    assert result.decision == "fail"
