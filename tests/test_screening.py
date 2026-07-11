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
