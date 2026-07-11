import pytest

from moldockpipe.services.molscrub import MolScrubService


def test_non_enumerating_generation_preserves_a_single_state() -> None:
    pytest.importorskip("rdkit")
    pytest.importorskip("molscrub")
    states, truncated = MolScrubService().generate("CCO", enumerate_states=False, max_states=4)
    assert len(states) == 1
    assert not truncated
    assert states[0].state_isomeric_smiles
    assert states[0].structure_hash
