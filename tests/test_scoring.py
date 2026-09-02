"""Programme objectives against chemistry that is known independently of this code."""

from __future__ import annotations

import pytest

from cosmepep.objectives import PROGRAMME_BY_KEY
from cosmepep.proteases import site_density
from cosmepep.reference import KNOWN_PEPTIDES


def copper(sequence: str) -> float:
    programme = PROGRAMME_BY_KEY["cu_carrier"]
    term = next(t for t in programme.terms if t.name == "copper_motif")
    return term.score(sequence)


def test_atcun_beats_ghk_type_beats_stray_histidine() -> None:
    assert copper("GGH") == 1.0          # Xaa-Xaa-His, four-nitrogen site
    assert copper("GHK") == 0.9          # Xaa-His, three-nitrogen site
    assert copper("GKH") == 1.0          # still position three
    assert copper("GKAH") == 0.15        # histidine present but out of position
    assert copper("GKAA") == 0.0         # no histidine at all


def test_proline_at_position_two_abolishes_the_atcun_site() -> None:
    """It replaces the backbone NH that has to deprotonate; this is not a soft penalty."""
    assert copper("GPH") == 0.0
    assert copper("PGH") == 0.75  # secondary alpha-amine: weaker, not absent


def test_site_density_counts_internal_bonds_only() -> None:
    """A rule matching at the C-terminus describes a bond the molecule does not have."""
    assert site_density("K") == 0.0
    assert 0.0 <= site_density("LSVD") <= 1.0


@pytest.mark.parametrize(
    "sequence", [k.sequence for k in KNOWN_PEPTIDES if k.mechanism_class == "matrix_signal" and k.evidence == "clinical"]
)
def test_clinically_studied_matrikines_are_not_rejected(sequence: str) -> None:
    """KTTKS, GPKG and LSVD went through human studies. An objective that zeroes them is
    measuring the author's expectations, which is how the first version of this scorer failed."""
    assert PROGRAMME_BY_KEY["matrix_signal"].score(sequence)["composite"] > 0.5


def test_order_weight_is_reported_and_below_one() -> None:
    for programme in PROGRAMME_BY_KEY.values():
        assert 0.0 < programme.order_weight < 1.0
