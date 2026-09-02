"""The skin extracellular matrix proteins that supply the design space.

Cosmetic peptides are, almost without exception, fragments of this proteome or close
analogues of them: KTTKS is a stretch of the collagen I alpha-1 C-propeptide, GHK is read
off collagen alpha-2(I), VGVAPG is the elastin hexapeptide repeat. The design premise here
is that the same proteome contains many more such fragments that nobody has looked at, and
that they can be enumerated rather than guessed.

Sequences are fetched from UniProt once and cached on disk, so a run is reproducible
without the network and the cache is the record of exactly which release was used.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import requests

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"

DEFAULT_CACHE = Path(__file__).resolve().parents[2] / "data" / "cache" / "ecm_proteins.json"


@dataclass(frozen=True)
class EcmProtein:
    """One matrix protein and the compartment of skin it belongs to."""

    accession: str
    gene: str
    name: str
    compartment: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)


#: Curated panel. Chosen for skin relevance rather than coverage: the fibrillar collagens
#: and elastic fibre that carry the known matrikines, the dermal-epidermal junction that
#: fails first in photoageing, the small leucine-rich proteoglycans that decorate collagen,
#: and filaggrin, whose proteolysis produces the natural moisturising factor.
ECM_PANEL: tuple[tuple[str, str, str, str], ...] = (
    ("P02452", "COL1A1", "Collagen alpha-1(I) chain", "dermis/fibrillar"),
    ("P08123", "COL1A2", "Collagen alpha-2(I) chain", "dermis/fibrillar"),
    ("P02461", "COL3A1", "Collagen alpha-1(III) chain", "dermis/fibrillar"),
    ("P20908", "COL5A1", "Collagen alpha-1(V) chain", "dermis/fibrillar"),
    ("P12111", "COL6A3", "Collagen alpha-3(VI) chain", "dermis/microfibrillar"),
    ("Q02388", "COL7A1", "Collagen alpha-1(VII) chain", "anchoring fibril/DEJ"),
    ("Q9UMD9", "COL17A1", "Collagen alpha-1(XVII) chain", "hemidesmosome/DEJ"),
    ("P39060", "COL18A1", "Collagen alpha-1(XVIII) chain", "basement membrane"),
    ("Q01955", "COL4A3", "Collagen alpha-3(IV) chain", "basement membrane"),
    ("P15502", "ELN", "Elastin", "dermis/elastic fibre"),
    ("P35555", "FBN1", "Fibrillin-1", "dermis/elastic fibre"),
    ("P02751", "FN1", "Fibronectin", "dermis/provisional matrix"),
    ("Q16787", "LAMA3", "Laminin subunit alpha-3", "DEJ"),
    ("Q13751", "LAMB3", "Laminin subunit beta-3", "DEJ"),
    ("P14543", "NID1", "Nidogen-1", "basement membrane"),
    ("P98160", "HSPG2", "Basement membrane proteoglycan (perlecan)", "basement membrane"),
    ("P07585", "DCN", "Decorin", "dermis/SLRP"),
    ("P51884", "LUM", "Lumican", "dermis/SLRP"),
    ("P20930", "FLG", "Filaggrin", "epidermis/NMF precursor"),
)


def fetch_panel(cache: Path = DEFAULT_CACHE, refresh: bool = False) -> list[EcmProtein]:
    """Return the ECM panel, downloading and caching sequences on first use."""
    if cache.exists() and not refresh:
        payload = json.loads(cache.read_text())
        return [EcmProtein(**record) for record in payload]

    proteins: list[EcmProtein] = []
    for accession, gene, name, compartment in ECM_PANEL:
        response = requests.get(UNIPROT_FASTA.format(accession=accession), timeout=60)
        response.raise_for_status()
        sequence = "".join(
            line.strip() for line in response.text.splitlines() if not line.startswith(">")
        )
        # Selenocysteine and the ambiguity codes break every downstream descriptor; the
        # panel has none today, but a UniProt revision could introduce one silently.
        sequence = "".join(residue for residue in sequence if residue in "ACDEFGHIKLMNPQRSTVWY")
        proteins.append(EcmProtein(accession, gene, name, compartment, sequence))

    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps([protein.__dict__ for protein in proteins], indent=2))
    return proteins


def residue_frequencies(proteins: list[EcmProtein]) -> dict[str, float]:
    """Pooled amino-acid frequency of the panel.

    This is the *right* naive null for matrix-derived peptides. The UniProt-wide background
    is not: the ECM is a third glycine and proline in places, so a sequence drawn from
    generic background frequencies is distinguishable from a matrix fragment on composition
    alone, and any filter would clear that bar without knowing anything.
    """
    counts: dict[str, int] = dict.fromkeys("ACDEFGHIKLMNPQRSTVWY", 0)
    for protein in proteins:
        for residue in protein.sequence:
            counts[residue] += 1
    total = sum(counts.values())
    return {residue: count / total for residue, count in counts.items()}
