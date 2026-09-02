"""Molecular-graph properties, and the delivery problem that defines cosmetic peptides.

A cosmetic peptide has to cross the stratum corneum, which is where almost all of them
fail. KTTKS on its own is a 563 Da molecule with a computed logP near -4.6 and a topological
polar surface area close to 300 A^2; it is not going anywhere. Palmitoylation is not a
formulation nicety, it is the entire reason the molecule works as a product. So every
candidate here is evaluated twice -- free acid, and N-terminally palmitoylated -- and the
comparison between the two is a design output in its own right.

The permeation number is Potts-Guy. It was fitted to small non-electrolytes and a
hexadecapeptide is far outside its applicability domain; it is used here to *rank* candidates
that are all outside that domain in the same way, never as an absolute flux prediction.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

#: The alpha-amino group specifically: a free amine on a carbon that also carries a
#: carbonyl. Written this way because the naive "primary amine" match also hits every
#: lysine epsilon-amine, and acylating one of those is a different molecule.
_ALPHA_AMINE = "[NX3;H1,H2;!$(NC=O);!$(NC=N):1][CX4:2][CX3:3](=[OX1:4])"

_PALMITOYLATE = AllChem.ReactionFromSmarts(
    f"{_ALPHA_AMINE}>>C(=O)(CCCCCCCCCCCCCCC)[N:1][C:2][C:3]=[O:4]"
)
_ACETYLATE = AllChem.ReactionFromSmarts(f"{_ALPHA_AMINE}>>C(=O)(C)[N:1][C:2][C:3]=[O:4]")
#: C-terminal carboxyl only: the carboxyl carbon must be bonded to the alpha carbon that
#: bears the backbone nitrogen, which excludes the Asp and Glu side chains.
_AMIDATE = AllChem.ReactionFromSmarts(
    "[NX3:1][CX4:2][CX3:3](=[OX1:4])[OX2H1]>>[N:1][C:2][C:3](=[O:4])N"
)


@dataclass(frozen=True)
class MolecularProfile:
    """RDKit descriptors for one chemical form of a peptide."""

    form: str
    molecular_weight: float
    clogp: float
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    heavy_atoms: int
    #: Potts-Guy skin permeability coefficient, log10(cm/s).
    log_kp: float

    @property
    def under_500_daltons(self) -> bool:
        """The Bos-Meinardi rule of thumb for what crosses intact skin at all."""
        return self.molecular_weight <= 500.0


def _apply(reaction: AllChem.ChemicalReaction, mol: Chem.Mol) -> Chem.Mol | None:
    products = reaction.RunReactants((mol,))
    if not products:
        return None
    product = products[0][0]
    try:
        Chem.SanitizeMol(product)
    except (Chem.KekulizeException, ValueError):
        return None
    return product


def build(sequence: str, form: str = "free") -> Chem.Mol | None:
    """Construct the molecule for a sequence in one of its formulated forms.

    Forms: ``free`` (H-peptide-OH), ``palmitoyl`` (the KTTKS/GHK product form),
    ``acetyl_amide`` (the Argireline form, Ac-peptide-NH2).
    """
    mol = Chem.MolFromSequence(sequence)
    if mol is None:
        return None
    if form == "free":
        return mol
    if form == "palmitoyl":
        return _apply(_PALMITOYLATE, mol)
    if form == "acetyl_amide":
        acetylated = _apply(_ACETYLATE, mol)
        return None if acetylated is None else _apply(_AMIDATE, acetylated)
    raise ValueError(f"unknown form: {form}")


def potts_guy(clogp: float, molecular_weight: float) -> float:
    """log Kp = -2.7 + 0.71 logP - 0.0061 MW, in log10 cm/s."""
    return -2.7 + 0.71 * clogp - 0.0061 * molecular_weight


@lru_cache(maxsize=100_000)
def profile(sequence: str, form: str = "free") -> MolecularProfile | None:
    """Descriptor set for a sequence in a given form, or None if it cannot be built."""
    mol = build(sequence, form)
    if mol is None:
        return None
    molecular_weight = Descriptors.MolWt(mol)
    clogp = Crippen.MolLogP(mol)
    return MolecularProfile(
        form=form,
        molecular_weight=molecular_weight,
        clogp=clogp,
        tpsa=rdMolDescriptors.CalcTPSA(mol),
        hbd=rdMolDescriptors.CalcNumHBD(mol),
        hba=rdMolDescriptors.CalcNumHBA(mol),
        rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(mol),
        heavy_atoms=mol.GetNumHeavyAtoms(),
        log_kp=potts_guy(clogp, molecular_weight),
    )


def smiles(sequence: str, form: str = "free") -> str | None:
    mol = build(sequence, form)
    return None if mol is None else Chem.MolToSmiles(mol)


#: Palmitoylation adds 238.4 Da and, because Crippen logP is a sum of atom contributions,
#: a fixed logP increment of 5.64. Potts-Guy is linear in both, so the permeability gain
#: from lipidation is the same 2.55 log units for *every* sequence. It looks like a
#: scoring term and is arithmetically a constant, so it is excluded from the objective.
LIPIDATION_LOG_KP_GAIN = 2.55


def delivery_gain(sequence: str) -> float | None:
    """Permeability gained by palmitoylation -- analytically constant, kept as a check.

    This returns ``LIPIDATION_LOG_KP_GAIN`` for every buildable sequence, by construction of
    the two additive models involved. It is retained so the constancy is testable rather
    than assumed, and as a standing warning: a lipidation term in a ranking function that
    is not compared against its own null will silently contribute nothing while appearing
    to encode delivery chemistry.
    """
    free = profile(sequence, "free")
    lipidated = profile(sequence, "palmitoyl")
    if free is None or lipidated is None:
        return None
    return lipidated.log_kp - free.log_kp
