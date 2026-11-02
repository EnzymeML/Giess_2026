"""Internal-standard math.

Concentrations in this case study are derived relative to mesitylene (species ``s5``),
which is present at a known prepared concentration and whose aromatic singlet is
integrated alongside the educt (species ``s1``) signals. Scaling the educt integrals
against the mesitylene integral, weighted by the number of protons giving rise to each
signal, yields the educt concentration.

The calculation does not depend on how the integrals were obtained, which is why it
lives here rather than in :mod:`modules.numeric_integration`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyenzyme import EnzymeMLDocument

# Protons giving rise to each integrated signal. The educt triplet is a single proton,
# each educt doublet is two, and the mesitylene reference singlet is three.
PROTONS_TRIPLET = 1
PROTONS_DOUBLET1 = 2
PROTONS_DOUBLET2 = 2
PROTONS_MESITYLENE = 3


def sum_integrals(*integrals) -> float:
    """Sum the integrals of peaks.

    Used to combine the per-line integrals of a multiplet (e.g. the two lines of a
    doublet) into a single integral for the whole multiplet.

    Args:
        *integrals: Individual peak integrals to add together.

    Returns:
        The sum of the given integrals.
    """
    return sum(integrals)


def calculate_c_from_integrals(
    c_mesitylene: float,
    mesitylene_integral: float,
    educt_triplet_integral: float = None,
    educt_doublet1_integral: float = None,
    educt_doublet2_integral: float = None,
) -> float:
    """Calculate the amount of educt from its integrals against the
    mesitylene integral.

    Any combination of the three educt signals may be supplied; only those actually
    given contribute to both the summed integral and the summed proton count, so the
    proton weighting stays consistent whether one, two or all three are available.

    Args:
        c_mesitylene: The concentration of mesitylene in mol/l.
        mesitylene_integral: The integral of the mesitylene peak.
        educt_triplet_integral: The integral of the educt triplet peak.
        educt_doublet1_integral: The integral of the educt first doublet peak.
        educt_doublet2_integral: The integral of the educt second doublet peak.

    Returns:
        The amount of educt in mol/l.

    Raises:
        ValueError: If none of the triplet, doublet1 or doublet2 integrals is given.
    """
    if educt_triplet_integral is None and educt_doublet1_integral is None and educt_doublet2_integral is None:
        raise ValueError("At least one of the triplet, doublet1, or doublet2 integrals must be provided.")

    sum_educt_integrals = 0
    sum_educt_protons = 0

    if educt_triplet_integral is not None:
        sum_educt_integrals += educt_triplet_integral
        sum_educt_protons += PROTONS_TRIPLET
    if educt_doublet1_integral is not None:
        sum_educt_integrals += educt_doublet1_integral
        sum_educt_protons += PROTONS_DOUBLET1
    if educt_doublet2_integral is not None:
        sum_educt_integrals += educt_doublet2_integral
        sum_educt_protons += PROTONS_DOUBLET2

    c_educt = ((sum_educt_integrals*PROTONS_MESITYLENE) / (mesitylene_integral*sum_educt_protons)) * c_mesitylene

    return c_educt


def standard_concentration(
    enzml_doc: "EnzymeMLDocument",
    measurement_id: str,
    species_id: str = "s5",
) -> float:
    """Look up the prepared concentration of the internal standard.

    Reads the ``prepared`` concentration recorded for the standard in the EnzymeML
    document, which is the reference the educt integrals are scaled against.

    Args:
        enzml_doc: The EnzymeML document holding the measurement.
        measurement_id: ID of the measurement to read the standard from.
        species_id: ID of the internal standard species. Defaults to ``"s5"``
            (mesitylene).

    Returns:
        The prepared concentration of the internal standard in mol/l.

    Raises:
        ValueError: If the measurement holds no species data for `species_id`.
    """
    for datum in enzml_doc.filter_measurements(id=measurement_id)[0].species_data:
        if datum.species_id == species_id:
            return datum.prepared
    raise ValueError(f"No species data for {species_id!r} in measurement {measurement_id!r}.")
