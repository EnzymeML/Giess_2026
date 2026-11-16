"""Non-interactive assignment of EnzymeML species to picked peaks.

:meth:`nmrpy.data_objects.FidArray.assign_peaks` and its per-FID counterpart only
instantiate a dropdown widget; NMRPy currently offers no scripted path, which would
leave a human gate in the middle of an otherwise reproducible notebook.

For this case study the assignment is entirely positional. ``Fid.ranges`` is picked in
ppm order as ``[mesitylene singlet, educt doublet 1, educt doublet 2]``, and NMRPy
builds ``Fid.fid_object.peaks`` by walking the ranges in that order and numbering the
peaks it finds with a running ``peak_index``. Peak 0 is therefore always the internal
standard and every remaining peak belongs to the educt, so the assignment the widget
would produce can be reproduced exactly without it.

That also covers the late 33 eq time points, where only the mesitylene range is picked:
there is a single peak, it is peak 0, and it is the standard.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nmrpy.data_objects import Fid, FidArray


def assign_species_to_fid(
    fid: "Fid",
    standard_species: str = "s5",
    educt_species: str = "s1",
) -> list[str]:
    """Assign species to one FID's peaks by position.

    Peak 0 — the peak picked inside the first (mesitylene) range — is assigned the
    internal standard; every other peak is assigned the educt. Writes both
    ``Fid.species`` and the ``species_id`` of each ``Peak`` in the NMRPy data model, and
    marks the FID as assigned, which is exactly what the assignment widget does.

    Calls ``Fid._setup_peak_objects()`` first, exactly as the widget does. That step is
    what creates the ``Peak`` objects in the NMRPy data model from ``Fid.peaks`` and
    ``Fid.ranges``, and it is otherwise only reached through deconvolution or the
    picking widgets — neither of which this workflow uses, since the peaks and ranges
    are assigned directly and the integration is numeric. Without it the data model
    would carry no peaks at all, and ``record_numeric_integration()`` would have nothing
    to attach its quantification records to.

    Args:
        fid: The FID to assign. Modified in place.
        standard_species: EnzymeML species id of the internal standard.
        educt_species: EnzymeML species id of the educt.

    Returns:
        The assigned species ids, in peak order.

    Raises:
        ValueError: If the FID has no peaks picked, or if its data model holds a
            different number of peaks than ``Fid.peaks``, which would mean the two have
            drifted out of step and positional assignment is unsafe.
    """
    if fid.peaks is None or len(fid.peaks) == 0:
        raise ValueError(f"{fid.id} has no peaks picked; pick peaks before assigning.")

    fid._setup_peak_objects()
    model_peaks = fid.fid_object.peaks
    if len(model_peaks) != len(fid.peaks):
        raise ValueError(
            f"{fid.id} has {len(fid.peaks)} peaks but {len(model_peaks)} in its data "
            f"model; refusing to assign by position."
        )

    species = [
        standard_species if index == 0 else educt_species
        for index in range(len(fid.peaks))
    ]
    fid.species = species
    for peak in model_peaks:
        peak.species_id = species[peak.peak_index]
    fid._flags["assigned"] = True
    return species


def assign_species(
    fid_array: "FidArray",
    standard_species: str = "s5",
    educt_species: str = "s1",
    verbose: bool = False,
) -> dict[str, list[str]]:
    """Assign species to every FID of an array by position.

    Scripted stand-in for :meth:`nmrpy.data_objects.FidArray.assign_peaks`. See
    :func:`assign_species_to_fid` for the rule and the module docstring for why it holds.

    Args:
        fid_array: The FID array to assign. Modified in place.
        standard_species: EnzymeML species id of the internal standard.
        educt_species: EnzymeML species id of the educt.
        verbose: Whether to print the assignment made for each FID.

    Returns:
        A mapping of FID id to that FID's assigned species ids, in peak order.

    Raises:
        ValueError: If any FID has no peaks picked, or if its data model and
            ``Fid.peaks`` disagree on how many peaks there are.
    """
    assigned = {}
    for fid in fid_array.get_fids():
        assigned[fid.id] = assign_species_to_fid(
            fid, standard_species=standard_species, educt_species=educt_species
        )
        if verbose:
            pairs = ", ".join(
                f"{position} ppm -> {species}"
                for position, species in zip(fid.peaks, assigned[fid.id])
            )
            print(f"{fid.id}: {pairs}")
    return assigned
