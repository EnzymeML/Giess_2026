"""Recording a numeric integration in the NMRPy data model.

The integrals produced by :mod:`modules.numeric_integration` are plain numbers; on their
own they do not say how they were obtained. The NMRPy data model has fields for exactly
that information, which is what makes the saved document a reusable FAIR artifact rather
than a bare table: every ``Peak`` carries a ``peak_quantification`` describing **how** its
area was obtained — a ``Quantification`` with ``method=NUMERIC_INTEGRATION``, the
``Quadrature`` rule applied and the number of spectral points that entered it, and the
local ``Baseline`` that was subtracted first — while its ``peak_range`` carries the
``proton_count`` used to weight that area when calculating concentrations.

One integral is produced per range (per whole multiplet), whereas the data model holds one
``Peak`` per picked line, several of which share the same range. The leftmost (highest-ppm)
``Peak`` of each range therefore carries that range's ``Quantification`` and the remaining
``Peak`` objects of the range are left without one, so summing ``peak_area`` over a FID's
peaks counts each multiplet exactly once. ``proton_count`` is set on every peak of the
range, since it describes the range itself.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

import numpy as np
from nmrpy.nmrpy_model import (
    Baseline,
    BaselineMethods,
    Quadrature,
    QuadratureRules,
    Quantification,
    QuantificationMethods,
)

from modules.numeric_integration import EDGE_POINTS, iter_fid_ranges

if TYPE_CHECKING:
    from nmrpy.data_objects import Fid, FidArray

# Re-exported so notebooks can name the quadrature rules without a second import.
__all__ = [
    "QuadratureRules",
    "local_linear_baseline_record",
    "proton_count_for_range",
    "peaks_of_range",
    "record_numeric_integration",
]


def local_linear_baseline_record(x: np.ndarray, y: np.ndarray, edge_points: int = EDGE_POINTS) -> Baseline:
    """`Baseline` record of the straight line that `local_linear_baseline()` subtracts.

    Since `x` runs from high to low ppm, `x[:edge_points]` are the anchor points at the
    left (higher-ppm) edge of the range and `x[-edge_points:]` those at the right
    (lower-ppm) edge. The fitted line is reported in ppm, as `slope * ppm + intercept`.

    Args:
        x: Sample positions in ppm inside the range.
        y: Spectral intensities at `x`.
        edge_points: Number of samples at each edge averaged to anchor the line. Must
            match what :func:`modules.numeric_integration.local_linear_baseline` used.

    Returns:
        A ``Baseline`` describing the line that was subtracted.
    """
    left = float(y[:edge_points].mean())
    right = float(y[-edge_points:].mean())
    slope = (right - left) / float(x[-1] - x[0])
    return Baseline(
        method=BaselineMethods.LOCAL_LINEAR,
        anchor_points_start=[float(v) for v in x[:edge_points]],
        anchor_points_end=[float(v) for v in x[-edge_points:]],
        slope=slope,
        intercept=left - slope * float(x[0]),
    )


def proton_count_for_range(lower: float, upper: float) -> int:
    """Number of protons giving rise to the multiplet covered by a range.

    Keyed by where the range sits in the spectrum, following the peak assignment
    documented in the notebooks and the proton counts used by
    :func:`modules.concentrations.calculate_c_from_integrals`:

    - ~6.71 ppm       mesitylene/BINOL singlet, 3 protons
    - ~5.08-5.12 ppm  educt triplet, 1 proton (triplet-doublet-doublet assignment only)
    - ~3.0-3.2 ppm    educt doublets, 2 protons each

    Args:
        lower: Lower (lower-ppm) bound of the range.
        upper: Upper (higher-ppm) bound of the range.

    Returns:
        The proton count for the multiplet in that range.
    """
    centre = 0.5 * (lower + upper)
    if centre > 6.0:
        return 3
    if centre > 4.0:
        return 1
    return 2


def peaks_of_range(fid: "Fid", lower: float, upper: float) -> list:
    """Peak objects of `fid` belonging to the given range, highest ppm first.

    Args:
        fid: The FID whose data-model peaks are searched.
        lower: Lower (lower-ppm) bound of the range.
        upper: Upper (higher-ppm) bound of the range.

    Returns:
        The ``Peak`` objects whose ``peak_range`` matches the given bounds, sorted by
        descending ppm so that element 0 is the leftmost peak of the range. Empty if the
        range has no peaks in the data model.
    """
    peaks = [
        peak
        for peak in fid.fid_object.peaks
        if peak.peak_range is not None
        and np.allclose(sorted((peak.peak_range.start, peak.peak_range.end)), (lower, upper))
    ]
    return sorted(peaks, key=lambda peak: peak.peak_position, reverse=True)


def record_numeric_integration(
    fid_array: "FidArray",
    integrals_dict: dict[str, list[float]],
    rule: QuadratureRules,
    baseline_record: Callable[[np.ndarray, np.ndarray], Baseline] = None,
) -> None:
    """Write the numeric integrals into the Peak objects of the NMRPy data model.

    For every FID, entry i of `integrals_dict[fid.id]` is the integral of range i of
    `fid.ranges`. That area, the `Quadrature` rule and point count, and the local
    `Baseline` are attached to the leftmost `Peak` of the range; any `Quantification`
    left on the other peaks of the range by an earlier call is cleared. The range's
    `proton_count` is set on all of its peaks.

    Calling this repeatedly overwrites the previous recording, so the last call made is
    what a subsequently saved data model holds.

    Args:
        fid_array: The FID array whose data model is written to.
        integrals_dict: Mapping of FID id to per-range integrals, as returned by
            :func:`modules.numeric_integration.integrate_fid_array`.
        rule: The quadrature rule these integrals were computed with.
        baseline_record: Returns the ``Baseline`` that was subtracted before integrating,
            e.g. :func:`local_linear_baseline_record`. Pass ``None`` for the raw variants,
            which are recorded with ``BaselineMethods.NONE``.

    Returns:
        None. The data model of `fid_array` is modified in place.
    """
    for fid in fid_array.get_fids():
        integrals = integrals_dict[fid.id]
        for (lower, upper, x, y), area in zip(iter_fid_ranges(fid), integrals):
            peaks = peaks_of_range(fid, lower, upper)
            if not peaks:
                continue
            for peak in peaks:
                peak.peak_range.proton_count = proton_count_for_range(lower, upper)
                peak.peak_quantification = None
            peaks[0].peak_quantification = Quantification(
                method=QuantificationMethods.NUMERIC_INTEGRATION,
                peak_area=float(area),
                local_baseline=(
                    baseline_record(x, y)
                    if baseline_record is not None
                    else Baseline(method=BaselineMethods.NONE)
                ),
                quadrature=Quadrature(rule=rule, n_points=len(x)),
            )
