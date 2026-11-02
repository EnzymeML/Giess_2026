"""Numeric quadrature of whole peak ranges, and the concentrations derived from it.

This is the path used for the final results of this case study, replacing NMRPy's
Lorentzian/Gaussian deconvolution. That deconvolution becomes poorly conditioned as peak
heights drop towards the noise floor — exactly the concentration regime the time-course
needs to resolve — and has no concept of a local baseline pedestal, so it converges on
inflated peak areas at low concentration. Numeric integration fits nothing; it sums
whatever is actually in the spectrum.

Each entry of ``fid.ranges`` covers one whole multiplet (e.g. the full educt doublet)
rather than an individual line, so **one integral is produced per range**. Nothing is
lost by that: the per-line deconvolution integrals were summed back together for the
concentration calculation anyway.

Three quadrature rules are provided. They are numerically interchangeable on data this
densely sampled (≤0.1% apart once baseline-corrected); the accuracy comes from
:func:`local_linear_baseline`, not from the rule. Simpson's rule is the method selected
for the final results. See ``integration_summary.md`` and ``quadrature_rule_spread.md``
for the full methodology and the comparison numbers.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

import numpy as np
from scipy.integrate import simpson

from modules.concentrations import calculate_c_from_integrals, standard_concentration

if TYPE_CHECKING:
    from nmrpy.data_objects import Fid, FidArray
    from pyenzyme import EnzymeMLDocument

# Number of samples at each edge of a range that anchor the local baseline.
EDGE_POINTS = 3


def integrate_squares(x: np.ndarray, y: np.ndarray) -> float:
    """Rectangle-rule (Riemann sum) integral of y over x, assuming uniform spacing in x.

    Args:
        x: Sample positions in ppm, assumed uniformly spaced.
        y: Spectral intensities at `x`.

    Returns:
        The integral of `y` over `x`.
    """
    dx = abs(x[1] - x[0])
    return float(np.sum(y) * dx)


def integrate_trapezoidal(x: np.ndarray, y: np.ndarray) -> float:
    """Composite trapezoidal-rule integral of y over x.

    Args:
        x: Sample positions in ppm.
        y: Spectral intensities at `x`.

    Returns:
        The absolute integral of `y` over `x`. The absolute value is taken because
        `x` runs from high to low ppm, which would otherwise flip the sign.
    """
    return float(abs(np.trapezoid(y, x)))


def integrate_simpson(x: np.ndarray, y: np.ndarray) -> float:
    """Composite Simpson's-rule integral of y over x.

    Uses :func:`scipy.integrate.simpson`, which handles the composite rule's
    even/odd-interval-count correction automatically.

    Args:
        x: Sample positions in ppm.
        y: Spectral intensities at `x`.

    Returns:
        The absolute integral of `y` over `x`. The absolute value is taken because
        `x` runs from high to low ppm, which would otherwise flip the sign.
    """
    return float(abs(simpson(y, x=x)))


def local_linear_baseline(x: np.ndarray, y: np.ndarray, edge_points: int = EDGE_POINTS) -> np.ndarray:
    """Subtract a straight-line baseline anchored on the mean of `edge_points`
    samples at each edge of the region.

    Global baseline correction (``FidArray.baseline_correct_fids()``) removes broad roll
    and offset across the whole spectrum, but does not guarantee the baseline is exactly
    zero under every individual multiplet — residual pedestal from apodisation ringing or
    overlapping peak tails can remain locally. This is what standard NMR integration
    tools (MestReNova, TopSpin, ...) correct for by default whenever a region is
    integrated.

    Anchoring on the mean of several samples per edge rather than on single points
    reduces sensitivity to noise. Applying it identically to every range and every FID
    (including the mesitylene reference range) is what keeps this a consistently-defined
    processing correction rather than a per-point adjustment.

    Args:
        x: Sample positions in ppm. Unused, but present so that the signature matches
            the `baseline` callable expected by :func:`integrate_fid_ranges`.
        y: Spectral intensities at `x`.
        edge_points: Number of samples at each edge to average when anchoring the line.

    Returns:
        `y` with the fitted straight line subtracted.
    """
    left = y[:edge_points].mean()
    right = y[-edge_points:].mean()
    baseline = np.linspace(left, right, len(y))
    return y - baseline


def iter_fid_ranges(fid: "Fid"):
    """Yield (lower, upper, x, y) for each of fid.ranges, sliced out of fid.data by ppm.

    fid._ppm runs from high to low ppm, so x[0] is the left (higher-ppm) edge of a
    range and x[-1] its right (lower-ppm) edge.

    Args:
        fid: The FID whose ``ranges`` are to be sliced out of its ``data``.

    Yields:
        One ``(lower, upper, x, y)`` tuple per range, in the order of ``fid.ranges``:
        the range bounds in ppm, the ppm positions inside the range, and the spectral
        intensities there.
    """
    ppm = fid._ppm
    for lower, upper in ((min(r), max(r)) for r in fid.ranges):
        mask = (ppm >= lower) & (ppm <= upper)
        yield lower, upper, ppm[mask], fid.data[mask]


def integrate_fid_ranges(
    fid: "Fid",
    method: Callable[[np.ndarray, np.ndarray], float] = integrate_squares,
    baseline: Callable[[np.ndarray, np.ndarray], np.ndarray] = None,
) -> list[float]:
    """Integrate fid.data over each of fid.ranges with the given numeric method.

    Returns one integral per range (i.e. per whole multiplet), in the same
    order as fid.ranges. If given, baseline(x, y) is applied to each range's
    data before integrating (see :func:`local_linear_baseline`).

    Args:
        fid: The FID to integrate.
        method: Quadrature rule to apply, e.g. :func:`integrate_simpson`.
        baseline: Optional local baseline correction applied to each range before
            integrating. ``None`` integrates the data as it is.

    Returns:
        One integral per range of `fid`, in the order of ``fid.ranges``.
    """
    integrals = []
    for _, _, x, y in iter_fid_ranges(fid):
        if baseline is not None:
            y = baseline(x, y)
        integrals.append(method(x, y))
    return integrals


def integrate_fid_array(
    fid_array: "FidArray",
    method: Callable[[np.ndarray, np.ndarray], float] = integrate_squares,
    baseline: Callable[[np.ndarray, np.ndarray], np.ndarray] = None,
) -> dict[str, list[float]]:
    """Integrate every FID of an array, keyed by FID id.

    Convenience wrapper around :func:`integrate_fid_ranges` for the whole time course.

    Args:
        fid_array: The FID array to integrate.
        method: Quadrature rule to apply, e.g. :func:`integrate_simpson`.
        baseline: Optional local baseline correction applied to each range before
            integrating. ``None`` integrates the data as it is.

    Returns:
        A mapping of FID id to that FID's per-range integrals.
    """
    return {fid.id: integrate_fid_ranges(fid, method, baseline=baseline) for fid in fid_array.get_fids()}


def calculate_educt_concentration_from_integrals(
    fid_array: "FidArray",
    integrals_dict: dict[str, list[float]],
    enzml_doc: "EnzymeMLDocument",
    measurement_id: str,
) -> dict[str, list[float]]:
    """Derive educt concentrations from per-range numeric integrals.

    Uses per-range numeric integrals (squares, trapezoidal or Simpson) rather than
    per-line deconvoluted integrals. It is generic over which integrals dict it is given,
    so the same function serves every quadrature rule and both baseline variants.

    A FID normally carries exactly three ranges, in ppm order: the mesitylene reference
    singlet, the first educt doublet, and the second educt doublet.

    A FID carrying **only** the mesitylene range is treated as a measured educt
    concentration of zero rather than an error. That is the situation of the late 33 eq
    time points, where the reaction has run to completion and no educt resonance is left
    to pick: there is nothing in those windows but noise, and picking something
    peak-shaped there would manufacture a false non-zero tail. Reporting zero records
    what was actually observed. Any other number of ranges is still an error, since it
    means the picking is incomplete rather than deliberately empty.

    Args:
        fid_array: The FID array the integrals were computed from; iterated to map FID
            ids onto time points in order.
        integrals_dict: Mapping of FID id to per-range integrals, as returned by
            :func:`integrate_fid_array`.
        enzml_doc: The EnzymeML document holding the prepared mesitylene concentration.
        measurement_id: ID of the measurement to read that concentration from.

    Returns:
        A dict with the educt concentrations in mol/l under ``"s1"``, in FID order.
        The mesitylene integrals are not returned: they are not concentrations, and
        the standard's concentration is its prepared value.

    Raises:
        ValueError: If a FID carries neither exactly three ranges nor the mesitylene
            range alone, or if the measurement holds no mesitylene species data.
    """
    c_mesitylene = standard_concentration(enzml_doc, measurement_id)

    concentrations = {"s1": []}
    for fid in fid_array.get_fids():
        integrals = integrals_dict[fid.id]
        if len(integrals) == 1:
            # Educt fully consumed: standard picked, no educt ranges to integrate.
            concentrations["s1"].append(0.0)
            continue
        if len(integrals) != 3:
            raise ValueError(
                f"{fid.id} carries {len(integrals)} ranges; expected three "
                f"(mesitylene, educt doublet 1, educt doublet 2) or the mesitylene "
                f"range alone where the educt is fully consumed."
            )
        mesitylene_integral, doublet1_integral, doublet2_integral = integrals
        concentrations["s1"].append(
            calculate_c_from_integrals(
                c_mesitylene=c_mesitylene,
                mesitylene_integral=mesitylene_integral,
                educt_doublet1_integral=doublet1_integral,
                educt_doublet2_integral=doublet2_integral,
            )
        )
    return concentrations
