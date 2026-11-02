"""Manual, per-FID phase correction for arrays that automatic phasing cannot handle.

:func:`nmrpy.data_objects.FidArray.phase_correct_fids` fits the zero- and first-order
phase of every FID automatically, and is what the 2 eq and 5 eq datasets use. It does
not converge usefully on the 33 eq array: the educt signals there fall close to the
noise floor within the first few time points, leaving the fit with too little signal to
key on, so each FID has to be phased by hand against angles read off
:meth:`nmrpy.data_objects.Fid.phaser`.

This module holds the helper that applies such hand-determined angles, so the notebooks
carry the angles as data and not the machinery that applies them.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nmrpy.data_objects import FidArray


def manual_phaser(
    fid_array: "FidArray",
    p0: float = 0.0,
    p1: float = 0.0,
    n_runs: int = 1,
    exclude_fid: list[str] | str | None = None,
    include_fid: list[str] | str | None = None,
    plot: bool = True,
    verbose: bool = False,
) -> None:
    """Apply a manual phase correction to selected FIDs of an array.

    Calls :meth:`nmrpy.data_objects.Fid.ps` with the given angles on each selected FID,
    optionally repeating the correction so that iterated small adjustments accumulate.

    Exactly one of `exclude_fid` and `include_fid` may be given. With neither, every FID
    of the array is phased.

    Args:
        fid_array: The FID array to phase. Modified in place.
        p0: Zero-order (frequency-independent) phase correction, in degrees.
        p1: First-order (frequency-dependent) phase correction, in degrees.
        n_runs: Number of times the correction is applied to each selected FID.
        exclude_fid: FID id, or list of ids, to leave untouched.
        include_fid: FID id, or list of ids, to phase to the exclusion of all others.
        plot: Whether to plot the whole array once phasing has finished.
        verbose: Whether to report progress FID by FID.

    Raises:
        ValueError: If both `p0` and `p1` are zero, which would be a no-op and is
            almost always a sign the angles were not filled in.
        ValueError: If both `exclude_fid` and `include_fid` are given.
    """
    if p0 == 0 and p1 == 0:
        raise ValueError("p0 and p1 are both 0, which would not change the data.")
    if exclude_fid and include_fid:
        raise ValueError("Give either exclude_fid or include_fid, not both.")

    if exclude_fid:
        if isinstance(exclude_fid, str):
            exclude_fid = [exclude_fid]
        fids = [fid for fid in fid_array.get_fids() if fid.id not in exclude_fid]
        if verbose:
            print(f"Excluding {exclude_fid} from phasing.")
    elif include_fid:
        if isinstance(include_fid, str):
            include_fid = [include_fid]
        fids = [fid for fid in fid_array.get_fids() if fid.id in include_fid]
        if verbose:
            print(f"Phasing only {include_fid}.")
    else:
        fids = fid_array.get_fids()
        if verbose:
            print("Phasing all FIDs.")

    for run in range(n_runs):
        if verbose and n_runs > 1:
            print(f"Phasing iteration {run + 1} of {n_runs}.")
        for fid in fids:
            if verbose:
                print(f"Phasing {fid.id} with p0={p0}, p1={p1}.")
            fid.ps(p0=p0, p1=p1)

    if plot:
        fid_array.plot_array()


def apply_phase_angles(
    fid_array: "FidArray",
    p0_values: list[float],
    p1_values: list[float],
    plot: bool = True,
    verbose: bool = False,
) -> None:
    """Phase every FID of an array with its own pair of hand-determined angles.

    Thin wrapper over :func:`manual_phaser` that walks the array and its angle lists
    together, which is how the 33 eq dataset is phased.

    Args:
        fid_array: The FID array to phase. Modified in place.
        p0_values: One zero-order angle per FID, in array order.
        p1_values: One first-order angle per FID, in array order.
        plot: Whether to plot the whole array once every FID has been phased.
        verbose: Whether to report progress FID by FID.

    Raises:
        ValueError: If the number of angles does not match the number of FIDs.
    """
    fids = fid_array.get_fids()
    if not len(p0_values) == len(p1_values) == len(fids):
        raise ValueError(
            f"Need one p0 and one p1 per FID: got {len(p0_values)} p0, "
            f"{len(p1_values)} p1 for {len(fids)} FIDs."
        )

    for fid, p0, p1 in zip(fids, p0_values, p1_values):
        manual_phaser(
            fid_array,
            p0=p0,
            p1=p1,
            include_fid=fid.id,
            plot=False,
            verbose=verbose,
        )

    if plot:
        fid_array.plot_array()
