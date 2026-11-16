"""Writing the generated artifacts, and displaying models inside a notebook.

The notebooks produce two kinds of file at the end of a run — the NMRPy data model as
JSON and a CSV of the derived concentrations — and repeatedly render EnzymeML and NMRPy
models to the notebook output with ``rich``. Both are collected here so the notebooks
carry call sites rather than definitions.

On JSON-LD: the data model serializes with ``ld_id``/``ld_type``/``ld_context`` keys that
are part of what makes the saved document a FAIR artifact, so :func:`save_data_model`
keeps them by default. Stripping them is a presentation convenience — for slides and
figures, where the context is noise — and has to be asked for explicitly.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd
from rich import print as rich_print
from rich.pretty import Pretty
from rich.table import Table

if TYPE_CHECKING:
    from nmrpy.data_objects import FidArray
    from pyenzyme import EnzymeMLDocument


def _measurement_name(fid_array: "FidArray", measurement_id: str) -> str:
    """Filename stem for a measurement's artifacts.

    Args:
        fid_array: The FID array holding the EnzymeML document.
        measurement_id: ID of the measurement to name.

    Returns:
        The measurement's own ``id``, as looked up in the EnzymeML document.
    """
    return fid_array.enzymeml_document.filter_measurements(id=measurement_id)[0].id


def strip_json_ld(obj: Any) -> Any:
    """Recursively remove JSON-LD (`ld_id`/`ld_type`/`ld_context`) fields from a nested dict/list.

    Args:
        obj: A nested structure of dicts and lists, e.g. a parsed data-model document.

    Returns:
        The same object, modified in place, with the JSON-LD keys removed at every level.
    """
    if isinstance(obj, dict):
        for key in ("ld_id", "ld_type", "ld_context"):
            obj.pop(key, None)
        for value in obj.values():
            strip_json_ld(value)
    elif isinstance(obj, list):
        for item in obj:
            strip_json_ld(item)
    return obj


def save_data_model(
    fid_array: "FidArray",
    measurement_id: str,
    directory: str = ".",
    strip_ld: bool = False,
) -> Path:
    """Write the NMRPy data model of a FID array to JSON.

    The file is named after the measurement, as ``<measurement id>_data_model.json``.

    Args:
        fid_array: The FID array whose ``data_model`` is written.
        measurement_id: ID of the measurement the data model belongs to.
        directory: Directory to write into.
        strip_ld: Whether to remove the JSON-LD ``ld_id``/``ld_type``/``ld_context``
            fields first. Defaults to ``False``, which keeps the context — that is what
            the notebooks want. Set it to ``True`` for a slimmed-down document
            to read in a presentation or figure.

    Returns:
        The path written to.
    """
    path = Path(directory) / f"{_measurement_name(fid_array, measurement_id)}_data_model.json"

    # Writing the file is the modification NMRPy's own save_data_model() records;
    # reading fid_array.data_model no longer stamps it.
    now = datetime.now().isoformat()
    fid_array.data_model.datetime_modified = now
    if not fid_array.data_model.datetime_created:
        fid_array.data_model.datetime_created = now

    if strip_ld:
        data_model_dict = strip_json_ld(json.loads(fid_array.data_model.model_dump_json()))
        with open(path, "w") as f:
            json.dump(data_model_dict, f, indent=2)
    else:
        with open(path, "w") as f:
            f.write(fid_array.data_model.model_dump_json(indent=2))

    return path


def write_concentrations_csv(
    fid_array: "FidArray",
    measurement_id: str,
    directory: str = ".",
    concentrations: dict[str, list[float]] = None,
) -> Path:
    """Write the time course and derived concentrations to CSV.

    The file is named after the measurement, as ``<measurement id>_concentrations.csv``.

    Args:
        fid_array: The FID array supplying the time array ``t``.
        measurement_id: ID of the measurement the concentrations belong to.
        directory: Directory to write into.
        concentrations: The concentrations to write. Defaults to ``None``, which reads
            them from ``fid_array.concentrations``.

    Returns:
        The path written to.

    Raises:
        KeyError: If the concentrations hold no ``"s1"`` entry.
    """
    if concentrations is None:
        concentrations = fid_array.concentrations

    path = Path(directory) / f"{_measurement_name(fid_array, measurement_id)}_concentrations.csv"

    df = pd.DataFrame()
    df["t / h"] = fid_array.t
    df["[s1] / M"] = concentrations["s1"]
    df.to_csv(path, index=False)

    return path


def print_model(model: Any, max_length: int = 3) -> None:
    """Render a model to the notebook output with `rich`.

    Args:
        model: The object to render, e.g. an ``EnzymeMLDocument`` or an NMRPy data model.
        max_length: Maximum number of elements shown per collection before truncating.

    Returns:
        None. The rendering is printed.
    """
    rich_print(Pretty(model, max_length=max_length))


def print_measurement(
    enzml_doc: "EnzymeMLDocument",
    measurement_id: str,
    max_length: int = 5,
) -> None:
    """Render one measurement of an EnzymeML document to the notebook output.

    Args:
        enzml_doc: The EnzymeML document to look the measurement up in.
        measurement_id: ID of the measurement to render.
        max_length: Maximum number of elements shown per collection before truncating.

    Returns:
        None. The rendering is printed.
    """
    print_model(enzml_doc.filter_measurements(id=measurement_id)[0], max_length=max_length)


def summarize_measurement(
    enzml_doc: "EnzymeMLDocument",
    measurement_id: str,
) -> None:
    """Render a one-line-per-species summary of a measurement.

    :func:`print_measurement` hands the whole measurement to ``rich``, which falls back
    to the pydantic ``__repr__`` for these models and fills the notebook with several
    thousand characters of nested unit definitions and JSON-LD ids. This prints only
    what is worth reading after a time-course has been written back: which species carry
    data, how many points each has, and the span of the values.

    Args:
        enzml_doc: The EnzymeML document to look the measurement up in.
        measurement_id: ID of the measurement to summarize.

    Returns:
        None. The summary is printed.

    Raises:
        ValueError: If the document holds no measurement with `measurement_id`.
    """
    matches = enzml_doc.filter_measurements(id=measurement_id)
    if not matches:
        raise ValueError(f"No measurement {measurement_id!r} in the document.")
    measurement = matches[0]

    table = Table(title=f"{measurement.id} — {measurement.name}")
    for column in ("species", "prepared", "initial", "points", "first", "last"):
        table.add_column(column, justify="right" if column != "species" else "left")

    for datum in measurement.species_data:
        values = list(datum.data or [])
        table.add_row(
            datum.species_id,
            "—" if datum.prepared is None else f"{datum.prepared:.5g}",
            "—" if datum.initial is None else f"{datum.initial:.5g}",
            str(len(values)),
            f"{values[0]:.6g}" if values else "—",
            f"{values[-1]:.6g}" if values else "—",
        )

    rich_print(table)
    times = next((list(d.time or []) for d in measurement.species_data if d.time), [])
    if times:
        rich_print(f"time / h: {[round(t, 4) for t in times]}")


def set_concentrations(
    fid_array: "FidArray",
    concentrations: dict[str, list[float]],
) -> None:
    """Attach calculated concentrations to a FID array.

    Normally just ``fid_array.concentrations = concentrations``. That setter, however,
    validates the species keys against ``fid.species`` where ``fid`` is a variable left
    over from the preceding loop — so it checks the **last** FID only. Where the last
    spectrum carries fewer species than the first, the educt key is rejected even though
    every value is well formed. That is exactly the 33 eq case: the reaction is complete
    before the third time point, so the late spectra have only the mesitylene peak
    picked and their assigned species is ``["s5"]`` alone.

    This performs the setter's checks against *every* FID rather than the last one, and
    on success writes the backing attribute directly.

    Args:
        fid_array: The FID array to attach the concentrations to. Modified in place.
        concentrations: Mapping of EnzymeML species id to one value per time point.

    Returns:
        None.

    Raises:
        TypeError: If `concentrations` is not a dict.
        ValueError: If a FID has no assigned species, if a species key is assigned
            nowhere in the array, if a series length does not match the time axis, or if
            a value is not numeric.
    """
    if not isinstance(concentrations, dict):
        raise TypeError("concentrations must be a dictionary.")

    assigned = set()
    for fid in fid_array.get_fids():
        if fid.species is None or not len(fid.species):
            raise ValueError(f"{fid.id} has no species assigned to its peaks.")
        assigned.update(s for s in fid.species if s is not None)

    unknown = set(concentrations) - assigned
    if unknown:
        raise ValueError(f"Species assigned nowhere in the array: {sorted(unknown)}")

    n_points = len(fid_array.t)
    for species, values in concentrations.items():
        if len(values) != n_points:
            raise ValueError(
                f"{species} has {len(values)} values for {n_points} time points."
            )
        if not all(isinstance(value, (int, float)) for value in values):
            raise ValueError(f"{species} holds a non-numeric value.")

    fid_array._FidArray__concentrations = concentrations
