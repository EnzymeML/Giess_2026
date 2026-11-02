"""Helper modules for the Pd(0)/BINOL-catalyzed enyne cycloisomerization (Scenario 1).

The notebooks in this directory hand state to each other through files in ``data/``
rather than through Python imports. This package holds the helper functions the
analysis notebook uses, so that the notebook reads as a documented pipeline instead
of carrying its own copy of the implementation.

Modules:
    assignment: Scripted species assignment, in place of the interactive widget.
    concentrations: Internal-standard math.
    numeric_integration: Numeric quadrature of whole peak ranges, and the
        concentrations derived from it.
    phasing: Applying recorded per-FID phase angles.
    quantification: Recording a numeric integration in the NMRPy data model.
    io: Writing the generated artifacts, and displaying models inside a notebook.
"""
