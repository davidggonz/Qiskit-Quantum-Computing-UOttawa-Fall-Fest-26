"""Error metrics vs the ED reference.  (Owner: Ririsha)

Convention: the previous report's "%" in Table 2 is 100 x the ABSOLUTE error (verified by
reproducing the table, see experiments/table2_trotter.py). Use `report_deviation` when
comparing with the report; `percent_deviation` is the true relative error.
"""
import numpy as np


def report_deviation(measured, exact):
    """100 * |measured - exact|, the previous report's convention."""
    return 100.0 * np.abs(np.asarray(measured, float) - np.asarray(exact, float))


def max_report_deviation(measured, exact):
    return float(np.max(report_deviation(measured, exact)))


def percent_deviation(measured, exact):
    """Relative deviation in % (guards tiny denominators). Not the report's convention."""
    measured, exact = np.asarray(measured, float), np.asarray(exact, float)
    return 100.0 * np.abs(measured - exact) / np.maximum(np.abs(exact), 1e-12)


def max_percent_deviation(measured, exact):
    return float(np.max(percent_deviation(measured, exact)))


def rmse(measured, exact):
    d = np.asarray(measured, float) - np.asarray(exact, float)
    return float(np.sqrt(np.mean(d ** 2)))
