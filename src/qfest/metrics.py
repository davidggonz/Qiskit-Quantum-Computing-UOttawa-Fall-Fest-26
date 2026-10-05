"""Error metrics vs the ED reference.  (Owner: Ririsha)

Convention: the previous report's "%" (Table 2, and its "<5% criterion") is
    100 * max_t |measured - exact| / max_t |exact|,
with t running over dt, 2dt, ..., T (t = 0 excluded). Taken from the old project's
src/run_dt_convergence.py and verified to reproduce its Table 2 exactly; see
experiments/table2_trotter.py. Use `max_report_deviation` when comparing with the report;
`percent_deviation` is the pointwise relative error.
"""
import numpy as np


def max_report_deviation(measured, exact):
    """100 * max|measured - exact| / max|exact|, the previous report's convention.

    Pass curves sampled at t = dt..T (without t = 0) to match the report's numbers.
    """
    measured, exact = np.asarray(measured, float), np.asarray(exact, float)
    return float(100.0 * np.max(np.abs(measured - exact)) / np.max(np.abs(exact)))


def percent_deviation(measured, exact):
    """Relative deviation in % (guards tiny denominators). Not the report's convention."""
    measured, exact = np.asarray(measured, float), np.asarray(exact, float)
    return 100.0 * np.abs(measured - exact) / np.maximum(np.abs(exact), 1e-12)


def max_percent_deviation(measured, exact):
    return float(np.max(percent_deviation(measured, exact)))


def rmse(measured, exact):
    d = np.asarray(measured, float) - np.asarray(exact, float)
    return float(np.sqrt(np.mean(d ** 2)))
