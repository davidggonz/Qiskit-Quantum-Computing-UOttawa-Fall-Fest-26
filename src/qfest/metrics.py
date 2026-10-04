"""Error metrics vs the ED reference.  (Owner: Ririsha)"""
import numpy as np


def percent_deviation(measured, exact):
    """Max-relative-style deviation in %, as in the report's Table 2 (guards tiny denominators)."""
    measured, exact = np.asarray(measured, float), np.asarray(exact, float)
    return 100.0 * np.abs(measured - exact) / np.maximum(np.abs(exact), 1e-12)


def max_percent_deviation(measured, exact):
    return float(np.max(percent_deviation(measured, exact)))


def rmse(measured, exact):
    d = np.asarray(measured, float) - np.asarray(exact, float)
    return float(np.sqrt(np.mean(d ** 2)))
