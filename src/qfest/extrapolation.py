"""ZNE extrapolators.  (Owner: Ririsha)

All take noise factors `lams` (e.g. [1, 3, 5]) and measured values `vals` and
return the zero-noise estimate.  The previous report used `linear` only.
"""
import numpy as np
from scipy.optimize import curve_fit


def linear(lams, vals):
    return float(np.polyval(np.polyfit(lams, vals, 1), 0.0))


def richardson(lams, vals):
    """Polynomial of degree len(lams)-1 through all points, evaluated at 0."""
    return float(np.polyval(np.polyfit(lams, vals, len(lams) - 1), 0.0))


def exponential(lams, vals, asymptote=None):
    """Fit a*exp(-b*lam) + c (c fixed to `asymptote` if given); evaluate at lam=0.
    Falls back to `linear` if the fit does not converge."""
    lams, vals = np.asarray(lams, float), np.asarray(vals, float)
    try:
        if asymptote is None:
            f = lambda x, a, b, c: a * np.exp(-b * x) + c
            p0 = [vals[0] - vals[-1], 0.3, vals[-1]]
            p, _ = curve_fit(f, lams, vals, p0=p0, maxfev=10000)
            return float(f(0.0, *p))
        f = lambda x, a, b: a * np.exp(-b * x) + asymptote
        p, _ = curve_fit(f, lams, vals, p0=[vals[0] - asymptote, 0.3], maxfev=10000)
        return float(f(0.0, *p))
    except (RuntimeError, ValueError):
        return linear(lams, vals)


EXTRAPOLATORS = {"linear": linear, "richardson": richardson, "exponential": exponential}
