"""Reproduce Table 2 of the previous report: max Mzz deviation of 2nd-order Trotter vs ED.

N = 6 periodic chain, quench from |0...0>, T = 20, deviation = 100 x max_t |dMzz| / max_t |Mzz_ED|
over t = dt..T (the report's convention, taken from the old project's run_dt_convergence.py). Run:  python experiments/table2_trotter.py
"""
import numpy as np

from qfest.ed import quench
from qfest.metrics import max_report_deviation, max_percent_deviation
from qfest.results import save
from qfest.tfim import trotter_curve

N, J, T = 6, 1.0, 20.0
HS = [0.5, 1.0, 2.0]
DTS = [0.2, 0.1, 0.05, 0.025, 0.0125]
REPORT = {  # previous report, Table 2
    0.5: [4.47, 1.16, 0.29, 0.073, 0.018],
    1.0: [50.53, 12.66, 3.12, 0.774, 0.193],
    2.0: [76.30, 31.75, 8.46, 2.100, 0.523],
}


def main():
    ours, rel = {}, {}
    for h in HS:
        ours[h], rel[h] = [], []
        for dt in DTS:
            steps = int(round(T / dt))
            tr = trotter_curve(N, J, h, dt, steps, order=2)
            ex = quench(N, tr["t"], J, h)
            # drop t = 0 like the old code, which sampled t = dt..T
            ours[h].append(max_report_deviation(tr["Mzz"][1:], ex["Mzz"][1:]))
            rel[h].append(max_percent_deviation(tr["Mzz"][1:], ex["Mzz"][1:]))

    print("Max Mzz deviation in % (report convention); report value in brackets")
    print("h/J  | " + " | ".join(f"dt={dt}" for dt in DTS))
    for h in HS:
        print(f"{h:<4} | " + " | ".join(f"{o:.3f} ({r})" for o, r in zip(ours[h], REPORT[h])))

    path = save({
        "experiment": "table2_trotter",
        "backend": "statevector",
        "params": {"n": N, "J": J, "T": T, "h": HS, "dt": DTS, "order": 2, "periodic": True},
        "method": "trotter2",
        "observables": {},
        "extra": {
            "max_report_dev_percent": {str(h): ours[h] for h in HS},
            "max_rel_dev_percent": {str(h): rel[h] for h in HS},
            "report_table2": {str(h): REPORT[h] for h in HS},
        },
    }, "table2_trotter")
    print(f"saved {path}")


if __name__ == "__main__":
    main()
