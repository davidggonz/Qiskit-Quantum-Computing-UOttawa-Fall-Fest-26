"""MPF vs plain 2nd-order Trotter at short times, at equal maximum circuit depth.

For each time t on a grid in (0, t_max], each method uses Trotter step size t/k:
  * Trotter:  one circuit with k = k_max steps
  * MPF:      circuits with k in `ks` (max(ks) = k_max), combined with Richardson weights
We report, per h/J:
  * noiseless error vs ED (statevector), and
  * error with a finite shot budget (same TOTAL shots for every method; MPF splits them
    across its circuits proportional to |x_j|), averaged over `reps` repetitions.
MPF coefficients: exact Richardson (= qiskit-addon-mpf exact LSE) and, with --l1, the addon's
well-conditioned sum-of-squares coefficients with ||x||_1 bounded (skipped when not binding).
Error metric: the previous report's convention, 100 * max_t |dMzz| / max_t |Mzz_ED|,
plus the full error-vs-time curves in the saved JSON.

Run:  python experiments/mpf_short_times.py            (defaults below, ~1 min)
      python experiments/mpf_short_times.py --n 8 --h 2 --kmax 10 --ks 4,6,8,10 --ks 6,10
"""
import argparse

import numpy as np

from qfest.ed import quench
from qfest.metrics import max_report_deviation
from qfest.mpf import richardson_coefficients, approx_coefficients, noise_amplification
from qfest.results import save
from qfest.shots import z_diagonals, sample_estimate, allocate_shots
from qfest.tfim import trotter_state


def parse():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n", type=int, default=6)
    p.add_argument("--J", type=float, default=1.0)
    p.add_argument("--h", type=float, nargs="+", default=[0.5, 1.0, 2.0])
    p.add_argument("--tmax", type=float, default=2.0)
    p.add_argument("--nt", type=int, default=20, help="number of time points in (0, tmax]")
    p.add_argument("--kmax", type=int, default=8, help="max Trotter steps (depth budget)")
    p.add_argument("--ks", action="append", default=None,
                   help="comma-separated MPF step counts; repeat flag for several sets")
    p.add_argument("--shots", type=int, nargs="+", default=[1000, 10000],
                   help="total shot budgets per time point")
    p.add_argument("--l1", type=float, nargs="*", default=[1.5, 2.0],
                   help="also run well-conditioned qiskit-addon-mpf coefficients with ||x||_1 <= each value")
    p.add_argument("--reps", type=int, default=30)
    p.add_argument("--seed", type=int, default=1234)
    a = p.parse_args()
    a.ks = [[int(k) for k in s.split(",")] for s in (a.ks or ["4,6,8", "6,8", "2,4,8"])]
    for ks in a.ks:
        if max(ks) != a.kmax:
            p.error(f"MPF set {ks} must have max(ks) == kmax ({a.kmax}) for an equal-depth comparison")
    return a


def run_h(a, h, rng):
    times = np.linspace(a.tmax / a.nt, a.tmax, a.nt)
    exact = quench(a.n, times, a.J, h)["Mzz"]
    diag = z_diagonals(a.n)["Mzz"]
    all_k = sorted({a.kmax, *[k for ks in a.ks for k in ks]})
    psi = {(i, k): trotter_state(a.n, a.J, h, t, k) for i, t in enumerate(times) for k in all_k}
    ideal = {(i, k): float(np.real(np.vdot(v, diag * v))) for (i, k), v in psi.items()}

    methods = {"trotter": {"ks": [a.kmax], "x": np.array([1.0])}}
    for ks in a.ks:
        tag = "mpf_" + "-".join(map(str, ks))
        methods[tag] = {"ks": ks, "x": richardson_coefficients(ks)}
        for l1 in a.l1:
            x = approx_coefficients(ks, l1)
            if noise_amplification(x) < noise_amplification(methods[tag]["x"]) - 1e-6:
                methods[f"{tag}_l1<={l1:g}"] = {"ks": ks, "x": x}

    out = {}
    for name, m in methods.items():
        curve = np.array([sum(x * ideal[(i, k)] for x, k in zip(m["x"], m["ks"]))
                          for i in range(len(times))])
        rec = {"ks": m["ks"], "coefficients": m["x"].tolist(),
               "noise_amplification": noise_amplification(m["x"]),
               "Mzz_noiseless": curve.tolist(),
               "noiseless_report_dev": max_report_deviation(curve, exact),
               "shots": {}}
        for S in a.shots:
            alloc = allocate_shots(m["x"], S)
            devs, abs_err = [], np.zeros(len(times))
            for _ in range(a.reps):
                est = np.array([sum(x * sample_estimate(psi[(i, k)], diag, s, rng)
                                    for x, k, s in zip(m["x"], m["ks"], alloc))
                                for i in range(len(times))])
                devs.append(max_report_deviation(est, exact))
                abs_err += (est - exact) ** 2
            rec["shots"][str(S)] = {
                "allocation": alloc.tolist(),
                "report_dev_mean": float(np.mean(devs)),
                "report_dev_std": float(np.std(devs)),
                "rms_error_vs_t": np.sqrt(abs_err / a.reps).tolist(),
            }
        out[name] = rec
    return times, exact, out


def main():
    a = parse()
    rng = np.random.default_rng(a.seed)
    for h in a.h:
        times, exact, out = run_h(a, h, rng)
        print(f"\nN={a.n}  h/J={h}  t in (0, {a.tmax}]  k_max={a.kmax}  "
              f"(error = report %, mean over {a.reps} reps)")
        hdr = f"{'method':<22}{'||x||1':>8}{'noiseless':>11}" + "".join(f"{'S=' + str(S):>16}" for S in a.shots)
        print(hdr)
        for name, r in out.items():
            cells = "".join(f"{r['shots'][str(S)]['report_dev_mean']:>9.3f} ±{r['shots'][str(S)]['report_dev_std']:<5.2f}"
                            for S in a.shots)
            print(f"{name:<22}{r['noise_amplification']:>8.2f}{r['noiseless_report_dev']:>11.4f}{cells}")
        path = save({
            "experiment": "mpf_short_times",
            "backend": "statevector+shot_sampling",
            "params": {"n": a.n, "J": a.J, "h": h, "tmax": a.tmax, "nt": a.nt, "kmax": a.kmax,
                       "shots": a.shots, "reps": a.reps, "seed": a.seed, "order": 2},
            "method": "trotter2_vs_mpf",
            "t": times,
            "observables": {"Mzz_ed": exact},
            "extra": {"methods": out},
        }, f"mpf_short_times_n{a.n}_h{h}")
        print(f"saved {path}")


if __name__ == "__main__":
    main()
