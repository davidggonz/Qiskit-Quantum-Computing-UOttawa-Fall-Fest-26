"""Noisy TFIM quench + ZNE on a (fake or real) IBM backend, compared against ED.

Pipeline (ported from Boutaina's hardware.py with the layout/folding/ISA fixes, see
src/qfest/hardware.py): logical circuit -> zero-SWAP layout -> ISA circuit -> folds 1,3,5 on the
same qubits -> counts (Aer noise or hardware) -> Zavg / Mzz / RMS Mz -> linear, Richardson and
exponential extrapolation -> JSON in results/.

Examples:
  python experiments/noisy_zne.py                                  # N=12 ring, FakeMarrakesh noise
  python experiments/noisy_zne.py --n 6 --open-chain --h 2 --t 2 --steps 8
  python experiments/noisy_zne.py --noise simple                   # depolarizing model
  python experiments/noisy_zne.py --hardware --backend ibm_marrakesh   # real QPU (uses QPU time!)
"""
import argparse

import numpy as np

from qfest import hardware as hw
from qfest.ed import quench_state
from qfest.extrapolation import EXTRAPOLATORS
from qfest.noise import simulator
from qfest.results import save
from qfest.shots import counts_expectations, state_expectations
from qfest.tfim import tfim_circuit, trotter_state

KEYS = ("Zavg", "Mzz", "Mz")


def parse():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n", type=int, default=12)
    p.add_argument("--open-chain", action="store_true")
    p.add_argument("--J", type=float, default=1.0)
    p.add_argument("--h", type=float, default=1.0)
    p.add_argument("--t", type=float, default=1.0, help="total evolution time")
    p.add_argument("--steps", type=int, default=4)
    p.add_argument("--scales", type=int, nargs="+", default=[1, 3, 5])
    p.add_argument("--shots", type=int, default=4000)
    p.add_argument("--opt-level", type=int, default=1)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--fake", default="marrakesh", help="fake backend for layout + noise")
    p.add_argument("--noise", choices=("backend", "simple", "ideal"), default="backend")
    p.add_argument("--hardware", action="store_true", help="run on a real IBM backend")
    p.add_argument("--backend", default=None, help="real backend name (with --hardware)")
    return p.parse_args()


def main():
    a = parse()
    periodic = not a.open_chain
    backend = hw.get_backend(name=a.backend) if a.hardware else hw.get_backend(fake=a.fake)
    layout = hw.find_layout(backend, a.n, periodic)
    logical = tfim_circuit(a.n, a.J, a.h, a.t / a.steps, a.steps, periodic=periodic)
    isa = hw.to_isa(logical, backend, layout, a.opt_level, a.seed)
    report = hw.transpile_report(logical, isa)
    circuits, factors = hw.zne_circuits(isa, a.n, backend, a.scales)
    ideal = state_expectations(trotter_state(a.n, a.J, a.h, a.t, a.steps, periodic=periodic), a.n, periodic)
    exact = state_expectations(quench_state(a.n, a.t, a.J, a.h, periodic), a.n, periodic)

    if a.hardware:
        counts, job_id = hw.run_batch(circuits, backend, a.shots)
        source = backend.name
    else:
        counts, job_id = hw.run_simulator(circuits, simulator(backend, a.noise), a.shots, a.seed), None
        source = f"aer_{a.noise}_{backend.name}"

    raw = [counts_expectations(c, a.n, periodic) for c in counts]
    mitigated = {k: {name: f(factors, [r[k] for r in raw]) for name, f in EXTRAPOLATORS.items()} for k in KEYS}

    print(f"{source}: N={a.n} {'ring' if periodic else 'chain'} h/J={a.h} t={a.t} steps={a.steps} "
          f"layout={layout}")
    print(f"ISA: depth {report['depth']}, 2q gates {report['two_qubit']} "
          f"(extra vs native: {report['extra_two_qubit']}), noise factors {factors}")
    print(f"{'obs':<5}{'ED':>8}{'ideal':>8}" + "".join(f"{'raw x' + str(s):>9}" for s in a.scales)
          + "".join(f"{name:>12}" for name in EXTRAPOLATORS))
    for k in KEYS:
        print(f"{k:<5}{exact[k]:>8.4f}{ideal[k]:>8.4f}" + "".join(f"{r[k]:>9.4f}" for r in raw)
              + "".join(f"{mitigated[k][name]:>12.4f}" for name in EXTRAPOLATORS))

    path = save({
        "experiment": "noisy_zne",
        "backend": source,
        "params": {"n": a.n, "periodic": periodic, "J": a.J, "h": a.h, "t": a.t, "steps": a.steps,
                   "scales": a.scales, "shots": a.shots, "opt_level": a.opt_level, "seed": a.seed},
        "method": "trotter2_zne",
        "t": [a.t],
        "observables": {"ed": exact, "ideal_trotter": ideal, "raw": raw, "mitigated": mitigated},
        "extra": {"layout": layout, "noise_factors": factors, "transpile": report, "job_id": job_id,
                  "counts": counts},
    }, f"noisy_zne_{source}_n{a.n}_{'ring' if periodic else 'chain'}_h{a.h}_t{a.t}_k{a.steps}")
    print(f"saved {path}")


if __name__ == "__main__":
    main()
