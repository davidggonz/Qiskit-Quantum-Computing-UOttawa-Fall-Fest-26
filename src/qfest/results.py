"""Shared results format.  Every experiment writes one JSON file under results/.

Schema (keep stable so Toto's plotting code works for everyone):
{
  "experiment": "trotter_vs_mpf",     # short id
  "backend": "aer_ideal | aer_noisy | ibm_<name>",
  "params": {"n": 6, "J": 1.0, "h": 0.5, "dt": 0.05, ...},
  "method": "trotter2 | mpf | zne_linear | zne_exp | iceberg | ...",
  "t": [...],
  "observables": {"Mz": [...], "Mx": [...], "Mzz": [...]},
  "extra": {"shots": 4000, "discard_rate": [...], "two_qubit_gates": 123, "depth": 45}
}
"""
import json
import pathlib

RESULTS_DIR = pathlib.Path(__file__).resolve().parents[2] / "results"


def save(record, name):
    RESULTS_DIR.mkdir(exist_ok=True)
    p = RESULTS_DIR / f"{name}.json"
    p.write_text(json.dumps(record, indent=2, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    return p


def load(name):
    return json.loads((RESULTS_DIR / f"{name}.json").read_text())
