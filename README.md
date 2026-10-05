# TFIM on IBM Quantum Hardware — hardware.py

Trotterized Transverse-Field Ising Model with transpilation comparison,
hardware-aware layout, ZNE mitigation, Aer + IBM hardware execution,
exact-diagonalization analysis. Qiskit 1.0+ (tested 2.5.2).

## Quick start (cloud — recommended for Fall Fest)

Local Windows Smart App Control can block `qiskit._accelerate.pyd`
(unsigned native lib). If `import qiskit` fails locally, run in Colab:

```python
!pip install "qiskit[visualization]==2.5.2" "qiskit-aer==0.17.2" \
  "qiskit-machine-learning==0.9.1" "qiskit-ibm-runtime==0.50.0" \
  "scipy==1.18.1" "scikit-learn==1.9.1" "matplotlib==3.11.2" "pandas==3.0.6"
```

Upload `hardware.py`, then open `demo.ipynb`.

## Local setup

```powershell
py -3.12 -m venv qff26
.\qff26\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install "qiskit[visualization]==2.5.2" "qiskit-aer==0.17.2" ^
  "qiskit-machine-learning==0.9.1" "qiskit-ibm-runtime==0.50.0" ^
  "scipy==1.18.1" "scikit-learn==1.9.1" "matplotlib==3.11.2" "pandas==3.0.6"
copy .env.example .env
# edit .env with your PINQ2 token (Fall Fest shared allocation, keep private)
```

## Usage

```powershell
# Aer only (no token needed)
python hardware.py --n 4 --steps 3 --shots 4096 --no-hardware --outdir results

# With IBM hardware (reads .env)
python hardware.py --n 4 --steps 3 --shots 4096 --outdir results
python hardware.py --backend ibm_marrakesh --n 5 --outdir results
```

## What it does

1. `get_service()` / `select_backend()` — `QiskitRuntimeService(channel, token)`,
   prefers `ibm_marrakesh`/`ibm_fez`; prints qubits, coupling map, basis gates,
   median 1q/2q/readout errors.
2. `build_tfim_circuit(N,J,h,dt,steps)` — RZZ(-2J·dt) + RX(-2h·dt) Trotter layers.
   `build_tfim_observables()` — energy + magnetization `SparsePauliOp`s.
3. `transpile_compare()` — levels 0–3: depth, gate counts, CX-like, est. error
   `1-Π(1-e)` from `backend.target`.
4. `select_best_qubits()` — greedy lowest-error connected chain from
   `backend.target` + `coupling_map`; `hardware_aware_transpile()` uses it as
   `initial_layout` (lv3) vs default (lv1).
5. ZNE — Runtime: `EstimatorV2(resilience_level=2, zne_mitigation=True)`;
   Aer: manual unitary folding of CX/CZ/ECR at scales 1/3/5 + linear/exponential
   extrapolation to zero noise.
6. Execution — Aer `AerSimulator`/`EstimatorV2`; hardware `SamplerV2`/`EstimatorV2`
   in a `Session`. Returns counts + expectation values.
7. Analysis — exact magnetization via `scipy.linalg.expm` vs Aer/hardware;
   plots `accuracy_vs_level.png`, `depth_vs_level.png`, `error_vs_mitigation.png`;
   results in `results/tfim_results.json`.

## Cost note

Real QPUs queue and consume QPU-seconds. Start Aer-only, then 1 short hardware
job (few thousand shots). Shared Fall Fest allocation: use only for event jobs.

## Files

- `hardware.py` — full module + CLI
- `demo.ipynb` — Aer demo + hardware snippet (skips cleanly without token)
- `.env.example` — token template
