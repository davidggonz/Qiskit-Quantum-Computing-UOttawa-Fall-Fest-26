# Quantum Hardware Simulation: Transverse-Field Ising Model (TFIM)

This project implements a hardware-side module for simulating the **Transverse-Field Ising Model (TFIM)** on IBM Quantum hardware using Trotterization. The goal is to evaluate the impact of different transpilation strategies and error mitigation techniques on simulation accuracy.

## 🚀 Features

- **Trotterized TFIM Circuit Generation**: Constructs the canonical periodic-ring TFIM circuit for even $N$, using symmetric second-order Trotterization and fused transverse-field rotations.
- **Transpilation Analysis**: Compares Qiskit's optimization levels (0 to 3) by measuring circuit depth, gate counts (specifically CNOTs), and estimating overall circuit error.
- **Hardware-Aware Optimization**: Implements a greedy qubit selection algorithm that picks physical qubits with the lowest readout and gate error rates, and maps the circuit to this optimal layout.
- **Error Mitigation (ZNE)**: Implements Zero-Noise Extrapolation (ZNE) via manual unitary folding to mitigate hardware noise.
- **Hybrid Execution**: Supports execution on `qiskit-aer` (simulators) and real IBM Quantum backends via `qiskit-ibm-runtime`.
- **Analysis Suite**: Calculates exact magnetization using numerical diagonalization as a baseline and generates accuracy and depth plots.

## 🛠️ Installation

### Prerequisites
- Python 3.10+
- An IBM Quantum account and API token.

### Dependencies
Install the required packages:
```bash
pip install "qiskit[visualization]" "qiskit-aer" "qiskit-ibm-runtime" scipy matplotlib pandas
```

## ⚙️ Configuration

1. Create a `.env` file in the root directory:
   ```bash
   cp .env.example .env
   ```
2. Edit the `.env` file and add your IBM Quantum token:
   ```env
   QISKIT_IBM_TOKEN=your_api_token_here
   ```

## 💻 Usage

### Running the Module
You can run the full simulation pipeline from the command line:

```bash
# Run with simulator only
python hardware.py --n 4 --steps 3 --shots 4096 --no-hardware

# Run on a specific IBM Quantum backend
python hardware.py --n 4 --steps 3 --shots 4096 --backend ibm_marrakesh
```

### Using the Notebook
For a step-by-step walkthrough, open `demo.ipynb` in Jupyter or Google Colab. The notebook demonstrates:
1. Backend setup and property inspection.
2. Circuit construction.
3. Transpilation level comparison.
4. Hardware-aware layout selection.
5. Execution and ZNE mitigation.

## 📊 Outputs & Visualization

The pipeline saves results to the root folder:
- `tfim_results.json`: Detailed metadata, transpilation stats, and expectation values.

### Results Preview

| Accuracy vs Optimization Level | Depth vs Optimization Level | Error Mitigation (ZNE) |
| :---: | :---: | :---: |
| ![Accuracy](accuracy_vs_level.png) | ![Depth](depth_vs_level.png) | ![ZNE](error_vs_mitigation.png) |

**N=6, periodic ring, 2nd-order Trotter, fake-backend ibm_marrakesh.**

- **Accuracy vs Level**: Compares noisy magnetization errors for fake-backend transpilation levels 0–3.
- **Depth vs Level**: Shows transpiled depths for optimization levels 0–3.
- **ZNE Plot**: Shows raw noisy samples and their zero-noise linear extrapolation.

### Reproducibility

Regenerate these plots with:

```bash
python hardware.py --n 6 --steps 20 --dt 0.05 --fake-backend ibm_marrakesh --noise simple --outdir results_canonical
```

The local `simple` noise model uses the project's 1% single-qubit and 3% two-qubit depolarizing rates.

## 📝 License
This project is developed for the Qiskit Fall Fest 2026.
