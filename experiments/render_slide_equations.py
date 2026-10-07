"""Render the slide equations as transparent LaTeX-style PNGs (docs/slide_equations/).

Run:  python experiments/render_slide_equations.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["mathtext.fontset"] = "cm"
OUT = "docs/slide_equations/"
NAVY = "#1e344c"
EQS = {
 "eq1_tfim_hamiltonian": r"$H = -J\,\sum_{i} Z_i Z_{i+1} - h\,\sum_{i} X_i$",
 "eq2_folding": r"$\mathrm{CZ} \rightarrow \mathrm{CZ}\,(\mathrm{CZ}^{\dagger}\mathrm{CZ})^{k}$",
 "eq3_iceberg_logical": r"$\bar{X}_i = X_t X_i, \quad \bar{Z}_i = Z_b Z_i$",
 "eq4_iceberg_checks": r"$S_X = X^{\otimes n}, \quad S_Z = Z^{\otimes n}$",
 "eq5_exp_model": r"$\langle O\rangle(\lambda) \approx a\,e^{-b\lambda} + c$",
 "eq6_noisy_limits": r"$M_z \to \sqrt{1/N} = 0.29, \quad M_x = M_{zz} \to 0$",
}
for name, tex in EQS.items():
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, tex, fontsize=40, color=NAVY)
    fig.savefig(OUT + name + ".png", dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(name)
