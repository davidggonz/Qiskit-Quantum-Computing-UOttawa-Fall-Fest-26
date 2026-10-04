
CSS = """
@page{size:Letter;margin:0.45in 0.55in}
body{font-family:Helvetica,Arial,sans-serif;font-size:9.8pt;line-height:1.32;color:#1a1a1a;margin:0}
h1{font-size:16pt;margin:0 0 1px;color:#0f2a4a}
.sub{color:#555;font-size:9pt;margin-bottom:6px;border-bottom:1.5px solid #0f2a4a;padding-bottom:3px}
h2{font-size:10.8pt;margin:6px 0 2px;color:#0f2a4a}
p{margin:2px 0} ul,ol{margin:2px 0 2px 15px;padding:0} li{margin:0}
table{border-collapse:collapse;width:100%;font-size:9pt;margin:3px 0}
th,td{border:1px solid #b9c4d2;padding:2px 4px;text-align:left;vertical-align:top}
th{background:#e6edf5}
"""

EN = dict(file="Proposal_EN.pdf", html="""
<h1>Hardware-Aware Quantum Materials Simulation on IBM Quantum</h1>
<div class="sub">Project proposal · UOttawa Qiskit Fall Fest 26 · 3–4 day build</div>
<h2>Goal</h2>
<p>Improve the documented limitations of our previous Quantinuum project (1D transverse-field Ising model and 2D Fermi–Hubbard; exact diagonalization, Trotter, VQE, ZNE, Iceberg error detection). We build a hardware-aware simulation and error-mitigation workflow on IBM Quantum hardware and test whether it is more accurate and efficient than the previous approach.</p>
<h2>Why Iceberg on heavy-hex is not our main path</h2>
<ul>
<li>Every logical gate uses a hub qubit (X̄ᵢ = X<sub>t</sub>Xᵢ, Z̄ᵢ = Z<sub>b</sub>Zᵢ) and both syndrome ancillas must reach all data qubits: a star graph. Heavy-hex qubits have at most 3 neighbours, so each Trotter step needs O(k)-depth SWAP networks.</li>
<li>The code has distance 2, and a failing SWAP is a correlated two-qubit error, so SWAPs create errors the code cannot detect.</li>
<li>Iceberg's advantage relies on all-to-all connectivity and native MS gates (trapped ions); even there the previous report saw 47–95% discard rates.</li>
</ul>
<p><b>Decision:</b> use hardware-native detection and mitigation on IBM hardware; run Iceberg only as a simulated cost-benefit study (all-to-all vs heavy-hex).</p>
<h2>Approach: limitation → fix → evidence</h2>
<table><tr><th>Previous limitation</th><th>Planned fix (IBM / Qiskit)</th><th>Evidence of improvement</th></tr>
<tr><td>Iceberg discard rate 47–95%</td><td>Ancilla-free symmetry post-selection (particle number for Hubbard, parity for TFIM)</td><td>Discard rate and error vs report</td></tr>
<tr><td>ZNE overcorrects after t ≈ 1.5</td><td>Runtime Estimator: Pauli twirling, dynamical decoupling, TREX readout, ZNE with exponential extrapolation</td><td>Error vs time vs linear 3-point ZNE</td></tr>
<tr><td>Trotter error at h/J = 2 (8.46% at Δt = 0.05)</td><td>Multi-product formulas (qiskit-addon-mpf)</td><td>&lt;5% at h/J = 2 with no deeper circuits</td></tr>
<tr><td>Depth and gate-count limits</td><td>Native RZZ gates, edge-coloured layers, 12-qubit ring mapped natively onto heavy-hex</td><td>2-qubit gate count and depth tables</td></tr>
<tr><td>Iceberg + ZNE cost-benefit left open</td><td>Aer with Heron-like noise; Iceberg transpiled to all-to-all vs heavy-hex</td><td>SWAP overhead and break-even depth plot</td></tr>
<tr><td>VQE does not converge</td><td>Deprioritized; adiabatic preparation + MPF (SQD as stretch goal)</td><td>n/a</td></tr></table>
<h2>3–4 day plan</h2>
<ol>
<li><b>Day 1:</b> Qiskit port of exact diagonalization and Trotter circuits, validated against the previous report's Tables 1–2; noise model, mitigation pipeline, MPF code.</li>
<li><b>Day 2:</b> Simulation studies (ZNE extrapolators, MPF vs 2nd-order Trotter, Iceberg topology study); freeze hardware circuits.</li>
<li><b>Day 3:</b> Batched hardware runs for the TFIM; 2×2 Hubbard with post-selection if budget allows.</li>
<li><b>Day 4:</b> Plots, comparison table vs the previous report, write-up.</li>
</ol>
<h2>Risks and early checks</h2>
<ul>
<li><b>QPU time:</b> Open Plan gives ~10 min/month; confirm the hackathon allocation and batch jobs.</li>
<li><b>Fractional gates:</b> verify twirling/PEA support with native RZZ; keep a CZ/ECR fallback.</li>
<li><b>2D Hubbard:</b> Jordan–Wigner strings cost many SWAPs on heavy-hex; limit to 2×2 or use a 1D chain.</li>
</ul>
<h2>Deliverables</h2>
<p>Reproducible Qiskit repository, hardware and simulation results, and a table comparing accuracy, discard rate and circuit cost against the previous Quantinuum results: lower Trotter error without deeper circuits, better mitigation than linear ZNE, and a quantified answer on when Iceberg pays off on heavy-hex.</p>
""")

ES = dict(file="Propuesta_ES.pdf", html="""
<h1>Simulación Cuántica de Materiales Adaptada al Hardware en IBM Quantum</h1>
<div class="sub">Propuesta de proyecto · UOttawa Qiskit Fall Fest 26 · desarrollo en 3–4 días</div>
<h2>Objetivo</h2>
<p>Mejorar las limitaciones documentadas de nuestro proyecto anterior en Quantinuum (modelo de Ising con campo transversal 1D y Fermi–Hubbard 2D; diagonalización exacta, Trotter, VQE, ZNE, detección de errores Iceberg). Construimos un flujo de simulación y mitigación de errores adaptado al hardware de IBM Quantum y comprobamos si es más preciso y eficiente que el enfoque previo.</p>
<h2>Por qué Iceberg en heavy-hex no es nuestra vía principal</h2>
<ul>
<li>Cada puerta lógica usa un qubit «hub» (X̄ᵢ = X<sub>t</sub>Xᵢ, Z̄ᵢ = Z<sub>b</sub>Zᵢ) y los dos ancillas de síndrome deben alcanzar todos los qubits de datos: un grafo en estrella. Los qubits de heavy-hex tienen como máximo 3 vecinos, así que cada paso de Trotter requiere redes de SWAP de profundidad O(k).</li>
<li>El código tiene distancia 2 y un SWAP defectuoso es un error correlacionado de dos qubits, por lo que los SWAP generan errores que el código no puede detectar.</li>
<li>La ventaja de Iceberg depende de conectividad total y puertas MS nativas (iones atrapados); incluso así, el informe previo observó descartes del 47–95 %.</li>
</ul>
<p><b>Decisión:</b> usar detección y mitigación nativas en el hardware de IBM; ejecutar Iceberg solo como estudio simulado de costo-beneficio (todos-con-todos vs heavy-hex).</p>
<h2>Enfoque: limitación → solución → evidencia</h2>
<table><tr><th>Limitación previa</th><th>Solución prevista (IBM / Qiskit)</th><th>Evidencia de mejora</th></tr>
<tr><td>Descarte de Iceberg 47–95 %</td><td>Postselección por simetría sin ancillas (número de partículas en Hubbard, paridad en TFIM)</td><td>Tasa de descarte y error vs informe</td></tr>
<tr><td>ZNE sobrecorrige tras t ≈ 1.5</td><td>Runtime Estimator: Pauli twirling, desacoplamiento dinámico, lectura TREX, ZNE con extrapolación exponencial</td><td>Error vs tiempo frente a ZNE lineal de 3 puntos</td></tr>
<tr><td>Error de Trotter en h/J = 2 (8.46 % con Δt = 0.05)</td><td>Fórmulas multiproducto (qiskit-addon-mpf)</td><td>&lt;5 % en h/J = 2 sin circuitos más profundos</td></tr>
<tr><td>Límites de profundidad y de puertas</td><td>Puertas RZZ nativas, capas coloreadas por aristas, anillo de 12 qubits mapeado de forma nativa en heavy-hex</td><td>Tablas de puertas de 2 qubits y profundidad</td></tr>
<tr><td>Costo-beneficio Iceberg + ZNE pendiente</td><td>Aer con ruido tipo Heron; Iceberg transpilado a todos-con-todos vs heavy-hex</td><td>Gráfica de sobrecosto de SWAP y profundidad de equilibrio</td></tr>
<tr><td>VQE no converge</td><td>Se desprioriza; preparación adiabática + MPF (SQD opcional)</td><td>n/a</td></tr></table>
<h2>Plan de 3–4 días</h2>
<ol>
<li><b>Día 1:</b> Port a Qiskit de la diagonalización exacta y los circuitos de Trotter, validado con las Tablas 1–2 del informe previo; modelo de ruido, pipeline de mitigación y código MPF.</li>
<li><b>Día 2:</b> Estudios de simulación (extrapoladores de ZNE, MPF vs Trotter de 2.º orden, estudio de topología de Iceberg); congelar los circuitos de hardware.</li>
<li><b>Día 3:</b> Ejecuciones por lotes en hardware para el TFIM; Hubbard 2×2 con postselección si el presupuesto lo permite.</li>
<li><b>Día 4:</b> Gráficas, tabla comparativa con el informe previo y redacción.</li>
</ol>
<h2>Riesgos y verificaciones tempranas</h2>
<ul>
<li><b>Tiempo de QPU:</b> el Open Plan da ~10 min/mes; confirmar la asignación del hackathon y agrupar trabajos.</li>
<li><b>Puertas fraccionarias:</b> verificar twirling/PEA con RZZ nativas; tener plan alterno con CZ/ECR.</li>
<li><b>Hubbard 2D:</b> las cadenas de Jordan–Wigner cuestan muchos SWAP en heavy-hex; limitar a 2×2 o usar una cadena 1D.</li>
</ul>
<h2>Entregables</h2>
<p>Repositorio reproducible en Qiskit, resultados de hardware y simulación, y una tabla que compare precisión, descarte y costo de circuito con los resultados previos de Quantinuum: menor error de Trotter sin circuitos más profundos, mejor mitigación que el ZNE lineal y una respuesta cuantificada sobre cuándo conviene Iceberg en heavy-hex.</p>
""")


for d in (EN,ES):
    open(d["file"].replace(".pdf",".html"),"w",encoding="utf-8").write(f"<html><head><meta charset='utf-8'><style>{CSS}</style></head><body>{d['html']}</body></html>")
