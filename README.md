# Heavy-flavor TEEC in proton-nucleus collisions

Numerical evaluation of the transverse energy-energy correlator (TEEC) associated
with forward heavy-quark pair production in the Color Glass Condensate framework.
The calculation retains the heavy-quark mass and full relative azimuthal-angle
dependence, including interference between gluon splitting before and after
scattering on the target shockwave.

The entry point is [heavy_quark_teec.ipynb](heavy_quark_teec.ipynb). Definitions,
physics inputs, integration, plotting and optional checks are separated into
notebook sections.

## Getting started

Use Python 3.10 or newer with NumPy, SciPy, pandas, Matplotlib, IPython and Jupyter.
For example, install the Python dependencies in your chosen environment with:

```bash
python -m pip install jupyterlab numpy scipy pandas matplotlib
```

Open a terminal in this `TEECS` directory, then launch:

```bash
jupyter lab heavy_quark_teec.ipynb
```

1. Select the Python kernel containing these dependencies.
2. Edit the choices in **Section 1** and any provider options in **Sections 4-5**.
3. Restart the kernel and run all cells in order after changing inputs.
4. Check the resolved configuration and output folder printed in **Section 6**.
5. Inspect the main plots in **Section 12** and the saved files under `results/`.

Keep the notebook's working directory here so that `data/` and `results/` resolve
as intended. A complete run performs numerical integrations; existing result
folders are not automatically loaded to skip them.

**CT18NLO additionally requires an installed LHAPDF C++ library, its data set and
a C++ compiler.** The notebook compiles and caches a small batch interface on
first use. LHAPDF Python bindings are not required. The notebook does not download
LHAPDF or PDF data. Choosing `pdf_name="toy"` avoids this dependency.

## Files

```text
TEECS/
  heavy_quark_teec.ipynb       Calculation and plotting notebook
  README.md                   This guide
  data/
    mv_bk_gamma1.npz           Bundled BK-evolved MV interpolation tables
    mv_bk_gamma1.json          Input parameters, source identities and checksums
    mv_bk_validation.json     Numerical validation and input-variant checks
  results/
    <run_id>/                 Figures, tables, replicas and provenance
```

## Observable and conventions

The angular variable is

$$
\tau=\frac{1+\cos\phi}{2}.
$$

Back-to-back pairs have $\tau=0$ and the near-side limit is $\tau=1$. Curves are
computed at fixed pseudorapidities, with the parent transverse momenta integrated.
The cross-branch pair is counted once and carries the partonic weight

$$
\frac{p_{T,1}p_{T,2}}{(p_{T,1}+p_{T,2})^2}.
$$

The fragmentation cancellation uses complete species and momentum sums in the
two heavy-quark fragmentation branches within the adopted collinear prescription.
Same-branch and self correlations are excluded. Exact final-state **parton**
kinematics does not replace that fragmentation approximation, and the interaction
with the target remains eikonal.

The plotted TEEC is normalized **per target nucleon**, in microbarns per unit
$\tau$, $\eta_1$ and $\eta_2$. Its nuclear modification factor is the ratio of
the per-nucleon nuclear TEEC to the proton TEEC; no further factor of $1/A$ is
needed when using the exported TEEC columns.

## Selecting a calculation

### Main choices: Section 1

| Setting | Meaning |
|---|---|
| `CHOICES["collider"]` | `"RHIC"` or `"LHC"` preset |
| `CHOICES["quark"]` | `"charm"` or `"bottom"` |
| `CHOICES["target_model"]` | `"GBW"`, `"toy_MV"`, `"MV"`, or a registered provider |
| `CHOICES["pdf_name"]` | `"toy"`, `"CT18NLO"`, another installed LHAPDF set, or a registered provider |
| `CHOICES["eta_frame"]` | `"CM"` or `"LAB"` for the supplied pseudorapidities |
| `OBSERVABLE["eta1"]`, `OBSERVABLE["eta2"]` | Fixed pseudorapidities; changing collider does not change them |
| `OBSERVABLE["kinematics"]` | `"exact"` massive parton kinematics or the forward `"eikonal"` approximation |
| `OBSERVABLE["xt_max"]` | Upper bound on the target momentum fraction |
| `OBSERVABLE["p_min"]`, `OBSERVABLE["p_max"]` | Optional parent transverse-momentum bounds, not detector track cuts |
| `C_VALUES` | Nuclear geometry choices; supplied evolved MV tables support `0.5` and `1.0` |
| `PROTON_AREA` | Proton transverse area in GeV$^{-2}$; `None` selects the provider default |

The presets are:

| Collider | Target | $A$ | $\sqrt{s_{NN}}$ | CM rapidity relative to LAB |
|---|---|---:|---:|---:|
| RHIC | Au | 197 | 200 GeV | 0 |
| LHC | Pb | 208 | 5020 GeV | 0.465 |

The positive longitudinal axis is proton-going. The numerical quark masses are
$m_c=1.5$ GeV and $m_b=4.75$ GeV. The paper benchmarks use CM pseudorapidities
$\eta_1=\eta_2=3.5$ for RHIC charm and $\eta_1=\eta_2=4$ for LHC charm/bottom.
RHIC bottom is available as a diagnostic, but its substantial large-$x_t$
sensitivity is why it is excluded from the paper's predictions.

`OVERRIDES` replaces resolved configuration values without changing the preset:

```python
OVERRIDES = {"sqrt_s": 8160.0}  # Example: change only the collision energy
```

For LAB inputs, the code converts the parton light-cone components to the NN-CM
frame before constructing momentum fractions. Set `cm_rapidity` to the value
appropriate to the actual beam configuration. Massive pseudorapidity does not
generally transform by a constant additive shift.

At the time this guide was written, Section 1 selected RHIC charm, `MV`, CT18NLO,
CM pseudorapidities 3.5, exact parton kinematics, `xt_max=1` and `paper=True`.
Treat the editable notebook and each run's provenance as the authoritative
settings.

### Target models: Section 4

| Model | Description |
|---|---|
| `GBW` | Analytic Gaussian momentum-space dipole, with a power-law prescription for $Q_s^2(x_t)$ |
| `toy_MV` | Unevolved MV-type coordinate-space dipole, with the GBW saturation-scale prescription and a cached numerical Hankel transform |
| `MV` | Interpolation of supplied fitted momentum-space dipoles after running-coupling BK evolution from a HERA-constrained MV-inspired initial condition |

`F(k,x)` is the Fourier transform of the dipole **S-matrix** $S(r,x)$, rather than
of $N=1-S$. The target input is per unit transverse area; the production integral
supplies the geometry factors.

Provider defaults are defined in their classes. Override them through, for example:

```python
TARGET_PROVIDER_OPTIONS["GBW"] = {"qs0_sq": 1.0, "lambda": 0.29}
TARGET_PROVIDER_OPTIONS["toy_MV"] = {"lambda_ir": 0.2}
TARGET_PROVIDER_OPTIONS["MV"] = {"table_variant": "fitted_rmin_1e-3"}
```

The bundled `MV` input was constructed from Meijian's shared numerical files.
The notebook loads the evolved tables; it does not rerun BK evolution during a
TEEC calculation. The initial condition has fixed $\gamma=1$, while the smoothing
fit to the evolved dipoles can have a varying effective $\gamma$.

The default `fitted_rmin_1e-3` variant uses the fit with a coordinate-space cutoff
$r\geq0.001006938631475812$ GeV$^{-1}$. `fitted_all_r` selects the alternative
supplied fit. Both variants are packaged in the NPZ file. The tables cover
approximately $1.125\times10^{-9}\leq x_t\leq0.01$ and
$0.01\leq k\leq100$ GeV, with separate proton, Au and Pb solutions for the
supported nuclear geometry choices.

Relevant `MV` options include:

| Option | Default | Alternatives / role |
|---|---|---|
| `table_variant` | `"fitted_rmin_1e-3"` | `"fitted_all_r"` |
| `interpolation` | `"linear"` | `"linear_logk"` for a momentum-interpolation check |
| `uv_tail` | `"matched_power"` | `"power4"`, `"zero"`, `"error"` |
| `x_above_x0` | `"freeze"` | `"error"` to reject $x_t>0.01$ |
| `table_path` | `"data/mv_bk_gamma1.npz"` | Alternative bundled-table path |

Below the momentum grid an even small-$k$ continuation is used. Values below the
minimum tabulated $x_t$ raise an error. With `x_above_x0="freeze"`, contributions
at $x_t>0.01$ use the initial condition: this is a diagnostic continuation, not a
BK prediction at large $x_t$. `xt_max=0.01` restricts the integral to the small-$x$
range; `xt_max=1` also includes that continuation wherever phase space allows it.

The nuclear density label is $\rho_A=cA^{1/3}$; the proton has $\rho_p=1$.
The geometry uses $S_{\perp,A}=A^{2/3}S_{\perp,p}/c$. With `PROTON_AREA=None`,
the evolved `MV` provider uses $S_{\perp,p}=36.211479\  \mathrm{GeV}^{-2}$, and `GBW` /
`toy_MV` use $58\ \mathrm{GeV}^{-2}$. For evolved MV, a density selects a separate evolved
input; it does not simply rescale the evolved proton distribution.

### PDF and hard coupling: Section 5

The PDF interface returns $xg(x,\mu_F^2)$. The toy model is
$Nx^{-a}(1-x)^b$, with its parameters in `PDF_PROVIDER_OPTIONS["toy"]`.

For LHAPDF, `member=0` is the central fit. Changing the member selects a different
PDF fit, not a different parton species; the gluon flavor ID is 21. A single
member run does not produce a PDF uncertainty band.

Normally `PDF_PROVIDER_OPTIONS["CT18NLO"] = {}` is sufficient. If discovery fails,
specify the installed locations:

```python
PDF_PROVIDER_OPTIONS["CT18NLO"] = {
    "member": 0,
    "lhapdf_config": "/path/to/bin/lhapdf-config",
    "data_path": "/path/to/share/LHAPDF",
}
```

The `pair_mt` prescription sets
$\mu_F^2=m_Q^2+(p_{T,1}^2+p_{T,2}^2)/2$. The `fixed` prescription uses
`pdf_mu_fixed`, which defaults to $2m_Q$ in GeV. These scales are separate from
the provider's PDF shape parameters.

`COUPLING_PARAMETERS` controls the hard-splitting coupling. With `alpha_s_mode="auto"`,
CT18NLO supplies its running coupling through LHAPDF at $\mu_R$; the toy PDF uses
the fixed `OBSERVABLE["alpha_s"]` value. By default `ren_scale="pdf_scale"` gives
$\mu_R=\mu_F$. Fixed or independent `pair_mt` renormalization scales and a
`ren_scale_factor` are available. The BK coordinate-space coupling is already
encoded in the evolved target tables and is independent of this hard coupling.

The default LHAPDF `range_policy="error"` reports out-of-grid physical points.
Selecting `"extrapolate"` explicitly delegates extrapolation to LHAPDF; points
are not silently clipped or discarded.

## Integration and precision

The full-angle yield and weighted diagnostics use a five-dimensional importance-
sampled integral. Continuous angular curves use a four-dimensional integral at
each fixed angle. Independent scrambled Sobol replicas estimate numerical errors.

| Setting | Role |
|---|---|
| `POWER`, `REPLICAS`, `SEED` | $2^{\mathrm{POWER}}$ points per replica and density for integrated quantities |
| `CURVE_POWER`, `CURVE_REPLICAS`, `CURVE_SEED` | Precision and seed for the pointwise angular curves |
| `N_TAU_POINTS` | Base angular-grid size; plotting-boundary nodes can be added |
| `TAU_PLOT_RANGE` | Displayed angular interval, normally `(0.10, 0.90)` |
| `SAMPLING_PARAMETERS["proposal"]` | `"auto"`, `"gaussian"` or `"powerlaw"`; changes sampling, not the target physics |
| `OVERRIDES["radial_scale"]` | Numerical transverse-momentum sampling scale, not a physical cut |

For exploration, the notebook suggests `POWER=15`, `REPLICAS=8` before increasing
precision. The same angular scrambles are reused to retain correlations in ratios.
Reported numerical errors do not include PDF, scale, target-model or geometry
uncertainties.

The angle coordinate $\beta=\pi-\phi$ is used to normalize curves without dividing
by a vanishing endpoint Jacobian. Every auxiliary unit-area shape is normalized
over **the full range $0\leq\tau\leq1$** in the LO model. `TAU_PLOT_RANGE` changes
the displayed interval, not that normalization. The paper uses the TEEC and
$R_{pA}$ panels; unit-area curves remain available as auxiliary outputs.

## Optional checks

Both switches are in Section 1 and can be enabled independently.

### `RUN_CONVERGENCE_CHECKS=True` — Section 13

Checks Sobol precision and changes in the importance proposal. For applicable
GBW configurations it also checks the dilute reference and three representative
angular points against independent Gaussian quadrature. Results are printed and
written to `convergence.csv` and, when applicable, `independent_quadrature.csv`.
The GBW analytic references are not used as references for MV targets.

### `RUN_APPROXIMATION_CHECKS=True` — Section 14

Produces two additional percentage-difference plots:

- $100[R_{pA}^{\mathrm{exact}}/R_{pA}^{\mathrm{eik}}-1]$, keeping the main $x_t$ cut fixed.
- $100[R_{pA}(x_{t,\max}=1)/R_{pA}(x_{t,\max}=0.01)-1]$, keeping the chosen parton kinematics fixed.

The second comparison always uses **1 and 0.01**, independently of the main
`xt_max`. For `MV`, it requires the `"freeze"` large-$x_t$ continuation. The section
also exports integrated comparison tables, including toy-PDF variations when the
toy PDF is selected. Matching angular configurations are reused within the check.

The TEEC-weighted mean target momentum fraction is already computed in
**Section 11**, as `mean_xt` in `weighted_kinematic_diagnostics.csv`. No extra
switch is needed. Interpret its `density` column using $1$ for the proton and
$cA^{1/3}$ for each nuclear target.

## Plotting modes

Set a Python boolean in Section 1, without quotation marks:

```python
paper = True
```

| Mode | Format |
|---|---|
| `paper=True` | 10-point Times-style text at 3.40-inch width, no model stamp, legends inside the axes |
| `paper=False` | Original 12-point, 3.15-inch plotting format, including the model stamp |

Paper-mode legends and kinematics were checked for the RHIC charm and LHC
charm/bottom figures in both PNG and PDF. The paper TEEC panels place legends
at the upper left and kinematics at the upper right. The TEEC and ratio panels
share a boundary with no vertical gap. All plots use the selected format,
including optional comparison plots.

The paper font size matches the manuscript when a standalone PDF is inserted
at about 3.4 inches wide. Rescaling the figure in LaTeX also rescales its text.
`PAPER_FONT_SIZE` / `PAPER_WIDTH_IN` and `PLOT_FONT_SIZE` / `PLOT_WIDTH_IN` control
the two modes separately. Local TeX fonts are used when available, with a
Matplotlib font fallback.

## Outputs and reproducibility

Run folders encode the collider, quark, optional paper label, target, PDF,
pseudorapidity frame, parton kinematics and an input hash:

```text
results/RHIC_charm_paper_MV_CT18NLO_CM_exact_<hash>/
results/RHIC_charm_MV_CT18NLO_CM_exact_<hash>/
```

Rerunning identical recorded inputs writes to the same folder. Changed physical
or recorded plotting settings can produce a different hash. Check the printed
`OUTPUT` path rather than selecting a folder by its modification time alone.

| File | Contents |
|---|---|
| `teec_saturation.pdf` / `.png` | Main TEEC and $R_{pA}$ panels |
| `teec_shape_ratio.pdf` / `.png` | Auxiliary unit-area shape ratio |
| `shape_and_dilute_control.pdf` / `.png` | Auxiliary combined overview |
| `teec_density_control.pdf` / `.png` | GBW dilute control, or an unavailable-reference message |
| `teec_tau_continuum.csv` | Pointwise TEEC per nucleon, numerical errors and normalized shapes |
| `nuclear_tau_ratios.csv` | Differential nuclear ratios and shape ratios |
| `density_scan.csv` | Integrated rates and density-scaling diagnostics |
| `weighted_kinematic_diagnostics.csv` | TEEC-weighted kinematics, including `mean_xt` |
| `nuclear_summary.csv` | Integrated nuclear modifications and angular diagnostics |
| `shape_normalization.csv` | Full-range normalization integrals and displayed-window fractions |
| `continuum_normalization_check.csv` | Closure between the pointwise angular integral and independent integrated estimator |
| `replica_results.npz` | Numerical replica arrays and integration grids |
| `provenance.json` | Resolved physical inputs, provider metadata, sampling settings and source identities |
| `publication_style.json` | Font sizes, dimensions and plotting conventions |

Optional sensitivity figures use `teec_kinematics_comparison` and
`teec_xt_max_comparison` stems, with corresponding CSV, JSON and replica files.
`approximation_and_pdf_checks.csv` contains their integrated diagnostic table.

The blue shaded band spans the selected nuclear geometry choices in `C_VALUES`.
It is not the numerical error band. Preserve the provenance and data files with
any figures used in the paper.

## Extending the notebook

Implement a target provider with `qs2(xt, density)` and
`log_F(q_sq, xt, density)`, or a PDF provider with `xg(x, mu2)` and, if needed,
`alphas(mu2)`. Register it in **Section 6** before initialization using
`register_target(...)` or `register_pdf(...)`, then select its name in `CHOICES`.
New providers require their own numerical validation.

The calculation is leading order in the hard pair-production process and does
not include Sudakov resummation or a matched treatment of the near-side region.
Restricting the displayed angular interval does not itself establish the size of
omitted radiative corrections. The notebook's paper-equation references follow
`In_jet_EECs-29.pdf`. Its opening equation map identifies the corresponding
kinematics, kernel, observable and target-model equations. Recheck these labels
if the manuscript is renumbered.
