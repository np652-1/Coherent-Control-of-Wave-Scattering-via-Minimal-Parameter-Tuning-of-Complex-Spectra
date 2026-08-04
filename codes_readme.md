# D-shaped reverberation cavity — scattering-matrix solvers, cross-checks, and CCON optimization

Example codes for a 2D scalar-Helmholtz **D-shaped reverberation cavity** with
four single-mode waveguide leads and internal ellipse scatterers. Two
independent solvers are provided — a boundary-integral-equation (BIE) solver
built on [chunkie](https://github.com/fastalgorithms/chunkie) and a
finite-element solver in **COMSOL Multiphysics** (via LiveLink for MATLAB) —
together with scripts that cross-check them and a small **coherent-control
(CCON) optimization** demo.

The cavity models a lossy/lossless reverberation chamber: a circular **arc**
(radius `R`, spanning 60°→300°) closed by a flat **chord wall** at `x = R/2`
that carries the four leads, plus `N` PEC **ellipse** scatterers inside. Walls
are either perfectly conducting (lossless) or carry a finite-conductivity
**Leontovich impedance** boundary (lossy).

---

## The codes

| # | file | what it does |
|---|------|--------------|
| 1 | `dcavity_chunkie_Smatrix_lossless.m` | 4×4 S-matrix of the **lossless** (PEC) cavity, chunkie BIE |
| 2 | `dcavity_comsol_Smatrix_lossless.m`  | 4×4 S-matrix of the **lossless** cavity, COMSOL FEM |
| 3 | `dcavity_compare_lossless.m`         | **lossless** cross-check: both solvers as subroutines, compares S **and** interior fields |
| 4 | `dcavity_chunkie_Smatrix_lossy.m`    | 4×4 S-matrix of the **lossy** (impedance) cavity, chunkie BIE |
| 5 | `dcavity_comsol_Smatrix_lossy.m`     | 4×4 S-matrix of the **lossy** cavity, COMSOL FEM |
| 6 | `dcavity_compare_lossy.m`            | **lossy** cross-check: both solvers as subroutines, compares S **and** interior fields |
| 7 | `dcavity_ccon_optimization/`         | CCON figure-of-merit **optimization** demo for 0% / 50% / 90% wall loss (~10 min) |

Each script is self-contained (helper functions are local to the file) and
writes its outputs into a sibling `*_outputs/` folder.

---

## Physics & conventions (shared by every code)

- **PDE / polarization.** 2D scalar Helmholtz, **TM** polarization: the scalar
  field is `H_z`. Time convention `exp(-iωt)` ⇒ outgoing waves are `H_0^{(1)}`.
- **Boundary conditions.**
  - PEC wall ⇒ **Neumann** `∂H_z/∂n = 0`.
  - Leontovich impedance wall ⇒ **Robin** `∂H_z/∂n + γ H_z = 0`, with
    `γ = i ω_phys ε₀ Z_s`, `Z_s = sqrt(i ω_phys μ₀ / (σ + i ω_phys ε₀))`.
  - The **arc + chord** carry the wall loss; the **lead walls + ellipses** are
    always PEC. Lossless = all PEC (`σ = ∞ ⇒ γ = 0`).
- **Ports.** Each lead supports one propagating **TEM** mode (`β = ω`). Driving
  each lead and projecting the field and its normal derivative onto the TEM mode
  on a cut-line gives incoming/outgoing amplitudes `a_in`, `a_out`, and
  `S = a_out / a_in`. Lossless ⇒ `S` unitary; lossy ⇒ sub-unitary (walls absorb).
- **BIE representation.** A plain **single-layer** density, mixed Robin/Neumann
  operator scaled to `(I + K)` form and assembled per-edge so chunkie's RCIP
  handles the wall/lead corners. This matches the production study and is
  validated here at a fixed, generic (non-resonant) frequency by a
  manufactured-solution gate and by the COMSOL cross-check. *(A plain
  single-layer carries spurious interior resonances; for complex-frequency
  root/pole finding switch to a combined-field representation.)*
- **COMSOL side.** `ElectromagneticWavesFrequencyDomain` (`ewfd`) with an
  in-plane E-field (⇒ scalar `H_z`), rectangular TEM ports, and an Impedance
  boundary (lossy) or the default PEC (lossless).

### S-matrix comparison caveat
The two solvers use a different per-port **phase reference** and the opposite
time convention (`e^{-iωt}` vs `e^{+jωt}`), under which the **magnitude** `|S|`
and the per-channel **absorption** `1 − Σ|S(:,m)|²` are invariant. The
comparison scripts therefore compare `|S|` and absorption (not the raw complex
`S`). Interior fields agree up to one complex **gauge scalar per channel**,
fit by complex least squares.

---

## How to run

### Getting chunkie (needed by codes 1, 3, 4, 6, 7)
The [chunkie](https://github.com/fastalgorithms/chunkie) BIE library is included
as a **git submodule** (pinned to the tested commit), and chunkie in turn uses
the pure-MATLAB **FLAM** library as a nested submodule — so initialise submodules
**recursively**:
```bash
git clone --recurse-submodules https://github.com/np652-1/Coherent-Control-of-Wave-Scattering-via-Minimal-Parameter-Tuning-of-Complex-Spectra.git
# or, in an already-cloned repo:
git submodule update --init --recursive
```
The scripts then auto-find the bundled `chunkie/` submodule. **No compilation is needed** —
these examples use *dense* chunkie, which is pure MATLAB (chunkie + FLAM). Only
chunkie's `fmm2d` sub-submodule would need compiling, and these dense examples do
not use it, so its absence is fine (chunkie's `startup` prints a harmless
"no fmm2d mex" warning). Code 7 additionally needs the MATLAB **Optimization
Toolbox** (`fmincon`).

### chunkie codes (1, 4, 7)
```matlab
dcavity_chunkie_Smatrix_lossless      % code 1
dcavity_chunkie_Smatrix_lossy         % code 4
cd dcavity_ccon_optimization; optimize_ccon_example   % code 7
```

### COMSOL codes (2, 3, 5, 6)
Need COMSOL with the **Wave Optics Module** (the `ewfd` "Electromagnetic Waves, Frequency Domain" interface) and **LiveLink for MATLAB**. Start a server in a terminal, then
run the script:

```bash
comsol mphserver -port 2036          # (add `-multi on` to keep it up for several runs)
```
```matlab
% set COMSOL_MLI (the COMSOL "mli" directory) at the top of each script, then:
dcavity_comsol_Smatrix_lossless       % code 2
dcavity_comsol_Smatrix_lossy          % code 5
dcavity_compare_lossless              % code 3  (also needs chunkie)
dcavity_compare_lossy                 % code 6  (also needs chunkie)
```
(You can equivalently launch `comsol mphserver matlab`, which starts a MATLAB
already linked to a server; the scripts detect an existing link.)

---

## Validation (verified end-to-end on a test machine)

- **Lossless `|S|` (codes 1 vs 2, cross-check in code 3):** chunkie and COMSOL
  agree to **max 1.6e-3 / Frobenius 2.1e-3** (COMSOL mesh-limited). Both `S`
  are unitary (`‖SᴴS − I‖`: chunkie 2.7e-9, COMSOL 2.5e-13) and reciprocal
  (`‖S − Sᵀ‖ ~ 2e-10`).
- **Lossy `|S|` (σ=50; codes 4 vs 5, cross-check in code 6):** agree to
  **max 1.7e-3 / Frobenius 2.7e-3**; per-channel absorption ≈ [0.67 0.59 0.68
  0.33] in both, agreeing to **7e-4**.
- **Interior fields (codes 3, 6):** the gauge-matched `H_z` fields agree
  pointwise to **~3e-3** across all four channels, both loss levels (the
  log₁₀-error panels are a flat ~−3 with no structure — see `example_outputs/`).
- **Manufactured-solution gate (codes 1/4):** rel. err ~1e-10.
- **CCON optimization (code 7):** for the headline word **NNDD** (2 ellipses),
  the FOM `det(CᴴC)` is driven from ~1e-1 to the numerical floor — **2.8e-17**
  (lossless), **2.2e-14** (50%), **1.9e-6** (90%) — a near-perfect
  coherent-control null except where heavy wall loss caps it. Total run ~7 min.

---

## The CCON optimization (code 7) in one paragraph

A **CCON word** assigns each of the 4 ports a role — **R** (reflectionless), **N**
(normal), **D** (dark), **T** (transmission) — which selects a sub-block
`C = S[rows,cols]`. The figure of merit `F(θ) = det(CᴴC)` is minimized over the
ellipse **rotation angles** `θ` (multi-start `fmincon` with an analytic adjoint
gradient + `fminsearch` polish). `optimize_ccon_example.m` runs one small fixed
realization per loss level; the full engine `run_ccon_study.m` also provides
`'verify'` (adjoint-vs-finite-difference gradient gate), `'smoke'` (one quick
optimization per CCON word), and the production `'run'` mode.
