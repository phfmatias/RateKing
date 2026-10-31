---
title: "RateKing — User Guide"
subtitle: "Rate constants from quantum-chemistry calculations · version 0.1.0"
author: "Pedro H. F. Matias — LEEDMOL, Universidade Federal de Goiás"
date: "October 2026"
lang: en-US
---

# Overview

RateKing is a Python library that computes rate constants of elementary reactions from **Gaussian 16**
and **ORCA** output files. You point it at the frequency calculations of the reactants, the
transition state (TS) and, optionally, the products; RateKing builds the partition functions,
computes the barriers and returns:

- the conventional Transition State Theory rate constant (TST, Eyring);
- the deformed TST rate constant (d-TST) of Carvalho-Silva, Coutinho and Aquilanti;
- five one-dimensional tunnelling transmission coefficients: Wigner, Bell 1935, Bell 1958,
  Bell two-term and Eckart;
- the reaction thermochemistry: forward and reverse Ea, ΔE, ΔH, ΔG, ΔH‡ and ΔG‡.

Files are parsed by MoleKing 1.6.0. Everything RateKing knows about a species comes from the log:
geometry, frequencies, rotational symmetry number, multiplicity and thermochemistry.

Every number in this guide was obtained by running RateKing 0.1.0 on the reaction
NH~3~ + OH → NH~2~ + H~2~O computed at the M06-2X level.

# Installation

RateKing is installed from source:

```sh
git clone https://github.com/phfmatias/RateKing.git
cd RateKing
pip install .            # or: make install
```

Requirements: Python ≥ 3.10, `MoleKing==1.6.0`, `numpy` and `scipy==1.15.3` (installed by `pip`).

For development use `make dev-install` (editable install). Note that a plain `pip install .`
**copies** the package into `site-packages`. If you edit `src/RateKing` and do not run
`make install` again, Python keeps using the old version.

Quick check:

```python
import RateKing
print(RateKing.__all__)
```

# Preparing the quantum-chemistry calculations

RateKing does not fix a bad calculation. Before using it, check:

1. **Same level of theory** for every species (functional, basis set, solvation).
2. **Optimisation + frequencies** for every species. Minima must have no imaginary frequency; the TS
   must have **exactly one**.
3. **The TS connects the right reactants to the right products** (check with an IRC).
4. **Correct multiplicity**: q~ele~ is taken as the spin multiplicity.
5. **Temperature of the log**: H, S and G are read from the log, at the temperature Gaussian or ORCA
   used (usually 298.15 K).

## Gaussian and `link`

A Gaussian log may hold several jobs chained with `--Link1--`. The `link` parameter selects which
one RateKing reads:

- counting starts at **1**;
- `-1` (default) is the **last** job in the file;
- `0` or a number larger than the number of jobs raises `RateKingError`.

The link **must** be the frequency job. If you point at a job without frequencies (opt, sp, IRC),
RateKing usually **does not raise an error**: it returns meaningless numbers and k comes out absurd.
In an opt + freq file the frequency job is the last one and `-1` works. In an opt → freq → sp file the
frequency job is link 2.

To count the jobs in a log:

```python
from RateKing import Molecule
print(Molecule.getAvailableLinks('ts.log'))   # number of jobs (Gaussian only)
```

## ORCA

ORCA files are detected automatically (by the `* O   R   C   A *` banner at the top of the file).
`link` is ignored: each output holds a single job. ORCA prints the entropy at 298.15 K only.

## Paths

`filePath` is relative to the directory Python was started in. Run the script in the folder that
holds the logs, or use absolute paths.

# How RateKing works

```
dicts {'name', 'filePath', 'link'}
        │
Reaction(name, reactants, products, ts, temperature=T)   # reads every log at temperature T
        │
rx.run()                                                 # barriers, ΔE/ΔH/ΔG, partition-function ratio
        │
TST(rx) / dTST(rx) / <tunnelling correction>(rx)         # read the attributes of rx
        │
k = kTST * kappa                                         # the product is up to you
```

Three ideas drive its use:

1. **`Reaction` reads the logs in its constructor**, at the given temperature. The partition
   functions are fixed there.
2. **`run()` is mandatory** before any rate calculation: it stores the barriers and the
   partition-function ratio.
3. **The calculation classes read attributes of the `Reaction`**. They do not recompute anything about
   the molecules.

# Species: `Molecule` and `TransitionState`

You usually do not create these by hand; `Reaction` does. They are handy to inspect a species.

```python
from RateKing import Molecule, TransitionState

nh3 = Molecule('NH3', 'nh3_m062x.out', link=-1, temperature=298.15)
print(nh3.zeroPointEnergy)      # -56.518173 Hartree (E + ZPE)
print(nh3.symmetryNumber)       # 3
print(nh3.qRotational)          # 70.9
print(nh3.qTotal)               # 1.998e8

ts = TransitionState('TS', 'nh3oh_ts_m062x.out')
print(ts.imaginaryFrequency)    # -947.5274 cm-1
```

| Attribute | Meaning | Unit |
|---|---|---|
| `zeroPointEnergy` | **E + ZPE** (not the ZPE alone) | Hartree |
| `enthalpy`, `gibbs` | H and G from the log | Hartree |
| `entropy` | S from the log | cal mol⁻¹ K⁻¹ |
| `symmetryNumber` | rotational symmetry number σ | — |
| `qElectronic` | spin multiplicity | — |
| `qVibrational`, `qRotational`, `qTranslational`, `qTotal` | partition functions at `temperature` | — |
| `imaginaryFrequency` (TS only) | imaginary frequency | cm⁻¹, negative |

Notes:

- `getZPE()`, `getEnthalpy()` and the other getters **return `None`**: they only store the attribute,
  and they already run in the constructor. Read the attribute.
- q~vib~ is counted from the v = 0 level and q~trans~ uses the 1 bar standard volume.
- `TransitionState` requires `name` and fails with `ExtractionError` if the file has no imaginary
  frequency (or more than one).

# Reaction: `Reaction`

```python
Reaction(name, reactantsPath, productsPath, transitionStatePath=None, temperature=298.15)
```

Each species is a **dictionary** with the three mandatory keys `name`, `filePath` and `link`:

```python
from RateKing import Reaction

reactants = [{'name': 'NH3', 'filePath': 'nh3_m062x.out', 'link': -1},
             {'name': 'OH',  'filePath': 'oh_m062x.out',  'link': -1}]
products  = [{'name': 'NH2', 'filePath': 'nh2_m062x.out', 'link': -1},
             {'name': 'H2O', 'filePath': 'h2o_m062x.out', 'link': -1}]
ts        =  {'name': 'TS',  'filePath': 'nh3oh_ts_m062x.out', 'link': -1}

rx = Reaction('NH3+OH', reactants, products, ts, temperature=298.15)
rx.run()
```

- **No products**: pass `[]`, never `None` (which raises `TypeError`). Without products, `deltaE`,
  `deltaH`, `deltaG` and `reverseActivationEnergy` become the string `'N/A'`, and Eckart cannot be
  used.
- **The TS is mandatory**: without it, `RateKingError`.
- Pass dictionaries, not `Molecule` objects.

## Attributes after `run()`

| Attribute | Definition | Unit | NH~3~ + OH |
|---|---|---|---|
| `activationEnergy` | (E+ZPE)~TS~ − Σ(E+ZPE)~R~ | kcal/mol | 3.30 |
| `reverseActivationEnergy` | (E+ZPE)~TS~ − Σ(E+ZPE)~P~ | kcal/mol | 16.41 |
| `activationEnergy_cal` | the same, in cal/mol | cal/mol | 3297 |
| `deltaE` | Σ(E+ZPE)~P~ − Σ(E+ZPE)~R~ | kcal/mol | −13.11 |
| `deltaH`, `deltaG` | ΣH~P~ − ΣH~R~, ΣG~P~ − ΣG~R~ | kcal/mol | −12.83; −13.73 |
| `enthalpyBarrier` | H~TS~ − ΣH~R~ (ΔH‡) | **Hartree** | 0.003158 |
| `gibbsBarrier` | G~TS~ − ΣG~R~ (ΔG‡) | **Hartree** | 0.015631 |
| `reactionTotalPartitionFunction` | q~TS~ / Π q~R~ | — | — |
| `transitionStateImaginaryFrequency` | ν‡ | cm⁻¹ | −947.5 |
| `reactionSymmetryNumber` | Πσ~R~ / σ~TS~ (informative) | — | 3 |

`enthalpyBarrier` and `gibbsBarrier` stay in Hartree. For kcal/mol, multiply by 627.509474
(or `rx.Constants.hartree2kcal`): ΔG‡ = 9.81 kcal/mol.

The barriers used in the rate constants are the **ZPE-corrected** ones (`activationEnergy`),
consistent with q~vib~ measured from the v = 0 level.

# Rate constant

## Conventional TST

```python
from RateKing import TST

k = TST(rx).calculateKTST('bi')     # 1.630e-14 cm3 molecule-1 s-1
```

For bimolecular reactions the implemented formula is Eyring's, in the 1 bar standard state:

$$
k^{\mathrm{TST}}(T) = \frac{k_B T}{h}\,\frac{Q_{\mathrm{TS}}}{Q_A Q_B}\,
\exp\!\left(-\frac{E_a}{RT}\right)\frac{RT}{P^\circ N_A}
$$

with the result in cm³ molecule⁻¹ s⁻¹. For M⁻¹ s⁻¹, multiply by N~A~/1000.

The default argument is `'bi'`. Any string other than `'uni'` (including `'Uni'`) is treated as
`'bi'`, silently.

## Unimolecular reactions

With a single reactant, use `'uni'`:

```python
k = TST(rx).calculateKTST('uni')
```

Do not use `'bi'` with a single reactant: the number is meaningless.

## Symmetry and reaction-path degeneracy

The q~rot~ of each species is already divided by its σ. The ratio Q~TS~/ΠQ~R~ therefore already
carries Πσ~R~/σ~TS~, and the symmetry-based path degeneracy enters k on its own.

**Do not multiply k by `reactionSymmetryNumber`**: that counts σ twice. The attribute is
informative only. Factors that σ does not cover (chiral TS with an enantiomer, several conformers)
must be applied by hand.

## Deformed TST (d-TST)

```python
from RateKing import dTST

d = dTST(rx)
k = d.calculateKDTST('bi')        # 3.364e-14 cm3 molecule-1 s-1
d.getdParameter()
print(d.dParameter)               # -0.05637
```

d-TST replaces the Boltzmann factor with a deformed exponential,

$$
k^{d\text{-TST}} = \frac{k_B T}{h}\,\frac{Q_{\mathrm{TS}}}{\prod Q_R}
\left(1 - d\,\frac{E_a}{RT}\right)^{1/d},
\qquad
d = -\frac{1}{3}\left(\frac{h\nu^\ddagger}{2E_a}\right)^{2},
$$

(in the bimolecular case multiplied by the same RT/(P°N~A~) factor as TST) and reduces to TST when
d → 0. It needs E~a~ > 0: with E~a~ = 0 it divides by zero; with E~a~ < 0 the result is meaningless.

# Tunnelling corrections

Each correction returns a transmission coefficient κ(T). **The library does not multiply k by κ**;
do it yourself:

```python
from RateKing import (TST, wignerCorrection, bell35Correction,
                      bell58Correction, bell2TCorrection)

kTST = TST(rx).calculateKTST('bi')
kW   = kTST * wignerCorrection(rx).calculateWignerCorrection()      # κ = 1.871
kB35 = kTST * bell35Correction(rx).calculateBell35Correction()      # κ = 3.340
kB58 = kTST * bell58Correction(rx).calculateBell58Correction()      # κ = 3.029
kB2T = kTST * bell2TCorrection(rx).calculateBell2TCorrection()      # κ = 2.696
```

| Class | Model | Needs |
|---|---|---|
| `wignerCorrection` | 1 + u²/24 | ν‡ |
| `bell35Correction` | truncated parabola (1935) | ν‡, E~a~ > 0 |
| `bell58Correction` | (u/2) / sin(u/2) (1958) | ν‡ |
| `bell2TCorrection` | Bell two-term | ν‡, E~a~ > 0 |
| `eckartCorrection` | asymmetric Eckart barrier | ν‡, E~a~, E~r~, products |

Here u = h c |ν‡| / (k~B~ T). The **crossover temperature** is

$$
T_c = \frac{h c\,|\nu^\ddagger|}{2\pi k_B}
$$

(217 K for NH~3~ + OH). Below T~c~ tunnelling dominates and the simple models break down:

- **Wigner** is only valid well above T~c~; below it, it underestimates badly.
- **Bell 1958** turns negative below T~c~: discard the value.
- **Bell two-term** removes the Bell 1958 divergence near T~c~.

## Negative barrier

If the TS lies below the reactants (E~a~ ≤ 0, common in OH additions), **no κ is valid** and d-TST does
not apply; the library gives no warning.
In those cases treat the pre-reactive complex and compute the barrier from it, which is outside
RateKing.

# Temperature scans

Partition functions are computed **when the `Reaction` is built**. Changing `rx.temperature`
afterwards recomputes nothing, and k comes out wrong without an error (at 500 K: 4.33×10⁻¹³ against the
correct 1.70×10⁻¹³). Build one `Reaction` per temperature; each takes about 0.01 s.

Full example, with TST, d-TST, Wigner, Bell two-term and an Arrhenius fit:

```python
import numpy as np
from RateKing import Reaction, TST, dTST, wignerCorrection, bell2TCorrection


reactants = [{'name': 'NH3', 'filePath': 'nh3_m062x.out', 'link': -1},
             {'name': 'OH',  'filePath': 'oh_m062x.out',  'link': -1}]
products  = [{'name': 'NH2', 'filePath': 'nh2_m062x.out', 'link': -1},
             {'name': 'H2O', 'filePath': 'h2o_m062x.out', 'link': -1}]
ts        =  {'name': 'TS',  'filePath': 'nh3oh_ts_m062x.out', 'link': -1}

temperatures = [250, 298.15, 400, 600, 1000]
rows = []

for T in temperatures:
    rx = Reaction('NH3+OH', reactants, products, ts, temperature=T)   # one Reaction per T
    rx.run()
    kTST = TST(rx).calculateKTST('bi')
    kdTST = dTST(rx).calculateKDTST('bi')
    kW = kTST * wignerCorrection(rx).calculateWignerCorrection()
    kB2T = kTST * bell2TCorrection(rx).calculateBell2TCorrection()
    rows.append((T, kTST, kdTST, kW, kB2T))

print('{:>8s} {:>11s} {:>11s} {:>11s} {:>11s}'.format('T (K)', 'kTST', 'kdTST', 'kTST/W', 'kTST/B2T'))
for T, kTST, kdTST, kW, kB2T in rows:
    print('{:8.2f} {:11.3e} {:11.3e} {:11.3e} {:11.3e}'.format(T, kTST, kdTST, kW, kB2T))

# Arrhenius fit: ln k = ln A - Ea/(R T)
T = np.array([r[0] for r in rows])
k = np.array([r[4] for r in rows])
slope, intercept = np.polyfit(1 / T, np.log(k), 1)
print('A  = {:.3e} cm3 molecule-1 s-1'.format(np.exp(intercept)))
print('Ea = {:.2f} kcal/mol'.format(-slope * 1.9872036e-3))
```

Output:

```
   T (K)        kTST       kdTST      kTST/W    kTST/B2T
  250.00   5.940e-15   1.613e-14   1.330e-14   2.596e-14
  298.15   1.630e-14   3.364e-14   3.050e-14   4.394e-14
  400.00   6.744e-14   1.027e-13   1.001e-13   1.135e-13
  600.00   3.390e-13   4.122e-13   4.119e-13   4.232e-13
 1000.00   1.954e-12   2.102e-12   2.105e-12   2.112e-12
A  = 5.965e-12 cm3 molecule-1 s-1
Ea = 2.86 kcal/mol
```

An Arrhenius plot with tunnelling is curved; the linear fit above is only illustrative. Over a wide
T range, fit a modified or deformed Arrhenius form.

What does **not** change with RateKing's `temperature`: H, S and G (hence ΔH, ΔG, ΔH‡ and ΔG‡) are
those of the log. For ΔG‡ at another temperature, rerun the frequency calculation with
`Temperature=` in the Gaussian route.

# Units

| Quantity | Unit |
|---|---|
| `zeroPointEnergy` (E+ZPE), `enthalpy`, `gibbs` | Hartree |
| `entropy` | cal mol⁻¹ K⁻¹ |
| `activationEnergy`, `reverseActivationEnergy`, `deltaE`, `deltaH`, `deltaG` | kcal/mol |
| `activationEnergy_cal`, `reverseActivationEnergy_cal` | cal/mol |
| `enthalpyBarrier`, `gibbsBarrier`, energy sums | **Hartree** |
| `imaginaryFrequency`, `transitionStateImaginaryFrequency` | cm⁻¹ (negative) |
| `temperature` | K |
| bimolecular k (`'bi'`) | cm³ molecule⁻¹ s⁻¹ (1 bar standard state) |
| κ | dimensionless |

Physical constants come from MoleKing (CODATA 2010) and live in `rx.Constants`: `k_B`, `h`, `hbar`,
`c` (in **cm/s**, so frequencies in cm⁻¹ go in directly), `R` (cal mol⁻¹ K⁻¹), `R_J`, `P0` (10⁵ Pa),
`NA`, `hartree2kcal` (627.509474) and `kcal2kj` (4.184).

# Errors

Exceptions live in `RateKing.exceptions`:

```
RateKingError
├── ExtractionError        # failed to read a property from the log
└── CalculationError       # failed calculation (e.g. Eckart without products)
```

```python
from RateKing.exceptions import RateKingError

try:
    rx = Reaction('NH3+OH', reactants, products, ts, temperature=298.15)
except RateKingError as err:
    print(err)
```

Inside `Reaction`, any parsing error arrives as `RateKingError` (the message carries the original
error under `Details:`). Catch `RateKingError`, not `ExtractionError`.

| Symptom | Likely cause |
|---|---|
| `RateKingError: Could not detect file format` | missing file or wrong path |
| `RateKingError: Invalid log number` | `link` = 0 or larger than the number of jobs |
| `RateKingError: No transition state provided` | `Reaction` without a TS |
| `ExtractionError` on the TS | no imaginary frequency, wrong link, or more than one imaginary |
| `TypeError` when building `Reaction` | `productsPath=None`; use `[]` |
| `AttributeError: ... reactionTotalPartitionFunction` | missing `rx.run()` |
| k around 10⁻³⁰⁰ or 10³⁰⁰, no error | link without frequencies |

# Checklist before trusting the number

1. All species at the same level of theory; minima with no imaginary frequency, TS with exactly one.
2. `link` points at the **frequency** job of every file.
3. One `Reaction` per temperature.
4. `'bi'` with two reactants, `'uni'` with one.
5. Corrected k = kTST × κ, without multiplying by `reactionSymmetryNumber`.
6. Bell 1958 discarded below T~c~; Wigner only well above T~c~.
7. E~a~ > 0; otherwise no tunnelling correction is valid.
8. `enthalpyBarrier`/`gibbsBarrier` converted from Hartree to kcal/mol.
9. Large molecules: if the Gaussian log shows `********` in the moments of inertia, MoleKing sets
   q~rot~ = 1 and k is wrong.
10. Atoms and radicals with low-lying states (Cl, OH ²Π): q~ele~ is only the multiplicity, with no
    spin–orbit coupling.

# API quick reference

| Call | Returns |
|---|---|
| `Molecule(name, filePath, link=-1, temperature=298.15)` | species with thermochemistry and qs |
| `Molecule.getAvailableLinks(filePath)` | number of jobs in a Gaussian log |
| `TransitionState(name, filePath, link=-1, temperature=298.15)` | `Molecule` + `imaginaryFrequency` |
| `Reaction(name, reacs, prods, ts, temperature=298.15)` | reaction; call `run()` |
| `Reaction.run()` | stores barriers, ΔE/ΔH/ΔG and Q‡/ΠQ~R~ |
| `TST(rx).calculateKTST('bi' or 'uni')` | TST k |
| `dTST(rx).calculateKDTST('bi' or 'uni')` | d-TST k |
| `dTST(rx).getdParameter()` | stores `dParameter` |
| `wignerCorrection(rx).calculateWignerCorrection()` | Wigner κ |
| `bell35Correction(rx).calculateBell35Correction()` | Bell 1935 κ |
| `bell58Correction(rx).calculateBell58Correction()` | Bell 1958 κ |
| `bell58Correction(rx).getBell58factor()` | u/2 |
| `bell2TCorrection(rx).calculateBell2TCorrection()` | Bell two-term κ |
| `eckartCorrection(rx).calculateEckartCorrection()` | Eckart κ (needs products) |
| `CONSTANTS()` | physical constants (also in `rx.Constants`) |

# References

- H. Eyring, *J. Chem. Phys.* **3**, 107 (1935) — Transition State Theory.
- V. H. C. Silva, V. Aquilanti, H. C. B. de Oliveira, K. C. Mundim, *Chem. Phys. Lett.* **590**, 201
  (2013) — deformed TST.
- V. H. Carvalho-Silva, N. D. Coutinho, V. Aquilanti, *Front. Chem.* **7**, 380 (2019) — d-TST and
  transitivity.
- E. Wigner, *Z. Phys. Chem. B* **19**, 203 (1932) — Wigner correction.
- R. P. Bell, *Proc. R. Soc. London A* **148**, 241 (1935); *Trans. Faraday Soc.* **55**, 1 (1959);
  *The Tunnel Effect in Chemistry*, Chapman & Hall (1980) — Bell corrections.
- C. Eckart, *Phys. Rev.* **35**, 1303 (1930); H. S. Johnston, J. Heicklen, *J. Phys. Chem.* **66**,
  532 (1962) — Eckart barrier.
