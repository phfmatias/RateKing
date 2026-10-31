# RateKing

**Rate constants from quantum-chemistry output files.** RateKing reads Gaussian and ORCA frequency
calculations for the reactants, transition state and products of a reaction, builds the partition
functions, and returns rate constants from conventional Transition State Theory (TST), deformed TST
(d-TST) and five one-dimensional tunnelling corrections.

Developed at LEEDMOL (Universidade Federal de Goiás) by Pedro H. F. Matias.

- Version: 0.1.0 (development)
- License: MIT
- Python: ≥ 3.10
- Parsing: MoleKing 1.6.0

## Features

| Area | What RateKing gives you |
|---|---|
| Input | Gaussian 16 `.log`/`.out` (any `Link1` job) and ORCA `.out`, detected automatically |
| Thermochemistry | E+ZPE, H, S, G per species; q_trans, q_rot, q_vib, q_ele, q_total at any temperature |
| Energetics | Forward and reverse barriers (ZPE-corrected), ΔE, ΔH, ΔG of reaction, ΔH‡, ΔG‡ |
| Rate theory | Conventional TST (Eyring), deformed TST (Carvalho-Silva, Coutinho & Aquilanti) |
| Tunnelling | Wigner, Bell 1935, Bell 1958, Bell two-term, Eckart |
| Symmetry | Rotational symmetry numbers taken from the log; reaction-path degeneracy enters through q_rot |

## Installation

RateKing is installed from source:

```sh
git clone https://github.com/phfmatias/RateKing.git
cd RateKing
pip install .          # or: make install
```

For development use `make dev-install` (editable install). The dependencies (`MoleKing==1.6.0`,
`numpy`, `scipy==1.15.3`) are pulled in by `pip`.

> A regular `pip install .` copies the package into `site-packages`. After editing `src/`, run
> `make install` again or the changes will not be seen.

## Quick start

```python
from RateKing import Reaction, TST, wignerCorrection

reactants = [{'name': 'NH3', 'filePath': 'nh3_m062x.out', 'link': -1},
             {'name': 'OH',  'filePath': 'oh_m062x.out',  'link': -1}]
products  = [{'name': 'NH2', 'filePath': 'nh2_m062x.out', 'link': -1},
             {'name': 'H2O', 'filePath': 'h2o_m062x.out', 'link': -1}]
ts        =  {'name': 'TS',  'filePath': 'nh3oh_ts_m062x.out', 'link': -1}

rx = Reaction('NH3+OH', reactants, products, ts, temperature=298.15)
rx.run()                                            # required before any rate calculation

kTST  = TST(rx).calculateKTST('bi')                 # cm3 molecule-1 s-1
kappa = wignerCorrection(rx).calculateWignerCorrection()

print(f'Ea      = {rx.activationEnergy:.2f} kcal/mol')
print(f'kTST    = {kTST:.3e} cm3 molecule-1 s-1')
print(f'k(TST/W)= {kTST * kappa:.3e} cm3 molecule-1 s-1')
```

Output (M06-2X logs):

```
Ea      = 3.30 kcal/mol
kTST    = 1.630e-14 cm3 molecule-1 s-1
k(TST/W)= 3.050e-14 cm3 molecule-1 s-1
```

## Package layout

```
src/RateKing/
├── core/
│   ├── molecule.py              Molecule, TransitionState
│   └── reaction.py              Reaction
├── calc/
│   ├── tst.py                   TST
│   ├── dtst.py                  dTST
│   └── tunneling_corrections.py wignerCorrection, bell35Correction, bell58Correction,
│                                bell2TCorrection, eckartCorrection
├── utils/constants.py           CONSTANTS (CODATA 2010, from MoleKing)
└── exceptions.py                RateKingError, ExtractionError, CalculationError
```

## Workflow

```
dicts {'name', 'filePath', 'link'}
        │
Reaction(name, reactants, products, ts, temperature=T)   # reads every log at T
        │
rx.run()                                                 # barriers, ΔE/ΔH/ΔG, partition-function ratio
        │
TST(rx) / dTST(rx) / <tunnelling>(rx)                    # read the attributes of rx
        │
k = kTST * kappa                                         # the product is up to you
```

Rules that matter:

1. `Reaction` takes **dicts**, not `Molecule` objects. With no products pass `[]`, never `None`.
2. `link` is 1-based and `-1` means the last `Link1` job. It must point at the **frequency** job.
3. Partition functions are computed when the `Reaction` is built. For a temperature scan, build one
   `Reaction` per temperature; changing `rx.temperature` afterwards does not update them.
4. The library does not multiply k by κ, and it must not be multiplied by `reactionSymmetryNumber`
   (σ is already inside q_rot).

## Units

| Quantity | Unit |
|---|---|
| `activationEnergy`, `reverseActivationEnergy`, `deltaE`, `deltaH`, `deltaG` | kcal/mol |
| `activationEnergy_cal`, `reverseActivationEnergy_cal` | cal/mol |
| `enthalpyBarrier`, `gibbsBarrier`, species energies | **Hartree** |
| `transitionStateImaginaryFrequency` | cm⁻¹ (negative) |
| Bimolecular k (`'bi'`) | cm³ molecule⁻¹ s⁻¹, standard state 1 bar |

## Documentation

- User guide (Portuguese): [`docs/guia_usuario_pt-BR.md`](docs/guia_usuario_pt-BR.md)
- User guide (English): [`docs/user_guide_en-US.md`](docs/user_guide_en-US.md)

## Citing

If RateKing is useful in your work, please cite the repository
(<https://github.com/phfmatias/RateKing>) and the original papers of the methods you use
(Eyring TST, d-TST, Wigner, Bell, Eckart).

## License

MIT. See [`LICENSE`](LICENSE).
