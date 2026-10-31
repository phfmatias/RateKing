---
title: "RateKing — Guia do Usuário"
subtitle: "Constantes de velocidade a partir de cálculos de química quântica · versão 0.1.0"
author: "Pedro H. F. Matias — LEEDMOL, Universidade Federal de Goiás"
date: "Outubro de 2026"
lang: pt-BR
---

# Visão geral

O RateKing é uma biblioteca Python que calcula constantes de velocidade de reações elementares a
partir dos arquivos de saída do **Gaussian 16** e do **ORCA**. Você aponta os cálculos de frequência
dos reagentes, do estado de transição (TS) e, se quiser, dos produtos; o RateKing monta as funções
de partição, calcula as barreiras e devolve:

- a constante de velocidade da Teoria do Estado de Transição convencional (TST, Eyring);
- a constante da TST deformada (d-TST), de Carvalho-Silva, Coutinho e Aquilanti;
- cinco coeficientes de transmissão de tunelamento unidimensional: Wigner, Bell 1935, Bell 1958,
  Bell de dois termos e Eckart;
- a termoquímica da reação: Ea direta e reversa, ΔE, ΔH, ΔG, ΔH‡ e ΔG‡.

A leitura dos arquivos é feita pelo MoleKing 1.6.0. Tudo o que o RateKing sabe sobre uma espécie
vem do log: geometria, frequências, número de simetria rotacional, multiplicidade e termoquímica.

Todos os números deste guia foram obtidos rodando o RateKing 0.1.0 na reação
NH~3~ + OH → NH~2~ + H~2~O, calculada em M06-2X.

# Instalação

O RateKing é instalado a partir do código-fonte:

```sh
git clone https://github.com/phfmatias/RateKing.git
cd RateKing
pip install .            # ou: make install
```

Requisitos: Python ≥ 3.10, `MoleKing==1.6.0`, `numpy` e `scipy==1.15.3` (instalados pelo `pip`).

Para desenvolver, use `make dev-install` (instalação editável). Atenção: o `pip install .` comum
**copia** o pacote para o `site-packages`. Se você editar `src/RateKing` e não rodar `make install` de
novo, o Python continua usando a versão antiga.

Teste rápido:

```python
import RateKing
print(RateKing.__all__)
```

# Preparando os cálculos de química quântica

O RateKing não corrige um cálculo ruim. Antes de usá-lo, confira:

1. **Mesmo nível de teoria** para todas as espécies (funcional, base, solvatação).
2. **Otimização + frequência** em todas as espécies. Os mínimos não podem ter frequência imaginária;
   o TS precisa ter **exatamente uma**.
3. **O TS liga os reagentes certos aos produtos certos** (confira com IRC).
4. **Multiplicidade correta**: o q~ele~ é tomado como a multiplicidade de spin.
5. **Temperatura do log**: H, S e G são lidos do log, na temperatura em que o Gaussian ou o ORCA
   rodou (em geral 298,15 K).

## Gaussian e o `link`

Um log do Gaussian pode ter vários jobs encadeados com `--Link1--`. O parâmetro `link` escolhe qual
deles o RateKing lê:

- a contagem começa em **1**;
- `-1` (padrão) é o **último** job do arquivo;
- `0` ou um número maior que o total de jobs gera `RateKingError`.

O link **precisa** ser o do cálculo de frequência. Se você apontar um job sem frequência (opt, sp,
IRC), o RateKing em geral **não dá erro**: devolve números sem sentido e o k sai absurdo. Num arquivo
opt + freq, o freq é o último e `-1` serve. Num arquivo opt → freq → sp, o freq é o link 2.

Para contar os jobs de um log:

```python
from RateKing import Molecule
print(Molecule.getAvailableLinks('ts.log'))   # número de jobs (só Gaussian)
```

## ORCA

Arquivos do ORCA são reconhecidos sozinhos (pela faixa `* O   R   C   A *` no início do arquivo). O
`link` é ignorado: cada saída tem um único job. O ORCA imprime a entropia apenas a 298,15 K.

## Caminhos

O `filePath` é relativo à pasta em que o Python foi aberto. Rode o script na pasta dos logs ou use
caminhos absolutos.

# Como o RateKing funciona

```
dicts {'name', 'filePath', 'link'}
        │
Reaction(nome, reagentes, produtos, ts, temperature=T)   # lê todos os logs na temperatura T
        │
rx.run()                                                 # barreiras, ΔE/ΔH/ΔG, razão de funções de partição
        │
TST(rx) / dTST(rx) / <correção de túnel>(rx)             # leem os atributos de rx
        │
k = kTST * kappa                                         # o produto fica por sua conta
```

Três ideias guiam o uso:

1. **A `Reaction` lê os logs no construtor**, na temperatura dada. As funções de partição ficam fixas
   ali.
2. **`run()` é obrigatório** antes de qualquer cálculo de velocidade: é ele que grava as barreiras e
   a razão de funções de partição.
3. **As classes de cálculo leem atributos da `Reaction`**. Elas não recalculam nada da molécula.

# Espécies: `Molecule` e `TransitionState`

Em geral você não cria essas classes à mão; a `Reaction` cria. Elas são úteis para inspecionar uma
espécie.

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

| Atributo | Significado | Unidade |
|---|---|---|
| `zeroPointEnergy` | **E + ZPE** (não a ZPE sozinha) | Hartree |
| `enthalpy`, `gibbs` | H e G do log | Hartree |
| `entropy` | S do log | cal mol⁻¹ K⁻¹ |
| `symmetryNumber` | número de simetria rotacional σ | — |
| `qElectronic` | multiplicidade de spin | — |
| `qVibrational`, `qRotational`, `qTranslational`, `qTotal` | funções de partição na `temperature` | — |
| `imaginaryFrequency` (só TS) | frequência imaginária | cm⁻¹, negativa |

Observações:

- Os métodos `getZPE()`, `getEnthalpy()` etc. **devolvem `None`**: eles apenas gravam o atributo e já
  rodam no construtor. Leia o atributo.
- q~vib~ é contada a partir do nível v = 0 e q~trans~ usa o volume padrão de 1 bar.
- `TransitionState` exige `name` e falha com `ExtractionError` se o arquivo não tiver uma frequência
  imaginária (ou se tiver mais de uma).

# Reação: `Reaction`

```python
Reaction(name, reactantsPath, productsPath, transitionStatePath=None, temperature=298.15)
```

Cada espécie é um **dicionário** com as três chaves obrigatórias `name`, `filePath` e `link`:

```python
from RateKing import Reaction

reagentes = [{'name': 'NH3', 'filePath': 'nh3_m062x.out', 'link': -1},
             {'name': 'OH',  'filePath': 'oh_m062x.out',  'link': -1}]
produtos  = [{'name': 'NH2', 'filePath': 'nh2_m062x.out', 'link': -1},
             {'name': 'H2O', 'filePath': 'h2o_m062x.out', 'link': -1}]
ts        =  {'name': 'TS',  'filePath': 'nh3oh_ts_m062x.out', 'link': -1}

rx = Reaction('NH3+OH', reagentes, produtos, ts, temperature=298.15)
rx.run()
```

- **Sem produtos**: passe `[]`, nunca `None` (que dá `TypeError`). Sem produtos, `deltaE`, `deltaH`,
  `deltaG` e `reverseActivationEnergy` viram a string `'N/A'`, e o Eckart não pode ser usado.
- **O TS é obrigatório**: sem ele, `RateKingError`.
- Passe dicionários, não objetos `Molecule`.

## Atributos depois de `run()`

| Atributo | Definição | Unidade | NH~3~ + OH |
|---|---|---|---|
| `activationEnergy` | (E+ZPE)~TS~ − Σ(E+ZPE)~R~ | kcal/mol | 3,30 |
| `reverseActivationEnergy` | (E+ZPE)~TS~ − Σ(E+ZPE)~P~ | kcal/mol | 16,41 |
| `activationEnergy_cal` | o mesmo, em cal/mol | cal/mol | 3297 |
| `deltaE` | Σ(E+ZPE)~P~ − Σ(E+ZPE)~R~ | kcal/mol | −13,11 |
| `deltaH`, `deltaG` | ΣH~P~ − ΣH~R~, ΣG~P~ − ΣG~R~ | kcal/mol | −12,83; −13,73 |
| `enthalpyBarrier` | H~TS~ − ΣH~R~ (ΔH‡) | **Hartree** | 0,003158 |
| `gibbsBarrier` | G~TS~ − ΣG~R~ (ΔG‡) | **Hartree** | 0,015631 |
| `reactionTotalPartitionFunction` | q~TS~ / Π q~R~ | — | — |
| `transitionStateImaginaryFrequency` | ν‡ | cm⁻¹ | −947,5 |
| `reactionSymmetryNumber` | Πσ~R~ / σ~TS~ (informativo) | — | 3 |

`enthalpyBarrier` e `gibbsBarrier` ficam em Hartree. Para kcal/mol, multiplique por 627,509474
(ou `rx.Constants.hartree2kcal`): ΔG‡ = 9,81 kcal/mol.

As barreiras usadas nas constantes de velocidade são as **com ZPE** (`activationEnergy`), coerentes
com q~vib~ medida a partir do nível v = 0.

# Constante de velocidade

## TST convencional

```python
from RateKing import TST

k = TST(rx).calculateKTST('bi')     # 1.630e-14 cm3 molecule-1 s-1
```

Para reações bimoleculares, a fórmula implementada é a de Eyring no estado padrão de 1 bar:

$$
k^{\mathrm{TST}}(T) = \frac{k_B T}{h}\,\frac{Q_{\mathrm{TS}}}{Q_A Q_B}\,
\exp\!\left(-\frac{E_a}{RT}\right)\frac{RT}{P^\circ N_A}
$$

com o resultado em cm³ molécula⁻¹ s⁻¹. Para M⁻¹ s⁻¹, multiplique por N~A~/1000.

O argumento padrão é `'bi'`. Qualquer string diferente de `'uni'` (inclusive `'Uni'`) é tratada como
`'bi'`, sem aviso.

## Reações unimoleculares

Com um único reagente, use `'uni'`:

```python
k = TST(rx).calculateKTST('uni')
```

Não use `'bi'` com um único reagente: o número não tem sentido.

## Simetria e degenerescência do caminho

O q~rot~ de cada espécie já vem dividido pelo seu σ. Assim, a razão Q~TS~/ΠQ~R~ já carrega
Πσ~R~/σ~TS~, e a degenerescência de caminho por simetria entra sozinha no k.

**Não multiplique o k por `reactionSymmetryNumber`**: isso conta σ duas vezes. Esse atributo é só
informativo. Fatores que o σ não cobre (TS quiral com enantiômero, vários confôrmeros) você aplica à
mão.

## TST deformada (d-TST)

```python
from RateKing import dTST

d = dTST(rx)
k = d.calculateKDTST('bi')        # 3.364e-14 cm3 molecule-1 s-1
d.getdParameter()
print(d.dParameter)               # -0.05637
```

A d-TST troca o fator de Boltzmann por uma exponencial deformada,

$$
k^{d\text{-TST}} = \frac{k_B T}{h}\,\frac{Q_{\mathrm{TS}}}{\prod Q_R}
\left(1 - d\,\frac{E_a}{RT}\right)^{1/d},
\qquad
d = -\frac{1}{3}\left(\frac{h\nu^\ddagger}{2E_a}\right)^{2},
$$

(no caso bimolecular, multiplicada pelo mesmo fator RT/(P°N~A~) da TST) e volta à TST quando
d → 0. Ela precisa de E~a~ > 0: com E~a~ = 0 há divisão por zero; com E~a~ < 0 o
resultado não tem sentido.

# Correções de tunelamento

Cada correção devolve um coeficiente de transmissão κ(T). **A biblioteca não multiplica k por κ**;
faça você:

```python
from RateKing import (TST, wignerCorrection, bell35Correction,
                      bell58Correction, bell2TCorrection)

kTST = TST(rx).calculateKTST('bi')
kW   = kTST * wignerCorrection(rx).calculateWignerCorrection()      # κ = 1.871
kB35 = kTST * bell35Correction(rx).calculateBell35Correction()      # κ = 3.340
kB58 = kTST * bell58Correction(rx).calculateBell58Correction()      # κ = 3.029
kB2T = kTST * bell2TCorrection(rx).calculateBell2TCorrection()      # κ = 2.696
```

| Classe | Modelo | Precisa de |
|---|---|---|
| `wignerCorrection` | 1 + u²/24 | ν‡ |
| `bell35Correction` | parábola truncada (1935) | ν‡, E~a~ > 0 |
| `bell58Correction` | (u/2) / sen(u/2) (1958) | ν‡ |
| `bell2TCorrection` | Bell com dois termos | ν‡, E~a~ > 0 |
| `eckartCorrection` | barreira de Eckart assimétrica | ν‡, E~a~, E~r~, produtos |

Aqui u = h c |ν‡| / (k~B~ T). A **temperatura de crossover** é

$$
T_c = \frac{h c\,|\nu^\ddagger|}{2\pi k_B}
$$

(217 K para NH~3~ + OH). Abaixo de T~c~ o tunelamento domina e os modelos simples divergem:

- **Wigner** só vale bem acima de T~c~; abaixo, subestima muito.
- **Bell 1958** fica negativo abaixo de T~c~: descarte o valor.
- **Bell de dois termos** corrige a divergência do Bell 1958 perto de T~c~.

## Barreira negativa

Se o TS ficar abaixo dos reagentes (E~a~ ≤ 0, comum em adição de OH), **nenhum κ vale** e a d-TST não
se aplica; a biblioteca não avisa. Nesses
casos trate o complexo pré-reativo e calcule a barreira a partir dele, o que fica fora do RateKing.

# Varrendo a temperatura

As funções de partição são calculadas **quando a `Reaction` é construída**. Mudar `rx.temperature`
depois não recalcula nada, e o k sai errado sem erro (a 500 K: 4,33·10⁻¹³ contra 1,70·10⁻¹³ do
correto). Construa uma `Reaction` por temperatura; cada uma leva cerca de 0,01 s.

Exemplo completo, com TST, d-TST, Wigner, Bell de dois termos e ajuste de Arrhenius:

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

Saída:

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

O gráfico de Arrhenius com tunelamento é curvo; o ajuste linear acima é só ilustrativo. Para uma
faixa larga de T, ajuste uma forma de Arrhenius modificada ou deformada.

A quantidade que **não** muda com a `temperature` do RateKing: H, S e G (logo ΔH, ΔG, ΔH‡ e ΔG‡) são
os do log. Para ΔG‡ em outra temperatura, refaça o cálculo de frequência com `Temperature=` no route
do Gaussian.

# Unidades

| Grandeza | Unidade |
|---|---|
| `zeroPointEnergy` (E+ZPE), `enthalpy`, `gibbs` | Hartree |
| `entropy` | cal mol⁻¹ K⁻¹ |
| `activationEnergy`, `reverseActivationEnergy`, `deltaE`, `deltaH`, `deltaG` | kcal/mol |
| `activationEnergy_cal`, `reverseActivationEnergy_cal` | cal/mol |
| `enthalpyBarrier`, `gibbsBarrier`, somas de energia | **Hartree** |
| `imaginaryFrequency`, `transitionStateImaginaryFrequency` | cm⁻¹ (negativa) |
| `temperature` | K |
| k bimolecular (`'bi'`) | cm³ molécula⁻¹ s⁻¹ (estado padrão 1 bar) |
| κ | adimensional |

As constantes físicas vêm do MoleKing (CODATA 2010) e ficam em `rx.Constants`: `k_B`, `h`, `hbar`,
`c` (em **cm/s**, para usar frequência em cm⁻¹ direto), `R` (cal mol⁻¹ K⁻¹), `R_J`, `P0` (10⁵ Pa),
`NA`, `hartree2kcal` (627,509474) e `kcal2kj` (4,184).

# Erros

As exceções ficam em `RateKing.exceptions`:

```
RateKingError
├── ExtractionError        # falha ao ler uma propriedade do log
└── CalculationError       # falha num cálculo (ex.: Eckart sem produtos)
```

```python
from RateKing.exceptions import RateKingError

try:
    rx = Reaction('NH3+OH', reagentes, produtos, ts, temperature=298.15)
except RateKingError as err:
    print(err)
```

Dentro da `Reaction`, qualquer erro de leitura chega como `RateKingError` (a mensagem traz o erro
original em `Details:`). Capture `RateKingError`, não `ExtractionError`.

| Sintoma | Causa provável |
|---|---|
| `RateKingError: Could not detect file format` | arquivo inexistente ou caminho errado |
| `RateKingError: Invalid log number` | `link` = 0 ou maior que o número de jobs |
| `RateKingError: No transition state provided` | `Reaction` sem TS |
| `ExtractionError` no TS | arquivo sem frequência imaginária, link errado ou mais de uma imaginária |
| `TypeError` ao criar a `Reaction` | `productsPath=None`; use `[]` |
| `AttributeError: ... reactionTotalPartitionFunction` | faltou `rx.run()` |
| k da ordem de 10⁻³⁰⁰ ou 10³⁰⁰, sem erro | link sem frequência |

# Checklist antes de confiar no número

1. Todas as espécies no mesmo nível de teoria; mínimos sem imaginária, TS com uma só.
2. `link` aponta o job de **frequência** de cada arquivo.
3. Uma `Reaction` por temperatura.
4. `'bi'` com dois reagentes, `'uni'` com um.
5. k corrigido = kTST × κ, sem multiplicar por `reactionSymmetryNumber`.
6. Bell 1958 descartado abaixo de T~c~; Wigner só bem acima de T~c~.
7. E~a~ > 0; se não, nenhuma correção de túnel vale.
8. `enthalpyBarrier`/`gibbsBarrier` convertidos de Hartree para kcal/mol.
9. Molécula grande: se o log do Gaussian mostra `********` nos momentos de inércia, o MoleKing põe
   q~rot~ = 1 e o k sai errado.
10. Átomos e radicais com estados de baixa energia (Cl, OH ²Π): o q~ele~ é só a multiplicidade, sem
    acoplamento spin-órbita.

# Referência rápida da API

| Chamada | Devolve |
|---|---|
| `Molecule(name, filePath, link=-1, temperature=298.15)` | espécie com termoquímica e qs |
| `Molecule.getAvailableLinks(filePath)` | número de jobs do log Gaussian |
| `TransitionState(name, filePath, link=-1, temperature=298.15)` | `Molecule` + `imaginaryFrequency` |
| `Reaction(name, reacs, prods, ts, temperature=298.15)` | reação; chamar `run()` |
| `Reaction.run()` | grava barreiras, ΔE/ΔH/ΔG e Q‡/ΠQ~R~ |
| `TST(rx).calculateKTST('bi' ou 'uni')` | k TST |
| `dTST(rx).calculateKDTST('bi' ou 'uni')` | k d-TST |
| `dTST(rx).getdParameter()` | grava `dParameter` |
| `wignerCorrection(rx).calculateWignerCorrection()` | κ Wigner |
| `bell35Correction(rx).calculateBell35Correction()` | κ Bell 1935 |
| `bell58Correction(rx).calculateBell58Correction()` | κ Bell 1958 |
| `bell58Correction(rx).getBell58factor()` | u/2 |
| `bell2TCorrection(rx).calculateBell2TCorrection()` | κ Bell de dois termos |
| `eckartCorrection(rx).calculateEckartCorrection()` | κ Eckart (precisa dos produtos) |
| `CONSTANTS()` | constantes físicas (também em `rx.Constants`) |

# Referências

- H. Eyring, *J. Chem. Phys.* **3**, 107 (1935) — Teoria do Estado de Transição.
- V. H. C. Silva, V. Aquilanti, H. C. B. de Oliveira, K. C. Mundim, *Chem. Phys. Lett.* **590**, 201
  (2013) — TST deformada.
- V. H. Carvalho-Silva, N. D. Coutinho, V. Aquilanti, *Front. Chem.* **7**, 380 (2019) — d-TST e
  transitividade.
- E. Wigner, *Z. Phys. Chem. B* **19**, 203 (1932) — correção de Wigner.
- R. P. Bell, *Proc. R. Soc. London A* **148**, 241 (1935); *Trans. Faraday Soc.* **55**, 1 (1959);
  *The Tunnel Effect in Chemistry*, Chapman & Hall (1980) — correções de Bell.
- C. Eckart, *Phys. Rev.* **35**, 1303 (1930); H. S. Johnston, J. Heicklen, *J. Phys. Chem.* **66**,
  532 (1962) — barreira de Eckart.
