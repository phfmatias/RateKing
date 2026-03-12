##   PYTHON FILE HEADER #
##
##   File:         [__init__.py]
##
##   Author(s):    ['Pedro H.F Matias']
##   Site(s):      ['https://github.com/phfmatias']
##   Email(s):     ['phfmatias@discente.ufg.br']
##   Credits:      ['Copyright © 2026 LEEDMOL. All rights reserved.']
##   Date:         ['10.02.2026']
##   Version:      ['0.1.0']
##   Status:       ['Development']
##   Language:     ['Python']
##   Description:  ['This is the top-level package for RateKing, which provides tools for calculating reaction rate constants using various methods including Transition State Theory (TST), deformed Transition State Theory (DTST), and tunneling corrections.']

import sys
from .utils.constants import CONSTANTS
from .core.molecule import Molecule, TransitionState
from .core.reaction import Reaction
from .calc.tst import TST
from .calc.dtst import dTST
from .calc import tunneling_corrections
from .calc.tunneling_corrections import eckartCorrection, bell35Correction, bell58Correction, bell2TCorrection, wignerCorrection

# Make tunneling_corrections accessible as RateKing.tunneling_corrections
sys.modules['RateKing.tunneling_corrections'] = tunneling_corrections

__all__ = [
    "CONSTANTS",
    "Molecule",
    "TransitionState",
    "Reaction",
    "TST",
    "dTST",
    "tunneling_corrections",
    "eckartCorrection",
    "bell35Correction",
    "bell58Correction",
    "bell2TCorrection",
    "wignerCorrection"
]
