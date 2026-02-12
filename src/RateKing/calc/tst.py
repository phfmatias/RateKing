##   PYTHON FILE HEADER #
##
##   File:         [tst.py]
##
##   Author(s):    ['Pedro H.F Matias']
##   Site(s):      ['https://github.com/phfmatias']
##   Email(s):     ['phfmatias@discente.ufg.br']
##   Credits:      ['Copyright © 2026 LEEDMOL. All rights reserved.']
##   Date:         ['10.02.2026']
##   Version:      ['0.1.0']
##   Status:       ['Development']
##   Language:     ['Python']
##   Description:  ['This module calculates the rate constant using the Transition State Theory (TST) for a given reaction.']

import numpy as np

from ..utils.constants import CONSTANTS

class TST:
    """Class to calculate the rate constant using Transition State Theory (TST) for a given reaction."""

    def __init__(self, reaction):
        self.reaction = reaction
        self.Constants = CONSTANTS()

    def calculateKTST(self, bi_or_uni='bi'):
        """Method to calculate the rate constant using TST."""

        f1 = ((self.Constants.k_B * self.reaction.temperature) / self.Constants.h)
        f2 = self.reaction.reactionTotalPartitionFunction
        f3 =  np.exp(-self.reaction.activationEnergy_cal/(self.Constants.R * self.reaction.temperature))        
        vol_correction = (self.Constants.R_J * self.reaction.temperature) / self.Constants.P0
        ktst = f1 * f2 * f3 * 1e6 * vol_correction
        
        if bi_or_uni == 'uni':
            return ktst / 1e6
        else:
            return ktst / self.reaction.Constants.NA


