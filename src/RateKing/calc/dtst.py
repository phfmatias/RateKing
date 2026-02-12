##   PYTHON FILE HEADER #
##
##   File:         [dtst.py]
##
##   Author(s):    ['Pedro H.F Matias']
##   Site(s):      ['https://github.com/phfmatias']
##   Email(s):     ['phfmatias@discente.ufg.br']
##   Credits:      ['Copyright © 2026 LEEDMOL. All rights reserved.']
##   Date:         ['10.02.2026']
##   Version:      ['0.1.0']
##   Status:       ['Development']
##   Language:     ['Python']
##   Description:  ['This module implements the deformed Transition State Theory (DTST) for rate constant calculations in RateKing.']

from ..utils.constants import CONSTANTS
import numpy as np

class dTST:
    """Class to calculate the rate constant using deformed Transition State Theory (DTST) for a given reaction."""

    def __init__(self, reaction):
        self.reaction = reaction
        self.Constants = CONSTANTS()

    def getdParameter(self):
        """Method to calculate the deformation parameter (d) for DTST."""

        self.dParameter = (-0.004769426830E0) * ((self.Constants.c * self.Constants.h * self.Constants.NA * self.reaction.transitionStateImaginaryFrequency) / (self.reaction.activationEnergy * 1000.0)) ** 2
        self.inverse_dParameter = 1 / self.dParameter

    def calculateKDTST(self, bi_or_uni='bi'):
        """Method to calculate the rate constant using DTST."""

        self.getdParameter()

        f1 = ((self.Constants.k_B*self.reaction.temperature)/self.Constants.h)
        f2 = self.reaction.reactionTotalPartitionFunction
        f3 = np.exp(self.inverse_dParameter * np.log(1 - self.dParameter * self.reaction.activationEnergy_cal / (self.Constants.R * self.reaction.temperature)))
        vol_correction = (self.Constants.R_J*self.reaction.temperature)/self.Constants.P0
        kdtst = f1 * f2 * f3 * 1e6 * vol_correction
        if bi_or_uni == 'uni':
                return kdtst / 1e6
        else:   
            return kdtst / self.Constants.NA    

       

# self.kdtst = self._kdtst(self.QTot, self.di, self.d, self.eacal, self.bi_or_uni)