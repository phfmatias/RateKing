##   PYTHON FILE HEADER #
##
##   File:         [tunneling_corrections.py]
##
##   Author(s):    ['Pedro H.F Matias']
##   Site(s):      ['https://github.com/phfmatias']
##   Email(s):     ['phfmatias@discente.ufg.br']
##   Credits:      ['Copyright © 2026 LEEDMOL. All rights reserved.']
##   Date:         ['10.02.2026']
##   Version:      ['0.1.0']
##   Status:       ['Development']
##   Language:     ['Python']
##   Description:  ['This module implements tunneling corrections for rate constant calculations in RateKing.']

from ..utils.constants import CONSTANTS
from ..exceptions import CalculationError
from scipy.integrate import quad
import numpy as np

class eckartCorrection:
    """Class to calculate the Eckart tunneling correction for a given reaction."""

    def __init__(self, reaction):
        self.reaction = reaction
        self.Constants = CONSTANTS()

        if len(self.reaction.products) == 0:
            raise CalculationError("Error! Eckart correction cannot be calculated for reactions without products.")

    def alphaEckart(self, E):
        """Method to calculate the alpha parameter for the Eckart correction."""
        return (self.Constants.kcal2kj / (self.Constants.c * self.Constants.NA)) * np.sqrt(2.0 * self.barrierWidth**2 * E / self.Constants.h**2)
    
    def betaEckart(self, E):
        """Method to calculate the beta parameter for the Eckart correction."""
        return (self.Constants.kcal2kj / (self.Constants.c * self.Constants.NA)) * np.sqrt(2.0 * self.barrierWidth**2 * (E - (self.reaction.activationEnergy_cal - self.reaction.reverseActivationEnergy_cal)) / self.Constants.h**2)
    
    def deltaEckart(self, E):
        """Method to calculate the delta parameter for the Eckart correction."""
        return (self.Constants.kcal2kj / (self.Constants.c * self.Constants.NA)) * np.sqrt(4.0 * self.reaction.activationEnergy_cal * self.reaction.reverseActivationEnergy_cal / (self.Constants.h * abs(self.reaction.transitionStateImaginaryFrequency))**2 - 0.25)
    
    def P_Eckart(self, E):
        a = self.alphaEckart(E)
        b = self.betaEckart(E)
        dlt = self.deltaEckart(E)

        return 1 - (np.cosh(2.0 * np.pi * (a + b)) - np.cosh(2.0 * np.pi * (a - b))) / (
            np.cosh(2.0 * np.pi * (a + b)) + np.cosh(2.0 * np.pi * dlt))
    
    def integradumEckart(self, E):
        return self.P_Eckart(E) * np.exp(-(E - self.reaction.activationEnergy_cal) / (self.Constants.R * self.reaction.temperature))
    
    def getLowerIntegrationLimit(self):
        self.lowerIntegrationLimit = max([0, self.reaction.activationEnergy_cal - self.reaction.reverseActivationEnergy_cal])

    def getBarrierWidth(self):
        self.barrierWidth = (self.reaction.activationEnergy_cal ** (-0.5) + self.reaction.reverseActivationEnergy_cal ** (-0.5)) ** (-1) * np.sqrt(2.0) / abs(self.reaction.transitionStateImaginaryFrequency)

    def calculateEckartCorrection(self):
        self.getLowerIntegrationLimit()
        self.getBarrierWidth()

        try:
            self.integralEckart, self.errorEckart = quad(self.integradumEckart, self.lowerIntegrationLimit, np.inf)
        except Exception as e:
            raise CalculationError(f"Error! Eckart correction calculation failed. Details: {str(e)}")
        
        self.kappaEckart = self.integralEckart / (self.Constants.R * self.reaction.temperature)

        return self.kappaEckart
    
class bell35Correction:
    """Class to calculate the Bell-35 tunneling correction for a given reaction."""

    def __init__(self, reaction):
        self.reaction = reaction
        self.Constants = CONSTANTS()

    def calculateBell35Correction(self):
        f2 = (np.exp(-self.reaction.activationEnergy_cal / (self.Constants.R * self.reaction.temperature)))

        self.kappaBell35 = (self.Constants.k_B * self.reaction.temperature - 1.0E0 * self.Constants.hbar * abs(self.reaction.transitionStateImaginaryFrequency) * self.Constants.c * np.exp(-(2.0E0 * self.Constants.kcal2kj) * np.pi * self.reaction.activationEnergy_cal / (self.Constants.c * self.Constants.h * abs(self.reaction.transitionStateImaginaryFrequency) * self.Constants.NA)) / f2) / (self.Constants.k_B * self.reaction.temperature - self.Constants.hbar * abs(self.reaction.transitionStateImaginaryFrequency) * self.Constants.c)

        return self.kappaBell35

class bell58Correction:
    """Class to calculate the Bell-58 tunneling correction for a given reaction."""

    def __init__(self, reaction):
        self.reaction = reaction
        self.Constants = CONSTANTS()

    def getBell58factor(self):
        self.fBell58 = (0.50E0 * self.Constants.h * abs(self.reaction.transitionStateImaginaryFrequency) * self.Constants.c) / (self.Constants.k_B * self.reaction.temperature)

        self.sin_fBell58 = np.sin(self.fBell58)

        if abs(np.sin(self.fBell58)) < 1e-10:
            raise CalculationError("Division by zero in Bell58 calculation. Check temperature and imaginary frequency.")
        
        return self.fBell58
    
    def calculateBell58Correction(self):
        self.getBell58factor()

        self.kappaBell58 = self.fBell58 / self.sin_fBell58

        return self.kappaBell58

class bell2TCorrection:
    """Class to calculate the Bell-2TC tunneling correction for a given reaction."""

    def __init__(self, reaction):
        self.reaction = reaction
        self.Constants = CONSTANTS()

        self.bell58Factor = bell58Correction(self.reaction).getBell58factor()
        self.kappaBell58 = bell58Correction(self.reaction).calculateBell58Correction()

    def calculateBell2TCorrection(self):
        f2 = (np.exp(-self.reaction.activationEnergy_cal / (self.Constants.R * self.reaction.temperature)))
        
        self.kappaBell2T = self.kappaBell58 - 1.0E0 * np.exp(-(2.0E0 * self.Constants.kcal2kj) * np.pi * self.reaction.activationEnergy_cal / (self.Constants.c * self.Constants.h * abs(self.reaction.transitionStateImaginaryFrequency) * self.Constants.NA)) / (f2 * (1.0E0 * np.pi / self.bell58Factor - 1.0E0))

        return self.kappaBell2T
    
class wignerCorrection:
    """Class to calculate the Wigner tunneling correction for a given reaction."""

    def __init__(self, reaction):
        self.reaction = reaction
        self.Constants = CONSTANTS()

    def calculateWignerCorrection(self):
        self.kappaWigner = 1 + (1/24) * ((self.Constants.h * abs(self.reaction.transitionStateImaginaryFrequency) * self.Constants.c) / (self.Constants.k_B * self.reaction.temperature)) ** 2

        return self.kappaWigner