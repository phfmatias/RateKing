##   PYTHON FILE HEADER #
##
##   File:         [constants.py]
##
##   Author(s):    ['Pedro H.F Matias']
##   Site(s):      ['https://github.com/phfmatias']
##   Email(s):     ['phfmatias@discente.ufg.br']
##   Credits:      ['Copyright © 2026 LEEDMOL. All rights reserved.']
##   Date:         ['09.02.2026']
##   Version:      ['0.1.0']
##   Status:       ['Development']
##   Language:     ['Python']
##   Description:  ['This module contains constants used throughout the RateKing package.']

## IMPORTS

from MoleKing import PeriodicTable

## CLASSES

class CONSTANTS:
    """Class to hold constants used in RateKing."""

    def __init__(self):
        self.PT = PeriodicTable()
        self.load_constants()

    def load_constants(self):
        """Method to load constants into the class."""

        self.k_B = self.PT.getConstant("k_B")                           # Boltzmann Constant (J/K)     
        self.h = self.PT.getConstant("h")                               # Planck Constant (J*s)       
        self.hbar = self.PT.getConstant("hbar")                         # Reduced Planck Constant (J*s) 
        self.c = self.PT.getConstant("C")                               # Speed of Light (cm/s)       
        self.R = self.PT.getConstant("R_cal")                           # Gas Constant (cal/(mol*K))
        self.R_J = self.PT.getConstant("R_J")                           # Gas Constant (J/(mol*K)) 
        self.PI = self.PT.getConstant("PI")                             # Pi                 
        self.P0  = self.PT.getConstant("P_0")                           # Pressure (Pa)                    
        self.NA = self.PT.getConstant("NA")                             # Avogadro's Number      
        self.AMU_2_KG = self.PT.getConversion("amu_to_kg")              # Conversion from AMU to Kg        
        self.AMU_2_KG_M2 = self.PT.getConversion("amu_to_kg_m2")        # Conversion from AMU to Kg/m^2    
        self.hartree2kcal = self.PT.getConversion("Hartree_to_Kcal")    # Conversion from Hartree to kcal/mol 
        self.kcal2kj = self.PT.getConversion("Kcal_to_KJ")              # Conversion from kcal to kJ        
