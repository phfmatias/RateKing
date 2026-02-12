##   PYTHON FILE HEADER #
##
##   File:         [reaction.py]
##
##   Author(s):    ['Pedro H.F Matias']
##   Site(s):      ['https://github.com/phfmatias']
##   Email(s):     ['phfmatias@discente.ufg.br']
##   Credits:      ['Copyright © 2026 LEEDMOL. All rights reserved.']
##   Date:         ['10.02.2026']
##   Version:      ['0.1.0']
##   Status:       ['Development']
##   Language:     ['Python']
##   Description:  ['Module to define the Reaction class for RateKing.']

from ..utils.constants import CONSTANTS
from ..exceptions import RateKingError
from .molecule import Molecule, TransitionState 

class Reaction(object):
    """Class to represent a chemical reaction in RateKing."""

    def __init__(self, name, reactantsPath, productsPath, transitionStatePath=None, temperature=298.15):
        self.Constants = CONSTANTS()
        self.name = name
        self.reactantsPath = reactantsPath  # List of dictionaries with keys: name, filePath, link
        self.productsPath = productsPath    # List of dictionaries with keys: name, filePath, link (optional)
        self.transitionStatePath = transitionStatePath  # Dictionary with keys: name, filePath, link (optional)
        self.reactants = []  # List of Molecule objects 
        self.products = []    # List of Molecule objects (optional)
        self.transitionState = None  # TransitionState object 
        self.temperature = temperature

        # Loading reaction data and initializing properties
        self.loadReaction(temperature)
        
    def loadReactants(self, temperature=298.15):
        """Method to load the reactants' properties."""
        if len(self.reactantsPath) == 0:
            raise RateKingError("Error! No reactants provided for the reaction.")
        try:
            for reactant in range(len(self.reactantsPath)):
                # reactant is a dictionary with keys: name, filePath, link
                name = self.reactantsPath[reactant]["name"]
                filePath = self.reactantsPath[reactant]["filePath"]
                link = self.reactantsPath[reactant]["link"]  
                self.reactants.append(Molecule(name, filePath, link=link, temperature=temperature))
        except Exception as e:
            raise RateKingError(f"Error! Failed to load reactants' properties. Details: {str(e)}")

    def loadProducts(self, temperature=298.15):
        """Method to load the products' properties."""
        try:
            for product in range(len(self.productsPath)):
                # product is a dictionary with keys: name, filePath, link
                name = self.productsPath[product]["name"]
                filePath = self.productsPath[product]["filePath"]
                link = self.productsPath[product]["link"]
                self.products.append(Molecule(name, filePath, link=link, temperature=temperature))
        except Exception as e:
            raise RateKingError(f"Error! Failed to load products' properties. Details: {str(e)}")

    def loadTransitionState(self, temperature=298.15):
        """Method to load the transition state's properties."""
        if self.transitionStatePath is None:
            raise RateKingError("Error! No transition state provided for the reaction.")
        try:
            name = self.transitionStatePath["name"]
            filePath = self.transitionStatePath["filePath"]
            link = self.transitionStatePath["link"]
            self.transitionState = TransitionState(name, filePath, link=link, temperature=temperature)
        except Exception as e:
            raise RateKingError(f"Error! Failed to load transition state's properties. Details: {str(e)}")

    def loadReaction(self, temperature=298.15):
        """Method to load all properties of the reaction."""
        self.loadReactants(temperature)
        self.loadTransitionState(temperature)
        if len(self.productsPath) > 0:
            self.loadProducts(temperature)

    def validateThermochemistryData(self):
        """Method to validate that all reactants, products, and transition state (if present) have valid thermochemical data."""
        if len(self.productsPath) == 0:
            self.productsEnergy = "N/A"
            self.productsEnthalpy = "N/A"
            self.productsGibbs = "N/A"
            self.deltaE = "N/A"
            self.deltaH = "N/A"
            self.deltaG = "N/A"
            self.reverseActivationEnergy = "N/A"
            return False
        return True

    def run(self):
        """Method to run the reaction calculations."""

        self.getReactantsEnergy()
        self.getReactantsEnthalpy()
        self.getReactantsGibbs()
        self.getActivationEnergy()
        self.getEnthalpyBarrier()
        self.getGibbsBarrier()
        self.getReactionSymmetryNumber()
        self.getReactionTotalPartitionFunction()
        self.getTransitionStateImaginaryFrequency()
        
        if self.validateThermochemistryData():
            self.getProductsEnergy()
            self.getProductsEnthalpy()
            self.getProductsGibbs()
            self.getDeltaE()
            self.getDeltaH()
            self.getDeltaG()
            self.getReverseActivationEnergy()

    def getReactantsEnergy(self):
        """Method to calculate the total energy of the reactants."""
        self.reactantsEnergy = 0
        
        for reactant in self.reactants:
            self.reactantsEnergy += reactant.zeroPointEnergy
  
    def getProductsEnergy(self):
        """Method to calculate the total energy of the products."""
        self.productsEnergy = 0
        
        for product in self.products:
            self.productsEnergy += product.zeroPointEnergy

    def getProductsEnthalpy(self):
        """Method to calculate the total enthalpy of the products."""
        self.productsEnthalpy = 0
        
        for product in self.products:
            self.productsEnthalpy += product.enthalpy
        
    def getReactantsEnthalpy(self):
        """Method to calculate the total enthalpy of the reactants."""
        self.reactantsEnthalpy = 0
        
        for reactant in self.reactants:
            self.reactantsEnthalpy += reactant.enthalpy

    def getProductsGibbs(self):
        """Method to calculate the total Gibbs free energy of the products."""
        self.productsGibbs = 0
        
        for product in self.products:
            self.productsGibbs += product.gibbs

    def getReactantsGibbs(self):
        """Method to calculate the total Gibbs free energy of the reactants."""
        self.reactantsGibbs = 0
        
        for reactant in self.reactants:
            self.reactantsGibbs += reactant.gibbs

    def getReactantsSymmetry(self):
        self.reactantsSymmetry = 1
        
        for reactant in self.reactants:
            self.reactantsSymmetry *= reactant.symmetryNumber

    def getTransitionStateSymmetry(self):
        self.transitionStateSymmetry = self.transitionState.symmetryNumber

    def getTransitionStateImaginaryFrequency(self):
        self.transitionStateImaginaryFrequency = self.transitionState.imaginaryFrequency
    
    def getReactantsTotalPartitionFunction(self):
        """Method to calculate the total partition function of the reactants."""
        self.reactantsTotalPartitionFunction = 1
        
        for reactant in self.reactants:
            self.reactantsTotalPartitionFunction *= reactant.qTotal
    
    def getTransitionStateTotalPartitionFunction(self):
        """Method to calculate the total partition function of the transition state."""
        self.transitionStateTotalPartitionFunction = self.transitionState.qTotal

    def getReactionTotalPartitionFunction(self):
        """Method to calculate the total partition function of the reaction."""
        self.getReactantsTotalPartitionFunction()
        self.getTransitionStateTotalPartitionFunction()

        self.reactionTotalPartitionFunction = self.transitionStateTotalPartitionFunction / self.reactantsTotalPartitionFunction
    
    # DeltaE = E_products - E_reactants
    def getDeltaE(self):
        """Method to calculate the energy change of the reaction."""
        self.getReactantsEnergy()
        self.getProductsEnergy()
        self.deltaE = (self.productsEnergy - self.reactantsEnergy) * self.Constants.hartree2kcal # kcal/mol

    # DeltaH = H_products - H_reactants
    def getDeltaH(self):
        """Method to calculate the enthalpy change of the reaction."""
        self.getReactantsEnthalpy()
        self.getProductsEnthalpy()
        self.deltaH = (self.productsEnthalpy - self.reactantsEnthalpy) * self.Constants.hartree2kcal # kcal/mol

    # DeltaG = G_products - G_reactants
    def getDeltaG(self):
        """Method to calculate the Gibbs free energy change of the reaction."""
        self.getReactantsGibbs()
        self.getProductsGibbs()
        self.deltaG = (self.productsGibbs - self.reactantsGibbs) * self.Constants.hartree2kcal # kcal/mol

    # Ea = E_TS - E_reactants
    def getActivationEnergy(self):
        """Method to calculate the activation energy of the reaction."""
        self.activationEnergy = (self.transitionState.zeroPointEnergy - self.reactantsEnergy)*self.Constants.hartree2kcal # kcal/mol
        self.activationEnergy_cal = self.activationEnergy * 1000 # cal/mol

    # Ea_reverse = E_TS - E_products
    def getReverseActivationEnergy(self):
        """Method to calculate the reverse activation energy of the reaction."""
        self.reverseActivationEnergy = (self.transitionState.zeroPointEnergy - self.productsEnergy) * self.Constants.hartree2kcal # kcal/mol
        self.reverseActivationEnergy_cal = self.reverseActivationEnergy * 1000 # cal/mol

    # DeltaH‡ = H_TS - H_reactants
    def getEnthalpyBarrier(self):
        """Method to calculate the enthalpy barrier of the reaction."""
        self.enthalpyBarrier = (self.transitionState.enthalpy - self.reactantsEnthalpy)

    # DeltaG‡ = G_TS - G_reactants
    def getGibbsBarrier(self):
        """Method to calculate the Gibbs free energy barrier of the reaction."""
        self.gibbsBarrier = (self.transitionState.gibbs - self.reactantsGibbs)

    def getReactionSymmetryNumber(self):
        """Method to calculate the reaction symmetry number."""
        self.getReactantsSymmetry()
        self.getTransitionStateSymmetry()
        self.reactionSymmetryNumber = (self.reactantsSymmetry / self.transitionStateSymmetry)

