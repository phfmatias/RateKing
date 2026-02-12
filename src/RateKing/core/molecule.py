##   PYTHON FILE HEADER #
##
##   File:         [molecule.py]
##
##   Author(s):    ['Pedro H.F Matias']
##   Site(s):      ['https://github.com/phfmatias']
##   Email(s):     ['phfmatias@discente.ufg.br']
##   Credits:      ['Copyright © 2026 LEEDMOL. All rights reserved.']
##   Date:         ['10.02.2026']
##   Version:      ['0.1.0']
##   Status:       ['Development']
##   Language:     ['Python']
##   Description:  ['Module to define the Molecule class for RateKing.']

from MoleKing import G16LOGfile
from ..exceptions import RateKingError, ExtractionError

class Molecule(object):
    """Class to represent a molecule in RateKing."""

    def __init__(self, name='Molecule', filePath=None, zeroPointEnergy=0, enthalpy=0, entropy=0, gibbs=0, qElectronic=0, qVibrational=0, qRotational=0, qTranslational=0, link=-1, temperature=298.15):
        self.name = name
        self.filePath = filePath
        self.zeroPointEnergy = zeroPointEnergy
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.gibbs = gibbs
        self.qElectronic = qElectronic
        self.qVibrational = qVibrational
        self.qRotational = qRotational
        self.qTranslational = qTranslational
        self.mkMolecule = None
        self.link = link
        self.temperature = temperature

        # Loading molecule and initializing properties
        self.loadMoleculeMoleking()
        self.initializeAllProperties()

    @staticmethod
    def getAvailableLinks(filePath):
        """Static method to get the number of available links in a file without full initialization.
        
        Use this in your web application to validate user input before loading.
        Returns the total number of links available in the file.
        """
        try:
            temp_molecule = G16LOGfile(filePath, link=-1, thermoAsw=True)
            total_links = temp_molecule.getLinkSize()
            return int(total_links)
        except Exception as e:
            raise RateKingError(f"Error! Could not read link information from file: {filePath}. Details: {str(e)}")

    def loadMoleculeMoleking(self):
        """Method to load the molecule data using MoleKing."""

        try:
            self.mkMolecule = G16LOGfile(self.filePath, link=self.link, thermoAsw=True)
        except Exception as e:
            raise RateKingError(f"Error! MoleKing failed to load the molecule from the file: {self.filePath}. Details: {str(e)}")
    
    def getZPE(self):
        """Method to get the zero-point energy of the molecule."""
        try:
            self.zeroPointEnergy = self.mkMolecule.getZPE()
        except Exception as e:
            raise ExtractionError(f"Error! Failed to extract Zero-Point Energy from the molecule file: {self.filePath}. Details: {str(e)}")
            
        
    def getEnthalpy(self):
        """Method to get the enthalpy of the molecule."""
        try:
            self.enthalpy = self.mkMolecule.getH()
        except Exception as e:
            raise ExtractionError(f"Error! Failed to extract Enthalpy from the molecule file: {self.filePath}. Details: {str(e)}")
        
    def getEntropy(self):
        """Method to get the entropy of the molecule."""
        try:
            self.entropy = self.mkMolecule.getS()
        except Exception as e:
            raise ExtractionError(f"Error! Failed to extract Entropy from the molecule file: {self.filePath}. Details: {str(e)}")
        
    def getGibbs(self):
        """Method to get the Gibbs free energy of the molecule."""
        try:
            self.gibbs = self.mkMolecule.getG()
        except Exception as e:
            raise ExtractionError(f"Error! Failed to extract Gibbs Free Energy from the molecule file: {self.filePath}. Details: {str(e)}")
        
    def getSymmetryNumber(self):
        """Method to get the symmetry number of the molecule."""
        try:
            self.symmetryNumber = self.mkMolecule.get_sigmaR()
        except Exception as e:
            raise ExtractionError(f"Error! Failed to extract Symmetry Number from the molecule file: {self.filePath}. Details: {str(e)}")

    def getElectronicPartitionFunction(self):
        """Method to get the electronic partition function of the molecule."""
        try:
            self.qElectronic = self.mkMolecule.get_qEle()
        except Exception as e:
            raise ExtractionError(f"Error! Failed to extract Electronic Partition Function from the molecule file: {self.filePath}. Details: {str(e)}")
    
    def getVibrationalPartitionFunction(self):
        """Method to get the vibrational partition function of the molecule."""
        try:
            self.qVibrational = self.mkMolecule.get_qVib(Temperature=self.temperature)
        except Exception as e:
            raise ExtractionError(f"Error! Failed to extract Vibrational Partition Function from the molecule file: {self.filePath}. Details: {str(e)}")
    
    def getRotationalPartitionFunction(self):
        """Method to get the rotational partition function of the molecule."""
        try:
            self.qRotational = self.mkMolecule.get_qRot(Temperature=self.temperature)
        except Exception as e:
            raise ExtractionError(f"Error! Failed to extract Rotational Partition Function from the molecule file: {self.filePath}. Details: {str(e)}")
    
    def getTranslationalPartitionFunction(self):
        """Method to get the translational partition function of the molecule."""
        try:
            self.qTranslational = self.mkMolecule.get_qTrans(Temperature=self.temperature)
        except Exception as e:
            raise ExtractionError(f"Error! Failed to extract Translational Partition Function from the molecule file: {self.filePath}. Details: {str(e)}")
    
    def getTotalPartitionFunction(self):
        """Method to get the total partition function of the molecule."""
        try:
            self.qTotal = self.mkMolecule.get_qTot(Temperature=self.temperature)
        except Exception as e: 
            raise ExtractionError(f"Error! Failed to extract Total Partition Function from the molecule file: {self.filePath}. Details: {str(e)}")

    def initializeAllProperties(self):
        """Method to initialize all thermodynamic properties and partition functions of the molecule."""
        self.getZPE()
        self.getEnthalpy()
        self.getEntropy()
        self.getGibbs()
        self.getSymmetryNumber()
        self.getElectronicPartitionFunction()
        self.getVibrationalPartitionFunction()
        self.getRotationalPartitionFunction()
        self.getTranslationalPartitionFunction()
        self.getTotalPartitionFunction()

class TransitionState(Molecule):
    """Class to represent a transition state in RateKing, inheriting from Molecule."""
    def __init__(self, name, filePath, zeroPointEnergy=0, enthalpy=0, entropy=0, gibbs=0, qElectronic=0, qVibrational=0, qRotational=0, qTranslational=0, link=-1, temperature=298.15):

        super().__init__(name, filePath, zeroPointEnergy, enthalpy, entropy, gibbs, qElectronic, qVibrational, qRotational, qTranslational, link, temperature)
        
        self.getImaginaryFrequency()

        
    def getImaginaryFrequency(self):
        """Method to get the imaginary frequency of the transition state."""
        try:
            freqs = self.mkMolecule.getVibFrequencies()
            negativeFreqs = [freq for freq in freqs if freq < 0]
            
            if len(negativeFreqs) > 1:
                raise ExtractionError(f"Error! More than one imaginary frequency found in the transition state file: {self.filePath}. Details: Found frequencies: {negativeFreqs}")
            
            self.imaginaryFrequency = negativeFreqs[0]

        except Exception as e:
            raise ExtractionError(f"Error! Failed to extract Imaginary Frequency from the transition state file: {self.filePath}. Details: {str(e)}")
    