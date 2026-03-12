
from RateKing import CONSTANTS
from RateKing import Molecule, TransitionState, Reaction
from RateKing import TST, dTST
from RateKing.tunneling_corrections import eckartCorrection, bell35Correction, bell58Correction, bell2TCorrection, wignerCorrection

from MoleKing import G16LOGfile

from pandas import DataFrame

class DEBUG_ERROS:
    """Class to hold debug errors for RateKing."""

    def __init__(self):
        pass

    def eval_RK_Error(self):
        """Method to evaluate the RK_Error for debugging."""
        self.debug.RK_Error()

    def eval_File_Validation_Error(self):
        """Method to evaluate the File_Validation_Error for debugging."""
        self.debug.File_Validation_Error()

    def eval_Calculation_Error(self):
        """Method to evaluate the Calculation_Error for debugging."""
        self.debug.Calculation_Error()

    def eval_Imaginary_Frequency_Not_Found_Error(self):
        """Method to evaluate the Imaginary_Frequency_Not_Found_Error for debugging."""
        self.debug.Imaginary_Frequency_Not_Found_Error()

class DEBUG_CONSTANTS:
    """Class to hold debug constants for RateKing."""

    def __init__(self):
        self.constants = CONSTANTS()
        self.constants.load_constants()

    def eval_constants(self):
        """Method to evaluate the constants for debugging."""
        print("Boltzmann Constant (k_B):", self.constants.k_B)
        print("Planck Constant (h):", self.constants.h)
        print("Reduced Planck Constant (hbar):", self.constants.hbar)
        print("Speed of Light (c):", self.constants.c)
        print("Gas Constant (R):", self.constants.R)
        print("Gas Constant (R_J):", self.constants.R_J)
        print("Pi (PI):", self.constants.PI)
        print("Pressure (P0):", self.constants.P0)
        print("Conversion from AMU to Kg (AMU_2_KG):", self.constants.AMU_2_KG)
        print("Conversion from AMU to Kg/m^2 (AMU_2_KG_M2):", self.constants.AMU_2_KG_M2)
        print("Conversion from Hartree to kcal/mol (hartree2kcal):", self.constants.hartree2kcal)
        print("Conversion from kcal to kJ (kcal2kj):", self.constants.kcal2kj)
        print("Avogadro's Number (NA):", self.constants.NA)

class MOLECULE_DEBUG:
    """Class to hold debug for Molecule class in RateKing."""

    def __init__(self):
        self.molecule = Molecule("TestMolecule", "../RateKingWeb/TEST_CASES/NH/nh3_m062x.out")

    def evalProperties(self):
        """Method to evaluate the Zero-Point Energy extraction."""
        print("Zero-Point Energy (ZPE):", self.molecule.zeroPointEnergy)
        print("Enthalpy (H):", self.molecule.enthalpy)
        print("Entropy (S):", self.molecule.entropy)
        print("Gibbs Free Energy (G):", self.molecule.gibbs)
        print("Electronic Partition Function (qElectronic):", self.molecule.qElectronic)
        print("Vibrational Partition Function (qVibrational):", self.molecule.qVibrational)
        print("Rotational Partition Function (qRotational):", self.molecule.qRotational)
        print("Translational Partition Function (qTranslational):", self.molecule.qTranslational)
        print("Total Partition Function (qTotal):", self.molecule.qTotal)

class TS_DEBUG:
    """Class to hold debug for TransitionState class in RateKing."""

    def __init__(self):
        self.ts = TransitionState("TestTS", "../RateKingWeb/TEST_CASES/NH/nh3oh_ts_m062x.out")

    def evalProperties(self):
        """Method to evaluate the properties of the transition state."""
        print("Zero-Point Energy (ZPE):", self.ts.zeroPointEnergy)
        print("Enthalpy (H):", self.ts.enthalpy)
        print("Entropy (S):", self.ts.entropy)
        print("Gibbs Free Energy (G):", self.ts.gibbs)
        print("Electronic Partition Function (qElectronic):", self.ts.qElectronic)
        print("Vibrational Partition Function (qVibrational):", self.ts.qVibrational)
        print("Rotational Partition Function (qRotational):", self.ts.qRotational)
        print("Translational Partition Function (qTranslational):", self.ts.qTranslational)
        print("Total Partition Function (qTotal):", self.ts.qTotal)
        print("Imaginary Frequency:", self.ts.imaginaryFrequency)

class REACTION_DEBUG:
    """Class to hold debug for Reaction class in RateKing."""

    def __init__(self):
        # self.reactant1 = Molecule("NH3", "../RateKingWeb/TEST_CASES/NH/nh3_m062x.out", link=2)
        # self.reactant2 = Molecule("OH", "../RateKingWeb/TEST_CASES/NH/oh_m062x.out", link=2)
        # self.ts = TransitionState("TestTS", "../RateKingWeb/TEST_CASES/NH/nh3oh_ts_m062x.out", link=2)
        # self.product1 = Molecule("NH2", "../RateKingWeb/TEST_CASES/NH/nh2_m062x.out", link=2)
        # self.product2 = Molecule("H2O", "../RateKingWeb/TEST_CASES/NH/h2o_m062x.out", link=2)
        # self.reaction = Reaction("TestReaction", [self.reactant1, self.reactant2], [self.product1, self.product2], self.ts)

        self.reactant1 = {"name": "NH3", "filePath": "../RateKingWeb/TEST_CASES/NH/nh3_m062x.out", "link": 2}
        self.reactant2 = {"name": "OH", "filePath": "../RateKingWeb/TEST_CASES/NH/oh_m062x.out", "link": 2}
        self.ts = {"name": "TestTS", "filePath": "../RateKingWeb/TEST_CASES/NH/nh3oh_ts_m062x.out", "link": 2}
        self.product1 = {"name": "NH2", "filePath": "../RateKingWeb/TEST_CASES/NH/nh2_m062x.out", "link": 2}
        self.product2 = {"name": "H2O", "filePath": "../RateKingWeb/TEST_CASES/NH/h2o_m062x.out", "link": 2}

        reactants = [self.reactant1, self.reactant2]
        products = [self.product1, self.product2]

        self.reaction = Reaction("TestReaction", reactants, products, self.ts, temperature=298)
        self.reaction.run()

    def evalProperties(self):
        """Method to evaluate the properties of the reaction."""
        print("Delta E:", self.reaction.deltaE)
        print("Delta H:", self.reaction.deltaH)
        print("Delta G:", self.reaction.deltaG)
        print("Energy Reactants:", self.reaction.reactantsEnergy)
        print("Energy Products:", self.reaction.productsEnergy)
        print("Activation Energy:", self.reaction.activationEnergy)
        print("Reverse Energy Profile:", self.reaction.reverseActivationEnergy)
        print("Enthalpy Barrier:", self.reaction.enthalpyBarrier)
        print("Gibbs Free Energy Barrier:", self.reaction.gibbsBarrier)
        print("Reaction Symmetry Number:", self.reaction.reactionSymmetryNumber)
        print("Reaction Total Partition Function:", self.reaction.reactionTotalPartitionFunction)

class KINECTICS_AND_TUNNELING_DEBUG:
    """Class to hold debug for TST class in RateKing."""

    def __init__(self):
        self.reactant1 = {"name": "NH3", "filePath": "../RateKingWeb/TEST_CASES/NH/nh3_m062x.out", "link": 2}
        self.reactant2 = {"name": "OH", "filePath": "../RateKingWeb/TEST_CASES/NH/oh_m062x.out", "link": 2}
        self.ts = {"name": "TestTS", "filePath": "../RateKingWeb/TEST_CASES/NH/nh3oh_ts_m062x.out", "link": 2}
        self.product1 = {"name": "NH2", "filePath": "../RateKingWeb/TEST_CASES/NH/nh2_m062x.out", "link": 2}
        self.product2 = {"name": "H2O", "filePath": "../RateKingWeb/TEST_CASES/NH/h2o_m062x.out", "link": 2}    
        self.reactants = [self.reactant1, self.reactant2]
        self.products = [self.product1, self.product2]

        # self.reactant1 = {"name": "et1", "filePath": "../RateKingWeb/TEST_CASES/ORCA/et1.out", 'link': 1}
        # self.reactant2 = {"name": "H", "filePath": "../RateKingWeb/TEST_CASES/ORCA/H.out", 'link': 1}
        # self.ts = {"name": "TSA_1", "filePath": "../RateKingWeb/TEST_CASES/ORCA/TSA_1.out", 'link': 1}

        # reactants = [self.reactant1, self.reactant2]
        # products = []

    def evalProperties(self):
        """Method to evaluate the rate constant calculated using TST."""
        initial = 200
        end = 1000
        step = 100
        fixed_temp = 298.15

        temperature_list = sorted([fixed_temp] + list(range(initial, end + step, step)))

        print("T (K) | kTST | kEckart | kBell35 | kBell58 | kBell2T | kWigner | kDTST")

        for i in temperature_list:
            self.reaction = Reaction("TestReaction", self.reactants, self.products, self.ts, temperature=i)
            self.reaction.run()
            
            tst = TST(self.reaction).calculateKTST()
            dtst = dTST(self.reaction).calculateKDTST()
            kEckart = tst * eckartCorrection(self.reaction).calculateEckartCorrection()
            kBell35 = tst * bell35Correction(self.reaction).calculateBell35Correction()
            kBell58 = tst * bell58Correction(self.reaction).calculateBell58Correction()
            kBell2T = tst * bell2TCorrection(self.reaction).calculateBell2TCorrection()
            kWigner = tst * wignerCorrection(self.reaction).calculateWignerCorrection()

            print(f"{i:.2f} | {tst:.2e} | {kEckart:.2e} | {kBell35:.2e} | {kBell58:.2e} | {kBell2T:.2e} | {kWigner:.2e} | {dtst:.2e}")\


if __name__ == "__main__":
    debug_errors = DEBUG_ERROS()
    debug_constants = DEBUG_CONSTANTS()
    molecule_debug = MOLECULE_DEBUG()
    ts_debug = TS_DEBUG()
    reaction_debug = REACTION_DEBUG()
    kinectics_and_tunneling_debug = KINECTICS_AND_TUNNELING_DEBUG()    

    #! Evaluate debug errors
    # debug_errors.eval_RK_Error()
    # debug_errors.eval_File_Validation_Error()
    # debug_errors.eval_Calculation_Error()
    # debug_errors.eval_Imaginary_Frequency_Not_Found_Error()

    #! Evaluate constants
    # debug_constants.eval_constants()

    #! Evaluate Molecule class debug
    # molecule_debug.evalProperties()

    #! Evaluate Transition State class debug
    #ts_debug.evalProperties()
   
    #! Evaluate Reaction class debug
    # reaction_debug.evalProperties()

    #! Evaluate TST class debug
    kinectics_and_tunneling_debug.evalProperties()