##   PYTHON FILE HEADER #
##
##   File:         [exceptions.py]
##
##   Author(s):    ['Pedro H.F Matias']
##   Site(s):      ['https://github.com/phfmatias']
##   Email(s):     ['phfmatias@discente.ufg.br']
##   Credits:      ['Copyright © 2026 LEEDMOL. All rights reserved.']
##   Date:         ['09.02.2026']
##   Version:      ['0.1.0']
##   Status:       ['Development']
##   Language:     ['Python']
##   Description:  ['This module contains custom exceptions for the RateKing package.']

class RateKingError(Exception):
    """Base exception class for RateKing package."""
    
    def __init__(self, message="An error occurred in RateKing package."):
        self.message = message
        super().__init__(self.message)

class ExtractionError(RateKingError):
    """Exception raised for errors in data extraction."""

    def __init__(self, message="Data extraction error occurred."):
        self.message = message
        super().__init__(self.message)

class FileValidationError(RateKingError):
    """Exception raised for errors in file validation."""

    def __init__(self, message="File validation error occurred."):
        self.message = message
        super().__init__(self.message)

class CalculationError(RateKingError):
    """Exception raised for errors in calculations."""

    def __init__(self, message="An error occurred during calculation."):
        self.message = message
        super().__init__(self.message)

class ImaginaryFrequencyNotFoundError(CalculationError):
    """Exception raised when no imaginary frequency is found in a transition state."""

    def __init__(self, message="No imaginary frequency found in the transition state."):
        self.message = message
        super().__init__(self.message)