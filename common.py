import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


def evaluate_encrypted(poly, input_encrypted, prime):
    '''Homomorphically evaluates a monomial.

    Args:
        poly: A SymPy monomial.
        input_encrypted: Monomial input encrypted by DLP.
        prime: Prime order of the finite field.
    '''

    # Get the coefficients of the monomial.
    coeffs = sp.Poly(poly).all_coeffs()
    # By default, coeffs will be in standard form, with coefficients in decending
    #   order of the power. Reverse the order.
    coeffs.reverse()

    # Initialize the accumulator for the monomial evaluation, E(p(s)).
    Eps = 1

    # Loop coefficients and evaluate the term of the monomial
    for idx, coeff in enumerate(coeffs):
        # `coeff ` is a SymPy type. Need to convert it to an integer.
        coeff = int(coeff)
        # Exponentiate the encrypted input by the coefficient. E.g.
        #   (g^s)^c_i = g^(s^c_i) and multiply that by the current value of the 
        #   accumulator `Eps`. 
        Eps *= pow(input_encrypted[idx], coeff, prime)
    
    # Ensure the accumulator is an element of the field.
    Eps = Eps % prime

    return(Eps)