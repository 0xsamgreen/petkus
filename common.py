import numpy as np
import sympy as sp


def evaluateEncrypted(poly, input_encrypted, prime):
    """Homomorphically evaluates a monomial.

    Args:
        poly: A SymPy monomial.
        input_encrypted: Monomial input encrypted by DLP.
        prime: Prime order of the finite field.
    """

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

    return Eps


def getCoefficients(y):
    """Solves interpolation problem using "set of equations with unknowns" method.

    Args:
        y: Vector of y-axis values that the resulting polynomial must pass through.

    Returns:
        Coefficients of the interpolated polynomial. Listed in increasing order.
    """

    # Let f(x) refer to the polynomial learned through interpolation.

    # A d-degree polynomial can fit d+1 points. Our interpolated polynomial
    #   must fit all points in the y vector.
    d = len(y) - 1

    # Need to eval f(x) at len(y) arbitrary points. Choosing x = 1, 2, ..., len(y).
    all_x = np.arange(1, len(y) + 1)

    # `A` matrix will contain powers of x. One row for each x in all_x.
    A = np.empty(shape=(len(all_x), len(y)))
    for row, x in enumerate(all_x):
        for term in range(d + 1):
            A[row, term] = x ** term

    #     | x_0^0      x_0^1 ...      x0^(d+1)       |
    # A = | x_1^0      x_1^1 ...      x1^(d+1)       |
    #     |                  ...                     |
    #     | x_len(y)^0 x_len(y)^1 ... x_len(y)^(d+1) |

    # Solve A*coeff = y for the `coeff` vector. The `coeff` vector
    #   will become the coefficients for f(x).
    A_inv = np.linalg.inv(A)
    coeff = A_inv @ y

    # Optional: test that `coeff` is the correct solution.
    assert (A @ coeff == y).all(), "Error with solution."

    return coeff
