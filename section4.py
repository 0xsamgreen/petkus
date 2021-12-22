from common import *

# Reserving x for use in building SymPy polynomials.
x = sp.symbols("x")

##############################################################################
# Example 4.1
##############################################################################

# 1. Verifier makes it publicly known that they would like a prover to take the
#   product of the following sequence.
vec = [2, 1, 3, 2]

# 2. Prover commits the following polynomial
t_x = sp.Poly((x - 1) * (x - 2) * (x - 3))

# 3. Prover makes polynomials out of the following system of equations:
#
#   2 x 1 = 2
#   2 x 3 = 6
#   6 x 2 = 12
#
# Note that the method by which polynomials are created from that
#   system is very usual from what you may have seen before.

# left-hand side (the green numbers in the document)
l_vec = [vec[0]]

# right-hand side (the blue numbers in the document)
r_vec = vec[1:]

# outputs (the red numbers)
o_vec = []

# Fill in the LHS and output vectors
for i in range(len(vec[:-1])):
    print(i, l_vec, r_vec)
    o_vec.append(l_vec[i] * r_vec[i])

    if i < len(vec) - 2:
        l_vec.append(o_vec[i])

# Each vector (l_vec, r_vec, and o_vec) now has n-1 elements.
#  Will now treat those vectors as the outputs for the l(x),
#   r(x), and o(x) polynomials respectively.

# The next step is to determine the coefficients for those
#   polynomials using interpolation.
l_coeff = getCoefficients(l_vec)
print(f"l(x) coeff.", l_coeff)
r_coeff = getCoefficients(r_vec)
print(f"r(x) coeff.", r_coeff)
o_coeff = getCoefficients(o_vec)
print(f"o(x) coeff.", o_coeff)

# Convert the coefficients to SymPy polynomials.

# Left operand polynomial. `l_x` => l(x)
l_x = sp.Poly.from_list(reversed(l_coeff), gens=x)
print("l(x)", l_x)

# Right operand polynomial. `r_x` => r(x)
r_x = sp.Poly.from_list(reversed(r_coeff), gens=x)
print("r(x)", r_x)

# Output polynomial. `o_x` => o(x)
o_x = sp.Poly.from_list(reversed(o_coeff), gens=x)
print("o(x)", o_x)

# 3. Prover makes the "operation polynomial", p(x):
p_x = l_x * r_x - o_x

# 4. Prover calculates a cofactor of the operation polynomial, using public
#   commitment t(x).
h_x = sp.Poly(sp.simplify(p_x / t_x))

# Prover sends h(x) and p(x) to Verifier.

# 5. Verifier checks that t(x)*h(x) = p(x)
assert h_x * t_x == p_x, "Verification failed."

# Note: this not an example of a secure proof, but it will be made more complete
#   in the sections ahead.


##############################################################################
# Example 4.2

# Proving
# l(x) * r(x) = o(x)
# a    * b    = r1
# a    * c    = r2
##############################################################################

# 1. Construct operand polynomial l(x)
l_vec = [3, 3, 1]
l_coeff = getCoefficients(l_vec)
print(l_coeff)

# 2. Sample random alpha and s
alpha = np.random.randint(2, 2 ** 8)
s = np.random.randint(2, 2 ** 8)

# 3. Set the proving key
