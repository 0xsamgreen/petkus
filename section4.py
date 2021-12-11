
from common import *

##############################################################################
# Algorithm 3.1
##############################################################################

# Setup
x = sp.symbols('x')

# Private to P
p = x**3 - 3*x**2 + 2*x

# The commitment which is public
t = x**2 - 3*x + 2

# 3.1.a
r = 23 #np.random.randint(1,10)
print('r', r)

t_r = t.evalf(subs={x:r})
print('t(r)', t_r)

# V sends r to P

# 3.1.b
h = sp.simplify(p/t)
print('h', h)

h_r = h.evalf(subs={x:r})
print(f'h({r}) {h_r}')

p_r = p.evalf(subs={x:r})
print(f'p({r}) {p_r}')

# Prover sends h and p to Verifier

# 3.1.c
assert p_r==t_r*h_r, 'Proof failed'
print(f'Proof passed! {p_r} == {t_r*h_r}')


# Using the "Set of equations with unknowns" interpolation method
#   to prove a computation of the form x1*x2*x3*...*xn.
vec = [2,1,3,2]

# left-hand side (the green numbers in the document)
lhs = [vec[0]]

# right-hand side (the blue numbers in the document)
rhs = vec[1:]

# intermediate products (the red numbers)
intermediates = []

# Fill in the LHS, RHS, and output vectors
for i in range(len(vec[:-1])):
    intermediates.append(lhs[i]*rhs[i])

    if i < len(vec) - 2:
        lhs.append(intermediates[i])

# Each vector (LHS, RHS, and intermediates) now has n-1 elements.
# Will now treat those vectors as the outputs for the l(x),
#   r(x), and o(x) polynomials respectively.
# The next step is to determine the coefficients for those
#   polynomials using interpolation.

# A degree d polynomial is needed to fit d+1 points.
# Recall that `lhs`, `rhs`, and `intermediates` all have the same
#   length.
d = len(lhs)-1

# Need to eval l(x) at len(lhs) points. Choosing x = 1, 2,  and 3.
all_x = [1,2,3]

A = np.empty(shape=(len(all_x), d+1))

for row, x in enumerate(all_x):
    for term in range(d+1):
        A[row, term] = x**term

# Solve A*coeff = lhs for the `coeff` vector. The `coeff` vector
#   will become the coefficients for the lhs(x) polynomial.
A_inv = np.linalg.inv(A)
coeff = A_inv@lhs

# Optional: test that `coeff` is the correct solution.
assert (A@coeff == lhs).all(), 'Error with solution.'

def getCoefficients(y_vals):
    '''asdf
    '''