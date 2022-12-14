from common import *
import numpy as np

##############################################################################
# Section 3.1, pages 9 and 10
##############################################################################

# Setup
x = sp.symbols("x")

# Private to P
p = x ** 3 - 3 * x ** 2 + 2 * x

# The commitment which is public
t = x ** 2 - 3 * x + 2

# V samples random number
r = 23  # np.random.randint(1,10)
print("r", r)

t_r = t.evalf(subs={x: r})
print("t(r)", t_r)

# V sends r to P

# P calculates h(x)=p(x)/t(x)
h = sp.simplify(p / t)
print("h", h)

h_r = h.evalf(subs={x: r})
print(f"h({r}) {h_r}")

p_r = p.evalf(subs={x: r})
print(f"p({r}) {p_r}")

# P sends h and p to VS

# P checks that p = t*h
assert p_r == t_r * h_r, "Proof failed"
print(f"Proof passed! {p_r} == {t_r*h_r}")


##############################################################################
# Section 3.1, page 10, see paragraph that starts with "On the contrary..."
##############################################################################

# Prover
# This is not p:
p_prime = 2 * x ** 3 - 3 * x ** 2 + 2 * x
p_prime_r = p_prime.evalf(subs={x: r})
print(f"p'({r}) {p_prime_r}")

h_prime = sp.simplify(p_prime / t)
print("h'", h_prime)

h_prime_r = h_prime.evalf(subs={x: r})
print(f"h'({r}) {h_prime_r}")

# Check that p' is an integer
assert p_prime_r - int(p_prime_r) == 0

# Check that h'(r) is an integer (this should fail)
assert h_prime_r - int(h_prime_r) != 0, f"h'({r}) = {h_prime_r} is not an integer"


##############################################################################
# Remark 3.1, page 10. How to attack V
##############################################################################

# Assume we are using the same commitment t(x) and challenge r. P then
#   creates bogus h' and p' and returns it to V.

# This is supposed to be h(x) = p(x)/t(x) evaluated at challenge r. Instead we
#   just makes something up that will pass the test.
h_prime_r = np.random.randint(10)

# P generates a forged p' for the attack and returns p'(r) and h'(r) to V.
p_prime_r = t_r * h_prime_r

# Looking back to Example 3.1, we see this will pass.

##############################################################################
# Section 3.3.4, instantiating the protocol given on page 15
##############################################################################

from sympy.ntheory.residue_ntheory import primitive_root

# Select field size by hand. This is public.
prime = 2 ** 7 - 1
print("Using prime", prime)

# Select the smallest primitive root of GF(prime). This is public.
g = primitive_root(prime)
print("Using base g", g)

# Select secret
s = np.random.randint(2, prime)
print("Selecting secret s", s)

# Get the coefficients of the polyomial p
p_coeffs = sp.Poly(p).all_coeffs()

# V calculates E(s^0), E(s^1), ..., E(s^p_deg) and gives them to P.
s_encs = [pow(g, s ** i, prime) for i in range(len(p_coeffs))]
print("s_encs", s_encs)

# Using the same `p` and `t` as used in the `Section 3.1` example code.
print(f"p {p}")
print(f"t {t}")

# P evaluates the polynomial in the ct domain. This is E(p(s)) = g^p.
gp = evaluateEncrypted(p, s_encs, prime)
print("E(p(s)) = g^p", gp)

# P calculates p(x)/t(x). This is g^h.
h = sp.simplify(p / t)
print("h", h)

# E(h(x)) = g^h
hs = evaluateEncrypted(h, s_encs, prime)
print("E(h(s))", hs)

# P sends g^p and g^h back to V.

# V checks that (g^h)^t(s) = g^p.
ts = int(t.evalf(subs={x: s}))
assert hs ** ts % prime == gp
print("Verification passed")


##############################################################################
# Section 3.4, instantiating KEA protocol, page 16
##############################################################################

# Alice setup

# Public parameters
n = 2 ** 13 - 1
a = primitive_root(n)
print("a", a)

# Secret
alpha = np.random.randint(2, n)

# encrypt alpha
a_prime = a ** alpha % n

# Alice sends (a, a') to Bob

# Bob's steps

c = 11

b = a ** c % n
b_prime = a_prime ** c % n

# Send (b, b') to Alice

# Alice checks Bob's work

assert b ** alpha % n == b_prime, f"b**alpha {b**alpha} != b_prime {b_prime}"


##############################################################################
# Section 3.4, instantiating protocol started at the bottom of page 17
#
# Note: builds on results for `Section 3.3.4``
##############################################################################

# Verifier setup

# Select secret shift parameter
alpha = np.random.randint(2, prime)

# V calculates E(alpha*s^0), E(alpha*s^1), ..., E(alpha*s^p_deg) and gives them to P.
s_shift_encs = [pow(s_encs[i], alpha, prime) for i in range(len(s_encs))]
print("s_shift_encs", s_shift_encs)

# Send E(s), E(s^2), ..., E(s^d) and E(alpha*s), E(alpha*s^2), ..., E(alpha*s^d) to Prover.

# Prover works

# P evaluates the polynomial in the ct domain. This is E(p(s)) = g^p.
gpprime = evaluateEncrypted(p, s_shift_encs, prime)
print("E(alpha*p(s))", gpprime)

# Verifier checks

assert gp ** alpha % prime == gpprime

##############################################################################
# Sections 3.5 and 3.6 Zero-Knowledge
#
# Note: Now that we are using bilinear maps, we are creating new variables and
#   redefining `g`. In fact, the g used in the whitepaper is an acceptable
#   abuse of notation to keep things simple for the reader.
##############################################################################

##############################################################################
# 3.6.2 Trusted Party Setup, page 21
##############################################################################

import bplib as bp
from bplib import *

# Known to Prover
# p(x) = 0 + 2x - 3x^2 + x^3
p_coeffs = [0, 2, -3, 1]
# h(x) = x
h_coeffs = [0, 1]

# The commitment which is public
# t(x) = 2 - 3x + x^2
t_coeffs = [2, -3, 1]

# Notice that h(x)*t(x) = p(x)

# Create a set of bilinear EC groups
G = bp.BpGroup()
g1 = G.gen1()
g2 = G.gen2()

# Verifier selects random s and alpha
s = np.random.randint(100)
alpha = np.random.randint(100)

# Verifier creates the CRS proving key and sends to Prover
#   g^(s^i) = `s_encs`
s_encs1 = [g1 * (s ** i) for i in range(len(p_coeffs))]
s_encs2 = [g2 * (s ** i) for i in range(len(p_coeffs))]
#   g^(as^i) = `s_shift_encs`
s_shift_encs1 = [g1 * (alpha * s ** i) for i in range(len(p_coeffs))]
s_shift_encs2 = [g2 * (alpha * s ** i) for i in range(len(p_coeffs))]

def evaluate(s, coeffs):
    """Evaluates a polynomial on encrypted inputs

    Args:
    s - List of encrypted inputs.
    coeffs - List of coefficients in increasing order.
    """
    result = s[0] * coeffs[0]
    for i in range(1, len(coeffs)):
        result += s[i] * coeffs[i]

    return result


# Prover evaluates the polynomial in the ciphertext domain. This is E(p(s)) = g^p.
gp = evaluate(s_encs1, p_coeffs)

# Prover evaluates g^h in ciphertext domain.
gh = evaluate(s_encs2, h_coeffs)

# Prover sends g^p and g^h to Verifier

# Verifier evaluates g^t(s)
gt = evaluate(s_encs1, t_coeffs)

# Verifier checks that e(g^p, g^1) == e(g^t, g^h)
G.pair(gp, g2) == G.pair(gt, gh)

##############################################################################
# Section 3.7 zk-SNARKOP protocol, page 24
##############################################################################

# Verifier setup

# The commitment which is public
# t(x) = 2 - 3x + x^2
t_coeffs = [2, -3, 1]

# Sample random values
s = np.random.randint(100)
alpha = np.random.randint(100)

# Proving key
s_encs1 = [g1 * (s ** i) for i in range(len(p_coeffs))]
s_encs2 = [g2 * (s ** i) for i in range(len(p_coeffs))]
s_shift_encs1 = [g1 * (alpha * s ** i) for i in range(len(p_coeffs))]
s_shift_encs2 = [g2 * (alpha * s ** i) for i in range(len(p_coeffs))]

# Verification key
g_alpha = g2 * alpha
gt = evaluate(s_encs1, t_coeffs)

# Proof

# Assign coefficients
# p(x) = 0 + 2x - 3x^2 + x^3
p_coeffs = [0, 2, -3, 1]

# Calculate h(x) = p(x)/t(x)
# h(x) = x
h_coeffs = [0, 1]

# Evaluate g^p(s) and g^h(s) using encrypted powers of s
gp = evaluate(s_encs1, p_coeffs)
gh = evaluate(s_encs2, h_coeffs)

# Evaluate encrypted shifted polynomial g^(alpha*p(s)) using encrypted shifted
#   powers of s
gp_shift = evaluate(s_shift_encs1, p_coeffs)

# Sample random delta
delta = np.random.randint(100)

# Set randomized proof
gdp = delta * gp
gdh = delta * gh
gdap = delta * gp_shift

# Verification
gp = gdp
gh = gdh
gp_prime = gdap

# Check polynomial restriction
G.pair(gp_prime, g2) == G.pair(gp, g_alpha)

# Check polynomial cofactors
G.pair(gp, g2) == G.pair(gt, gh)
