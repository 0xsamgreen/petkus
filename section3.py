from common import *

##############################################################################
# Example 3.1
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


##############################################################################
# Example 3.2
##############################################################################

# Prover
# This is not p:
p_prime = 2*x**3 - 3*x**2 + 2*x
p_prime_r = p_prime.evalf(subs={x:r})
print(f'p\'({r}) {p_prime_r}')

h_prime = sp.simplify(p_prime/t)
print('h\'', h_prime)

h_prime_r = h_prime.evalf(subs={x:r})
print(f'h\'({r}) {h_prime_r}')

# Check that p' is an integer
assert p_prime_r - int(p_prime_r) == 0

# Check that h'(r) is an integer (this should fail)
assert h_prime_r - int(h_prime_r) != 0, f'h\'({r}) = {h_prime_r} is not an integer'


##############################################################################
# Remark 3.1
##############################################################################

# Assume we are using the same commitment t(x) and challenge r. P then
#   creates bogus h' and p' and returns it to V.

# This is supposed to be h(x) = p(x)/t(x) evaluated at challenge r. Instead we
#   just makes something up that will pass the test. 
h_prime_r = np.random.randint(10)

# P generates a forged p' for the attack and returns p'(r) and h'(r) to V.
p_prime_r = t_r * h_prime_r

# Looking back to Example 3.1, we see this will pass.
# Looking back to Example 3.1, we see this will pass.


##############################################################################
# Example 3.3
##############################################################################

# 3.3.a 

# Select field size by hand. This is public.
prime = 17
print('Using prime', prime)

# Select a random base in the field. This is public.
g = np.random.randint(2, prime)
print('Using base g', g)

# Select secret
# using a smaller s for testing
s = np.random.randint(2, prime)
print('Selecting secret s', s)

# Get the coefficients of the polyomial p
p_coeffs = sp.Poly(p).all_coeffs()

# V calculates E(s^0), E(s^1), ..., E(s^p_deg) and gives them to P.
s_encs = [pow(g, s**i, prime) for i in range(len(p_coeffs))]
print('s_encs', s_encs)

# 3.3.b
# Using the same `p` and `t` as used in 3.1
print(f'p {p}')
print(f't {t}')

# P evaluates the polynomial in the ct domain. This is E(p(s)) = g^p.
ps = evaluateEncrypted(p, s_encs, prime)
print('E(p(s))', ps)

# P calculates p(x)/t(x). This is g^h.
h = sp.simplify(p/t)
print('h', h)

# E(h(x)) = g^h
hs = evaluateEncrypted(h, s_encs, prime)
print('E(h(s))', hs)

# P sends g^p and g^h back to V.

# 3.3.c. V checks that (g^h)^t(s) = g^p.
ts = int(t.evalf(subs={x:s}))
assert hs**ts % prime == ps
print('Verification passed')