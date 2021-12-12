# This Program is distributed under the terms of the MIT license.
# Author: Rui Sun
# Date: 2021-12-12
# Copyright: Rui Sun from SAA, SJTU
# Description: Task 2 - Numerical-Integration
# Environment: Python 3.9.5 64-bit
#              Numpy 1.20.3

import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.polynomial import roots

EPS = 1e-7
MAXSIZE = 200

# Romberg Integration
# \int_a^b{f(x)}\mathrm{d}x
def romberg(f, a, b):
    h = b - a
    R = np.zeros(MAXSIZE)
    R[0] = (h / 2) * (f(a) + f(b))
    for k in range(1, MAXSIZE):
        h /= 2
        R[k] = R[k - 1] / 2
        for j in range(1, 2**(k - 1) + 1):
            R[k] += h * f(a + (2 * j - 1) * h)
        if abs(R[k] - R[k - 1]) < EPS:
            return R[k]
    print("Error: Romberg Integration failed.")
    return 0


# Gauss-Legendre Integration
# \int_a^b{f(x)}\mathrm{d}x
def legendre(n, x):
    if n == 0:
        return 1.0
    elif n == 1:
        return x
    else:
        return ((2.0 * n - 1.0) * x * legendre(n - 1, x) - (n - 1.0) * legendre(n - 2, x)) / n

def dlegendre(n, x):
    if n == 0:
        return 0.0
    elif n == 1:
        return 1.0
    else:
        return (n / (x**2 - 1.0)) * (x * legendre(n, x) - legendre(n - 1, x))

def legendre_roots(polyorder):
    if polyorder < 2:
        print("Error: Polynomial order must be greater than 1.")
        return None
    else:
        roots = []
        for i in range(1, int(polyorder / 2) + 1):
            x = np.cos(np.pi * (i - 0.25) / (polyorder + 0.5))
            error = 10 * EPS**2
            while error > EPS**2:
                dx = -legendre(polyorder, x) / dlegendre(polyorder, x)
                x += dx
                error = abs(dx)
            roots.append(x)
        roots = np.array(roots)
        if polyorder % 2 == 0:
            roots = np.concatenate((-roots, roots[::-1]))
        else:
            roots = np.concatenate((-roots, [0], roots[::-1]))
        return roots

def gauss_legendre_weights(polyorder):
    w = []
    xs = legendre_roots(polyorder)
    w = 2.0 / ((1 - xs**2) * dlegendre(polyorder, xs)**2)
    return w, xs

def gauss_legendre(f, a, b):
    delta = 99999
    I = 0
    n = 1
    while delta > EPS:
        n += 1
        w, xs = gauss_legendre_weights(n)
        I_old = I
        I = (b - a) / 2 * np.sum(w * f((b - a) / 2 * xs + (b + a) / 2))
        delta = abs(I - I_old)
        if n > 1000:
            print("Error: Gauss-Legendre Integration failed.")
            return 0
    return I

def fa(x):
    return (x**2) * np.exp(x)

def fb(x):
    return np.exp(x) * np.sin(x)

def fc(x):
    return 4 / (1 + x**2)

def fd(x):
    return 1 / (1 + x)

# main
print('Romberg Integration')
print('f(x) = x^2 e^x')
print('a = 0, b = 1')
print('I = ', romberg(fa, 0, 1))
print('\n')
print('f(x) = e^x sin(x)')
print('a = 1, b = 3')
print('I = ', romberg(fb, 1, 3))
print('\n')
print('f(x) = 4/(1+x^2)')
print('a = 0, b = 1')
print('I = ', romberg(fc, 0, 1))
print('\n')
print('f(x) = 1/(1+x)')
print('a = 0, b = 1')
print('I = ', romberg(fd, 0, 1))
print('\n')

print('Gauss-Legendre Integration')
print('f(x) = x^2 e^x')
print('a = 0, b = 1')
print('I = ', gauss_legendre(fa, 0, 1))
print('\n')
print('f(x) = e^x sin(x)')
print('a = 1, b = 3')
print('I = ', gauss_legendre(fb, 1, 3))
print('\n')
print('f(x) = 4/(1+x^2)')
print('a = 0, b = 1')
print('I = ', gauss_legendre(fc, 0, 1))
print('\n')
print('f(x) = 1/(1+x)')
print('a = 0, b = 1')
print('I = ', gauss_legendre(fd, 0, 1))
print('\n')