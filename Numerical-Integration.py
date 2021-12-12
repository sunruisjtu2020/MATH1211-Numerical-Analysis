# This Program is distributed under the terms of the MIT license.
# Author: Rui Sun
# Date: 2021-12-12
# Copyright: Rui Sun from SAA, SJTU
# Description: Task 2 - Numerical-Integration
# Environment: Python 3.9.5 64-bit
#              Numpy 1.20.3
#              Matplotlib 3.4.2

import numpy as np
import matplotlib.pyplot as plt

EPS = 1e-7
MAXSIZE = 200

# Romberg Integration
# \int_a^b{f(x)}\mathrm{d}x
# FIXME: Error Results
def romberg(f, a, b):
    n = 1
    I = 0
    delta = 9999999
    while delta > EPS:
        h = (b - a) / n
        R = np.zeros((n, n))
        R[0, 0] = (f(a) + f(b)) * h / 2
        for i in range(1, n):
            for j in range(i + 1):
                R[i, j] = R[i - 1, j] + (R[i - 1, j] - R[i - 1, j - 1]) / (4 ** j - 1)
        delta = np.abs(R[n - 1, n - 1] - R[n - 2, n - 2])
        I = R[n - 1, n - 1]
        n += 1
        if n > MAXSIZE:
            print("Error: Romberg Integration failed.")
            return 0
    return I

# Gauss-Legendre Integration
# \int_a^b{f(x)}\mathrm{d}x
# FIXME: No Results Output
def gauss_legendre(f, a, b):
    n = 1
    I = 0
    delta = 9999999
    while delta > EPS:
        h = (b - a) / n
        x = np.zeros(n)
        w = np.zeros(n)
        for i in range(n):
            x[i] = a + (i + 0.5) * h
            w[i] = h / 2
        I_ = 0
        for i in range(n):
            I_ += w[i] * f(x[i])
        delta = np.abs(I - I_)
        n += 1
        if n > MAXSIZE:
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