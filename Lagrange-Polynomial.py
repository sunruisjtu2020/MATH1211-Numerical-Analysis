# This Program is distributed under the terms of the MIT license.
# Author: Rui Sun
# Date: 2021-12-12
# Copyright: Rui Sun from SAA, SJTU
# Description: Task 1 - Lagrange-Polynomial
# Environment: Python 3.9.5 64-bit
#              Numpy 1.20.3
#              Matplotlib 3.4.2

import numpy as np
import matplotlib.pyplot as plt

# Lagrange-Polynomial
def L(n, x, f):
    x_ = np.linspace(-1, 1, n + 1)
    y_ = f(x_)
    l_ = np.zeros(x_.size)
    res = 0.0
    for i in range(n + 1):
        deno, nume = 1.0, 1.0
        for j in range(n + 1):
            if i != j:
                deno *= (x_[i] - x_[j])
                nume *= (x - x_[j])
        l_[i] = nume / deno
        res += l_[i] * y_[i]
    return res

def LLL(n, x, f):
    res = np.zeros(200)
    for i in range(200):
        res[i] = L(n, x[i], f)
    return res

def redefinedL(n, x, f):
    x_ = np.zeros(n + 1)
    for i in range(n + 1):
        x_[i] = np.cos((2 * i + 1) * np.pi / (2 * (n + 1)))
    y_ = f(x_)
    l_ = np.zeros(x_.size)
    res = 0.0
    for i in range(n + 1):
        deno, nume = 1.0, 1.0
        for j in range(n + 1):
            if i != j:
                deno *= (x_[i] - x_[j])
                nume *= (x - x_[j])
        l_[i] = nume / deno
        res += l_[i] * y_[i]
    return res

def redefinedLLL(n, x, f):
    res = np.zeros(200)
    for i in range(200):
        res[i] = redefinedL(n, x[i], f)
    return res

def f(x):
    return 1.0 / (1.0 + 25.0 * x**2)

def g(x):
    return np.exp(x)

# main
n = np.array([5, 10, 20])
x = np.array([-0.95, -0.47, 0.1, 0.37, 0.93])

# function f(x) = 1 / (1 + 25x^2)
for i in n:
    print(f'n = {i}:')
    for j in x:
        print(f'f({j})-p_n({j}) = {f(j) - L(i, j, f)}')
    plt.title(r'Lagrange-Polynomial of $f(x)=\dfrac{1}{25+x^2}$')
    plt.grid('on')
    xx = np.linspace(-1, 1, 200)
    yy = f(xx)
    zz = LLL(i, xx, f)
    plt.plot(xx, yy, 'r-', label='$f(x)=\dfrac{1}{25+x^2}$')
    plt.plot(xx, zz, 'g:', label=f'$p{i}(x)$')
    plt.legend()
    plt.show()

# function g(x) = e^x
for i in n:
    print(f'n = {i}:')
    for j in x:
        print(f'g({j})-p_n({j}) = {g(j) - L(i, j, g)}')
    plt.title(r'Lagrange-Polynomial of $g(x)=e^x$')
    plt.grid('on')
    xx = np.linspace(-1, 1, 200)
    yy = g(xx)
    zz = LLL(i, xx, g)
    plt.plot(xx, yy, 'r-', label='$g(x)=e^x$')
    plt.plot(xx, zz, 'g:', label=f'$p{i}(x)$')
    plt.legend()
    plt.show()


# Redefined Lagrange-Polynomial
# function f(x) = 1 / (1 + 25x^2)
for i in n:
    print(f'n = {i}:')
    for j in x:
        print(f'f({j})-p_n({j}) = {f(j) - redefinedL(i, j, f)}')
    plt.title(r'Lagrange-Polynomial of $f(x)=\dfrac{1}{25+x^2}$')
    plt.grid('on')
    xx = np.linspace(-1, 1, 200)
    yy = f(xx)
    zz = redefinedLLL(i, xx, f)
    plt.plot(xx, yy, 'r-', label='$f(x)=\dfrac{1}{25+x^2}$')
    plt.plot(xx, zz, 'g:', label=f'$p{i}(x)$')
    plt.legend()
    plt.show()