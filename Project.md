# Project

## Lagrange polynomial

### (a)

Let $f(x)=\dfrac{1}{1+25x^2}$ be a function in $[-1,1]$, and let $p_n(x)$ be the polynomial of degree at most $n$ that interpolates the function $f$ at $n+1$ equally-spaced points $x_0,x_1,\cdots,x_n$ in the interval $[-1,1]$. Thus one gets
$$
h=\frac{2.0}{n}, x_k=-1.0+kh, k=0,1,\cdots,n
$$
.

$p_n(x)$ is an approximation of $f(x)$ at these points. Choose $n=5,10,20$ and compute $p_n(x)$ at
$$
x=-0.95,-0.47,0.1,0.37,0.93
$$
and compare them with $f(x)$ at these points.

### (b)

Let $f(x)=e^x$ be a function in $[-1,1]$ and $p_n(x)$ be the interpolating polynomia of $n+1$ equally-spaced points. Choose $n=5,10,20$ and compute $p_n(x)$ at $x=-0.95,-0.47, 0.1,0.37,0.93$ and compare them with $f(x)$ at these points.

Discribe your observation. Based on your observation, do you think the higher the interpolating polynomial degree, the better?

### (c)

Let $f(x)=\dfrac{1}{25+x^2}$ be a function in $[-1,1]$ and consider the interpolating polynomial with respect to thr distinct points $x_i=\cos{\dfrac{(2k+1)\pi}{2(n+1)}},\, k=0,1,\cdots,n$. Choose $n=5,10,20$ and compare them with $f(x)$ at these points.

Describe the observation. Read reference [1], p.315-323 and explain this phenomenon.

## Numerical Integration

### (a)

To compute the integral $\displaystyle\int_{a}^{b}{f(x)}\mathrm{d}x$ numerically, one can use the Romberg integration formula. To be more specific, let $R(n, 0)$ denote the trapezoid estimate with $2^n$ subintervals, we have the Romberg formula:

$$
\begin{cases}
R(0,0)=\dfrac{1}{2}(b-a)[f(a)+f(b)],\\
\displaystyle{R(n,0)=\dfrac{1}{2}R(n-1,0)+h_n\sum_{i=1}^{2^{n-1}}{f(a+(2i-1)h_n)}}
\end{cases}
$$
where 
$$h_0=b-a, h_n=\dfrac{h_{n-1}}{2}$$

Use Romberg integration to compute the following integrals with the tolerance of accuracy to be $\epsilon=10^{-7}$, which means if $[R(n+1,0)-R(n,0)]\le\epsilon$, then we use $R(n,0)$ as the desired numerical integration.

+ (a) $\displaystyle{\int_{0}^{1}{x^2e^x}\mathrm{d}x}$
+ (b) $\displaystyle{\int_{1}^{3}{e^x\sin{x}}\mathrm{d}x}$
+ (c) $\displaystyle{\int_{0}^{1}{\dfrac{4}{1+x^2}}\mathrm{d}x}$
+ (d) $\displaystyle{\int_{0}^{1}{\dfrac{1}{x+1}}\mathrm{d}x}$

### (b)

Use Gauss-Legendre quadrature rule ( the code to generate the weights and Gauss nodes will be given! ) to compute the above integrals to $\epsilon=10^{-7}$ and compare the number of quadrature points used.

### (c)

Compute $\displaystyle{\int_{-\infty}^{+\infty}}{e^{-x^2}\mathrm{d}x}$ to the tolerance of accuracy $\epsilon =10^{-7}$. ( Hint: although the integral is in unbounded domain, the integrand itself decays fast. )

## References

+ [1] D. Kincaid and W. Cheney. Numerical Analysis: Mathematics of Scientific Computing, 3rd edition. American Mathematical Society, 1996.