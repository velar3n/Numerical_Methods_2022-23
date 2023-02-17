import cmath
import numpy as np

polynomial1 = np.poly1d([243, -486, 783, -990, 558, -28, -72, 16])      # kolejne wielomiany z zadania
polynomial2 = np.poly1d([1, 0 + 1j, -1, 0 - 1j, 1])
polynomial3 = np.poly1d([1, 1, 3, 2, -1, -3, -11, -8, -12, -4, -4])

def new_polynomial(f_n, x0):     # deflacja / obniżenie stopnia wielomianu
    n = f_n.order
    zeroarr = np.full(n - 1, -x0, dtype=complex)
    coefs_fn = f_n.coef
    coefs_fn = np.delete(coefs_fn, -1)
    A = np.diag(np.ones((n,), dtype=int)) + np.diag(zeroarr, -1)
    return np.linalg.solve(A, coefs_fn)

def laguerre(polynomial, x0):
    n = polynomial.order                 # stopień wielomianu
    first_der = np.polyder(polynomial)   # pierwsza pochodna
    second_der = np.polyder(first_der)   # druga pochodna
    xk = x0
    while abs(polynomial(xk)) > 0.00000001:     # iterujemy, aż różnica między miejscami zerowymi będzie mniejsza niż zadany błąd
        l = n * polynomial(xk)
        root = cmath.sqrt((n - 1) * ((n - 1) * first_der(xk) ** 2 - n * polynomial(xk) * second_der(xk)))
        m = max([first_der(xk) + root, first_der(xk) - root], key=abs)
        a = l / m
        xk -= a
    return xk

def calculate_roots(polynomial):     # obliczanie miejsc zerowych
    print(polynomial)
    x_roots = np.array([])
    x = laguerre(polynomial, 0)
    x_roots = np.append(x_roots, x)
    order = polynomial.order

    f_p = polynomial
    while (order > 1):      # obliczamy dopóki wielomian nie zostanie obniżony do stopnia pierszego
        f = np.poly1d(new_polynomial(f_p, x))
        x_ = laguerre(f, 0)
        x = laguerre(polynomial, x_)
        x_roots = np.append(x_roots, x)
        f_p = f
        order -= 1

    for i in x_roots:       # wypisanie wyników
        print(i)
    print("\n")


calculate_roots(polynomial1)
calculate_roots(polynomial2)
calculate_roots(polynomial3)
