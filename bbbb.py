import math
import cmath

# 1. Bisection Method
def bisection(f, a, b, eps=1e-6, Nmax=100):
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        return "Error: f(a) and f(b) must have opposite signs"

    for _ in range(Nmax):
        c = (a + b) / 2
        fc = f(c)

        if abs(fc) < eps or (b - a) / 2 < eps:
            return c

        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc

    return "Failed to converge"

# 2. Fixed-Point Method
def fixed_point(g, x0, eps=1e-6, Nmax=100):
    x = x0
    for _ in range(Nmax):
        x_new = g(x)
        if abs(x_new - x) < eps:
            return x_new
        x = x_new
    return "Failed to converge"

# 3. Newton–Raphson Method
def newton_raphson(f, df, x0, eps=1e-6, Nmax=100):
    x = x0
    for _ in range(Nmax):
        fx = f(x)
        dfx = df(x)

        if dfx == 0:
            return "Error: derivative is zero"

        x_new = x - fx / dfx
        if abs(x_new - x) < eps:
            return x_new
        x = x_new

    return "Failed to converge"

# 4. Secant Method
def secant(f, x0, x1, eps=1e-6, Nmax=100):
    f0, f1 = f(x0), f(x1)

    for _ in range(Nmax):
        if f1 - f0 == 0:
            return "Error: division by zero"

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if abs(f(x2)) < eps:
            return x2

        x0, x1 = x1, x2
        f0, f1 = f1, f(x1)

    return "Failed to converge"

# 5. False Position Method
def false_position(f, a, b, eps=1e-6, Nmax=100):
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        return "Error: f(a) and f(b) must have opposite signs"

    for _ in range(Nmax):
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)

        if abs(fc) < eps:
            return c

        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc

    return "Failed to converge"

# 6. Muller's Method
def muller(f, x0, x1, x2, eps=1e-6, Nmax=100):
    for _ in range(Nmax):
        y0, y1, y2 = f(x0), f(x1), f(x2)

        h0 = x1 - x0
        h1 = x2 - x1
        d0 = (y1 - y0) / h0
        d1 = (y2 - y1) / h1

        a = (d1 - d0) / (h1 + h0)
        b = a * h1 + d1
        c = y2

        D = cmath.sqrt(b**2 - 4*a*c)
        denom = b + D if abs(b + D) > abs(b - D) else b - D

        if denom == 0:
            return "Error: division by zero"

        x_new = x2 - (2 * c) / denom

        if abs(f(x_new)) < eps:
            return x_new

        x0, x1, x2 = x1, x2, x_new

    return "Failed to converge"

# Example Test
if __name__ == "__main__":
    f  = lambda x: x**3 - 2*x - 5
    df = lambda x: 3*x**2 - 2
    g  = lambda x: math.cos(x)

    print("Bisection:", bisection(f, 2, 3))
    print("Fixed-point:", fixed_point(g, 0.5))
    print("Newton-Raphson:", newton_raphson(f, df, 2))
    print("Secant:", secant(f, 2, 3))
    print("False Position:", false_position(f, 2, 3))
    print("Muller:", muller(f, 2, 2.5, 3))
