import numpy as np

def f(x):
    return np.tan(x) - x

def phi(x):
    return f(x) + x

def psi(x):
    return compl_polinom(x) + x

def df(x):
    return 1 / (np.cos(x) ** 2) - 1

def newton_next(x):
    return x - f(x)/df(x)

def secants_next(x_prev, x):
    return x - ((x - x_prev) * f(x))/(f(x) - f(x_prev))

def compl_polinom(z):
    return z**3 - 1

def real_equation(a):
    return a**3 - 3*a - 1

def dreal(a):
    return 3*(a**2) - 3

def imaginary_equation(a, b):
    return 3*(a**2) - b**3

# dimaginary == -3 (by b of course

def real_newton_next(a):
    return a - real_equation(a)/dreal(a)

def imaginary_newton_next(a, b):
    return b + imaginary_equation(a, b)/3

# works only if we f(bounds) have different signs:
#   otherwise we can't guarantee there will be no false negative or no infinite cicle,
#   therefore in those cases "no guarantee solutions" will be returned
def divide_by_two(a, b, accuracy):
    print("DIVIDE BY TWO METHOD")
    left = a
    right = b
    left_f = f(a)
    right_f = f(b)

    if (f(a) * f(b) > 0):
        print(f"No guarantee solution on interval [{a}, {b}]\n")
        return

    while(np.abs(left_f) > accuracy) or (np.abs(right_f) > accuracy):
        new_bound = right - (right - left) / 2
        new_bound_f = f(new_bound)
        if (new_bound_f * right_f > 0):
            right = new_bound
            right_f= new_bound_f
        elif (new_bound_f * right_f < 0):
            left = new_bound
            left_f = new_bound_f
        else:
            right = new_bound
            right_f= new_bound_f
            break

    if (np.abs(left_f) <= accuracy):
        print(f"One of solutions on interval [{a}, {b}] is x={left}, f(x)={f(left)}\n")
        return
    elif (np.abs(right_f) <= accuracy):
        print(f"One of solutions on interval [{a}, {b}] is x={right}, f(x)={f(right)}\n")
        return

def simple_iterations(x0, accurancy, max_iterations=1000):
    print("SIMPLE ITERATIONS METHOD")
    _helper(x0, accurancy, max_iterations, 0)

def newton(x0, accurancy, max_iterations=1000):
    print("NEWTON METHOD")
    _helper(x0, accurancy, max_iterations, 1)

def secants(x0, x1, accurancy, max_iterations=1000):
    print("SECANTS METHOD")
    x_prev = x0
    x = x1
    for _ in range(max_iterations):
        x_next = secants_next(x_prev, x)
        if f(x_next) <= accurancy:
            print(f"One of solutions with accurancy={accurancy} is x={x_next}, f(x)={f(x_next)}\n")
            return
        x_prev = x
        x = x_next
    print(f"Not enough iterations for accurancy={accurancy} or there is no convergence, last x={x}, f(x)={f(x)}\n")
    return


# 0 for si, 1 for newton
def _helper(x0, accurancy, max_iterations, indicator):
    x = x0
    next_f = f(x)
    for _ in range(max_iterations):
        if indicator == 0:
            x_next = phi(x)
            next_f = f(x_next)
        elif indicator == 1:
            x_next = newton_next(x)
            next_f = f(x_next)
        if np.abs(next_f) <= accurancy:
            print(f"One of solutions with accurancy={accurancy} is x={x_next}, f(x)={next_f}\n")
            return
        x = x_next
    print(f"Not enough iterations for accurancy={accurancy} or there is no convergence, last x={x}, f(x)={next_f}\n")
    return

# As long as we want to get compl solution, let it be z = a + i*b
def newton_for_polinom(a0, b0, accurancy_a, accurancy_b, max_iterations=100000):
    print("NEWTON FOR POLINOM z**3 - 1 = 0")
    success_a = False
    success_b = False
    a = a0
    b = b0

    for _ in range(max_iterations):
        a = real_newton_next(a)
        if np.abs(a) <= accurancy_a:
            success_a = True
            break
    for _ in range(max_iterations):
        b = imaginary_newton_next(a, b)
        if np.abs(b) <= accurancy_b:
            success_b = True
            break
    if (success_a and success_b):
        print(f"One of solutions is z={a}+i{b}\n")
    else:
        print(f"Not enough iterations, last z={a}+i{b}\n")



divide_by_two(-2, 1, 0.001)
simple_iterations(1.5,0.001)
newton(1.5, 0.001, 10)
secants(1, 1.5, 0.001, 10)
newton_for_polinom(-0.35, 0.75, 0.001, 0.001)