import numpy as np
import matplotlib.pyplot as plt


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

# z0 = 1
# z12 = -0.5 +- i(0.75)^0.5
def compl_polinom(z):
    return z**3 - 1


def real_equation(a):
    return a**3 - 3*a - 1

def dreal(a):
    return 3*(a**2) - 3

def imaginary_equation(a, b):
    return 3*(a**2) - b**3

# dimaginary == -3 (by b of course

def newton_next(z):
    # epsilon = 1e-8
    # max_attempts = 10
    # derivative = dreal(a)
    #
    # attempt = 0
    # while abs(derivative) < 1e-12 and attempt < max_attempts:
    #     print(f"Function real newton next: derivative close to 0, slide point. Attempt {attempt}")
    #     a += epsilon
    #     derivative = dreal(a)
    #     attempt += 1
    #
    # if abs(derivative) < 1e-12:
    #     return a
    return (2 * z**3 + 1) / (3 * z**2)

# def imaginary_newton_next(a, b):
#     return b + imaginary_equation(a, b)/3

# works only if f(bounds) have different signs:
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
def newton_for_polinom(starts, accurancy=0.001, max_iterations=100):
    print("\nNEWTON FOR POLINOM z**3 - 1 = 0")

    found = []

    for z0 in starts:
        z = z0
        for i in range(max_iterations):
            if np.abs(z) < 1e-12:
                print(f"Warning: Reached zero at point {z}. Skipping.")
                break
            z_next = newton_next(z)
            if np.abs(z_next - z) <= accurancy:
                found.append(z_next)
                break
            if abs(z_next) > 1e6:
                print(f"Diverged (too large value) from start {z0} on {i} step.")
                break
            z = z_next


    unique = []
    for z in found:
        if not any(np.abs(z - w) < accurancy for w in unique):
            unique.append(z)
    print(f"Found solutions: {unique}\n")


def find_all_roots():
    roots = [
        1,
        np.exp(2j * np.pi / 3),
        np.exp(4j * np.pi / 3)
    ]
    print("All solutions for z^3 = 1:")
    for i, root in enumerate(roots):
        print(f"  z_{i} = {root}")
    return roots


def plot_basins(xmin=-2, xmax=2, ymin=-2, ymax=2, resolution=5000, max_iter=50):
    print(f"\nBuilding up Newton basins... \n")

    x_vals = np.linspace(xmin, xmax, resolution)
    y_vals = np.linspace(ymin, ymax, resolution)

    X, Y = np.meshgrid(x_vals, y_vals)
    Z = X + 1j * Y

    roots = [
        (1, 0),
        (np.exp(2j * np.pi / 3), 1),
        (np.exp(4j * np.pi / 3), 2)
    ]

    basins_indices = np.zeros((resolution, resolution, 3))

    colors = {
        0: [1.0, 0.0, 0.0],
        1: [0.0, 1.0, 0.0],
        2: [0.0, 0.0, 1.0]
    }

    for i in range(resolution):
        for j in range(resolution):
            z = Z[i, j]

            if abs(z) < 1e-10:
                continue

            basins_indices[i, j] = [0.0, 0.0, 0.0]

            for _ in range(max_iter):
                z_next = newton_next(z)

                if abs(z_next) > 1e6:
                    break


                for root_val, idx in roots:
                    if abs(z_next - root_val) < 1e-3:
                        basins_indices[i, j] = colors[idx]
                        break
                else:
                    z = z_next
                    continue
                break

    plt.figure(figsize=(10, 10))

    #cmap = plt.cm.get_cmap('hsv', len(roots))

    plt.imshow(basins_indices, extent=(xmin, xmax, ymin, ymax), origin='lower')

    plt.scatter([1, -0.5, -0.5], [0, np.sqrt(3) / 2, -np.sqrt(3) / 2], color='black', s=150, marker='x', linewidths=2)

    plt.title("Newton basins of attraction for $z^3 - 1 = 0$")
    plt.xlabel("Re(z)")
    plt.ylabel("Im(z)")

    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

starts = [
    complex(1.0, 0.0),
    complex(-0.5, 0.9),
    complex(-0.5, -0.9)
]

#divide_by_two(-2, 2, 0.001)
#simple_iterations(1.5,0.001)
#newton(1.5, 0.001, 1000)
#secants(1, 1.5, 0.001, 10)
find_all_roots()
newton_for_polinom(starts)
plot_basins()