import matplotlib.pyplot as plt
import numpy as np
import qade


x = np.linspace(0, 1, 100)
t = np.linspace(0, 1, 100)
phi = qade.function(n_in=2, n_out=1)  # The function to be solved for

# The wave equation (evaluated at the "bulk" (x, t) sample points)
xt = np.array([[x_elem, t_elem] for x_elem in x for t_elem in t])
eq = qade.equation(phi[2, 0] - phi[0, 2], xt)

# Initial conditions for phi and dphi/dt at t = 0
xt_init = np.array([[x_elem, 0] for x_elem in x])
ic_1 = qade.equation(phi[0, 0] - np.sin(2 * np.pi * x), xt_init)
ic_2 = qade.equation(phi[0, 1] - np.pi * np.cos(2 * np.pi * x), xt_init)

# Boundary conditions for phi at x = 0 and x = 1
xt_left = np.array([[0, t_elem] for t_elem in t])
bc_left = qade.equation(phi[0, 0] - np.sin(2 * np.pi * t), xt_left)
xt_right = np.array([[1, t_elem] for t_elem in t])
bc_right = qade.equation(phi[0, 0] - np.sin(2 * np.pi * t), xt_right)

# Solve the equation using a basis of Fourier functions
phi_sol = qade.solve(
    [eq, ic_1, ic_2, bc_left, bc_right], qade.basis("fourier", 3), n_spins=2
)
print(f"loss = {phi_sol.loss:.3}, weights = {np.around(phi_sol.weights, 2)}")


def phi_true(xt):  # Analytical solution, for comparison purposes
    x, t = xt[:, 0], xt[:, 1]
    tau, xi = 2 * np.pi * t, 2 * np.pi * x
    return np.cos(tau) * np.sin(xi) + 0.5 * np.sin(tau) * np.cos(xi)


# Show the results
for t_slice in [0, 0.3, 0.6, 0.9]:
    xt_slice = np.array([[x_elem, t_slice] for x_elem in x])
    plt.plot(x, phi_sol(xt_slice), label=f"t = {t_slice}", linewidth=5)
    plt.plot(x, phi_true(xt_slice), linestyle="dashed", color="black")

plt.legend()
plt.xlabel("$x$")
plt.ylabel("$\\phi(x, t)$")
plt.show()
