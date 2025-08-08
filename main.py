import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
beta = 0.1
g = 0.1
alpha = 0.5
l = 0.9
Delta = 0.3
sigma = 1
k = 1
t0 = 20
tau_min = -t0
tau_max = t0

def xi_func(tau, OmegaA):
    return Delta * np.tanh(tau) + OmegaA

def beta_A_func(xi_val):
    return alpha + beta - 1 - xi_val**2

def beta_Z_func(xi_val):
    numerator = 2 * alpha**2 * (xi_val**4 + 2*g*xi_val**3 + 2*g**2*xi_val**2 - 5*xi_val**2 - 6*g*xi_val + 3)
    denom1 = xi_val**2 - 1
    denom2 = xi_val**4 - 6*xi_val**2 - 4*g*xi_val + 3
    denom = denom1 * denom2
    if np.abs(denom) < 1e-12:
        return np.nan
    return beta + 2*alpha + numerator / denom

def compute_P(tau, OmegaA):
    xi_val = xi_func(tau, OmegaA)
    beta_A_val = beta_A_func(xi_val)
    beta_Z_val = beta_Z_func(xi_val)
    if np.isnan(beta_Z_val):
        return np.nan
    denom = (1 - l) * beta_Z_val + l * beta_A_val
    if np.abs(denom) < 1e-12:
        return np.nan
    P_val = (beta_A_val * beta_Z_val) / denom
    return P_val

def compute_Q(tau, OmegaA):
    xi_val = xi_func(tau, OmegaA)
    beta_A_val = beta_A_func(xi_val)
    return (k**2 / sigma**2) * beta_A_val

def compute_asymptotic_params(OmegaA, sign):
    if sign == -1:
        xi_asympt = OmegaA - Delta
    else:
        xi_asympt = OmegaA + Delta
    beta_A_asympt = beta_A_func(xi_asympt)
    if np.isnan(beta_A_asympt):
        return np.nan, np.nan
    beta_Z_asympt = beta_Z_func(xi_asympt)
    if np.isnan(beta_Z_asympt):
        return np.nan, np.nan
    denom = (1 - l) * beta_Z_asympt + l * beta_A_asympt
    if np.abs(denom) < 1e-12:
        return np.nan, np.nan
    P_asympt = (beta_A_asympt * beta_Z_asympt) / denom
    Q_asympt = (k**2 / sigma**2) * beta_A_asympt
    return P_asympt, Q_asympt

def compute_lambda(OmegaA, sign):
    P_asympt, Q_asympt = compute_asymptotic_params(OmegaA, sign)
    if np.isnan(P_asympt) or np.isnan(Q_asympt) or np.abs(P_asympt) < 1e-12:
        return np.nan
    ratio = Q_asympt / P_asympt
    lam = np.sqrt(ratio)
    if np.real(lam) <= 0:
        return np.nan
    return lam

def compute_error_for_OmegaA(OmegaA, tau_min, tau_max):
    lambda_minus = compute_lambda(OmegaA, -1)
    if np.isnan(lambda_minus):
        return np.nan
    lambda_plus = compute_lambda(OmegaA, 1)
    if np.isnan(lambda_plus):
        return np.nan
    P_min_exact = compute_P(tau_min, OmegaA)
    if np.isnan(P_min_exact):
        return np.nan
    P_max_exact = compute_P(tau_max, OmegaA)
    if np.isnan(P_max_exact):
        return np.nan
    y1_0 = 1.0 + 0.0j
    y2_0 = P_min_exact * lambda_minus
    Y0 = np.array([y1_0, y2_0], dtype=np.complex128)
    def ode_system(tau, Y):
        y1, y2 = Y
        P_val = compute_P(tau, OmegaA)
        if np.isnan(P_val):
            return [0.0, 0.0]
        Q_val = compute_Q(tau, OmegaA)
        if np.isnan(Q_val):
            return [0.0, 0.0]
        dy1dtau = y2 / P_val
        dy2dtau = Q_val * y1
        return [dy1dtau, dy2dtau]
    sol = solve_ivp(ode_system, [tau_min, tau_max], Y0, method='RK45', rtol=1e-6, atol=1e-8)
    if not sol.success:
        return np.nan
    y1_end = sol.y[0, -1]
    y2_end = sol.y[1, -1]
    residual = y2_end + P_max_exact * lambda_plus * y1_end
    return np.abs(residual)

def find_and_plot_solution(OmegaA_opt, tau_min, tau_max):
    lambda_minus = compute_lambda(OmegaA_opt, -1)
    P_min_exact = compute_P(tau_min, OmegaA_opt)
    y1_0 = 1.0 + 0.0j
    y2_0 = P_min_exact * lambda_minus
    Y0 = np.array([y1_0, y2_0], dtype=np.complex128)
    def ode_system(tau, Y):
        y1, y2 = Y
        P_val = compute_P(tau, OmegaA_opt)
        Q_val = compute_Q(tau, OmegaA_opt)
        dy1dtau = y2 / P_val
        dy2dtau = Q_val * y1
        return [dy1dtau, dy2dtau]
    t_eval = np.linspace(tau_min, tau_max, 1000)
    sol = solve_ivp(ode_system, [tau_min, tau_max], Y0, t_eval=t_eval, method='RK45', rtol=1e-6, atol=1e-8)
    y_sol = sol.y[0]
    plt.figure(figsize=(10, 6))
    plt.plot(sol.t, np.real(y_sol), label='Real part')
    plt.plot(sol.t, np.imag(y_sol), label='Imaginary part')
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$y(\tau)$')
    plt.title(f'Solution for $\Omega_A = {OmegaA_opt}$')
    plt.legend()
    plt.grid(True)
    plt.show()

# Grid setup
n_real = 501
n_imag = 501
real_grid = np.linspace(-1, 1, n_real)
imag_grid = np.linspace(0, 1, n_imag)
error_grid = np.full((n_real, n_imag), np.nan)

# Grid scan
for i, re in enumerate(real_grid):
    for j, im in enumerate(imag_grid):
        OmegaA = re + 1j * im
        error = compute_error_for_OmegaA(OmegaA, tau_min, tau_max)
        if not np.isnan(error):
            error_grid[i, j] = error

# Find optimal OmegaA
min_error = np.inf
optimal_OmegaA = None
for i, re in enumerate(real_grid):
    for j, im in enumerate(imag_grid):
        if not np.isnan(error_grid[i, j]) and error_grid[i, j] < min_error:
            min_error = error_grid[i, j]
            optimal_OmegaA = real_grid[i] + 1j * imag_grid[j]

# Plot error distribution
plt.figure(figsize=(10, 8))
extent = [imag_grid[0], imag_grid[-1], real_grid[0], real_grid[-1]]
plt.imshow(np.log10(error_grid.T), extent=extent, origin='lower', aspect='auto', cmap='viridis')
plt.colorbar(label='log10(error)')
plt.xlabel('Imag($\Omega_A$)')
plt.ylabel('Real($\Omega_A$)')
plt.title('Error Distribution over Complex $\Omega_A$')
plt.show()

# Plot solution for optimal OmegaA
if optimal_OmegaA is not None:
    print(f"Optimal OmegaA: {optimal_OmegaA}, Error: {min_error}")
    find_and_plot_solution(optimal_OmegaA, tau_min, tau_max)
else:
    print("No valid OmegaA found.")
