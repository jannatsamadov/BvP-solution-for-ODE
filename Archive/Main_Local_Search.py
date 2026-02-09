#hələlik əsas kod budur. main code
#bu kodda metod əlavəsi də var
#bir iki əlavə etməyi düşünürəm, ikinci sərhədlə ədədi hesabın fərqinin aşağı olduğu oblasstlarda mesh götürüb yenidən hesablamaq.


import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sympy import symbols, tanh, diff, simplify, latex, Function

# Create plots directory
os.makedirs("plots", exist_ok=True)

# Global constants from the problem
beta = 0.1
g = 0.1
alpha = 0.5
l = 0.9
Delta = 0.3
sigma = 1.0
k = 1.0
t0 = 20.0  # Approximate infinity for tau range, e.g., tau from -t0 to +t0

# Def 1: xi
def xi(tau, Omega_A):
    return Delta * np.tanh(tau) + Omega_A

# Def 2: beta_A
def beta_A(tau, Omega_A):
    xi_val = xi(tau, Omega_A)
    return alpha + beta - 1 - xi_val**2

# Def 3: beta_Z
def beta_Z(tau, Omega_A):
    xi_val = xi(tau, Omega_A)
    num = 2 * alpha**2 * (xi_val**4 + 2 * g * xi_val**3 + 2 * g**2 * xi_val**2 - 5 * xi_val**2 - 6 * g * xi_val + 3)
    den = (xi_val**2 - 1) * (xi_val**4 - 6 * xi_val**2 - 4 * g * xi_val + 3)
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12  # Safeguard
    return beta + 2 * alpha + num / den

# Def 4: P
def P(tau, Omega_A):
    bA = beta_A(tau, Omega_A)
    bZ = beta_Z(tau, Omega_A)
    return bA * bZ / ((1 - l) * bZ + l * bA)

# Def 5: Q
def Q(tau, Omega_A):
    return k**2 * beta_A(tau, Omega_A)

# Def 6: xi_minus
def xi_minus(Omega_A):
    return Omega_A - Delta

# Def 7: xi_plus
def xi_plus(Omega_A):
    return Omega_A + Delta

# Def 8: beta_A_limit
def beta_A_limit(xi_lim):
    return alpha + beta - 1 - xi_lim**2

# Def 9: beta_Z_limit
def beta_Z_limit(xi_lim):
    num = 2 * alpha**2 * (xi_lim**4 + 2 * g * xi_lim**3 + 2 * g**2 * xi_lim**2 - 5 * xi_lim**2 - 6 * g * xi_lim + 3)
    den = (xi_lim**2 - 1) * (xi_lim**4 - 6 * xi_lim**2 - 4 * g * xi_lim + 3)
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12
    return beta + 2 * alpha + num / den

# Def 10: P_limit
def P_limit(xi_lim):
    bA = beta_A_limit(xi_lim)
    bZ = beta_Z_limit(xi_lim)
    den = (1 - l) * bZ + l * bA
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12
    return bA * bZ / den

# Def 11: Q_limit
def Q_limit(xi_lim):
    return k**2 * beta_A_limit(xi_lim)

# Def 12: lambda_limit
def lambda_limit(xi_lim):
    p_lim = P_limit(xi_lim)
    if np.abs(p_lim) < 1e-12:
        p_lim = 1e-12 + 1j * 1e-12
    return np.sqrt(Q_limit(xi_lim) / p_lim)

# Def 13: lambda_minus
def lambda_minus(Omega_A):
    return lambda_limit(Omega_A - Delta)

# Def 14: lambda_plus
def lambda_plus(Omega_A):
    return lambda_limit(Omega_A + Delta)

# Def 15: dxi_dtau
def dxi_dtau(tau):
    return Delta * (1 - np.tanh(tau)**2)

# Def 16: dbeta_A_dtau
def dbeta_A_dtau(tau, Omega_A):
    return -2 * xi(tau, Omega_A) * dxi_dtau(tau)

# Def 17: dbeta_Z_dtau
def dbeta_Z_dtau(tau, Omega_A):
    xi_val = xi(tau, Omega_A)
    dxi = dxi_dtau(tau)
    num = 2 * alpha**2 * (xi_val**4 + 2*g*xi_val**3 + 2*g**2*xi_val**2 -5*xi_val**2 -6*g*xi_val +3)
    dnum = 2 * alpha**2 * dxi * (4*xi_val**3 + 6*g*xi_val**2 + 4*g**2*xi_val -10*xi_val -6*g)
    a = xi_val**2 - 1
    da = 2 * xi_val * dxi
    b = xi_val**4 -6*xi_val**2 -4*g*xi_val +3
    db = dxi * (4*xi_val**3 -12*xi_val -4*g)
    den = a * b
    dden = da * b + a * db
    if np.abs(den) < 1e-12:
        return 0.0 + 0j
    return (dnum / den) - (num * dden / den**2)

# Def 18: dP_dtau
def dP_dtau(tau, Omega_A):
    bA = beta_A(tau, Omega_A)
    bZ = beta_Z(tau, Omega_A)
    dbA = dbeta_A_dtau(tau, Omega_A)
    dbZ = dbeta_Z_dtau(tau, Omega_A)
    num = bA * bZ
    dnum = dbA * bZ + bA * dbZ
    den = (1 - l) * bZ + l * bA
    dden = (1 - l) * dbZ + l * dbA
    if np.abs(den) < 1e-12:
        return 0.0 + 0j
    return (dnum * den - num * dden) / den**2

# Def 19: ode_fun_real
def ode_fun_real(tau, z, Omega_A):
    yr, yi, ypr, ypi = z
    y = yr + 1j * yi
    yp = ypr + 1j * ypi
    p = P(tau, Omega_A)
    dp = dP_dtau(tau, Omega_A)
    q = Q(tau, Omega_A)
    if np.abs(p) < 1e-12:
        p = 1e-12 + 1j * 1e-12
    ypp = (q * y - dp * yp) / p
    return np.array([ypr, ypi, np.real(ypp), np.imag(ypp)])

# Execution Controls
methods = ['RK45', 'DOP853']  # istəsən 'BS5' və ya 'Radau' da əlavə edə bilərsən
top_results_overall = []

def process_grid(LReal, LImag, method):
    mismatch_grid = np.full_like(LReal, np.nan, dtype=float)
    grid_results = []
    
    for i in range(LReal.shape[0]):
        for j in range(LReal.shape[1]):
            Omega_A = LReal[i, j] + 1j * LImag[i, j]
            
            lam_m = lambda_minus(Omega_A)
            z0 = [1.0, 0.0, np.real(lam_m), np.imag(lam_m)]
            
            if not np.all(np.isfinite(z0)):
                continue
                
            tau_span = (-t0, t0)
            sol = solve_ivp(ode_fun_real, tau_span, z0, args=(Omega_A,), method=method,
                            rtol=1e-8, atol=1e-8, dense_output=True)
            
            if not sol.success:
                continue
            
            z_right = sol.sol(t0)
            yr, yi, ypr, ypi = z_right
            y_right = yr + 1j * yi
            yp_right = ypr + 1j * ypi
            
            lam_p = lambda_plus(Omega_A)
            mismatch = yp_right + lam_p * y_right
            abs_mismatch = np.abs(mismatch)
            
            mismatch_grid[i, j] = abs_mismatch
            grid_results.append((Omega_A, abs_mismatch, sol, method))
            
    return mismatch_grid, grid_results

def plot_heatmap(Real_Grid, Imag_Grid, Mismatch_Grid, title, filename):
    plt.figure(figsize=(10, 6))
    plt.pcolormesh(Real_Grid, Imag_Grid, np.log10(Mismatch_Grid + 1e-10), shading='auto', cmap='viridis_r')
    plt.colorbar(label='log10(|mismatch|)')
    plt.xlabel('Re(Omega_A)')
    plt.ylabel('Im(Omega_A)')
    plt.title(title)
    plt.savefig(filename)
    plt.close()

# -----------------
# TWO-LEVEL SEARCH
# -----------------

dt_real = 2.0 / 19  # Coarse step width
dt_imag = 1.0 / 19  # Coarse step width

for method in methods:
    print(f"Running global search for method {method}...")
    
    # Grid as specified (Coarse)
    real_parts = np.linspace(-1, 1, 20)
    imag_parts = np.linspace(0, 1, 20)
    Real, Imag = np.meshgrid(real_parts, imag_parts)
    
    global_mismatch_grid, global_results = process_grid(Real, Imag, method)
    
    # Plot global heatmap
    plot_heatmap(Real, Imag, global_mismatch_grid, f'Global Heatmap of Mismatch ({method})', f'plots/Global_Heatmap_{method}.png')
    
    # Sort and get top 5 global
    global_results.sort(key=lambda x: x[1])
    top_5_global = global_results[:5]
    
    print(f"Top 5 global points for {method}:")
    for omega, mm, _, _ in top_5_global:
        print(f"  Omega_A = {omega:.4f}, |mismatch| = {mm:.4e}")
        
    # Local Refinement phase
    for r_idx, (global_omega, global_mm, _, _) in enumerate(top_5_global):
        print(f"  Running local search around {global_omega:.4f} for {method}...")
        
        # Build local fine grid
        local_real_parts = np.linspace(global_omega.real - dt_real, global_omega.real + dt_real, 20)
        local_imag_parts = np.linspace(global_omega.imag - dt_imag, global_omega.imag + dt_imag, 20)
        LReal, LImag = np.meshgrid(local_real_parts, local_imag_parts)
        
        local_mismatch_grid, local_results = process_grid(LReal, LImag, method)
        
        # Plot local heatmap
        plot_heatmap(LReal, LImag, local_mismatch_grid, f'Local Heatmap around {global_omega:.4f} ({method})', f'plots/Local_Heatmap_{method}_{r_idx+1}.png')
        
        # Find local best to keep for final plotting
        if local_results:
            local_results.sort(key=lambda x: x[1])
            best_local = local_results[0]
            top_results_overall.append(best_local)
            print(f"    Best local: {best_local[0]:.4f}, |mismatch| = {best_local[1]:.4e}")

# -----------------
# PLOTS AND ERRORS
# -----------------

print("\nFinalizing and plotting overall absolute best trajectory maps...")
plots = []
errors = []

for omega, mm, sol, method in top_results_overall:
    tau_eval = np.linspace(-t0, t0, 500)
    z_eval = sol.sol(tau_eval)
    y_eval = z_eval[0] + 1j * z_eval[1]
    plots.append((omega, tau_eval, y_eval, method))
    
    # Residual computation
    residual = []
    for tt in tau_eval:
        z_tt = sol.sol(tt)
        y_tt = z_tt[0] + 1j * z_tt[1]
        yp_tt = z_tt[2] + 1j * z_tt[3]
        
        p_tt = P(tt, omega)
        dp_tt = dP_dtau(tt, omega)
        q_tt = Q(tt, omega)
        ypp_calc = (q_tt * y_tt - dp_tt * yp_tt) / p_tt if np.abs(p_tt) > 1e-12 else 0j
        
        dt = 1e-5
        yp_tt_dt = sol.sol(tt + dt)[2] + 1j * sol.sol(tt + dt)[3]
        ypp_num = (yp_tt_dt - yp_tt) / dt
        
        res = np.abs(ypp_calc - ypp_num)
        residual.append(res)
    
    errors.append((omega, tau_eval, np.array(residual), method))

# Analitik approx funksiya
def y_anal_left(tau, Omega_A):
    lam_m = lambda_minus(Omega_A)
    return np.exp(lam_m * (tau + t0))  # Normalized to 1 at -t0

def y_anal_right(tau, Omega_A, C_r):
    lam_p = lambda_plus(Omega_A)
    return C_r * np.exp(-lam_p * tau)

# Plots hissəsi
for omega, tau, y, method in plots:
    # C_r for right: match at τ=0
    idx0 = np.argmin(np.abs(tau))
    C_r = y[idx0] / np.exp(-lambda_plus(omega) * tau[idx0])
    
    # Left region: τ < -15
    mask_left = tau < -18
    y_anal_l = y_anal_left(tau[mask_left], omega)
    
    # Right region: τ > 18
    mask_right = tau > 18
    y_anal_r = y_anal_right(tau[mask_right], omega, C_r)
    
    # 1. 3D
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(tau, np.real(y), np.imag(y), label='Numerical')
    ax.plot(tau[mask_left], np.real(y_anal_l), np.imag(y_anal_l), '--')
    ax.plot(tau[mask_right], np.real(y_anal_r), np.imag(y_anal_r), '--')
    # Boundaries
    ax.scatter([-t0, t0], [np.real(y[0]), np.real(y[-1])], [np.imag(y[0]), np.imag(y[-1])], color='red', s=50)
    ax.scatter([-t0, t0], [np.real(y_anal_left(-t0, omega)), np.real(y_anal_right(t0, omega, C_r))], 
               [np.imag(y_anal_left(-t0, omega)), np.imag(y_anal_right(t0, omega, C_r))], color='green', s=50)
    ax.legend()
    ax.set_title(f'3D y(tau) for Omega_A = {omega:.4f}, method = {method}')
    plt.savefig(f'plots/3D_y_OmegaA_{np.real(omega):.2f}_{np.imag(omega):.2f}j_{method}.png')
    plt.close()

    # 2. Im(y)
    plt.figure(figsize=(8, 5))
    plt.plot(tau, np.imag(y), label='Numerical')
    plt.plot(tau[mask_left], np.imag(y_anal_l), '--')
    plt.plot(tau[mask_right], np.imag(y_anal_r), '--')
    plt.scatter([-t0, t0], [np.imag(y[0]), np.imag(y[-1])], color='red', s=50)
    plt.scatter([-t0, t0], [np.imag(y_anal_left(-t0, omega)), np.imag(y_anal_right(t0, omega, C_r))], color='green', s=50)
    plt.legend()
    plt.title(f'Im(y) vs tau for Omega_A = {omega:.4f}, method = {method}')
    plt.savefig(f'plots/Im_y_OmegaA_{np.real(omega):.2f}_{np.imag(omega):.2f}j_{method}.png')
    plt.close()

# Yeni Difference Plot
for omega, tau, y, method in plots:
    y_anal_full = np.zeros_like(y, dtype=complex)
    y_anal_full[mask_left] = y_anal_left(tau[mask_left], omega)
    y_anal_full[mask_right] = y_anal_right(tau[mask_right], omega, C_r)
    diff = np.abs(y - y_anal_full)
    
    plt.figure(figsize=(8, 5))
    plt.plot(tau, diff)
    plt.yscale('log')
    plt.xlabel('tau')
    plt.ylabel('|Num - Analytic|')
    plt.title(f'Difference Num vs Analytic y for Omega_A = {omega:.4f} ({method})')
    plt.savefig(f'plots/Difference_OmegaA_{np.real(omega):.2f}_{np.imag(omega):.2f}j_{method}.png')
    plt.close()