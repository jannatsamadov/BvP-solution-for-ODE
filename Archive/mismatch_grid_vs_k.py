"""
mismatch_grid_vs_k.py — K Parametri üzrə Mismatch Heatmap Animasiyası

Delta = 0.3 sabit tutularaq, k = 0.1..2.0 aralığında 20 addımla dəyişir.
Hər k üçün 100x50 Omega grid üzərində mismatch hesablanır, heatmap çəkilir.
Nəticədə bütün kadrlar birləşdirilərək 6 looplu GIF yaradılır.
Data: data/k_sweep/*.npz | Kadrlar: plots/k_sweep/*.png

Müstəqil skriptdir, physics.py-dan asılı deyil (öz fizika funksiyalarını saxlayır).
"""
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Fallback import for imageio
try:
    import imageio.v2 as imageio
except ImportError:
    print("Warning: imageio module not found, GIF will not be generated automatically.")

# Create plots directory
os.makedirs("plots/k_sweep", exist_ok=True)
os.makedirs("data/k_sweep", exist_ok=True)

# Global constants from the problem
beta = 0.1
g = 0.1
alpha = 0.5
l = 0.9
sigma = 1.0
t0 = 20.0  # Approximate infinity for tau range
M = 5

# FIXED Delta
DELTA = 0.3

# Def 1: xi
def xi(tau, Omega_A, Delta):
    return M * (Delta * np.tanh(tau) + Omega_A)

# Def 2: beta_A
def beta_A(tau, Omega_A, Delta, k_val):
    xi_val = xi(tau, Omega_A, Delta)
    return alpha + beta - 1 - xi_val**2

# Def 3: beta_Z
def beta_Z(tau, Omega_A, Delta):
    xi_val = xi(tau, Omega_A, Delta)
    num = 2 * alpha**2 * (xi_val**4 + 2 * g * xi_val**3 + 2 * g**2 * xi_val**2 - 5 * xi_val**2 - 6 * g * xi_val + 3)
    den = (xi_val**2 - 1) * (xi_val**4 - 6 * xi_val**2 - 4 * g * xi_val + 3)
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12
    return beta + 2 * alpha + num / den

# Def 4: P
def P(tau, Omega_A, Delta):
    bA = beta_A(tau, Omega_A, Delta, 0)  # k_val not used in beta_A
    bZ = beta_Z(tau, Omega_A, Delta)
    den = (1 - l) * bZ + l * bA
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12
    return bA * bZ / den

# Def 5: Q (k-dependent)
def Q(tau, Omega_A, Delta, k_val):
    return k_val**2 * beta_A(tau, Omega_A, Delta, k_val)

# Asymptotic limit funksiyaları
def beta_A_limit(xi_lim):
    return alpha + beta - 1 - xi_lim**2

def beta_Z_limit(xi_lim):
    num = 2 * alpha**2 * (xi_lim**4 + 2 * g * xi_lim**3 + 2 * g**2 * xi_lim**2 - 5 * xi_lim**2 - 6 * g * xi_lim + 3)
    den = (xi_lim**2 - 1) * (xi_lim**4 - 6 * xi_lim**2 - 4 * g * xi_lim + 3)
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12
    return beta + 2 * alpha + num / den

def P_limit(xi_lim):
    bA = beta_A_limit(xi_lim)
    bZ = beta_Z_limit(xi_lim)
    den = (1 - l) * bZ + l * bA
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12
    return bA * bZ / den

def Q_limit(xi_lim, k_val):
    return k_val**2 * beta_A_limit(xi_lim)

def lambda_limit(xi_lim, k_val):
    p_lim = P_limit(xi_lim)
    if np.abs(p_lim) < 1e-12:
        p_lim = 1e-12 + 1j * 1e-12
    return np.sqrt((Q_limit(xi_lim, k_val) / p_lim) + 0j)

def lambda_minus(Omega_A, Delta, k_val):
    return lambda_limit(Omega_A - Delta, k_val)

def lambda_plus(Omega_A, Delta, k_val):
    return lambda_limit(Omega_A + Delta, k_val)

# Törəmələr
def dxi_dtau(tau, Delta):
    return Delta * (1 - np.tanh(tau)**2)

def dbeta_A_dtau(tau, Omega_A, Delta):
    return -2 * xi(tau, Omega_A, Delta) * dxi_dtau(tau, Delta)

def dbeta_Z_dtau(tau, Omega_A, Delta):
    xi_val = xi(tau, Omega_A, Delta)
    dxi = dxi_dtau(tau, Delta)
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

def dP_dtau(tau, Omega_A, Delta):
    bA = beta_A(tau, Omega_A, Delta, 0)
    bZ = beta_Z(tau, Omega_A, Delta)
    dbA = dbeta_A_dtau(tau, Omega_A, Delta)
    dbZ = dbeta_Z_dtau(tau, Omega_A, Delta)
    num = bA * bZ
    dnum = dbA * bZ + bA * dbZ
    den = (1 - l) * bZ + l * bA
    dden = (1 - l) * dbZ + l * dbA
    if np.abs(den) < 1e-12:
        return 0.0 + 0j
    return (dnum * den - num * dden) / den**2

# ODE with k_val passed through
def ode_fun_real(tau, z, Omega_A, Delta, k_val):
    yr, yi, ypr, ypi = z
    y = yr + 1j * yi
    yp = ypr + 1j * ypi
    p = P(tau, Omega_A, Delta)
    dp = dP_dtau(tau, Omega_A, Delta)
    q = Q(tau, Omega_A, Delta, k_val)
    if np.abs(p) < 1e-12:
        p = 1e-12 + 1j * 1e-12
    ypp = (q * y - dp * yp) / p
    return np.array([ypr, ypi, np.real(ypp), np.imag(ypp)])

# Worker function — now k_val aware
def compute_single_point(Omega_A, Delta, k_val):
    lam_m = lambda_minus(Omega_A, Delta, k_val)
    z0 = [1.0, 0.0, np.real(lam_m), np.imag(lam_m)]
    
    if not np.all(np.isfinite(z0)):
        return np.nan
        
    tau_span = (-t0, t0)
    sol = solve_ivp(ode_fun_real, tau_span, z0, args=(Omega_A, Delta, k_val), method='RK45',
                    rtol=1e-8, atol=1e-8, dense_output=False)
    
    if not sol.success:
        return np.nan
    
    z_right = sol.y[:, -1]
    yr, yi, ypr, ypi = z_right
    y_right = yr + 1j * yi
    yp_right = ypr + 1j * ypi
    
    lam_p = lambda_plus(Omega_A, Delta, k_val)
    mismatch = yp_right + lam_p * y_right
    return np.abs(mismatch)

def plot_heatmap(Real_Grid, Imag_Grid, Mismatch_Grid, k_val, frame_idx):
    plt.figure(figsize=(10, 6))
    plt.pcolormesh(Real_Grid, Imag_Grid, np.log10(Mismatch_Grid + 1e-10), shading='auto', cmap='viridis', vmax=6)
    plt.colorbar(label='log10(|mismatch|)')
    plt.xlabel('Re(Omega_A)')
    plt.ylabel('Im(Omega_A)')
    plt.title(f'Mismatch Heatmap (Delta = {DELTA:.3f}, k = {k_val:.3f})')
    
    frame_filename = f'plots/k_sweep/frame_{frame_idx:02d}.png'
    plt.savefig(frame_filename)
    plt.close()
    return frame_filename

if __name__ == '__main__':
    print("====================================")
    print(f"K-Sweep for Fixed Delta = {DELTA}")
    print("====================================")
    
    # Omega grid
    real_parts = np.linspace(-1, 1, 100)
    imag_parts = np.linspace(0, 1, 50)
    Real, Imag = np.meshgrid(real_parts, imag_parts)
    
    # K values: 0.1 to 2.0, 20 steps
    k_values = np.linspace(0.1, 2.0, 20)
    frame_filenames = []
    
    for idx, k_val in enumerate(k_values):
        print(f"\nProcessing Frame {idx+1}/{len(k_values)} [k = {k_val:.3f}]...")
        
        mismatch_grid = np.full_like(Real, np.nan, dtype=float)
        
        total = Real.shape[0] * Real.shape[1]
        count = 0
        for i in range(Real.shape[0]):
            for j in range(Real.shape[1]):
                Omega_A = Real[i, j] + 1j * Imag[i, j]
                mismatch_grid[i, j] = compute_single_point(Omega_A, DELTA, k_val)
                count += 1
            # Print progress per row
            if (i + 1) % 10 == 0:
                print(f"  Row {i+1}/{Real.shape[0]} done ({count}/{total} points)")
        
        # Save Heatmap as frame
        fname = plot_heatmap(Real, Imag, mismatch_grid, k_val, idx)
        frame_filenames.append(fname)
        print(f"-> Saved {fname}")
        
        # Save raw data as compressed NPZ
        npz_filename = f'data/k_sweep/mismatch_k_{k_val:.3f}.npz'
        np.savez_compressed(npz_filename, Real=Real, Imag=Imag, Mismatch=mismatch_grid, k=k_val, Delta=DELTA)
        print(f"-> Data saved to {npz_filename}")
        
    print("\nRendering Animation GIF...")
    try:
        gif_path = 'plots/K_Sweep_Evolution.gif'
        images = []
        for filename in frame_filenames:
            images.append(imageio.imread(filename))
        # 6 loops: repeat frames 6 times
        images_looped = images * 6
        imageio.mimsave(gif_path, images_looped, fps=3)
        print(f"Successfully generated {gif_path} (6 loops)!")
    except Exception as e:
        print(f"Could not construct final animated GIF automatically: {e}")
        print("Your frame PNGs are available in the plots/k_sweep/ folder.")

    print("DONE.")
