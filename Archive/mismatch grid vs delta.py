#Bu kod digər dəyərlərin fix edilmiş qiymətində deltaya uyğun olaraq omega mismatch dəyərləri gridinin necə dəyişdiyini müəyyənləşdirmək üçündür.
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import multiprocessing

# Fallback import for imageio
try:
    import imageio.v2 as imageio
except ImportError:
    print("Warning: imageio module not found, GIF will not be generated automatically.")

# Create plots directory
os.makedirs("plots", exist_ok=True)

# Global constants from the problem (except Delta)
beta = 0.1
g = 0.1
alpha = 0.5
l = 0.9
sigma = 1.0
k = 10.0
t0 = 20.0 # Approximate infinity for tau range
M = 5

# Def 1: xi
def xi(tau, Omega_A, Delta):
    return M * (Delta * np.tanh(tau) + Omega_A)

# Def 2: beta_A
def beta_A(tau, Omega_A, Delta):
    xi_val = xi(tau, Omega_A, Delta)
    return alpha + beta - 1 - xi_val**2

# Def 3: beta_Z
def beta_Z(tau, Omega_A, Delta):
    xi_val = xi(tau, Omega_A, Delta)
    num = 2 * alpha**2 * (xi_val**4 + 2 * g * xi_val**3 + 2 * g**2 * xi_val**2 - 5 * xi_val**2 - 6 * g * xi_val + 3)
    den = (xi_val**2 - 1) * (xi_val**4 - 6 * xi_val**2 - 4 * g * xi_val + 3)
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12  # Safeguard
    return beta + 2 * alpha + num / den

# Def 4: P
def P(tau, Omega_A, Delta):
    bA = beta_A(tau, Omega_A, Delta)
    bZ = beta_Z(tau, Omega_A, Delta)
    den = (1 - l) * bZ + l * bA
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12
    return bA * bZ / den

# Def 5: Q
def Q(tau, Omega_A, Delta):
    return k**2 * beta_A(tau, Omega_A, Delta)

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

def Q_limit(xi_lim):
    return k**2 * beta_A_limit(xi_lim)

def lambda_limit(xi_lim):
    p_lim = P_limit(xi_lim)
    if np.abs(p_lim) < 1e-12:
        p_lim = 1e-12 + 1j * 1e-12
    return np.sqrt(Q_limit(xi_lim) / p_lim)

def lambda_minus(Omega_A, Delta):
    return lambda_limit(Omega_A - Delta)

def lambda_plus(Omega_A, Delta):
    return lambda_limit(Omega_A + Delta)

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
    bA = beta_A(tau, Omega_A, Delta)
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

# Def 19: ode_fun_real
def ode_fun_real(tau, z, Omega_A, Delta):
    yr, yi, ypr, ypi = z
    y = yr + 1j * yi
    yp = ypr + 1j * ypi
    p = P(tau, Omega_A, Delta)
    dp = dP_dtau(tau, Omega_A, Delta)
    q = Q(tau, Omega_A, Delta)
    if np.abs(p) < 1e-12:
        p = 1e-12 + 1j * 1e-12
    ypp = (q * y - dp * yp) / p
    return np.array([ypr, ypi, np.real(ypp), np.imag(ypp)])

# Worker function for multiprocessing
def compute_single_point(args):
    Omega_A, Delta = args
    lam_m = lambda_minus(Omega_A, Delta)
    z0 = [1.0, 0.0, np.real(lam_m), np.imag(lam_m)]
    
    if not np.all(np.isfinite(z0)):
        return np.nan
        
    tau_span = (-t0, t0)
    sol = solve_ivp(ode_fun_real, tau_span, z0, args=(Omega_A, Delta), method='RK45',
                    rtol=1e-8, atol=1e-8, dense_output=False)
    
    if not sol.success:
        return np.nan
    
    # We only care about the right boundary at t0
    z_right = sol.y[:, -1]
    yr, yi, ypr, ypi = z_right
    y_right = yr + 1j * yi
    yp_right = ypr + 1j * ypi
    
    lam_p = lambda_plus(Omega_A, Delta)
    mismatch = yp_right + lam_p * y_right
    return np.abs(mismatch)

def plot_heatmap(Real_Grid, Imag_Grid, Mismatch_Grid, Delta, frame_idx):
    plt.figure(figsize=(10, 6))
    plt.pcolormesh(Real_Grid, Imag_Grid, np.log10(Mismatch_Grid + 1e-10), shading='auto', cmap='viridis_r', vmax=6)
    plt.colorbar(label='log10(|mismatch|)')
    plt.xlabel('Re(Omega_A)')
    plt.ylabel('Im(Omega_A)')
    plt.title(f'Mismatch Heatmap (Delta = {Delta:.3f})')
    
    frame_filename = f'plots/frame_{frame_idx:02d}.png'
    plt.savefig(frame_filename)
    plt.close()
    return frame_filename

if __name__ == '__main__':
    print("====================================")
    print("Starting High-Resolution Grid Evolution")
    print("====================================")
    
    # Grid as specified (100 Im, 200 Re -> 200 col, 100 row?) 
    # Usually linspace maps to axis. We want 200 horizontal (Re), 100 vertical (Im).
    real_parts = np.linspace(-1.2, 1.2, 60)
    imag_parts = np.linspace(0, 1.2, 30)
    Real, Imag = np.meshgrid(real_parts, imag_parts)
    Omega_grid = Real + 1j * Imag
    Omega_flat = Omega_grid.flatten() # 20,000 points
    
    # Delta cycle
    delta_values = np.linspace(0, 1, 10)
    frame_filenames = []
    
    for idx, Delta in enumerate(delta_values):
        print(f"\nProcessing Frame {idx+1}/{len(delta_values)} [Delta = {Delta:.3f}]...")
        
        mismatch_grid = np.full_like(Omega_grid, np.nan, dtype=float)
        
        print(f"Executing {len(Omega_flat)} integrations sequentially...")
        for i in range(Real.shape[0]):
            for j in range(Real.shape[1]):
                Omega_A = Real[i, j] + 1j * Imag[i, j]
                args = (Omega_A, Delta)
                mismatch_grid[i, j] = compute_single_point(args)
        
        # Save Heatmap as frame
        fname = plot_heatmap(Real, Imag, mismatch_grid, Delta, idx)
        frame_filenames.append(fname)
        print(f"-> Saved {fname}")
        
        # Save raw data to CSV (X, Y, Z/Mismatch format)
        csv_filename = f'plots/data_frame_{idx:02d}_delta_{Delta:.3f}.csv'
        data_to_save = np.column_stack((Real.flatten(), Imag.flatten(), mismatch_grid.flatten()))
        np.savetxt(csv_filename, data_to_save, delimiter=',', header='Re(Omega),Im(Omega),Mismatch', comments='')
        
    print("\nRendering Animation Video...")
    try:
        gif_path = 'Delta_Evolution.gif'
        images = []
        
        for filename in frame_filenames:
            images.append(imageio.imread(filename))
            
        # 6 loops: repeat frames 6 times
        images_looped = images * 6
        imageio.mimsave(gif_path, images_looped, fps=3)
        print(f"Successfully generated {gif_path} in the project folder!")
    except Exception as e:
        print(f"Could not construct final animated GIF automatically: {e}")
        print("Your frame PNGs are available in the plots/ folder to manually stitch together.")

    print("DONE.")
