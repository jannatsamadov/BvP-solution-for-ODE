"""
solvers.py — ODE Həlledici və Optimizasiya Mühərriki

Midpoint Matching metodu: Dalğanı -T0→0 və +T0→0 istiqamətlərindən atır,
Wronskian determinantı ilə mismatch hesablayır (partlama yox edilir).
Hibrid axtarış: Kobud grid + scipy.optimize.root ilə dəqiq kök tapıcısı.
Çoxlu kökləri (Mode_0, Mode_1...) eyni anda aşkarlayır.
10x10 lokal grid ilə itmiş kökləri bərpa etmə qabiliyyəti.

K qlobal sabit deyil — bütün funksiyalara k_val arqumenti ilə ötürülür.
Asılılıqlar: physics.py
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
import scipy.ndimage as ndimage
import physics as ph

def ode_fun_real(tau, z, Omega_A, Delta, k_val):
    """
    Real-valued representation of the 2nd order complex ODE:
    P(tau)*y'' + P'(tau)*y' - Q(tau)*y = 0
    """
    yr, yi, ypr, ypi = z
    y = yr + 1j * yi
    yp = ypr + 1j * ypi
    
    p = ph.P(tau, Omega_A, Delta)
    dp = ph.dP_dtau(tau, Omega_A, Delta)
    q = ph.Q(tau, Omega_A, Delta, k_val)
    
    if np.abs(p) < 1e-12:
        p = 1e-12 + 1j * 1e-12
        
    ypp = (q * y - dp * yp) / p
    return np.array([ypr, ypi, np.real(ypp), np.imag(ypp)])

def calculate_mismatch(omega_complex, Delta, k_val, return_sol=False):
    """
    Midpoint Matching Method:
    Shoots from -T0 to 0 and from T0 to 0.
    Mismatch is the Wronskian at tau=0.
    """
    Omega_A = omega_complex
    
    # --- Left Integration (-T0 to 0) ---
    lam_m = ph.lambda_minus(Omega_A, Delta, k_val)
    start_val_L = 1e-3 + 1e-3 * 1j
    y_prime_start_L = lam_m * start_val_L
    z0_L = [np.real(start_val_L), np.imag(start_val_L), np.real(y_prime_start_L), np.imag(y_prime_start_L)]
    
    if not np.all(np.isfinite(z0_L)):
        return 1e10 if not return_sol else (1e10, None)
        
    sol_L = solve_ivp(ode_fun_real, (-ph.T0, 0), z0_L, args=(Omega_A, Delta, k_val), method='RK45',
                      rtol=1e-10, atol=1e-10, max_step=0.1, dense_output=True)
                      
    if not sol_L.success:
        return 1e10 if not return_sol else (1e10, None)
        
    z_L_0 = sol_L.y[:, -1]
    y_L_0 = z_L_0[0] + 1j * z_L_0[1]
    yp_L_0 = z_L_0[2] + 1j * z_L_0[3]
    
    # --- Right Integration (+T0 to 0) ---
    lam_p = ph.lambda_plus(Omega_A, Delta, k_val)
    start_val_R = 1e-3 + 1e-3 * 1j
    y_prime_start_R = -lam_p * start_val_R  # Negative sign for decaying mode at +infinity
    z0_R = [np.real(start_val_R), np.imag(start_val_R), np.real(y_prime_start_R), np.imag(y_prime_start_R)]
    
    if not np.all(np.isfinite(z0_R)):
        return 1e10 if not return_sol else (1e10, None)
        
    sol_R = solve_ivp(ode_fun_real, (ph.T0, 0), z0_R, args=(Omega_A, Delta, k_val), method='RK45',
                      rtol=1e-10, atol=1e-10, max_step=0.1, dense_output=True)
                      
    if not sol_R.success:
        return 1e10 if not return_sol else (1e10, None)
        
    z_R_0 = sol_R.y[:, -1]
    y_R_0 = z_R_0[0] + 1j * z_R_0[1]
    yp_R_0 = z_R_0[2] + 1j * z_R_0[3]
    
    # --- Wronskian Mismatch at tau=0 ---
    mismatch = y_L_0 * yp_R_0 - yp_L_0 * y_R_0
    
    if return_sol:
        if np.abs(y_R_0) < 1e-15:
            C = 1.0
        else:
            C = y_L_0 / y_R_0
            
        t_L = sol_L.t
        y_L_complex = sol_L.y[0, :] + 1j * sol_L.y[1, :]
        
        t_R = sol_R.t
        y_R_complex = sol_R.y[0, :] + 1j * sol_R.y[1, :]
        
        # Reverse right arrays if they were integrated backwards in time natively
        if len(t_R) > 1 and t_R[0] > t_R[-1]:
            t_R = t_R[::-1]
            y_R_complex = y_R_complex[::-1]
            
        y_R_scaled = C * y_R_complex
        
        # Merge exactly at the midpoint to avoid duplicate 0 points
        if t_R[0] == 0.0 and t_L[-1] == 0.0:
            t_full = np.concatenate((t_L[:-1], t_R))
            y_full = np.concatenate((y_L_complex[:-1], y_R_scaled))
        else:
            t_full = np.concatenate((t_L, t_R))
            y_full = np.concatenate((y_L_complex, y_R_scaled))
            
        class DummySol:
            pass
        dsol = DummySol()
        dsol.t = t_full
        dsol.y = np.array([np.real(y_full), np.imag(y_full)])
        return mismatch, dsol
        
    return mismatch

def objective_function(omega_vec, Delta, k_val):
    """
    Wrapper for calculate_mismatch to be used with scipy.optimize.root.
    omega_vec: [Re(Omega), Im(Omega)]
    Returns: [Re(mismatch), Im(mismatch)]
    """
    omega_complex = omega_vec[0] + 1j * omega_vec[1]
    mismatch = calculate_mismatch(omega_complex, Delta, k_val)
    if isinstance(mismatch, float) and mismatch == 1e10:
        return [1e10, 1e10]
    return [np.real(mismatch), np.imag(mismatch)]

def find_optimal_omega(Delta, k_val, initial_guess_complex, tol=1e-6):
    """
    Finds exactly one Omega_A using a specific initial guess.
    """
    x0 = [initial_guess_complex.real, initial_guess_complex.imag]
    res = root(objective_function, x0, args=(Delta, k_val), method='hybr', tol=tol)
    opt_omega = res.x[0] + 1j * res.x[1]
    final_mismatch = calculate_mismatch(opt_omega, Delta, k_val)
    err = 1e10 if isinstance(final_mismatch, float) else np.abs(final_mismatch)
    return opt_omega, err, res.success

# ============================================================
# Grid & Multi-Root Utilities
# ============================================================

def compute_mismatch_grid(Delta, k_val, re_bounds=(-1, 1), im_bounds=(0, 0.8), res=100):
    """
    Computes a 2D mismatch grid over the Omega complex plane.
    Returns: (Re_grid, Im_grid, Mismatch_grid)
    """
    re_vals = np.linspace(re_bounds[0], re_bounds[1], res)
    im_vals = np.linspace(im_bounds[0], im_bounds[1], res)
    Re, Im = np.meshgrid(re_vals, im_vals)
    
    mismatch_grid_abs = np.zeros_like(Re)
    mismatch_grid_complex = np.zeros_like(Re, dtype=complex)
    total = res * res
    count = 0
    for i in range(res):
        for j in range(res):
            omega_c = Re[i, j] + 1j * Im[i, j]
            m = calculate_mismatch(omega_c, Delta, k_val)
            mismatch_grid_complex[i, j] = m
            mismatch_grid_abs[i, j] = np.abs(m) if isinstance(m, complex) else np.abs(m)
            count += 1
            if count % 500 == 0:
                print(f"  Grid progress: {count}/{total}")
                
    return Re, Im, mismatch_grid_abs, mismatch_grid_complex

def has_phase_vortex(mismatch_grid_complex, r, c):
    """
    Checks if there is a topological phase vortex around the pixel (r, c).
    Sums the phase differences along a 3x3 loop (8 boundary pixels).
    """
    rows, cols = mismatch_grid_complex.shape
    if r == 0 or r == rows - 1 or c == 0 or c == cols - 1:
        return False # Edge pixel, cannot form a full loop
        
    # Anti-clockwise loop around (r, c)
    loop_indices = [
        (r-1, c-1), (r, c-1), (r+1, c-1), (r+1, c),
        (r+1, c+1), (r, c+1), (r-1, c+1), (r-1, c)
    ]
    
    phases = [np.angle(mismatch_grid_complex[ir, ic]) for ir, ic in loop_indices]
    
    total_winding = 0.0
    for i in range(len(phases)):
        p1 = phases[i]
        p2 = phases[(i + 1) % len(phases)]
        
        diff = p2 - p1
        # Wrap diff to [-pi, pi]
        diff = (diff + np.pi) % (2 * np.pi) - np.pi
        total_winding += diff
        
    # The winding number is total_winding / (2*pi).
    # If it is a vortex, winding number should be +1 or -1 (i.e. abs(total_winding) ~ 2*pi)
    return np.abs(total_winding) > np.pi

def extract_local_minima(Re, Im, mismatch_grid_abs, mismatch_grid_complex=None):
    """
    Finds local minima in the mismatch grid using a 3x3 neighborhood filter.
    Returns list of (omega_guess, mismatch_value) sorted by |omega| descending (furthest from 0 = Mode 0).
    """
    neighborhood = np.ones((3, 3))
    local_min = ndimage.minimum_filter(mismatch_grid_abs, footprint=neighborhood) == mismatch_grid_abs
    min_indices = np.argwhere(local_min)
    
    candidates = []
    for idx in min_indices:
        r_i, c_i = idx[0], idx[1]
        mval = mismatch_grid_abs[r_i, c_i]
        
        # 1. Şərt: Mismatch kifayət qədər kiçik olmalıdır (10^-3 tərtibində)
        if mval > 1e-3:
            continue
            
        # 2. Şərt: Topoloji Faza Vortex-i olmalıdır
        if mismatch_grid_complex is not None:
            if not has_phase_vortex(mismatch_grid_complex, r_i, c_i):
                continue
                
        omega_guess = Re[r_i, c_i] + 1j * Im[r_i, c_i]
        candidates.append((omega_guess, mval))
    
    # Sort by |Omega| descending (furthest from origin = Mode 0)
    candidates.sort(key=lambda x: np.abs(x[0]), reverse=True)
    return candidates

def find_all_optimal_omegas(Delta, k_val, re_bounds=(-1, 1), im_bounds=(0.0, 0.8), res=30, tol=1e-6):
    """
    Hybrid Strategy: 
    1. Fast Coarse Grid over parameter space to map the mismatch trenches.
    2. Local Minima extraction.
    3. Multi-Root finding starting from all local minima.
    """
    Re, Im, mismatch_grid_abs, mismatch_grid_complex = compute_mismatch_grid(Delta, k_val, re_bounds, im_bounds, res)
    candidates = extract_local_minima(Re, Im, mismatch_grid_abs, mismatch_grid_complex)
    
    valid_roots = []
    
    for guess_complex, _ in candidates:
        opt_omega, err, success = find_optimal_omega(Delta, k_val, guess_complex, tol=tol)
        
        if success and err <= tol:
            # Check for duplicates
            is_dup = False
            for v_omega, v_err in valid_roots:
                if np.abs(v_omega - opt_omega) < 1e-4:
                    is_dup = True
                    break
            if not is_dup:
                valid_roots.append((opt_omega, err))
                
    return valid_roots

def recover_root(Delta, k_val, last_omega, radius=0.1, local_res=10, tol=1e-6):
    """
    Attempts to recover a lost root by searching a 10x10 local grid
    centered around the last known position of the root.
    Returns: (opt_omega, err, success)
    """
    re_center = last_omega.real
    im_center = last_omega.imag
    re_bounds = (re_center - radius, re_center + radius)
    im_bounds = (max(0, im_center - radius), im_center + radius)
    
    Re, Im, mismatch_grid_abs, mismatch_grid_complex = compute_mismatch_grid(Delta, k_val, re_bounds, im_bounds, local_res)
    candidates = extract_local_minima(Re, Im, mismatch_grid_abs, mismatch_grid_complex)
    
    for guess_complex, _ in candidates:
        opt_omega, err, success = find_optimal_omega(Delta, k_val, guess_complex, tol=tol)
        if success and err <= tol:
            return opt_omega, err, True
            
    return last_omega, 1e10, False
