"""
physics.py — Fiziki Sabitlər və Analitik Funksiyalar

Fiziki parametrlər (BETA, G, ALPHA, L, SIGMA, T0, M) burada saxlanılır.
K (dalğa ədədi) qlobal sabit DEYİL — bütün K-asılı funksiyalara k_val arqumenti ilə ötürülür.
Əsas funksiyalar: xi, beta_A, beta_Z, P, Q və onların törəmələri (dP/dtau).
Asimptotik limit funksiyaları (lambda_minus, lambda_plus) sərhəd şərtlərini hesablayır.

Bu modul digər bütün skriptlər tərəfindən import edilir.
"""
import numpy as np

# Physical Constants
BETA = 0.1
G = 0.1
ALPHA = 0.5
L = 0.9
SIGMA = 1.0
T0 = 20.0  # Approximate infinity for tau range
M = 5.0    # Mach Number for growth rate calculation

# --- Analytical Functions ---

def xi(tau, Omega_A, Delta):
    return M * (Delta * np.tanh(tau) + Omega_A)

def beta_A(tau, Omega_A, Delta):
    xi_val = xi(tau, Omega_A, Delta)
    return ALPHA + BETA - 1 - xi_val**2

def beta_Z(tau, Omega_A, Delta):
    xi_val = xi(tau, Omega_A, Delta)
    num = 2 * ALPHA**2 * (xi_val**4 + 2 * G * xi_val**3 + 2 * G**2 * xi_val**2 - 5 * xi_val**2 - 6 * G * xi_val + 3)
    den = (xi_val**2 - 1) * (xi_val**4 - 6 * xi_val**2 - 4 * G * xi_val + 3)
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12
    return BETA + 2 * ALPHA + num / den

def P(tau, Omega_A, Delta):
    bA = beta_A(tau, Omega_A, Delta)
    bZ = beta_Z(tau, Omega_A, Delta)
    den = (1 - L) * bZ + L * bA
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12
    return bA * bZ / den

def Q(tau, Omega_A, Delta, k_val):
    return (k_val / SIGMA)**2 * beta_A(tau, Omega_A, Delta)

# --- Asymptotic Limits ---

def beta_A_limit(xi_lim):
    return ALPHA + BETA - 1 - xi_lim**2

def beta_Z_limit(xi_lim):
    num = 2 * ALPHA**2 * (xi_lim**4 + 2 * G * xi_lim**3 + 2 * G**2 * xi_lim**2 - 5 * xi_lim**2 - 6 * G * xi_lim + 3)
    den = (xi_lim**2 - 1) * (xi_lim**4 - 6 * xi_lim**2 - 4 * G * xi_lim + 3)
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12
    return BETA + 2 * ALPHA + num / den

def P_limit(xi_lim):
    bA = beta_A_limit(xi_lim)
    bZ = beta_Z_limit(xi_lim)
    den = (1 - L) * bZ + L * bA
    if np.abs(den) < 1e-12:
        den = 1e-12 + 1j * 1e-12
    return bA * bZ / den

def Q_limit(xi_lim, k_val):
    return (k_val / SIGMA)**2 * beta_A_limit(xi_lim)

def lambda_limit(xi_lim, k_val):
    p_lim = P_limit(xi_lim)
    if np.abs(p_lim) < 1e-12:
        p_lim = 1e-12 + 1j * 1e-12
    return np.sqrt((Q_limit(xi_lim, k_val) / p_lim) + 0j)

def lambda_minus(Omega_A, Delta, k_val):
    return lambda_limit(M * (Omega_A - Delta), k_val)

def lambda_plus(Omega_A, Delta, k_val):
    return lambda_limit(M * (Omega_A + Delta), k_val)

# --- Derivatives ---

def dxi_dtau(tau, Delta):
    return M * Delta * (1 - np.tanh(tau)**2)

def dbeta_A_dtau(tau, Omega_A, Delta):
    return -2 * xi(tau, Omega_A, Delta) * dxi_dtau(tau, Delta)

def dbeta_Z_dtau(tau, Omega_A, Delta):
    xi_val = xi(tau, Omega_A, Delta)
    dxi = dxi_dtau(tau, Delta)
    num = 2 * ALPHA**2 * (xi_val**4 + 2*G*xi_val**3 + 2*G**2*xi_val**2 -5*xi_val**2 -6*G*xi_val +3)
    dnum = 2 * ALPHA**2 * dxi * (4*xi_val**3 + 6*G*xi_val**2 + 4*G**2*xi_val -10*xi_val -6*G)
    a = xi_val**2 - 1
    da = 2 * xi_val * dxi
    b = xi_val**4 -6*xi_val**2 -4*G*xi_val +3
    db = dxi * (4*xi_val**3 -12*xi_val -4*G)
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
    den = (1 - L) * bZ + L * bA
    dden = (1 - L) * dbZ + L * dbA
    if np.abs(den) < 1e-12:
        return 0.0 + 0j
    return (dnum * den - num * dden) / den**2
