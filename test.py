"""
test.py — Manual Omega Test və Diaqnostika Skripti

Vizual olaraq seçilmiş Omega dəyərinin mismatch xətasını hesablayır.
Optimizatorun (scipy.root) həmin nöqtədən hara sürüşdüyünü müqayisə edir.
Verilmiş Omega üçün y(tau) eigenfunction qrafiki çəkir (plots/test_eigenfunction.png).

Asılılıqlar: solvers.py → physics.py
"""

import numpy as np
import matplotlib.pyplot as plt
import solvers
import os

Delta = 0.48
K_VAL = 0.01  # Dalğa ədədi
guess = -0.1579+0.0526j
#guess = -0.01274956402778681 - 7.86651924488455e-13j


# 1. Sırf sizin verdiyiniz dəyər üçün Mismatch hesablanması (Sürüşmədən)
err_complex, sol_guess = solvers.calculate_mismatch(guess, Delta, K_VAL, return_sol=True)
err_abs = np.abs(err_complex) if isinstance(err_complex, complex) else err_complex

print("\n=== SİZİN VERDİYİNİZ XÜSUSİ NÖQTƏ ===")
print(f"Kök (Omega): {guess}")
print(f"Mismatch (Xəta): {err_abs:.2e}")
print("======================================\n")

# 2. Optimizatorun hara sürüşdüyünü yoxlayaq:
opt_omega, err_opt, success = solvers.find_optimal_omega(Delta, K_VAL, guess)

print("=== OPTİMİZATORUN SÜRÜŞDÜYÜ NÖQTƏ ===")
print(f"Tapılan Kök: {opt_omega}")
print(f"Mismatch (Xəta): {err_opt:.2e}")
print(f"Uğurludurmu?: {success}")
print("======================================\n")

# Müqayisə nəticəsi:
if err_abs < err_opt:
    print("!!! SİZ HAQLISINIZ: Mühərrik daha pis yerə sürüşüb !!!\n")
else:
    print("Mühərrik xətanı daha da kiçiltməyə nail olub.\n")

# 3. Sizin verdiyiniz Dəyər üçün y vs tau qrafikinin qurulması
if sol_guess is not None:
    os.makedirs('plots', exist_ok=True)
    tau = sol_guess.t
    y = sol_guess.y[0, :] + 1j * sol_guess.y[1, :]
    
    plt.figure(figsize=(10, 6))
    plt.plot(tau, np.abs(y), 'k-', linewidth=2, label='|y(tau)| (Amplitude)')
    plt.plot(tau, np.real(y), 'b--', linewidth=1.5, label='Re(y)')
    plt.plot(tau, np.imag(y), 'r--', linewidth=1.5, label='Im(y)')
    
    plt.xlabel('$\\tau$ (Spatial grid)')
    plt.ylabel('Eigenfunction y($\\tau$)')
    plt.title(f'Eigenfunction for $\Delta$={Delta} | Təxmininiz: $\Omega$={guess}')
    plt.legend()
    plt.grid(True)
    
    # Sağlam görmək üçün limit
    plt.ylim(-0.002, 0.002)
    
    filename = "plots/test_eigenfunction.png"
    plt.savefig(filename)
    print(f"Qrafik uğurla yaddaşa verildi: {filename}")
else:
    print("ODE Solver çökdü, qrafik qurula bilmədi.")
