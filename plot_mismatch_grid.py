"""
plot_mismatch_grid.py — Tək Delta/K üçün Mismatch Heatmap Vizuallaşdırıcısı

Verilmiş bir (Delta, K) cütü üçün Omega_re x Omega_im müstəvisində
Wronskian Mismatch-in Log10 rəngli xəritəsini (heatmap) qurur.
Lokal minimumları (çuxurları) qırmızı ulduzlarla işarələyir və
hər çuxurun Omega dəyərini annotasiya edir.
Sürətli vizual diaqnostika və eksperiment üçün istifadə olunur.

Asılılıqlar: solvers.py → physics.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import solvers
import pandas as pd

def generate_and_plot_grid(Delta, k_val, re_bounds=(-1, 1), im_bounds=(0, 0.8), res=30):
    print(f"Generating {res}x{res} mismatch grid for Delta={Delta}, K={k_val}...")
    
    Re, Im, mismatch_grid, mismatch_grid_complex = solvers.compute_mismatch_grid(Delta, k_val, re_bounds, im_bounds, res)
    
    # Log10 for visualization
    safe_grid = np.where(mismatch_grid < 1e-15, 1e-15, mismatch_grid)
    log_mismatch = np.log10(safe_grid)
    
    # Extract local minima with Phase Vortex filtering
    candidates = solvers.extract_local_minima(Re, Im, mismatch_grid, mismatch_grid_complex)
    
    os.makedirs('plots', exist_ok=True)
    
    plt.figure(figsize=(10, 8))
    cp = plt.pcolormesh(Re, Im, log_mismatch, cmap='viridis', shading='auto')
    plt.colorbar(cp, label='Log10(Mismatch Error)')
    
    # Annotate local minima
    for mode_id, (omega_guess, mval) in enumerate(candidates):
        plt.plot(omega_guess.real, omega_guess.imag, 'r*', markersize=12)
        plt.annotate(f'Mode {mode_id}\n{omega_guess.real:.3f}+{omega_guess.imag:.3f}i',
                     xy=(omega_guess.real, omega_guess.imag),
                     xytext=(8, 8), textcoords='offset points',
                     fontsize=7, color='white',
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.7))
        
    plt.xlabel('Re($\\Omega$)')
    plt.ylabel('Im($\\Omega$)')
    plt.title(f'Wronskian Mismatch Heatmap | $\\Delta$ = {Delta}, K = {k_val} | Grid: {res}x{res}')
    
    filename = f'plots/mismatch_grid_delta_{Delta:.3f}_k_{k_val:.3f}.png'
    plt.savefig(filename, dpi=150)
    plt.close()
    
    # --- YENİ: Phase Vortex Heatmap ---
    plt.figure(figsize=(10, 8))
    phase_grid = np.angle(mismatch_grid_complex)
    cp2 = plt.pcolormesh(Re, Im, phase_grid, cmap='twilight', shading='auto')
    plt.colorbar(cp2, label='Phase of Mismatch (radians)')
    
    # Çuxurları fazada da göstərək
    for mode_id, (omega_guess, mval) in enumerate(candidates):
        plt.plot(omega_guess.real, omega_guess.imag, 'wo', markersize=4) # Ağ nöqtələr
        
    plt.xlabel('Re($\\Omega$)')
    plt.ylabel('Im($\\Omega$)')
    plt.title(f'Phase Vortex Map (Argument Principle) | $\\Delta$ = {Delta}, K = {k_val}')
    
    phase_filename = f'plots/phase_grid_delta_{Delta:.3f}_k_{k_val:.3f}.png'
    plt.savefig(phase_filename, dpi=150)
    plt.close()
    
    print(f"Grid saved to {filename}")
    print(f"Found {len(candidates)} local minima:")
    for i, (omega, mval) in enumerate(candidates):
        print(f"  Mode {i}: Omega ≈ {omega.real:+.4f}{omega.imag:+.4f}i  |  Mismatch ≈ {mval:.2e}")
        
    # --- YENİ: CSV Eksport ---
    os.makedirs('data', exist_ok=True)
    df = pd.DataFrame({
        'Omega_Re': Re.flatten(),
        'Omega_Im': Im.flatten(),
        'Mismatch_Re': mismatch_grid_complex.real.flatten(),
        'Mismatch_Im': mismatch_grid_complex.imag.flatten(),
        'Mismatch_Abs': mismatch_grid.flatten(),
        'Mismatch_Phase': np.angle(mismatch_grid_complex).flatten()
    })
    csv_filename = f'data/grid_complex_mismatch_Delta_{Delta:.3f}_K_{k_val:.3f}.csv'
    df.to_csv(csv_filename, index=False)
    print(f"Complex grid data exported to {csv_filename}")

if __name__ == '__main__':
    # Siz istədiyiniz vaxt bu parametrləri dəyişərək fərqli grid eksperimentləri edə bilərsiniz:
    DELTA_TO_TEST = 0.3
    K_TO_TEST = 0.1
    RESOLUTION = 240
    RE_BOUNDS = (-1.2, 1.2)
    IM_BOUNDS = (0, 1.2)
    
    generate_and_plot_grid(DELTA_TO_TEST, K_TO_TEST, RE_BOUNDS, IM_BOUNDS, RESOLUTION)
