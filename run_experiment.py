"""
run_experiment.py — İnteraktiv Parametrik Sweep Mühərriki

İki rejim dəstəkləyir:
  [1] K-Sweep: Delta sabit, K dəyişir (0.1 → 2.0)
  [2] Delta-Sweep: K sabit, Delta dəyişir (0.0 → 1.0)

Axın:
  1) İlk (Delta₀, K₀) üçün 100x100 grid hesablanır
  2) Heatmap çəkilir, çuxurlar tapılır və istifadəçiyə göstərilir
  3) İstifadəçi hansı modları izləmək istədiyini seçir (y/n)
  4) Sweep dövrü: Continuation Method + 10x10 Lokal Bərpa
  5) Nəticələr: data/sweep_results.csv + data/eigenfunctions/*.npz

Gamma = omega_im / (1 + omega_re)
Asılılıqlar: physics.py, solvers.py
"""
import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import physics as ph
import solvers

def plot_initial_heatmap(Re, Im, mismatch_grid, candidates, Delta, k_val):
    """Plots the initial 100x100 grid heatmap with annotated local minima."""
    os.makedirs('plots', exist_ok=True)
    
    safe_grid = np.where(mismatch_grid < 1e-15, 1e-15, mismatch_grid)
    log_mismatch = np.log10(safe_grid)
    
    plt.figure(figsize=(12, 8))
    cp = plt.pcolormesh(Re, Im, log_mismatch, cmap='viridis', shading='auto')
    plt.colorbar(cp, label='Log10(Wronskian Mismatch)')
    
    for mode_id, (omega_guess, mval) in enumerate(candidates):
        plt.plot(omega_guess.real, omega_guess.imag, 'r*', markersize=14)
        plt.annotate(f'Mode {mode_id}\n$\\Omega$ ≈ {omega_guess.real:.3f}+{omega_guess.imag:.3f}i',
                     xy=(omega_guess.real, omega_guess.imag),
                     xytext=(10, 10), textcoords='offset points',
                     fontsize=8, color='white',
                     bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.7))
    
    plt.xlabel('Re($\\Omega$)')
    plt.ylabel('Im($\\Omega$)')
    plt.title(f'Initial Mismatch Heatmap | $\\Delta$ = {Delta:.3f}, K = {k_val:.3f}')
    plt.savefig('plots/initial_grid_heatmap.png', dpi=150)
    plt.close()
    print("  Heatmap saved: plots/initial_grid_heatmap.png")

def run():
    print("=" * 60)
    print("  PARAMETRIK SWEEP MÜHƏRRİKİ")
    print("=" * 60)
    
    # --- Step 1: Rejim Seçimi ---
    print("\nRejim seçin:")
    print("  [1] K-Sweep  (Həndəsi / Faizlə artım)")
    print("  [2] K-Sweep  (Xətti addım)")
    print("  [3] Delta-Sweep  (K sabit, Delta dəyişir)")
    
    while True:
        choice = input("\nSeçiminiz (1, 2 və ya 3): ").strip()
        if choice in ('1', '2', '3'):
            break
        print("Xahiş olunur 1, 2 və ya 3 daxil edin.")
    
    if choice == '1':
        sweep_mode = 'k_sweep'
        fixed_delta = float(input("Sabit Delta dəyərini daxil edin [default=0.3]: ").strip() or "0.3")
        k_start = float(input("K başlanğıc [default=0.01]: ").strip() or "0.01")
        k_end = float(input("K son [default=100.0]: ").strip() or "100.0")
        growth_rate = float(input("Böyümə faizi (%) [default=1.0]: ").strip() or "1.0")
        
        sweep_values = []
        k_curr = k_start
        multiplier = 1.0 + (growth_rate / 100.0)
        while k_curr <= k_end:
            sweep_values.append(k_curr)
            k_curr *= multiplier
            
        sweep_values = np.array(sweep_values)
        initial_k = sweep_values[0]
        initial_delta = fixed_delta
        print(f"\n→ K-Sweep (Həndəsi): Delta={fixed_delta}, K={k_start}→{k_end} ({growth_rate}% artım), {len(sweep_values)} addım")
    elif choice == '2':
        sweep_mode = 'k_sweep'
        fixed_delta = float(input("Sabit Delta dəyərini daxil edin [default=0.3]: ").strip() or "0.3")
        k_start = float(input("K başlanğıc [default=0.1]: ").strip() or "0.1")
        k_end = float(input("K son [default=2.0]: ").strip() or "2.0")
        n_steps = int(input("Addım sayı [default=20]: ").strip() or "20")
        sweep_values = np.linspace(k_start, k_end, n_steps)
        initial_k = sweep_values[0]
        initial_delta = fixed_delta
        print(f"\n→ K-Sweep (Xətti): Delta={fixed_delta}, K={k_start}→{k_end}, {n_steps} addım")
    else:
        sweep_mode = 'delta_sweep'
        fixed_k = float(input("Sabit K dəyərini daxil edin [default=1.0]: ").strip() or "1.0")
        d_start = float(input("Delta başlanğıc [default=0.0]: ").strip() or "0.0")
        d_end = float(input("Delta son [default=1.0]: ").strip() or "1.0")
        n_steps = int(input("Addım sayı [default=21]: ").strip() or "21")
        sweep_values = np.linspace(d_start, d_end, n_steps)
        initial_k = fixed_k
        initial_delta = sweep_values[0]
        print(f"\n→ Delta-Sweep: K={fixed_k}, Delta={d_start}→{d_end}, {n_steps} addım")
    
    # --- Step 2: İlk 100x100 Grid ---
    print(f"\n{'=' * 60}")
    print(f"  100x100 Grid hesablanır (Delta={initial_delta:.3f}, K={initial_k:.3f})...")
    print(f"{'=' * 60}")
    
    Re, Im, mismatch_grid_abs, mismatch_grid_complex = solvers.compute_mismatch_grid(
        initial_delta, initial_k, 
        re_bounds=(-1.2, 1.2), im_bounds=(0.2, 1.2), res=100
    )
    
    # --- Step 3: Çuxurları Tap ---
    candidates = solvers.extract_local_minima(Re, Im, mismatch_grid_abs, mismatch_grid_complex)
    
    if not candidates:
        print("XƏTA: Heç bir lokal minimum tapılmadı! Grid hüdudlarını dəyişməyə çalışın.")
        return
    
    # Heatmap çək
    plot_initial_heatmap(Re, Im, mismatch_grid_abs, candidates, initial_delta, initial_k)
    
    print(f"\n{len(candidates)} potensial çuxur tapıldı:")
    print("-" * 50)
    for mode_id, (omega_guess, mval) in enumerate(candidates):
        print(f"  Mode {mode_id}: Omega ≈ {omega_guess.real:+.4f}{omega_guess.imag:+.4f}i  |  Grid Mismatch ≈ {mval:.2e}")
    print("-" * 50)
    
    # --- Step 4: İstifadəçi Mod Seçimi ---
    selected_modes = []
    for mode_id, (omega_guess, mval) in enumerate(candidates):
        answer = input(f"  Mode {mode_id} ({omega_guess.real:+.4f}{omega_guess.imag:+.4f}i) — Hesablansın? [y/n]: ").strip().lower()
        if answer == 'y':
            # Dəqiq kök tapılsın
            print(f"    → scipy.root ilə dəqiq kök axtarılır...")
            opt_omega, err, success = solvers.find_optimal_omega(initial_delta, initial_k, omega_guess)
            if success and err < 1e-6:
                print(f"    ✓ Tapıldı: Omega = {opt_omega.real:+.6f}{opt_omega.imag:+.6f}i | Mismatch = {err:.2e}")
                selected_modes.append({'mode_id': mode_id, 'omega': opt_omega, 'error': err})
            else:
                print(f"    ✗ Bu çuxurda dəqiq kök tapılmadı (mismatch = {err:.2e}). Atlanır.")
    
    if not selected_modes:
        print("\nHeç bir mod seçilmədi. Proqram bitir.")
        return
    
    print(f"\n{len(selected_modes)} mod seçildi. Sweep başlayır...")
    
    # --- Step 5: Data Hazırlığı ---
    os.makedirs('data', exist_ok=True)
    os.makedirs('data/eigenfunctions', exist_ok=True)
    
    csv_path = 'data/sweep_results.csv'
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Delta', 'K', 'Mode_ID', 'Omega_re', 'Omega_im', 'Gamma', 'Mismatch_Error'])
    
    # --- Step 6: Sweep Dövrü ---
    # Track current omega for each mode (for continuation)
    tracked_modes = []
    for sm in selected_modes:
        tracked_modes.append({
            'mode_id': sm['mode_id'],
            'current_omega': sm['omega'],
            'alive': True
        })
    
    for step_idx, sweep_val in enumerate(sweep_values):
        if sweep_mode == 'k_sweep':
            current_delta = fixed_delta
            current_k = sweep_val
            param_label = f"K = {current_k:.4f}"
        else:
            current_delta = sweep_val
            current_k = fixed_k
            param_label = f"Delta = {current_delta:.4f}"
        
        # Məlumat kirliliyinin qarşısını almaq üçün uğurlu addımları çap etmirik
        # print(f"\n--- Step {step_idx+1}/{len(sweep_values)}: {param_label} ---")
        
        for tm in tracked_modes:
            if not tm['alive']:
                continue
                
            mode_id = tm['mode_id']
            guess = tm['current_omega']
            
            # Continuation: use previous root as initial guess
            opt_omega, err, success = solvers.find_optimal_omega(current_delta, current_k, guess)
            
            if success and err < 1e-6:
                # Uğurlu!
                tm['current_omega'] = opt_omega
                gamma = opt_omega.imag / (1 + opt_omega.real)
                
                # Uğurlu tapıntını səssizcə yaddaşa yazırıq
                # print(f"  Mode {mode_id}: Omega = {opt_omega.real:+.6f}{opt_omega.imag:+.6f}i | Mismatch = {err:.2e} | Gamma = {gamma:.6f}")
                
                # CSV-yə yaz
                with open(csv_path, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([f"{current_delta:.4f}", f"{current_k:.4f}", mode_id, 
                                     opt_omega.real, opt_omega.imag, gamma, err])
                
                # Eigenfunction saxla
                mismatch_val, sol = solvers.calculate_mismatch(opt_omega, current_delta, current_k, return_sol=True)
                if sol is not None:
                    tau = sol.t
                    y_complex = sol.y[0, :] + 1j * sol.y[1, :]
                    npz_path = f'data/eigenfunctions/step_{step_idx:03d}_mode_{mode_id}.npz'
                    np.savez_compressed(npz_path, tau=tau, y=y_complex, 
                                        Delta=current_delta, K=current_k, 
                                        Omega=opt_omega, Mode=mode_id)
            else:
                # Kök itdi — 10x10 lokal bərpa cəhdi
                print(f"  Mode {mode_id}: ⚠ Kök itdi (mismatch = {err:.2e}). 10x10 lokal bərpa cəhdi...")
                rec_omega, rec_err, rec_success = solvers.recover_root(current_delta, current_k, guess)
                
                if rec_success:
                    tm['current_omega'] = rec_omega
                    gamma = rec_omega.imag / (1 + rec_omega.real)
                    print(f"  Mode {mode_id}: ✓ BƏRPA UĞURLU! Omega = {rec_omega.real:+.6f}{rec_omega.imag:+.6f}i | Mismatch = {rec_err:.2e}")
                    
                    with open(csv_path, 'a', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow([f"{current_delta:.4f}", f"{current_k:.4f}", mode_id,
                                         rec_omega.real, rec_omega.imag, gamma, rec_err])
                    
                    mismatch_val, sol = solvers.calculate_mismatch(rec_omega, current_delta, current_k, return_sol=True)
                    if sol is not None:
                        tau = sol.t
                        y_complex = sol.y[0, :] + 1j * sol.y[1, :]
                        npz_path = f'data/eigenfunctions/step_{step_idx:03d}_mode_{mode_id}.npz'
                        np.savez_compressed(npz_path, tau=tau, y=y_complex, 
                                            Delta=current_delta, K=current_k,
                                            Omega=rec_omega, Mode=mode_id)
                else:
                    tm['alive'] = False
                    print(f"  Mode {mode_id}: ✗ MOD ÖLDÜ ({param_label} dəyərində). Daha izlənilməyəcək.")
        
        # Əgər bütün modlar öldüsə
        if not any(tm['alive'] for tm in tracked_modes):
            print("\n⚠ BÜTÜN MODLAR ÖLDÜ. Sweep dayandırılır.")
            break
    
    print(f"\n{'=' * 60}")
    print(f"  Sweep tamamlandı!")
    print(f"  CSV: {csv_path}")
    print(f"  Eigenfunctions: data/eigenfunctions/")
    print(f"{'=' * 60}")

if __name__ == '__main__':
    run()
