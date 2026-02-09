"""
plotter.py — Qlobal Nəticə və Eigenfunction Vizuallaşdırıcısı

data/sweep_results.csv-dən Omega_re, Omega_im, Gamma vs (K və ya Delta) qrafiklərini qurur.
Hər bir Mode_ID (məs. Mode 0, Mode 1) fərqli rəngdə və etiketlə göstərilir.
Bütün qrafiklər plots/ qovluğuna saxlanılır.

Asılılıqlar: pandas, matplotlib
"""
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_global_results():
    csv_path = 'data/sweep_results.csv'
    if not os.path.exists(csv_path):
        print(f"XƏTA: {csv_path} tapılmadı. Əvvəlcə sweep hesablamasını tamamlayın.")
        return

    df = pd.read_csv(csv_path)
    
    # Detect sweep variable
    if df['Delta'].nunique() == 1:
        sweep_col = 'K'
        sweep_label = 'K (Wave Number)'
    elif df['K'].nunique() == 1:
        sweep_col = 'Delta'
        sweep_label = '$\\Delta$ (Delta)'
    else:
        sweep_col = 'Delta'
        sweep_label = '$\\Delta$'
    
    os.makedirs('plots', exist_ok=True)
    
    # Modlara görə qruplaşdıraq (Hər mod fərqli rəngdə olsun)
    modes = df['Mode_ID'].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(modes)))
    
    metrics = [
        ('Omega_re', 'Re($\\Omega$)', 'plots/Omega_re_vs_sweep.png', 'Real Part of $\\Omega$'),
        ('Omega_im', 'Im($\\Omega$)', 'plots/Omega_im_vs_sweep.png', 'Imaginary Part of $\\Omega$'),
        ('Gamma', '$\\Gamma$ (Growth Rate)', 'plots/Gamma_vs_sweep.png', 'Growth Rate $\\Gamma$')
    ]
    
    for col_name, y_label, save_path, title_prefix in metrics:
        plt.figure(figsize=(10, 6))
        
        for i, m_id in enumerate(sorted(modes)):
            subset = df[df['Mode_ID'] == m_id].sort_values(by=sweep_col)
            plt.plot(subset[sweep_col], subset[col_name], 'o-', 
                     color=colors[i], label=f'Mode {m_id}', markersize=5, linewidth=1.5)
        
        plt.xlabel(sweep_label)
        plt.ylabel(y_label)
        plt.title(f'{title_prefix} vs {sweep_label}')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.savefig(save_path, dpi=150)
        plt.close()
        print(f"  Qrafik hazır: {save_path}")

def plot_eigenfunctions():
    """Plots saved eigenfunctions from data/eigenfunctions/*.npz"""
    files = glob.glob('data/eigenfunctions/*.npz')
    if not files:
        return
        
    os.makedirs('plots/eigenfunctions', exist_ok=True)
    print("Eigenfunction qrafikləri çəkilir...")
    
    for f in files:
        data = np.load(f)
        tau = data['tau']
        y = data['y']
        delta = data['Delta']
        k_val = data['K']
        omega = data['Omega']
        mode = data['Mode']
        
        plt.figure(figsize=(10, 6))
        plt.plot(tau, np.abs(y), 'k-', label='|y| (Amplitude)')
        plt.plot(tau, np.real(y), 'b--', alpha=0.6, label='Re(y)')
        plt.plot(tau, np.imag(y), 'r--', alpha=0.6, label='Im(y)')
        
        plt.xlabel('$\\tau$')
        plt.ylabel('y($\\tau$)')
        plt.title(f'Mode {mode} | $\\Delta$={delta:.3f}, K={k_val:.3f}\n$\\Omega$ = {omega.real:.4f}+{omega.imag:.4f}i')
        plt.legend()
        plt.grid(True)
        
        # Basename for filename
        base = os.path.basename(f).replace('.npz', '.png')
        plt.savefig(f'plots/eigenfunctions/{base}')
        plt.close()

if __name__ == '__main__':
    plot_global_results()
    # plot_eigenfunctions() # Çox sayda fayl varsa, bunu manual açın


def plot_eigenfunctions():
    os.makedirs('plots/eigenfunctions', exist_ok=True)
    
    npz_files = glob.glob('data/eigenfunctions/*.npz')
    
    for file in npz_files:
        data = np.load(file)
        tau = data['tau']
        y = data['y']
        Delta = data['Delta']
        Omega = data['Omega']
        mode_id = data['Mode'] if 'Mode' in data else 0
        
        y_abs = np.abs(y)
        y_re = np.real(y)
        y_im = np.imag(y)
        
        plt.figure(figsize=(10, 6))
        
        plt.plot(tau, y_abs, 'k-', linewidth=2, label='|y(tau)| (Amplitude)')
        plt.plot(tau, y_re, 'b--', linewidth=1.5, label='Re(y)')
        plt.plot(tau, y_im, 'r--', linewidth=1.5, label='Im(y)')
        
        plt.xlabel('$\\tau$ (Spatial grid)')
        plt.ylabel('Eigenfunction y($\\tau$)')
        plt.title(f'Eigenfunction for $\Delta$ = {Delta:.3f} | Mode {mode_id} | $\Omega$ = {Omega:.4f}')
        plt.legend(loc='upper right')
        plt.grid(True)
        
        # Adjusted y-limits for the 1e-3 initial amplitude scale
        # This cuts off the numerical explosion at the right boundary so the true wave is visible
        plt.ylim(-0.005, 0.005)
        
        filename = os.path.basename(file).replace('.npz', '.png')
        plt.savefig(f'plots/eigenfunctions/{filename}')
        plt.close()
        
    print(f"Generated {len(npz_files)} eigenfunction plots in 'plots/eigenfunctions/'.")

if __name__ == '__main__':
    os.makedirs('plots', exist_ok=True)
    print("Generating visualizations...")
    plot_global_results()
    plot_eigenfunctions()
    print("Done!")
