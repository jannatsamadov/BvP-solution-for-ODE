import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# 1. CSV Faylının Oxunması
file_path = r"D:\solar wind\Numerical solutions\Main_antigravity\data\grid_complex_mismatch_Delta_1.500_K_0.100.csv"

try:
    # Datanı yükləyirik (Sütun adlarının Re(Omega), Im(Omega), Mismatch olduğunu fərz edərək)
    df = pd.read_csv(file_path)
    
    # Əgər sütun adlarında boşluqlar varsa təmizləyirik
    df.columns = df.columns.str.strip()
    
    # Sütunları ayırırıq
    re_cols = [c for c in df.columns if 'Re' in c or 'real' in c.lower()][0]
    im_cols = [c for c in df.columns if 'Im' in c or 'imag' in c.lower()][0]
    val_cols = [c for c in df.columns if 'Mismatch' in c or 'val' in c.lower()][0]
    
    real_vals = df[re_cols].values
    imag_vals = df[im_cols].values
    mismatch_vals = df[val_cols].values

    # Kompleks ədədlərin string kimi yazılma ehtimalına qarşı təmizləmə
    if mismatch_vals.dtype == object:
        mismatch_vals = mismatch_vals.astype(str)
        mismatch_vals = np.array([complex(s.replace(' ', '')) for s in mismatch_vals])

except Exception as e:
    print(f"Faylı oxuyarkən xəta baş verdi: {e}")
    exit()

# 2. 1D Datanın 2D Qridə Çevrilməsi (Pivot)
# Unik dəyərləri tapıb matris formasına salırıq
unique_real = np.unique(real_vals)
unique_imag = np.unique(imag_vals)

# Qridin ölçüləri
num_re = len(unique_real)
num_im = len(unique_imag)

# Datanın tam qrid formasına salınması
Real_Grid, Imag_Grid = np.meshgrid(unique_real, unique_imag)
Mismatch_Matrix = np.zeros(Real_Grid.shape, dtype=complex if np.iscomplexobj(mismatch_vals) else float)

# Datanı qridə doldururuq
for r, i, m in zip(real_vals, imag_vals, mismatch_vals):
    r_idx = np.where(unique_real == r)[0][0]
    i_idx = np.where(unique_imag == i)[0][0]
    Mismatch_Matrix[i_idx, r_idx] = m

# 3. Vizualların Hazırlanması (Elmi Keyfiyyətdə)
fig, axes = plt.subplots(1, 2, figsize=(16, 6), dpi=300)

# A Qrafiki: log10(|Mismatch|) - Amplituda xəritəsi
abs_matrix = np.abs(Mismatch_Matrix)
# Sıfıra getmə logaritmik xətanın qarşısını almaq üçün kiçik epsilon əlavə edirik
log_mismatch = np.log10(abs_matrix + 1e-15)

im1 = axes[0].pcolormesh(Real_Grid, Imag_Grid, log_mismatch, shading='auto', cmap='viridis_r')
axes[0].set_title(r'Amplituda Xəritəsi: $\log_{10}(|\text{Mismatch}|)$', fontsize=14, fontweight='bold')
axes[0].set_xlabel(r'$\text{Re}(\Omega_A)$', fontsize=12)
axes[0].set_ylabel(r'$\text{Im}(\Omega_A)$', fontsize=12)
axes[0].grid(True, linestyle='--', alpha=0.5)
fig.colorbar(im1, ax=axes[0], label=r'$\log_{10}(|\text{Mismatch}|)$')

# B Qrafiki: np.angle(Mismatch) - Faza Xəritəsi
# Faza burulğanlarını (vortices) görmək üçün ən ideal rəng palitrası 'twilight' və ya 'hsv'-dir
phase_matrix = np.angle(Mismatch_Matrix)

im2 = axes[1].pcolormesh(Real_Grid, Imag_Grid, phase_matrix, shading='auto', cmap='twilight')
axes[1].set_title(r'Faza Xəritəsi: $\text{Angle}(\text{Mismatch})$', fontsize=14, fontweight='bold')
axes[1].set_xlabel(r'$\text{Re}(\Omega_A)$', fontsize=12)
axes[1].set_ylabel(r'$\text{Im}(\Omega_A)$', fontsize=12)
axes[1].grid(True, linestyle='--', alpha=0.5)
fig.colorbar(im2, ax=axes[1], label='Faza (Radian)')

plt.suptitle(r'MHD Sərhəd Şərtləri Spektral Analizi ($\Delta=0.5, K=0.1$)', fontsize=16, y=1.02)
plt.tight_layout()

# Rəsmi yaddaşa vermək
output_filename = "Scientific_Visuals_Delta_0.5_K_0.1.png"
plt.savefig(output_filename, bbox_inches='tight')
print(f"Vizuallar uğurla yaradıldı və '{output_filename}' olaraq yaddaşa yazıldı.")