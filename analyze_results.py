import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV files
df1 = pd.read_csv('bs_fd_vs_complex_scenario1.csv')
df2 = pd.read_csv('bs_fd_vs_complex_scenario2.csv')

# Create figure with subplots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Black-Scholes Greeks: Error Analysis vs Step Size', fontsize=16, fontweight='bold')

# Scenario 1 - Delta errors
ax = axes[0, 0]
ax.loglog(df1['h_rel'], df1['err_D_fd'], 'o-', label='FD Delta', alpha=0.7)
ax.loglog(df1['h_rel'], df1['err_D_cs'], 's-', label='CS Delta', alpha=0.7)
ax.set_xlabel('Relative Step Size (h_rel)', fontsize=10)
ax.set_ylabel('Absolute Error', fontsize=10)
ax.set_title('Scenario 1 (ATM): Delta Errors', fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3, which='both')
ax.legend(fontsize=9)

# Scenario 1 - Gamma errors
ax = axes[0, 1]
ax.loglog(df1['h_rel'], df1['err_G_fd'], 'o-', label='FD Gamma', alpha=0.7)
ax.loglog(df1['h_rel'], df1['err_G_cs_real'], 's-', label='CS Gamma (real)', alpha=0.7)
ax.loglog(df1['h_rel'], df1['err_G_cs_45'], '^-', label='CS Gamma (45°)', alpha=0.7)
ax.set_xlabel('Relative Step Size (h_rel)', fontsize=10)
ax.set_ylabel('Absolute Error', fontsize=10)
ax.set_title('Scenario 1 (ATM): Gamma Errors', fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3, which='both')
ax.legend(fontsize=9)

# Scenario 2 - Delta errors
ax = axes[1, 0]
ax.loglog(df2['h_rel'], df2['err_D_fd'], 'o-', label='FD Delta', alpha=0.7)
ax.loglog(df2['h_rel'], df2['err_D_cs'], 's-', label='CS Delta', alpha=0.7)
ax.set_xlabel('Relative Step Size (h_rel)', fontsize=10)
ax.set_ylabel('Absolute Error', fontsize=10)
ax.set_title('Scenario 2 (Near-expiry, low-vol): Delta Errors', fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3, which='both')
ax.legend(fontsize=9)

# Scenario 2 - Gamma errors
ax = axes[1, 1]
ax.loglog(df2['h_rel'], df2['err_G_fd'], 'o-', label='FD Gamma', alpha=0.7)
ax.loglog(df2['h_rel'], df2['err_G_cs_real'], 's-', label='CS Gamma (real)', alpha=0.7)
ax.loglog(df2['h_rel'], df2['err_G_cs_45'], '^-', label='CS Gamma (45°)', alpha=0.7)
ax.set_xlabel('Relative Step Size (h_rel)', fontsize=10)
ax.set_ylabel('Absolute Error', fontsize=10)
ax.set_title('Scenario 2 (Near-expiry, low-vol): Gamma Errors', fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3, which='both')
ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig('greeks_error_analysis.png', dpi=300, bbox_inches='tight')
print("Plot saved as greeks_error_analysis.png")

# Statistical summary
print("\n" + "="*80)
print("STATISTICAL SUMMARY")
print("="*80)

for scenario_num, df in enumerate([df1, df2], 1):
    print(f"\nScenario {scenario_num}:")
    print("-" * 40)
    
    # Find optimal step size ranges
    print("\nDelta Methods:")
    fd_delta_min_idx = df['err_D_fd'].idxmin()
    cs_delta_min_idx = df['err_D_cs'].idxmin()
    print(f"  FD:  Min error = {df.loc[fd_delta_min_idx, 'err_D_fd']:.2e} at h_rel = {df.loc[fd_delta_min_idx, 'h_rel']:.2e}")
    print(f"  CS:  Min error = {df.loc[cs_delta_min_idx, 'err_D_cs']:.2e} at h_rel = {df.loc[cs_delta_min_idx, 'h_rel']:.2e}")
    
    print("\nGamma Methods:")
    fd_gamma_min_idx = df['err_G_fd'].idxmin()
    cs_real_gamma_min_idx = df['err_G_cs_real'].idxmin()
    cs_45_gamma_min_idx = df['err_G_cs_45'].idxmin()
    print(f"  FD:         Min error = {df.loc[fd_gamma_min_idx, 'err_G_fd']:.2e} at h_rel = {df.loc[fd_gamma_min_idx, 'h_rel']:.2e}")
    print(f"  CS (real):  Min error = {df.loc[cs_real_gamma_min_idx, 'err_G_cs_real']:.2e} at h_rel = {df.loc[cs_real_gamma_min_idx, 'h_rel']:.2e}")
    print(f"  CS (45°):   Min error = {df.loc[cs_45_gamma_min_idx, 'err_G_cs_45']:.2e} at h_rel = {df.loc[cs_45_gamma_min_idx, 'h_rel']:.2e}")
    
    # Best achievable accuracy across all step sizes
    print("\nBest Achievable Accuracy:")
    print(f"  Delta FD:     {df['err_D_fd'].min():.2e}")
    print(f"  Delta CS:     {df['err_D_cs'].min():.2e}")
    print(f"  Gamma FD:     {df['err_G_fd'].min():.2e}")
    print(f"  Gamma CS-real:{df['err_G_cs_real'].min():.2e}")
    print(f"  Gamma CS-45°: {df['err_G_cs_45'].min():.2e}")

print("\n" + "="*80)
