# Black‑Scholes Greeks: Build & Run Guide

This file focuses on how to build, run and inspect results for the repository.

## Prerequisites (macOS)
- C++17 toolchain (clang++ or g++)
- make
- Python 3 (optional, for plotting)
  - Python packages: numpy, pandas, matplotlib (install via pip if needed)

Install Python deps:
```bash
python3 -m pip install --user numpy pandas matplotlib
```

## Quick Build & Run (recommended)
From repository root (/Users/mananaggarwal/Downloads/files-3):

```bash
# Compile
make

# Run validation (produces CSVs and greeks_error_analysis.png)
make run

# Generate / refresh plots from CSV (if needed)
make analyze

# Clean build artifacts and outputs
make clean
```

make targets:
- `make` — compile binary (C++17, -O3)
- `make run` — run validation program and write CSV output:
  - bs_fd_vs_complex_scenario1.csv
  - bs_fd_vs_complex_scenario2.csv
  - greeks_error_analysis.png
- `make analyze` — run analyze_results.py to generate/refresh plots
- `make clean` — remove binaries and generated files

## Manual (no Makefile)
Compile and run manually:

```bash
# Compile (example using g++)
g++ -std=c++17 -O3 -o bs_greeks_validation bs_greeks_validation.cpp

# Run program
./bs_greeks_validation

# Generate plots with Python
python3 analyze_results.py
```

Open the generated plot on macOS:
```bash
open greeks_error_analysis.png
```

## Files produced
- bs_fd_vs_complex_scenario1.csv — scenario 1 validation data
- bs_fd_vs_complex_scenario2.csv — scenario 2 validation data
- greeks_error_analysis.png — combined error plots
- bs_greeks_validation (binary)

## Notes
- Ensure `bs_call_price.h` is present in the same directory when building.
- The Python script expects the CSVs produced by the run step.
- For reproducibility, use the provided Makefile; manual steps match the Makefile's behavior.

# Black‑Scholes Greeks — Consolidated Validation Study

This document consolidates the core findings, tables, short explanations, and plotting scripts from the Black‑Scholes Greeks validation comparing Analytic, Forward Finite Differences (FD) and Complex‑Step (CS) methods for Delta (Δ) and Gamma (Γ).

## 1. Overview
Three methods were validated:
- Analytic (closed form) — ground truth.
- Forward finite differences (FD) — classical numerical differentiation.
- Complex‑step differentiation (CS) — uses complex arithmetic to avoid catastrophic cancellation.  
Major conclusion: CS (including the 45° CS method for second derivatives) is substantially more accurate and robust than FD in all tested scenarios.

## 2. Test Scenarios
| Scenario | S | K | r | q | σ | T | Δ (analytic) | Γ (analytic) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 (ATM) | 100 | 100 | 0 | 0 | 0.20 | 1.0 | 0.5398 | 0.0198 |
| 2 (Near‑expiry, stress) | 100 | 100 | 0 | 0 | 0.01 | 1/365 | 0.5001 | 7.6218 |

Short note: Scenario 2 is ill‑conditioned (low vol + short time), amplifying Gamma and numerical sensitivity.

## 3. Best Achieved Accuracy (summary tables)

Delta (Δ) best absolute errors:
| Method | Scenario 1 | Scenario 2 | Winner |
|---|---:|---:|:---:|
| FD | 2.7×10⁻¹⁰ | 5.7×10⁻⁸ | |
| CS | ≈ 1×10⁻¹⁵ (machine eps) | 7.2×10⁻¹⁵ | CS ✓ |

Gamma (Γ) best absolute errors:
| Method | Scenario 1 | Scenario 2 | Winner |
|---|---:|---:|:---:|
| FD | 2.5×10⁻⁷ | 2.5×10⁻⁶ | |
| CS (real part) | 2.0×10⁻² | 6.59 | (poor) |
| CS (45°) | 1.9×10⁻¹² | 5.1×10⁻⁹ | CS (45°) ✓ |

Short explanation: CS for Delta extracts the imaginary part and avoids cancellation; the 45° CS variant for Gamma attains O(h⁴) truncation and far better accuracy than FD.

## 4. Recommended Settings (production)
Delta:
- Method: Complex‑step
- Step: h_rel = 1e‑8 (h = h_rel * S)
- Formula: Δ ≈ Im[C(S + i·h)] / h
- Expected accuracy: ~machine epsilon (~1e‑15)

Gamma:
- Method: Complex‑step 45° variant
- Step: h_rel = 1e‑6 (scenario dependent; 1e‑8 for stress)
- ω = (1 + i)/√2
- Formula: Γ ≈ Im[C(S + h·ω) + C(S − h·ω)] / h²
- Expected accuracy: 1e‑9 to 1e‑12

If complex arithmetic is unavailable, use FD with carefully tuned h_rel (≈1e‑8 for Δ, ≈1e‑6 for Γ); expect substantially worse accuracy.

## 5. CSV / Data format (used for analysis)
Each CSV contains 26 rows: header + 25 h values (log spaced).
Header columns:
h_rel,h,
Delta_analytic,Delta_fd,Delta_cs,err_D_fd,err_D_cs,
Gamma_analytic,Gamma_fd,Gamma_cs_real,Gamma_cs_45,
err_G_fd,err_G_cs_real,err_G_cs_45

Files:
- bs_fd_vs_complex_scenario1.csv
- bs_fd_vs_complex_scenario2.csv

## 6. Plotting (generated the requested plots using `analyze_results.py`)
Below is a compact Python script to load the CSVs and create the 4‑panel error plot used in the study. Save as `analyze_results.py` and run with Python 3 (requires numpy, pandas, matplotlib).
![Plot](https://github.com/MajorMask/CS_FD_Analytic_Greek_Comparision/blob/main/analyze_results.py)
