# Black-Scholes Greeks Validation - Quick Summary

## Test Scenarios

| Scenario  | S   | K   | r | q | σ    | T       | Δ (analytic) | Γ (analytic) | 
| 1 (ATM)   | 100 | 100 | 0 | 0 | 0.20 | 1.0     | 0.5398       | 0.0198       |
| 2 (Stress)| 100 | 100 | 0 | 0 | 0.01 | 1/365   | 0.5001       | 7.6218       |

## Best Accuracy Achieved

### Delta (Δ)

| Method | Scenario 1 | Scenario 2     | Winner  |
| FD     | 2.7×10⁻¹⁰  | 5.7×10⁻⁸       |         |
| CS     | **≈ 0**    | **7.2×10⁻¹⁵**  | **      |

### Gamma (Γ)

| Method      | Scenario 1   | Scenario 2   | Winner |
| FD          | 2.5×10⁻⁷     | 2.5×10⁻⁶     |        |
| CS (real)   | 2.0×10⁻²     |  6.59        |        |
| CS (45°)    | **1.9×10⁻¹²**| **5.1×10⁻⁹** | **     |

## Recommended Settings

### Production Use - Delta

```cpp
// Complex-step method (recommended)
double h = 1e-8 * S;  // h_rel = 1e-8
Delta = imag(C(S + i*h)) / h;

// Expected accuracy: ~10⁻¹⁵ (machine epsilon)
// Robust across wide range of h
```

### Production Use - Gamma

```cpp
// Complex-step 45° method (recommended)
double h = 1e-6 * S;  // h_rel = 1e-6 for scenario 1
// OR
double h = 1e-8 * S;  // h_rel = 1e-8 for scenario 2

complex<double> omega = (1.0 + i) / sqrt(2.0);
Gamma = imag(C(S + h*omega) + C(S - h*omega)) / (h*h);

// Expected accuracy: 10⁻⁹ to 10⁻¹²
```

## Key Takeaways

**Complex-step is superior** for both Delta and Gamma
**Delta**: Machine precision achievable with CS, only 10⁻⁸ to 10⁻¹⁰ with FD
**Gamma**: CS 45° method is 10³ to 10⁵ times more accurate than FD
**Robustness**: CS works across 8+ orders of magnitude of step sizes
**Near-expiry**: Both methods degrade but CS maintains advantage
**Avoid**: CS real-part method for Gamma (poor performance)

## Implementation Notes

1. **Templated pricing function** enables complex-step with single codebase
2. **Type-dispatched Φ(z)** handles both real and complex arguments
3. **Taylor expansion** for complex CDF avoids full complex error function
4. **Analytic Greeks** always preferred when available (exact, no approximation error)

## File Manifest

- `bs_greeks.cpp` - Complete implementation
- `bs_fd_vs_complex_scenario1.csv` - Scenario 1 validation data
- `bs_fd_vs_complex_scenario2.csv` - Scenario 2 validation data  
- `greeks_error_analysis.png` - Error plots
- `DESIGN.md` - Detailed analysis and methodology
- `SUMMARY.md` - This quick reference
- `analyze_results.py` - Visualization creator
