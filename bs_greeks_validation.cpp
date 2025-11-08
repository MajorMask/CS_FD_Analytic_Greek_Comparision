/**
 * @file bs_greeks_validation.cpp
 * @brief Black-Scholes Greeks validation: Analytic vs FD vs Complex-Step
 * 
 * Uses provided header files:
 * - bs_call_price.h: Black-Scholes pricing with Phi_real, phi, and bs_price_call
 * - InverseCumulativeNormal.h: (Available but not needed for this assignment)
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <algorithm>

// Include the provided headers
#include "bs_call_price.h"
// #include "InverseCumulativeNormal.h"  // Not needed for this assignment

using namespace std;

// Type-dispatched Φ_t for complex-step differentiation

// Overload for double - uses Phi_real from bs_call_price.h
inline double Phi_t(double z) {
    return Phi_real(z);
}

// Overload for complex - first-order Taylor expansion
// Φ(z_r + i*z_i) ≈ Φ(z_r) + i*z_i*φ(z_r)

inline complex<double> Phi_t(const complex<double>& z) {
    double z_real = z.real();
    double z_imag = z.imag();
    double phi_real = Phi_real(z_real);
    double phi_derivative = phi(z_real);  // Φ'(z) = φ(z)
    return complex<double>(phi_real, z_imag * phi_derivative);
}

// Templated Black-Scholes call price for complex-step
    

template<class T>
T bs_price_call_t(T S, T K, T r, T q, T sigma, T Tmat) {
    const T DF = exp(-r * Tmat);
    const T F = S * exp((r - q) * Tmat);
    const T sigmaT = sigma * sqrt(Tmat);
    
    // For complex type, assume valid positive inputs
    T ln_F_over_K = log(F / K);
    
    const T d1 = (ln_F_over_K + T(0.5) * sigma * sigma * Tmat) / sigmaT;
    const T d2 = d1 - sigmaT;

    return DF * (F * Phi_t(d1) - K * Phi_t(d2));
}

// TASK 1: Analytic Greeks


struct AnalyticGreeks {
    double delta;
    double gamma;
};

AnalyticGreeks compute_analytic_greeks(double S, double K, double r, double q, 
                                       double sigma, double T) {
    AnalyticGreeks greeks;
    
    const double F = S * exp((r - q) * T);
    const double sigmaT = sigma * sqrt(max(T, 0.0));
    
    if (sigmaT < 1e-15) {
        // At expiry or zero vol, Greeks are discontinuous
        greeks.delta = (F > K) ? exp(-q * T) : 0.0;
        greeks.gamma = 0.0;
        return greeks;
    }

    double ln_F_over_K;
    if (K > 0.0) {
        const double x = (F - K) / K;
        ln_F_over_K = (abs(x) <= 1e-12) ? log1p(x) : log(F / K);
    } else {
        ln_F_over_K = log(F / K);
    }
    
    const double d1 = (ln_F_over_K + 0.5 * sigma * sigma * T) / sigmaT;
    
    // Delta = e^(-qT) * Φ(d1)
    greeks.delta = exp(-q * T) * Phi_real(d1);
    
    // Gamma = e^(-qT) * φ(d1) / (S * σ * √T)
    // Use log-space computation to avoid underflow: log(φ(d1)) = -d1²/2 - log(√(2π))
    const double log_phi_d1 = -0.5 * d1 * d1 - 0.5 * log(2.0 * M_PI);
    const double phi_d1 = exp(log_phi_d1);
    
    greeks.gamma = exp(-q * T) * phi_d1 / (S * sigmaT);
    
    return greeks;
}

// TASK 2: Forward Finite Difference Greeks

struct FDGreeks {
    double delta;
    double gamma;
};

FDGreeks compute_fd_greeks(double S, double K, double r, double q, 
                           double sigma, double T, double h) {
    FDGreeks greeks;
    
    // Use bs_price_call from bs_call_price.h
    double C_S = bs_price_call(S, K, r, q, sigma, T);
    double C_Sph = bs_price_call(S + h, K, r, q, sigma, T);
    double C_Sp2h = bs_price_call(S + 2.0*h, K, r, q, sigma, T);
    
    // Delta_fwd = (C(S+h) - C(S)) / h
    greeks.delta = (C_Sph - C_S) / h;
    
    // Gamma_fwd = (C(S+2h) - 2*C(S+h) + C(S)) / h²
    greeks.gamma = (C_Sp2h - 2.0*C_Sph + C_S) / (h * h);
    
    return greeks;
}

// TASK 3 & 4: Complex-Step Differentiation Greeks


struct CSGreeks {
    double delta;
    double gamma_real;  // Using real-part method
    double gamma_45;    // Using 45-degree method
};

CSGreeks compute_cs_greeks(double S, double K, double r, double q, 
                           double sigma, double T, double h) {
    CSGreeks greeks;
    
    // Convert all parameters to complex type (with zero imaginary part)
    complex<double> K_c(K, 0.0);
    complex<double> r_c(r, 0.0);
    complex<double> q_c(q, 0.0);
    complex<double> sigma_c(sigma, 0.0);
    complex<double> T_c(T, 0.0);
    
    // DELTA: Δ_cs = Im[C(S + ih)] / h
    complex<double> S_plus_ih(S, h);
    complex<double> C_complex = bs_price_call_t(S_plus_ih, K_c, r_c, q_c, sigma_c, T_c);
    greeks.delta = C_complex.imag() / h;
    
    // GAMMA (real-part method): Γ = -2 * [Re(C(S+ih)) - C(S)] / h²
    double C_S = bs_price_call(S, K, r, q, sigma, T);
    greeks.gamma_real = -2.0 * (C_complex.real() - C_S) / (h * h);
    
    // GAMMA (45-degree method): Γ = Im[C(S+hω) + C(S-hω)] / h²
    // where ω = e^(iπ/4) = (1+i)/√2
    const double sqrt2 = sqrt(2.0);
    complex<double> omega(1.0/sqrt2, 1.0/sqrt2);
    
    complex<double> S_plus_homega = S + h * omega;
    complex<double> S_minus_homega = S - h * omega;
    
    complex<double> C_plus = bs_price_call_t(S_plus_homega, K_c, r_c, q_c, sigma_c, T_c);
    complex<double> C_minus = bs_price_call_t(S_minus_homega, K_c, r_c, q_c, sigma_c, T_c);
    
    greeks.gamma_45 = (C_plus + C_minus).imag() / (h * h);
    
    return greeks;
}

// Validation Sweep


struct Scenario {
    string name;
    double S, K, r, q, sigma, T;
};

void run_validation_sweep(const Scenario& scenario, const string& output_file) {
    cout << "\n=== Running validation for " << scenario.name << " ===" << endl;
    cout << "S=" << scenario.S << ", K=" << scenario.K 
         << ", r=" << scenario.r << ", q=" << scenario.q 
         << ", σ=" << scenario.sigma << ", T=" << scenario.T << endl;
    
    // Compute analytic Greeks (ground truth)
    AnalyticGreeks analytic = compute_analytic_greeks(
        scenario.S, scenario.K, scenario.r, scenario.q, scenario.sigma, scenario.T
    );
    
    cout << "Analytic Delta = " << setprecision(15) << analytic.delta << endl;
    cout << "Analytic Gamma = " << setprecision(15) << analytic.gamma << endl;
    
    // Create logarithmic grid: h_rel from 10^-16 to 10^-4
    // Using 24 intervals = 25 points
    vector<double> h_rel_values;
    for (int i = 0; i <= 24; ++i) {
        double log_h_rel = -16.0 + i * (12.0 / 24.0);  // From -16 to -4
        h_rel_values.push_back(pow(10.0, log_h_rel));
    }
    
    // Open CSV file for output
    ofstream csv(output_file);
    csv << setprecision(16) << scientific;
    
    // Write header (exact format required by assignment)
    csv << "h_rel,h,"
        << "Delta_analytic,Delta_fd,Delta_cs,err_D_fd,err_D_cs,"
        << "Gamma_analytic,Gamma_fd,Gamma_cs_real,Gamma_cs_45,"
        << "err_G_fd,err_G_cs_real,err_G_cs_45\n";
    
    // Sweep over step sizes
    for (double h_rel : h_rel_values) {
        double h = h_rel * scenario.S;
        
        // Compute finite difference Greeks
        FDGreeks fd = compute_fd_greeks(
            scenario.S, scenario.K, scenario.r, scenario.q, scenario.sigma, scenario.T, h
        );
        
        // Compute complex-step Greeks
        CSGreeks cs = compute_cs_greeks(
            scenario.S, scenario.K, scenario.r, scenario.q, scenario.sigma, scenario.T, h
        );
        
        // Compute errors
        double err_D_fd = abs(fd.delta - analytic.delta);
        double err_D_cs = abs(cs.delta - analytic.delta);
        double err_G_fd = abs(fd.gamma - analytic.gamma);
        double err_G_cs_real = abs(cs.gamma_real - analytic.gamma);
        double err_G_cs_45 = abs(cs.gamma_45 - analytic.gamma);
        
        // Write to CSV
        csv << h_rel << "," << h << ","
            << analytic.delta << "," << fd.delta << "," << cs.delta << ","
            << err_D_fd << "," << err_D_cs << ","
            << analytic.gamma << "," << fd.gamma << "," 
            << cs.gamma_real << "," << cs.gamma_45 << ","
            << err_G_fd << "," << err_G_cs_real << "," << err_G_cs_45 << "\n";
    }
    
    csv.close();
    cout << "Results written to " << output_file << endl;
}

// Main Program

int main() {
    cout << setprecision(15);
    
    // Scenario 1: ATM reference (happy path)
    Scenario scenario1 {
        "Scenario 1 (ATM reference)",
        100.0,  // S
        100.0,  // K
        0.0,    // r
        0.0,    // q
        0.20,   // sigma
        1.0     // T
    };
    
    // Scenario 2: Near-expiry, low-vol, ATM (stress test)
    Scenario scenario2 {
        "Scenario 2 (Near-expiry, low-vol, ATM)",
        100.0,              // S
        100.0,              // K
        0.0,                // r
        0.0,                // q
        0.01,               // sigma
        1.0 / 365.0         // T = 1/365
    };
    
    // Run validation sweeps for both scenarios
    run_validation_sweep(scenario1, "bs_fd_vs_complex_scenario1.csv");
    run_validation_sweep(scenario2, "bs_fd_vs_complex_scenario2.csv");
    
    cout << "\n=== Validation Complete ===" << endl;
    cout << "\nGenerated files:" << endl;
    cout << "  - bs_fd_vs_complex_scenario1.csv" << endl;
    cout << "  - bs_fd_vs_complex_scenario2.csv" << endl;
    cout << "\nNext steps:" << endl;
    cout << "  Run: python3 analyze_results.py" << endl;
    cout << "  to generate plots and statistical analysis." << endl;
    
    return 0;
}
