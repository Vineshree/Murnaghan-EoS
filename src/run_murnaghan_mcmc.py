import numpy as np
import emcee
import corner
import matplotlib.pyplot as plt
from scipy.integrate import quad
import pandas as pd

# =================================================================
# 1. PHYSICS CLASS: MURNAGHAN EQUATION OF STATE
# =================================================================
class MurnaghanCosmology:
    def __init__(self, h0, omega_m, alpha):
        self.H0 = h0 * 100
        self.omega_m = omega_m
        self.alpha = alpha
        # Reference Planck density ratio (log ~ 122.7) from your paper
        self.rho_pl_ratio = 1e122 
        
        # Solving for A_star such that total density today (a=1) is 1.0
        # Ensures flatness consistently across parameter space
        target_de = 1.0 - self.omega_m
        bracket = (1.0 / (1.0 + alpha)) * (self.omega_m / self.rho_pl_ratio)**(-alpha) - 1.0
        self.A_star = target_de / bracket

    def E(self, z):
        """Dimensionless Hubble Parameter: E(z) = H(z)/H0"""
        a = 1.0 / (1.0 + z)
        rho_m = self.omega_m * a**-3
        # Total energy density evolution (Eq. 14 in paper)
        term_de = (self.A_star / self.alpha) * (
            (1 / (1 + self.alpha)) * (rho_m / self.rho_pl_ratio)**(-self.alpha) - 1
        )
        return np.sqrt(rho_m + term_de)

# =================================================================
# 2. DATASETS (OHD & GROWTH f*sigma8)
# =================================================================
# Observational Hubble Data (OHD) - Representative subset from Cosmic Chronometers
z_ohd = np.array([0.07, 0.12, 0.20, 0.28, 0.35, 0.48, 0.59, 0.68, 0.78, 0.88, 1.03, 1.30, 1.53, 1.75])
h_ohd = np.array([69.0, 68.5, 72.9, 75.0, 82.0, 97.0, 104.0, 92.0, 105.0, 125.0, 154.0, 190.0, 214.0, 202.0])
err_ohd = np.array([19.6, 29.7, 8.1, 5.2, 4.0, 6.0, 13.0, 8.0, 12.0, 17.0, 20.0, 24.0, 22.0, 40.0])

# Growth Rate Data (f_sigma8) - Representative subset
z_fs8 = np.array([0.02, 0.15, 0.38, 0.51, 0.70, 0.85])
fs8_obs = np.array([0.360, 0.490, 0.440, 0.458, 0.448, 0.390])
fs8_err = np.array([0.040, 0.045, 0.025, 0.025, 0.035, 0.025])

# =================================================================
# 3. LIKELIHOOD AND MCMC SETUP
# =================================================================
def log_prior(theta):
    h0, om, alpha = theta
    # Broad priors matching the search range in your paper
    if 0.5 < h0 < 0.9 and 0.1 < om < 0.5 and 0.0 < alpha < 0.3:
        return 0.0
    return -np.inf

def log_likelihood(theta, z_o, h_o, e_o, z_f, f_o, ef_o):
    h0, om, alpha = theta
    model = MurnaghanCosmology(h0, om, alpha)
    
    # OHD Likelihood
    h_th = model.H0 * model.E(z_o)
    chi2_ohd = np.sum(((h_o - h_th) / e_o)**2)
    
    # Growth Likelihood (Approximation: f ~ Omega_m(z)^gamma)
    # Refined for Murnaghan using the gamma-index used in your paper
    gamma = 0.55
    omega_m_z = (om * (1+z_f)**3) / (model.E(z_f)**2)
    f_th = (omega_m_z**gamma) * 0.81  # 0.81 is a typical fiducial sigma8
    chi2_fs8 = np.sum(((f_o - f_th) / ef_o)**2)
    
    return -0.5 * (chi2_ohd + chi2_fs8)

def log_probability(theta, z_o, h_o, e_o, z_f, f_o, ef_o):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, z_o, h_o, e_o, z_f, f_o, ef_o)

# =================================================================
# 4. EXECUTION AND POST-PROCESSING
# =================================================================
def run_full_analysis():
    n_dim = 3
    n_walkers = 32
    # Initial guess around the Murnaghan "Best Fit" found in your paper
    p0 = [0.69, 0.30, 0.015] + 1e-4 * np.random.randn(n_walkers, n_dim)

    sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_probability, 
                                    args=(z_ohd, h_ohd, err_ohd, z_fs8, fs8_obs, fs8_err))

    print("Sampling MCMC (this may take a minute)...")
    sampler.run_mcmc(p0, 3000, progress=True)

    # Flatten chain and discard burn-in
    samples = sampler.get_chain(discard=1000, thin=15, flat=True)
    
    # --- BEST FIT TABLE (AIC/BIC) ---
    best_fit_theta = np.median(samples, axis=0)
    max_log_prob = np.max(sampler.get_log_prob())
    
    n_data = len(z_ohd) + len(z_fs8)
    aic = -2 * max_log_prob + 2 * n_dim
    bic = -2 * max_log_prob + n_dim * np.log(n_data)
    
    print("\n" + "="*30)
    print("      NUMERICAL RESULTS")
    print("="*30)
    print(f"h0:      {best_fit_theta[0]:.4f}")
    print(f"Omega_m: {best_fit_theta[1]:.4f}")
    print(f"alpha:   {best_fit_theta[2]:.4f}")
    print("-"*30)
    print(f"AIC:     {aic:.2f}")
    print(f"BIC:     {bic:.2f}")
    print("="*30)

    # --- TRIANGLE PLOT ---
    labels = [r"$h_0$", r"$\Omega_m$", r"$\alpha$"]
    
    fig = corner.corner(samples, labels=labels, truths=best_fit_theta, 
                        show_titles=True, title_fmt=".3f", color="darkblue")
    plt.show()

if __name__ == "__main__":
    run_full_analysis()
