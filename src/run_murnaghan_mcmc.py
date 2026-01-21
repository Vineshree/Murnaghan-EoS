import numpy as np
import emcee
import corner
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad

# ==========================================
# 1. PHYSICS: Murnaghan Cosmology
# ==========================================
class MurnaghanCosmology:
    def __init__(self, h0, omega_m, alpha):
        self.H0 = h0 * 100
        self.omega_m = omega_m
        self.alpha = alpha
        self.rho_pl_ratio = 1e122 
        
        # Solving for A_star based on Omega_m and alpha (Eq. 32c)
        target_de = 1.0 - self.omega_m
        bracket = (1.0 / (1.0 + alpha)) * (self.omega_m / self.rho_pl_ratio)**(-alpha) - 1.0
        self.A_star = target_de / bracket

    def E(self, z):
        """Dimensionless Hubble parameter E(z) = H(z)/H0"""
        a = 1.0 / (1.0 + z)
        rho_m = self.omega_m * a**-3
        term_de = (self.A_star / self.alpha) * (
            (1 / (1 + self.alpha)) * (rho_m / self.rho_pl_ratio)**(-self.alpha) - 1
        )
        # Avoid numerical issues with sqrt of negative numbers during MCMC exploration
        return np.sqrt(np.maximum(rho_m + term_de, 1e-10))

    def calculate_distance_modulus(self, z_array):
        c = 299792.458 # km/s
    
        # 1. Create a fine grid of redshifts up to the max z in your data
        z_max = np.max(z_array)
        z_grid = np.linspace(0, z_max, 100) # 100 points is plenty for smoothness
    
        # 2. Compute the integral for each point in the grid
        integrals = [quad(lambda zp: 1.0 / self.E(zp), 0, z)[0] for z in z_grid]
    
        # 3. Create an interpolator
        # This allows you to get the integral for any z instantly
        interp_func = interp1d(z_grid, integrals, kind='cubic')
    
        # 4. Apply to your actual data vector (instantaneous vectorization)
        integral_at_z = interp_func(z_array)
        dL = (1 + z_array) * (c / self.H0) * integral_at_z
    
        # 5. Calculate mu using the standard cosmology constant
        # Use np.maximum to avoid log of zero errors during random MCMC jumps
        mu = 5 * np.log10(np.maximum(dL, 1e-10)) + 25
    
        return mu

# ==========================================
# 2. DATA LOADING
# ==========================================
# Load OHD
# Columns: z, H, error
z_ohd, h_obs, h_err = np.loadtxt('OHD_3.txt', unpack=True)

# Load Pantheon+ using Pandas for speed and header handling
df_sn = pd.read_csv('Pantheon+SH0ES.txt', sep=' ')
df_sn = df_sn[df_sn['USED_IN_SH0ES_HF'] == 1] # Select Hubble Flow sample
z_sn = df_sn['zHD'].values
mu_obs = df_sn['MU_SH0ES'].values
mu_err = df_sn['MU_SH0ES_ERR_DIAG'].values

# ==========================================
# 3. STATISTICAL ENGINE (Likelihoods)
# ==========================================
def log_prior(theta):
    h0, om, alpha = theta
    if 0.6 < h0 < 0.8 and 0.2 < om < 0.4 and 0.0 < alpha < 0.3:
        return 0.0
    return -np.inf

def log_likelihood_ohd(theta):
    model = MurnaghanCosmology(*theta)
    h_th = model.H0 * model.E(z_ohd)
    return -0.5 * np.sum(((h_obs - h_th) / h_err)**2)

def log_likelihood_sn(theta):
    model = MurnaghanCosmology(*theta)
    mu_th = model.calculate_distance_modulus(z_sn)
    return -0.5 * np.sum(((mu_obs - mu_th) / mu_err)**2)

def log_probability(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    # Total Log Likelihood (Eq. 18 in your paper)
    return lp + log_likelihood_ohd(theta) + log_likelihood_sn(theta)

# ==========================================
# 4. RUNNER
# ==========================================
n_dim, n_walkers = 3, 32
p0 = [0.70, 0.30, 0.10] + 1e-4 * np.random.randn(n_walkers, n_dim)

sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_probability)
print("Starting MCMC Sampling...")
sampler.run_mcmc(p0, 1000, progress=True) # Start with 1000 for a test run

# Save results for your Notebook
samples = sampler.get_chain(discard=200, thin=15, flat=True)
np.save('murnaghan_combined_samples.npy', samples)

# Quick Triangle Plot

labels = [r"$h_0$", r"$\Omega_m$", r"$\alpha$"]
fig = corner.corner(samples, labels=labels, show_titles=True)
plt.show()

# Calculate Information Criteria
max_log_prob = np.max(sampler.get_log_prob())
n_data = len(z_ohd) + len(z_sn)
p = n_dim

aic = -2 * max_log_prob + 2 * p
bic = -2 * max_log_prob + p * np.log(n_data)

print(f"Final Statistics for Murnaghan Model:")
print(f"Maximum Log-Likelihood: {max_log_prob:.2f}")
print(f"AIC: {aic:.2f}")
print(f"BIC: {bic:.2f}")