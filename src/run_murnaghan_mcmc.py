import numpy as np
import emcee
import corner
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d

# --- Murnaghan Physics ---
class MurnaghanCosmology:
    def __init__(self, h0, omega_m, alpha):
        self.H0 = h0 * 100
        self.omega_m = omega_m
        self.alpha = alpha
        # We set a 'normalized' ratio to avoid the 10^122 issue
        # This keeps the math within Python's floating point limits
        self.ratio_norm = 1.0 
        
        # Solving for A_star (Eq. 32c)
        target_de = 1.0 - self.omega_m
        # Using a safer power law approach
        bracket = (1.0 / (1.0 + alpha)) * (self.omega_m)**(-alpha) - 1.0
        self.A_star = target_de / bracket

    def E(self, z):
        a = 1.0 / (1.0 + z)
        rho_m = self.omega_m * a**-3
        
        # Murnaghan Term (Stable version)
        # term = (A*/alpha) * [ (1/(1+alpha)) * (rho_m)^(-alpha) - 1 ]
        term_de = (self.A_star / self.alpha) * (
            (1.0 / (1.0 + self.alpha)) * (rho_m**(-self.alpha)) - 1.0
        )
        
        # If term_de becomes NaN or negative, the MCMC is in a "bad" zone
        return np.sqrt(np.maximum(rho_m + term_de, 1e-10))

    def calculate_distance_modulus(self, z_array):
        c = 299792.458 
        z_grid = np.linspace(0, np.max(z_array), 50)
        integrals = [quad(lambda zp: 1.0 / self.E(zp), 0, z)[0] for z in z_grid]
        interp_func = interp1d(z_grid, integrals, kind='cubic')
        dL = (1 + z_array) * (c / self.H0) * interp_func(z_array)
        return 5 * np.log10(np.maximum(dL, 1e-10)) + 25

# --- Likelihoods ---
def log_probability(theta, z_ohd, h_obs, h_err, z_sn, mu_obs, mu_err):
    h0, om, alpha = theta
    if not (0.6 < h0 < 0.8 and 0.2 < om < 0.4 and 0.0 < alpha < 0.3):
        return -np.inf
    
    model = MurnaghanCosmology(h0, om, alpha)
    # H(z) Likelihood
    lp_ohd = -0.5 * np.sum(((h_obs - model.H0 * model.E(z_ohd)) / h_err)**2)
    # Pantheon Likelihood
    lp_sn = -0.5 * np.sum(((mu_obs - model.calculate_distance_modulus(z_sn)) / mu_err)**2)
    return lp_ohd + lp_sn
def log_prior(theta):
    h0, om, alpha = theta
    # Keep alpha strictly positive and small to avoid the 1/alpha singularity
    if 0.6 < h0 < 0.8 and 0.1 < om < 0.5 and 0.001 < alpha < 0.3:
        return 0.0
    return -np.inf

# --- Execution ---
if __name__ == "__main__":
    # Load Data
    z_ohd, h_obs, h_err = np.loadtxt('data/OHD_31.txt', unpack=True)
    df_sn = pd.read_csv('data/Pantheon+SH0ES.txt', sep=' ')
    df_sn = df_sn[df_sn['USED_IN_SH0ES_HF'] == 1]
    z_sn, mu_obs, mu_err = df_sn['zHD'].values, df_sn['MU_SH0ES'].values, df_sn['MU_SH0ES_ERR_DIAG'].values

    n_dim, n_walkers = 3, 32
    # Spread starting points
    p0 = np.array([
        np.random.uniform(0.67, 0.73, n_walkers),
        np.random.uniform(0.25, 0.35, n_walkers),
        np.random.uniform(0.05, 0.15, n_walkers)
    ]).T

    sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_probability, 
                                    args=(z_ohd, h_obs, h_err, z_sn, mu_obs, mu_err))

    print("Running MCMC...")
    sampler.run_mcmc(p0, 2000, progress=True)
    # Use more points for the plot
    # Discard the first 500 steps (burn-in), but don't thin the rest
    samples = sampler.get_chain(discard=300, flat=True)
    print(f"Total points for plotting: {len(samples)}")
    
    # Generate Figure
    fig = corner.corner(samples, labels=[r"$h_0$", r"$\Omega_m$", r"$\alpha$"], 
                        show_titles=True, 
                        plot_datapoints=True, # Show the points if contours fail
                        levels=(0.68, 0.95))
    plt.savefig("Murnaghan_Contours.png")
    plt.show()