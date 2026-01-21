import numpy as np
from scipy.integrate import quad

class GCGModel:
    def __init__(self, h=0.6774, Omega_m=0.3089, alpha=0.2):
        self.Omega_m = Omega_m
        self.alpha = alpha
        # As is usually tuned so that rho_ude(a=1) = 1 - Omega_m
        self.As = 0.774 

    def rho_ude(self, a):
        term = self.As + (1 - self.As) * a**(-3 * (1 + self.alpha))
        return term**(1 / (1 + self.alpha))

    def w_ude(self, a):
        return -self.As / (self.As + (1 - self.As) * a**(-3 * (1 + self.alpha)))

    def cs2(self, a):
        return -self.alpha * self.w_ude(a)

class LogotropicModel:
    def __init__(self, h=0.6774, Omega_m=0.3089, B=0.01):
        self.Omega_m = Omega_m
        self.B = B
        self.Omega_de0 = 1.0 - Omega_m

    def rho_ude(self, a):
        return self.Omega_m * a**(-3) + self.Omega_de0 * (1 - 3 * self.B * np.log(a))

    def w_ude(self, a):
        rho = self.rho_ude(a)
        p = -self.Omega_de0 * (1 + self.B - 3 * self.B * np.log(a))
        return p / rho

    def cs2(self, a):
        # Based on your notebook: B / (rho/rho_crit - 1)
        return self.B / (self.rho_ude(a) - 1 + 1e-10)
    
class MurnaghanModel:
    def __init__(self, h=0.6774, Omega_m=0.3089, alpha=0.018):
        self.h = h
        self.Omega_m = Omega_m
        self.alpha = alpha
        self.rho_star_ratio = 1e120 # Normalized Planck density ratio
        
        # NORMALIZATION LOGIC:
        # We solve for A_star such that total density today (a=1) is exactly 1.0
        # epsilon_M(1) = Omega_m + A_star * [ (1/(1+alpha)) * (Omega_m/rho_star)^-alpha - 1 ] = 1
        target_de = 1.0 - Omega_m
        bracket = (1.0 / (1.0 + self.alpha)) * (self.Omega_m / self.rho_star_ratio)**(-self.alpha) - 1.0
        
        # This defines the attribute the functions are looking for!
        self.A_star = target_de / bracket 

    def rho_m(self, a):
        return self.Omega_m * a**(-3)

    def epsilon_M(self, a):
        """Total Energy Density"""
        rho = self.rho_m(a)
        term = self.A_star * ((1/(1+self.alpha)) * (rho / self.rho_star_ratio)**(-self.alpha) - 1)
        return rho + term

    def p_murn(self, a):
        """Pressure"""
        rho = self.rho_m(a)
        return -self.A_star * ((rho / self.rho_star_ratio)**(-self.alpha) - 1)

    def w_ude(self, a):
        """Equation of State"""
        return self.p_murn(a) / self.epsilon_M(a)

    def cs2(self, a):
        """Sound Speed Squared"""
        rho = self.rho_m(a)
        # s_M = (A / rho_M) * (rho_m / rho_pl)^-alpha
        # Note: In your text A = A_star * alpha
        return (self.alpha * self.A_star / self.epsilon_M(a)) * (rho / self.rho_star_ratio)**(-self.alpha)
    def get_kinematics(self):
        """Calculates q0 and j0 at a=1 (today)."""
        # We use a small step in e-folds (N) for numerical differentiation
        # N = ln(a), so today is N=0
        dN = 0.0001
        
        def ln_H(N):
            a = np.exp(N)
            # H^2 proportional to rho_M
            return 0.5 * np.log(self.epsilon_M(a))
            
        # First derivative: d(ln H)/dN
        f_plus = ln_H(dN)
        f_minus = ln_H(-dN)
        f_prime = (f_plus - f_minus) / (2 * dN)
        
        # Second derivative: d^2(ln H)/dN^2
        f_0 = ln_H(0)
        f_double_prime = (f_plus - 2*f_0 + f_minus) / (dN**2)
        
        # Physics definitions from your text:
        # q = -1 - d(ln H)/dN
        q0 = -1 - f_prime
        
        # j = q(2q + 1) + dq/dN
        # After some algebra, this relates to the second derivative:
        j0 = f_double_prime + (f_prime)**2 + 3*f_prime + 1
        
        return q0, j0
    
    def calculate_distance_modulus(self, z_array):
        c = 299792.458 # Speed of light in km/s
        mu_list = []
    
        for z in z_array:
            # Integral of 1/E(z)
            integral, _ = quad(lambda zp: 1.0 / self.E(zp), 0, z)
            dL = (1 + z) * (c / self.H0) * integral
        
            # Convert dL from Mpc to pc (multiply by 10^6) for the log
            mu = 5 * np.log10(dL * 1e6 / 10)
            mu_list.append(mu)
        
        return np.array(mu_list)