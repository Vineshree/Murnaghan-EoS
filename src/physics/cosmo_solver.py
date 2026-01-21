from scipy.integrate import odeint
import numpy as np

class CosmoSolver:
    def __init__(self, model):
        self.model = model

    def growth_equations(self, var, N, alpha):
        delta, d_delta = var
        a = np.exp(N)
        
        # Pull physics from the model
        cs2 = self.model.sound_speed_sq(a)
        w = self.model.w_chap(a)
        
        # Calculate derived quantities for the ODE
        # These reflect the equations in your eqs() function
        H_prime = -1.5 * (1 + w) # Simplified for UDE dominated era
        Omega_ude = 1.0 / (1 + (self.model.Omega_m * a**-3 / self.model.rho_chap(a)))
        f_term = 1.5 * (1 + 3*cs2) * Omega_ude
        
        # The growth ODE logic
        dd_delta = (-d_delta * (3*cs2 + 2 + H_prime) - 
                    delta * (0 + 3*cs2*(2 + H_prime) - f_term)) # s_prime approx 0 here
        return [d_delta, dd_delta]

    def solve_growth(self, N_array, alpha, initial_delta=0.0007):
        y0 = [initial_delta, 0.0]
        return odeint(self.growth_equations, y0, N_array, args=(alpha,))
    
    def growth_ode(self, y, N, *args):
        delta, d_delta = y
        a = np.exp(N)
        
        # Get background quantities from our model
        w = self.model.w_ude(a)
        cs2 = self.model.cs2(a)
        
        # Calculate H'/H and Omega_UDE
        # H'/H = -3/2 * (1 + w_eff)
        # In a UDE dominated universe, w_eff approx w_ude * (rho_ude/rho_tot)
        h_prime_h = -1.5 * (1 + w) 
        
        # Approximation for s_prime (derivative of cs2) 
        # For simplicity in comparison, we'll keep it as 0 or 
        # calculate numerically if your specific model requires it.
        cs2_prime = 0 
        
        # Growth term
        f_term = 1.5 * (1 + 3*cs2) 

        # The ODE system
        dd_delta = (-d_delta * (3*cs2 + 2 + h_prime_h) - 
                    delta * (3*cs2_prime + 3*cs2*(2 + h_prime_h) - f_term))
        
        return [d_delta, dd_delta]

    def solve(self, N_range, delta0=0.0007):
        y0 = [delta0, 0.0] # Initial delta and delta_prime
        sol = odeint(self.growth_ode, y0, N_range)
        return sol
