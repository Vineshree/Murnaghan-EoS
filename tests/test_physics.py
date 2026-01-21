import sys
import os

# Adds the parent directory to the python path
sys.path.append(os.path.abspath('../'))

import pytest
import numpy as np
from src.physics.ude_models import MurnaghanModel

def test_murnaghan_today_density():
    """Test if total density at a=1 is roughly 1.0."""
    model = MurnaghanModel(Omega_m=0.3)
    # Check your ude_models.py file—if the function is epsilon_M, use that here:
    assert np.isclose(model.epsilon_M(a=1.0), 1.0, rtol=1e-2)

def test_murnaghan_sound_speed_positivity():
    """Physical constraint: Sound speed squared must be positive (stability)."""
    model = MurnaghanModel(alpha=0.018)
    a_vals = np.logspace(-3, 0, 10)
    for a in a_vals:
        # Testing the sound speed method (cs2)
        assert model.cs2(a) >= 0

def test_murnaghan_early_dust_limit():
    """Test if EoS w -> 0 at early times (a << 1)."""
    model = MurnaghanModel(alpha=0.018)
    # Testing the Equation of State method (w_ude)
    # At high redshift, UDE should mimic dust (w=0)
    w_early = model.w_ude(a=1e-5)
    assert np.isclose(w_early, 0.0, atol=1e-3)