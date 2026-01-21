import pytest
import numpy as np
from src.physics.ude_models import MurnaghanModel

def test_murnaghan_today_density():
    """Test if total density at a=1 is roughly 1.0 (normalized)."""
    model = MurnaghanModel(Omega_m=0.3)
    # Today (a=1), total density should be close to 1.0
    assert np.isclose(model.rho_M(a=1.0), 1.0, rtol=1e-2)

def test_murnaghan_sound_speed_positivity():
    """Physical constraint: Sound speed squared must be positive (stability)."""
    model = MurnaghanModel(alpha=0.018)
    a_vals = np.logspace(-3, 0, 10)
    for a in a_vals:
        assert model.s_M(a) >= 0

def test_murnaghan_early_dust_limit():
    """Test if EoS w -> 0 at early times (a << 1)."""
    model = MurnaghanModel(alpha=0.018)
    # At very high redshift, it should behave like dust (w=0)
    w_early = model.w_total(a=1e-5)
    assert np.isclose(w_early, 0.0, atol=1e-3)
