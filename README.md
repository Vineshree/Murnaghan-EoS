# Double Polytropic Cosmic Acceleration from the Murnaghan EoS

T# Murnaghan-EoS: Unified Dark Fluid Cosmology

This repository contains the numerical implementation and perturbation analysis of the **Murnaghan Equation of State (EoS)** as a Unified Dark Fluid (UDF) model. The project compares the Murnaghan model against standard $\Lambda$CDM, Generalized Chaplygin Gas (GCG), and Logotropic models, specifically focusing on the growth of matter density perturbations.

## 🌌 Scientific Context
The Murnaghan EoS offers a unique approach to cosmic acceleration by representing the entire dark sector (Dark Matter and Dark Energy) as a single perfect fluid. Unlike other UDE models, the Murnaghan formulation:

* **Dust-like Early Phase:** Recovers a pressureless dust phase ($w \approx 0$, $c_s^2 \approx 0$) at early times ($V \ll V_0$), which is a critical requirement for healthy structure formation.
* **Late-time Acceleration:** Transitions into a Chaplygin-like acceleration phase as the volume increases.
* **Stability:** Maintains a positive-defined sound speed throughout cosmic evolution, avoiding the typical instabilities of unified models.

## 🛠 Features
* **Physics Engine:** Implementation of background evolution ($H(z)$, $w(z)$) and sound speed ($c_s^2$) for multiple UDE models.
* **Perturbation Solver:** Linear growth solver for the density contrast $\delta_m$ and the growth rate $f = d \ln \delta / d \ln a$.
* **Kinematic Analysis:** Functions to compute the current deceleration ($q_0$) and jerk ($j_0$) parameters.
* **Unit Testing:** A `pytest` suite ensuring physical normalization ($\Omega_{tot}=1$) and stability constraints.

## 📁 Repository Structure
```text
Murnaghan-EoS/
├── src/
│   ├── physics/
│   │   ├── ude_models.py    # Core Murnaghan, GCG, and Logotropic classes
│   │   └── cosmo_solver.py  # ODE solver for matter perturbations
├── tests/
│   └── test_physics.py      # Physical sanity checks (Normalization, w-limit)
├── Notebook_UDE_Comparison.ipynb # Main analysis and plotting
└── README.md

## Publication
For the theoretical background and detailed derivation, please refer to the full paper:
> Dunsby, P. K. S., Luongo, O., Muccino, M., & Pillay, V. (2024). 
> *Double polytropic cosmic acceleration from the Murnaghan equation of state*. 
> **Physics of the Dark Universe**, 46, 101563. 
> [DOI: 10.1016/j.dark.2024.101563](https://doi.org/10.1016/j.dark.2024.101563)
