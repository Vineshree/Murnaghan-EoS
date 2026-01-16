# Double Polytropic Cosmic Acceleration from the Murnaghan EoS

This repository contains the numerical implementation and data analysis for the research 
published in *Physics of the Dark Universe* (2024): 
**"Double polytropic cosmic acceleration from the Murnaghan equation of state."**

## Project Overview
This project explores a Unified Dark Energy (UDE) model based on a double polytropic fluid 
inspired by the Murnaghan equation of state from solid-state physics. The code solves 
the background cosmological evolution equations to test the model's viability against 
$\Lambda$ CDM and other dark energy paradigms.

## Key Computational Features
* **Equation of State (EoS) Solver**: Implements the Murnaghan-inspired double polytropic fluid 
    dynamics to compute the evolution of the density parameter ($\Omega$) and Hubble rate ($H$).
* **Dynamical Analysis**: Scripts to calculate the deceleration parameter ($q$) and 
    Statefinder parameters $\{r, s\}$ to characterize cosmic acceleration phases.
* **Numerical Integration**: Uses Scipy's integration suite to solve the coupled 
    cosmological differential equations across redshift ranges $z \in [0, 10]$.
* **Visualization**: High-quality plotting scripts for comparing UDE trajectories against 
    standard cosmological models.

## Repository Structure
* `Murnaghan_Model_Main.ipynb`: The primary Jupyter Notebook containing the model 
    derivations, numerical solvers, and results.
* `src/cosmo_utils.py`: Modular Python functions for cosmological distance and 
    density calculations.
* `plots/`: Generated figures showing the evolution of cosmological parameters.

## Publication
For the theoretical background and detailed derivation, please refer to the full paper:
> Dunsby, P. K. S., Luongo, O., Muccino, M., & Pillay, V. (2024). 
> *Double polytropic cosmic acceleration from the Murnaghan equation of state*. 
> **Physics of the Dark Universe**, 46, 101563. 
> [DOI: 10.1016/j.dark.2024.101563](https://doi.org/10.1016/j.dark.2024.101563)
