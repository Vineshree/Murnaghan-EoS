import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

def plot_ude_comparison(results_dict, alpha_vals, a_array):
    """
    results_dict: Dictionary containing arrays for rho, w, cs2, and delta.
    alpha_vals: The range of alpha parameters used.
    a_array: The scale factor array (x-axis).
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharex=True)
    cmap = plt.get_cmap('tab20b')
    norm = mpl.colors.Normalize(vmin=min(alpha_vals), vmax=max(alpha_vals))
    
    # --- Plot 1: Density Evolution ---
    ax = axes[0, 0]
    for i, alpha in enumerate(alpha_vals):
        ax.plot(a_array, results_dict['rho'][i], color=cmap(norm(alpha)))
    ax.set_yscale('log')
    ax.set_ylabel(r'$\rho/\rho_{cr}$')
    ax.set_title("Energy Density")

    # --- Plot 2: Equation of State ---
    ax = axes[0, 1]
    for i, alpha in enumerate(alpha_vals):
        ax.plot(a_array, results_dict['w'][i], color=cmap(norm(alpha)))
    ax.axhline(-1, color='k', linestyle=':')
    ax.set_ylabel(r'$w$')
    ax.set_title("Equation of State")

    # --- Plot 3: Sound Speed ---
    ax = axes[1, 0]
    for i, alpha in enumerate(alpha_vals):
        ax.plot(a_array, results_dict['cs2'][i], color=cmap(norm(alpha)))
    ax.set_ylabel(r'$c_s^2$')
    ax.set_title("Sound Speed Squared")

    # --- Plot 4: Growth Function (f) ---
    ax = axes[1, 1]
    for i, alpha in enumerate(alpha_vals):
        ax.plot(a_array, results_dict['f'][i], color=cmap(norm(alpha)))
    ax.set_ylabel(r'Growth Rate $f$')
    ax.set_title("Growth Rate")

    for ax in axes.flat:
        ax.set_xscale('log')
        ax.set_xlabel('$a$')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig
