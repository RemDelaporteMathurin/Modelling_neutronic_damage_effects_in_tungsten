from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import odeint

from neutron_trap_creation_models import (
    neutron_trap_creation_numerical,
)

"""
orginal data from A.Zaloznik et al, available at https://doi.org/10.1088/0031-8949/t167/1/014031
subsequently fitted by E.Hodille et al, available at https://doi.org/10.1088/1741-4326/aa5aa5
"""

test_temperatures = [298, 600, 800, 1000, 1200]
trap_3_densities =  [0.09, 0.08, 0.06, 0.00, 0.00]
trap_4_densities = [0.28, 0.23, 0.19, 0.15, 0.05]
annealing_time = 3600

# ##### standard variables ##### #
atom_density_W = 6.3e28
t_values = np.linspace(0, annealing_time, num=1000)
T_values = np.linspace(1, 1400, num=1000)
A_0_optimised = 6.1838e-03
E_A_optimised = 0.2792

# dummy values for damage parameters
phi = 0
K = 1
n_max = 1

# ##### post-processing ##### #

def evaluate_annealed_trap_density_fitting():
    """ Runs a numerical model to evaluate how the trap
    densities vary over the annealing time and takes the final
    value
    """

    n_0_trap_3 = trap_3_densities[0] * atom_density_W * 1e-02
    n_0_trap_4 = trap_4_densities[0] * atom_density_W * 1e-02

    annealed_trap_3_densities = []
    annealed_trap_4_densities = []

    for T in T_values:
        extra_args = (phi, K, n_max, A_0_optimised, E_A_optimised, T)
        # trap 3
        n_traps_annleaing_trap_3 = odeint(
            neutron_trap_creation_numerical, n_0_trap_3, t_values, args=extra_args
        )
        end_value_trap_3 = float(n_traps_annleaing_trap_3[-1])
        annealed_trap_3_densities.append(end_value_trap_3)

        # trap 4
        n_traps_annleaing_trap_4 = odeint(
            neutron_trap_creation_numerical, n_0_trap_4, t_values, args=extra_args
        )
        end_value_trap_4 = float(n_traps_annleaing_trap_4[-1])
        annealed_trap_4_densities.append(end_value_trap_4)

    return annealed_trap_3_densities, annealed_trap_4_densities


annealed_trap_3_densities, annealed_trap_4_densities = evaluate_annealed_trap_density_fitting()
trap_3_densities = (np.array(trap_3_densities) / 100) * 6.3e28
trap_4_densities = (np.array(trap_4_densities) / 100) * 6.3e28

# ##### Plotting ##### #

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

green_ryb = (117 / 255, 184 / 255, 42 / 255)
firebrick = (181 / 255, 24 / 255, 32 / 255)

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(5, 7.5))

# trap 3
plt.sca(axs[0])
# ref values
plt.scatter(
    test_temperatures,
    trap_3_densities,
    marker="x",
    color=green_ryb,
)
# fitting
plt.plot(
    T_values,
    annealed_trap_3_densities,
    color=green_ryb,
    label=r"Trap 3",
)

# trap 4
plt.sca(axs[1])
# ref values
plt.scatter(
    test_temperatures,
    trap_4_densities,
    marker="x",
    color=firebrick,
)
# fitting
plt.plot(
    T_values,
    annealed_trap_4_densities,
    color=firebrick,
    label=r"Trap 4",
)

for ax in [axs[0], axs[1]]:
    plt.sca(ax)
    plt.ylim(bottom=0)
    plt.xlim(0, 1400)
    plt.ylabel(r"Trap density, n$_{\mathrm{t}}$ (m$^{-3}$)")
    plt.legend(loc="lower left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
plt.xlabel(r"Annealing temperature (K)")
plt.subplots_adjust(wspace=0.112, hspace=0.2)
plt.tight_layout()

plt.show()
