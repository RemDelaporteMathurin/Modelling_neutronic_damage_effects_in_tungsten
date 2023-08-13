import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from neutron_trap_creation_models import neutron_trap_creation_numerical

"""
TDS data from T.Swartz-Selinger, currently unpublished
"""
dpa_values = [0, 0.001, 0.005, 0.023, 0.1, 0.23, 0.5, 2.5]

# trap density variations, fitted from TDS data
trap_D1_densities = [0, 3.5e24, 5.3e24, 2.0e25, 4.2e25, 4.5e25, 4.7e25, 5.3e25]
trap_D2_densities = [0, 5.0e23, 1.9e24, 9.5e24, 2.6e25, 3.4e25, 3.6e25, 4.5e25]
trap_D3_densities = [0, 5.0e23, 1.0e24, 6.0e24, 2.0e25, 2.7e25, 3.2e25, 3.9e25]
trap_D4_densities = [0, 1.0e24, 2.0e24, 1.7e25, 3.2e25, 3.4e25, 3.8e25, 4.2e25]

dpa_x_values = np.linspace(0, 3, num=100)
t_damage = 86400
A_0_optimised = 6.1838e-03
E_A_optimised = 0.2792
T = 800  # K


def obtain_damaged_trap_density_varience_fitting(dpa_x_values, t_damage, A_0, E_A, T):
    """
    Obtains trap density variation from a numerical model of the trap creation, taking
    the final value after a period of damage, for a number of dpa values

    Args:
        dpa_x_values (list, np.array) : dpa values to plot over
        t_damage (float, int) : damaging time (s)
        A_0 (float, int) : Trap annealing pre-exponential factor (s-1)
        E_A (float, int): annealing activation energy (eV).
        T (float, int): temperature during damage (K)

    return:
        damaged_trap_D1_densities (list): trap D1 concentrations (m-3)
        damaged_trap_D2_densities (list): trap D2 concentrations (m-3)
        damaged_trap_D3_densities (list): trap D3 concentrations (m-3)
        damaged_trap_D4_densities (list): trap D4 concentrations (m-3)
    """
    t = np.linspace(0, t_damage, t_damage)
    n_0 = 0

    trap_D1_K = 1.5e28
    trap_D1_n_max = 5.2e25
    trap_D2_K = 4e27
    trap_D2_n_max = 4.5e25
    trap_D3_K = 3e27
    trap_D3_n_max = 4e25
    trap_D4_K = 9e27
    trap_D4_n_max = 4.2e25

    trap_D1_fitting = []
    trap_D2_fitting = []
    trap_D3_fitting = []
    trap_D4_fitting = []

    for dpa in dpa_x_values:
        phi = dpa / t_damage
        # trap 1
        trap1_extra_args = (phi, trap_D1_K, trap_D1_n_max, A_0, E_A, T)
        n_trap1_damaged = odeint(
            neutron_trap_creation_numerical, n_0, t, args=trap1_extra_args
        )
        trap_D1_fitting.append(float(n_trap1_damaged[-1]))
        # trap 2
        trap2_extra_args = (phi, trap_D2_K, trap_D2_n_max, A_0, E_A, T)
        n_trap2_damaged = odeint(
            neutron_trap_creation_numerical, n_0, t, args=trap2_extra_args
        )
        trap_D2_fitting.append(float(n_trap2_damaged[-1]))
        # trap 3
        trap3_extra_args = (phi, trap_D3_K, trap_D3_n_max, A_0, E_A, T)
        n_trap3_damaged = odeint(
            neutron_trap_creation_numerical, n_0, t, args=trap3_extra_args
        )
        trap_D3_fitting.append(float(n_trap3_damaged[-1]))
        # trap 4
        trap4_extra_args = (phi, trap_D4_K, trap_D4_n_max, A_0, E_A, T)
        n_trap4_damaged = odeint(
            neutron_trap_creation_numerical, n_0, t, args=trap4_extra_args
        )
        trap_D4_fitting.append(float(n_trap4_damaged[-1]))

    return (trap_D1_fitting, trap_D2_fitting, trap_D3_fitting, trap_D4_fitting)


(
    trap_D1_fitting,
    trap_D2_fitting,
    trap_D3_fitting,
    trap_D4_fitting,
) = obtain_damaged_trap_density_varience_fitting(
    dpa_x_values=dpa_x_values,
    t_damage=t_damage,
    A_0=A_0_optimised,
    E_A=E_A_optimised,
    T=T,
)

# ##### plotting ##### #

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

green_ryb = (117 / 255, 184 / 255, 42 / 255)
firebrick = (181 / 255, 24 / 255, 32 / 255)
pewter_blue = (113 / 255, 162 / 255, 182 / 255)
electric_blue = (83 / 255, 244 / 255, 255 / 255)

fig = plt.figure(figsize=(4.5, 9))

# create one big plot to have a common y label
ax = fig.add_subplot(111)
ax.spines["top"].set_color("none")
ax.spines["bottom"].set_color("none")
ax.spines["left"].set_color("none")
ax.spines["right"].set_color("none")
ax.tick_params(labelcolor="w", top=False, bottom=False, left=False, right=False)
ax.set_ylabel(r"Trap density (m$^{-3}$)")

# highlight trap 1
ax1 = fig.add_subplot(411)
plt.sca(ax1)
plt.scatter(dpa_values, trap_D1_densities, marker="x", color=firebrick)
plt.plot(dpa_x_values, trap_D1_fitting, color=firebrick, label=r"Trap D1")
plt.plot(dpa_x_values, trap_D2_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D3_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D4_fitting, color="grey", alpha=0.2)
plt.legend(loc="lower right")

# highlight trap D2
ax2 = fig.add_subplot(412)
plt.sca(ax2)
plt.scatter(dpa_values, trap_D2_densities, color=pewter_blue, marker="x")
plt.plot(dpa_x_values, trap_D2_fitting, color=pewter_blue, label=r"Trap D2")
plt.plot(dpa_x_values, trap_D1_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D3_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D4_fitting, color="grey", alpha=0.2)
plt.legend(loc="lower right")

# highlight trap D3
ax3 = fig.add_subplot(413)
plt.sca(ax3)
plt.scatter(dpa_values, trap_D3_densities, color=electric_blue, marker="x")
plt.plot(dpa_x_values, trap_D3_fitting, color=electric_blue, label=r"Trap D3")
plt.plot(dpa_x_values, trap_D1_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D2_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D4_fitting, color="grey", alpha=0.2)
plt.legend(loc="lower right")

# highlight trap D4
ax4 = fig.add_subplot(414)
plt.sca(ax4)
plt.scatter(dpa_values, trap_D4_densities, color=green_ryb, marker="x")
plt.plot(dpa_x_values, trap_D4_fitting, color=green_ryb, label=r"Trap D4")
plt.plot(dpa_x_values, trap_D1_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D2_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D3_fitting, color="grey", alpha=0.2)
plt.legend(loc="lower right")
plt.xlabel(r"Damage (dpa)")

for ax in [ax1, ax2, ax3, ax4]:
    plt.sca(ax)
    ax.get_shared_x_axes().join(ax, ax4)
    plt.ylim(0, 6e25)
    plt.xlim(left=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# remove the xticks for top plots
ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.set_xticklabels([])

plt.tight_layout()

plt.show()
