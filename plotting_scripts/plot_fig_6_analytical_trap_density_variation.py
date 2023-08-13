import matplotlib.pyplot as plt
import numpy as np
from neutron_trap_creation_models import neutron_trap_creation_analytical_steady


A_0_optimised = 6.1838e-03
E_A_optimised = 0.2792
dpa_values = np.geomspace(1e-3, 1e03, num=1000)
fpy = 3600 * 24 * 365.25


def obtain_analytical_trap_density_variation(dpa_values, T_damage, t_damage, A_0, E_A):
    """
    Obtains steady-state trap density variation from an analytical model of the trap
    creation, at a given temperarture, T_damage, for a number of dpa values

    Args:
        dpa_values (list, np.array) : dpa values to plot over
        T_damage (float, int) : temperature during damage (K)
        t_damage (flot, int): damage exposure time (s)
        A_0 (float, int) : Trap annealing pre-exponential factor (s-1)
        E_A (float, int): annealing activation energy (eV).

    return:
        trap_D1_density_variation (list): trap D1 concentrations (m-3)
        trap_D2_density_variation (list): trap D2 concentrations (m-3)
        trap_D3_density_variation (list): trap D3 concentrations (m-3)
        trap_D4_density_variation (list): trap D4 concentrations (m-3)
    """
    trap_D1_K = 1.5e28
    trap_D1_n_max = 5.2e25
    trap_D2_K = 4e27
    trap_D2_n_max = 4.5e25
    trap_D3_K = 3e27
    trap_D3_n_max = 4e25
    trap_D4_K = 9e27
    trap_D4_n_max = 4.2e25

    trap_D1_density_variation = []
    trap_D2_density_variation = []
    trap_D3_density_variation = []
    trap_D4_density_variation = []

    for dpa in dpa_values:
        phi = dpa / t_damage
        trap_D1_density = neutron_trap_creation_analytical_steady(
            T=T_damage, phi=phi, K=trap_D1_K, n_max=trap_D1_n_max, A_0=A_0, E_A=E_A
        )
        trap_D1_density_variation.append(trap_D1_density)
        trap_D2_density = neutron_trap_creation_analytical_steady(
            T=T_damage, phi=phi, K=trap_D2_K, n_max=trap_D2_n_max, A_0=A_0, E_A=E_A
        )
        trap_D2_density_variation.append(trap_D2_density)
        trap_D3_density = neutron_trap_creation_analytical_steady(
            T=T_damage, phi=phi, K=trap_D3_K, n_max=trap_D3_n_max, A_0=A_0, E_A=E_A
        )
        trap_D3_density_variation.append(trap_D3_density)
        trap_D4_density = neutron_trap_creation_analytical_steady(
            T=T_damage, phi=phi, K=trap_D4_K, n_max=trap_D4_n_max, A_0=A_0, E_A=E_A
        )
        trap_D4_density_variation.append(trap_D4_density)

    return (
        trap_D1_density_variation,
        trap_D2_density_variation,
        trap_D3_density_variation,
        trap_D4_density_variation,
    )


(
    trap_D1_density_variation,
    trap_D2_density_variation,
    trap_D3_density_variation,
    trap_D4_density_variation,
) = obtain_analytical_trap_density_variation(
    dpa_values=dpa_values,
    T_damage=700,
    t_damage=fpy,
    A_0=A_0_optimised,
    E_A=E_A_optimised,
)

# ##### plotting ##### #
plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

green_ryb = (117 / 255, 184 / 255, 42 / 255)
firebrick = (181 / 255, 24 / 255, 32 / 255)
pewter_blue = (113 / 255, 162 / 255, 182 / 255)
electric_blue = (83 / 255, 244 / 255, 255 / 255)

plt.figure()
plt.plot(
    dpa_values,
    trap_D1_density_variation,
    label=r"Trap D1",
    color=firebrick,
)
plt.plot(
    dpa_values,
    trap_D2_density_variation,
    label=r"Trap D2",
    color=pewter_blue,
)
plt.plot(
    dpa_values,
    trap_D3_density_variation,
    label=r"Trap D3",
    color=electric_blue,
)
plt.plot(
    dpa_values,
    trap_D4_density_variation,
    label=r"Trap D4",
    color=green_ryb,
)
plt.xlabel(r"Damage rate (dpa/fpy)")
plt.ylabel(r"Trap density (m$^{-3}$)")
plt.xlim(dpa_values[0], dpa_values[-1])
plt.ylim(1e23, 1e26)
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.tight_layout()
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.show()
