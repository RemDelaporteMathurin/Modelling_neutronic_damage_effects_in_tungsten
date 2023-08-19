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

# exporting
np.savetxt(
    "data/damage_trap_D1_density_variation.txt",
    trap_D1_density_variation,
)
np.savetxt(
    "data/damage_trap_D2_density_variation.txt",
    trap_D2_density_variation,
)
np.savetxt(
    "data/damage_trap_D3_density_variation.txt",
    trap_D3_density_variation,
)
np.savetxt(
    "data/damage_trap_D4_density_variation.txt",
    trap_D4_density_variation,
)
