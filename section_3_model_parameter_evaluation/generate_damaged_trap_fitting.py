import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from neutron_trap_creation_models import neutron_trap_creation_numerical

"""
TDS data from T.Swartz-Selinger, currently unpublished
"""
dpa_values = [0, 0.001, 0.005, 0.023, 0.1, 0.23, 0.5, 2.5]

# trap density variations, fitted from TDS data
trap_D1_densities = [0, 4.5e24, 7.0e24, 2.4e25, 5.4e25, 5.8e25, 6.0e25, 6.8e25]
trap_D2_densities = [0, 1.0e24, 2.5e24, 1.4e24, 3.8e25, 4.4e25, 4.8e25, 6.1e25]
trap_D3_densities = [0, 5.0e23, 1.0e24, 6.0e24, 2.8e25, 3.5e25, 4.3e25, 5.0e25]
trap_D4_densities = [0, 1.0e24, 1.9e24, 2.1e25, 3.6e25, 4.0e25, 4.3e25, 5.0e25]
trap_D5_densities = [0, 2.0e23, 1.6e24, 6.0e24, 1.1e25, 1.4e25, 1.8e25, 2.0e25]

dpa_x_values = np.linspace(0, 3, num=100)
t_damage = 86400
A_0_1 = 6.1838e-03
E_A_1 = 0.26
E_A_2 = 0.34

# dummy values for no annealing
A_0_2 = 0
E_A_3 = 1

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

# exporting
np.savetxt("data/damage_trap_D1_fitting.txt", trap_D1_fitting)
np.savetxt("data/damage_trap_D2_fitting.txt", trap_D2_fitting)
np.savetxt("data/damage_trap_D3_fitting.txt", trap_D3_fitting)
np.savetxt("data/damage_trap_D4_fitting.txt", trap_D4_fitting)
