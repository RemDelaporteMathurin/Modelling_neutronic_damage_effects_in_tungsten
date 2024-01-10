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
t_damage = 86400 # s
T_damage = 800  # K

t = np.linspace(0, t_damage, t_damage)
n_0 = 0

# defect type I
A_0_1 = 6.1838e-03
E_A_1 = 0.24

trap_D1_K = 1.9e28
trap_D1_n_max = 6.8e25
trap_D2_K = 8.0e27
trap_D2_n_max = 6.2e25

# defect type II
A_0_2 = 6.1838e-03
E_A_2 = 0.30

trap_D3_K = 3.5e27
trap_D3_n_max = 5.3e25
trap_D4_K = 7.0e27
trap_D4_n_max = 4.9e25

# defect type III
A_0_3 = 0
E_A_3 = 1

trap_D5_K = 1.2e26
trap_D5_n_max = 2.0e25


trap_D1_fitting = []
trap_D2_fitting = []
trap_D3_fitting = []
trap_D4_fitting = []
trap_D5_fitting = []

for dpa in dpa_x_values:
    phi = dpa / t_damage
    # trap 1
    trap1_extra_args = (phi, trap_D1_K, trap_D1_n_max, A_0_1, E_A_1, T_damage)
    n_trap1_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap1_extra_args
    )
    trap_D1_fitting.append(float(n_trap1_damaged[-1]))
    # trap 2
    trap2_extra_args = (phi, trap_D2_K, trap_D2_n_max, A_0_1, E_A_1, T_damage)
    n_trap2_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap2_extra_args
    )
    trap_D2_fitting.append(float(n_trap2_damaged[-1]))
    # trap 3
    trap3_extra_args = (phi, trap_D3_K, trap_D3_n_max, A_0_2, E_A_2, T_damage)
    n_trap3_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap3_extra_args
    )
    trap_D3_fitting.append(float(n_trap3_damaged[-1]))
    # trap 4
    trap4_extra_args = (phi, trap_D4_K, trap_D4_n_max, A_0_2, E_A_2, T_damage)
    n_trap4_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap4_extra_args
    )
    trap_D4_fitting.append(float(n_trap4_damaged[-1]))
    # trap 5
    trap5_extra_args = (phi, trap_D5_K, trap_D5_n_max, A_0_3, E_A_3, T_damage)
    n_trap5_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap5_extra_args
    )
    trap_D5_fitting.append(float(n_trap5_damaged[-1]))


# exporting
np.savetxt("data/damage_trap_D1_fitting.txt", trap_D1_fitting)
np.savetxt("data/damage_trap_D2_fitting.txt", trap_D2_fitting)
np.savetxt("data/damage_trap_D3_fitting.txt", trap_D3_fitting)
np.savetxt("data/damage_trap_D4_fitting.txt", trap_D4_fitting)
np.savetxt("data/damage_trap_D5_fitting.txt", trap_D5_fitting)
