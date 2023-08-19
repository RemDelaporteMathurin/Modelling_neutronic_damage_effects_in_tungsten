import numpy as np
from scipy.integrate import odeint
from neutron_trap_creation_models import neutron_trap_creation_numerical

"""
orginal data from A.Zaloznik et al, available at https://doi.org/10.1088/0031-8949/t167/1/014031
subsequently fitted by E.Hodille et al, available at https://doi.org/10.1088/1741-4326/aa5aa5
"""

test_temperatures = [298, 600, 800, 1000, 1200]
trap_3_densities = [0.09, 0.08, 0.06, 0.00, 0.00]
trap_4_densities = [0.28, 0.23, 0.19, 0.15, 0.05]
annealing_time = 3600

# ##### standard variables ##### #
atom_density_W = 6.3e28
t_annealing = np.linspace(0, annealing_time, num=1000)
T_values = np.linspace(1, 1400, num=1000)

A_0_optimised = 6.1838e-03
E_A_optimised = 0.2792


def obtain_annealed_trap_densities(T_values, A_0, E_A, t_annealing):
    """_summary_

    Args:
        T_values (list, np.array): _description_
        A_0 (float, int): _description_
        E_A (float, int): _description_
        t_annealing (list, np.array): _description_

    Returns:
        list: _description_
    """
    # dummy values for damage parameters
    phi = 0
    K = 1
    n_max = 1

    # convert trap densities to m-3
    n_0_trap_3 = trap_3_densities[0] * atom_density_W * 1e-02
    n_0_trap_4 = trap_4_densities[0] * atom_density_W * 1e-02

    annealed_trap_3_densities, annealed_trap_4_densities = [], []

    for T in T_values:
        extra_args = (phi, K, n_max, A_0, E_A, T)
        # trap 3
        n_traps_annleaing_trap_3 = odeint(
            neutron_trap_creation_numerical, n_0_trap_3, t_annealing, args=extra_args
        )
        annealed_trap_3_densities.append(n_traps_annleaing_trap_3[-1])

        # trap 4
        n_traps_annleaing_trap_4 = odeint(
            neutron_trap_creation_numerical, n_0_trap_4, t_annealing, args=extra_args
        )
        annealed_trap_4_densities.append(n_traps_annleaing_trap_4[-1])

    return annealed_trap_3_densities, annealed_trap_4_densities


annealed_trap_3_densities, annealed_trap_4_densities = obtain_annealed_trap_densities(
    T_values=T_values, A_0=A_0_optimised, E_A=E_A_optimised, t_annealing=t_annealing
)

# exporting
np.savetxt("data/annealed_trap_3_densities.txt", annealed_trap_3_densities)
np.savetxt("data/annealed_trap_4_densities.txt", annealed_trap_4_densities)
