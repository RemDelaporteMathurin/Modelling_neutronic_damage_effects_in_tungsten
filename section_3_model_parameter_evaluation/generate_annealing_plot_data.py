import numpy as np
from scipy.integrate import odeint
from neutron_trap_creation_models import neutron_trap_creation_numerical

"""
orginal data from A.Zaloznik et al, available at https://doi.org/10.1088/0031-8949/t167/1/014031
subsequently fitted by E.Hodille et al, available at https://doi.org/10.1088/1741-4326/aa5aa5
"""

test_temperatures = [370, 400, 500, 600, 800]
defect_type_1_densities = [0.230, 0.230, 0.225, 0.153, 0.107]
defect_type_2_densities = [0.290, 0.290, 0.280, 0.280, 0.189]
defect_type_3_densities = [0.05, 0.05, 0.05, 0.05, 0.06]
annealing_time = 7200

# ##### standard variables ##### #
atom_density_W = 6.3e28
t_annealing = np.linspace(0, annealing_time, num=1000)
T_values = np.linspace(1, 800, num=1000)

A_0_1 = 6.1834e-03
E_A_1 = 0.24
E_A_2 = 0.30

# dummy values for no annealing
A_0_2 = 0
E_A_3 = 1

# dummy values for damage parameters
phi = 0
K = 1
n_max = 1

n_0_defect_1, n_0_defect_2, n_0_defect_3 = (
    defect_type_1_densities[0],
    defect_type_2_densities[0],
    defect_type_3_densities[0],
)

(
    annealed_defect_1_densities,
    annealed_defect_2_densities,
    annealed_defect_3_densities,
) = ([], [], [])

for T in T_values:
    # defect type 1
    defect_1_extra_args = (phi, K, n_max, A_0_1, E_A_1, T)
    n_traps_annleaing_defect_1 = odeint(
        neutron_trap_creation_numerical,
        n_0_defect_1,
        t_annealing,
        args=defect_1_extra_args,
    )
    annealed_defect_1_densities.append(n_traps_annleaing_defect_1[-1])

    # defect type 2
    defect_2_extra_args = (phi, K, n_max, A_0_1, E_A_2, T)
    n_traps_annleaing_defect_2 = odeint(
        neutron_trap_creation_numerical,
        n_0_defect_2,
        t_annealing,
        args=defect_2_extra_args,
    )
    annealed_defect_2_densities.append(n_traps_annleaing_defect_2[-1])

    # defect type 3
    defect_3_extra_args = (phi, K, n_max, A_0_2, E_A_3, T)
    n_traps_annleaing_defect_3 = odeint(
        neutron_trap_creation_numerical,
        n_0_defect_3,
        t_annealing,
        args=defect_3_extra_args,
    )
    annealed_defect_3_densities.append(n_traps_annleaing_defect_3[-1])


# exporting
np.savetxt("data/annealed_defect_1_densities.txt", annealed_defect_1_densities)
np.savetxt("data/annealed_defect_2_densities.txt", annealed_defect_2_densities)
np.savetxt("data/annealed_defect_3_densities.txt", annealed_defect_3_densities)
