import numpy as np
from neutron_trap_creation_models import neutron_trap_creation_analytical_steady


dpa_values = np.geomspace(1e-3, 1e03, num=1000)
fpy = 3600 * 24 * 365.25  # full power year
t_damage = fpy
T_damage = 800

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


trap_D1_density_variation = []
trap_D2_density_variation = []
trap_D3_density_variation = []
trap_D4_density_variation = []
trap_D5_density_variation = []

total_trap_density_800 = []
total_trap_density_295 = []

for dpa in dpa_values:
    phi = dpa / t_damage
    trap_D1_density = neutron_trap_creation_analytical_steady(
        T=T_damage, phi=phi, K=trap_D1_K, n_max=trap_D1_n_max, A_0=A_0_1, E_A=E_A_1
    )
    trap_D1_density_variation.append(trap_D1_density)
    trap_D2_density = neutron_trap_creation_analytical_steady(
        T=T_damage, phi=phi, K=trap_D2_K, n_max=trap_D2_n_max, A_0=A_0_1, E_A=E_A_1
    )
    trap_D2_density_variation.append(trap_D2_density)
    trap_D3_density = neutron_trap_creation_analytical_steady(
        T=T_damage, phi=phi, K=trap_D3_K, n_max=trap_D3_n_max, A_0=A_0_2, E_A=E_A_2
    )
    trap_D3_density_variation.append(trap_D3_density)
    trap_D4_density = neutron_trap_creation_analytical_steady(
        T=T_damage, phi=phi, K=trap_D4_K, n_max=trap_D4_n_max, A_0=A_0_2, E_A=E_A_2
    )
    trap_D4_density_variation.append(trap_D4_density)
    trap_D5_density = neutron_trap_creation_analytical_steady(
        T=T_damage, phi=phi, K=trap_D5_K, n_max=trap_D5_n_max, A_0=A_0_3, E_A=E_A_3
    )
    trap_D5_density_variation.append(trap_D5_density)
    
    total_trap_density = total_trap_density = (
        trap_D1_density + trap_D2_density + trap_D3_density + trap_D4_density + trap_D5_density
    )
    total_trap_density_800.append(total_trap_density)
    
    
    # alt for T_damage = 295 K
    trap_D1_density = neutron_trap_creation_analytical_steady(
        T=295, phi=phi, K=trap_D1_K, n_max=trap_D1_n_max, A_0=A_0_1, E_A=E_A_1
    )
    trap_D2_density = neutron_trap_creation_analytical_steady(
        T=295, phi=phi, K=trap_D2_K, n_max=trap_D2_n_max, A_0=A_0_1, E_A=E_A_1
    )
    trap_D3_density = neutron_trap_creation_analytical_steady(
        T=295, phi=phi, K=trap_D3_K, n_max=trap_D3_n_max, A_0=A_0_2, E_A=E_A_2
    )
    trap_D4_density = neutron_trap_creation_analytical_steady(
        T=295, phi=phi, K=trap_D4_K, n_max=trap_D4_n_max, A_0=A_0_2, E_A=E_A_2
    )
    trap_D5_density = neutron_trap_creation_analytical_steady(
        T=295, phi=phi, K=trap_D5_K, n_max=trap_D5_n_max, A_0=A_0_3, E_A=E_A_3
    )
    total_trap_density = (
        trap_D1_density + trap_D2_density + trap_D3_density + trap_D4_density
    )
    total_trap_density_295.append(total_trap_density)
    

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
np.savetxt(
    "data/total_trap_density_295.txt",
    total_trap_density_295,
)
np.savetxt(
    "data/total_trap_density_800.txt",
    total_trap_density_800,
)
