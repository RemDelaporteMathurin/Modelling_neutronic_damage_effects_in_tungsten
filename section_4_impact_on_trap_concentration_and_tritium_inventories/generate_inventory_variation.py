import numpy as np
from neutron_trap_creation_models import analytical_model

T_range = np.linspace(400, 1300, num=50)
dpa_range = np.geomspace(1e-5, 1e03, num=10)
T_range_contour = np.linspace(400, 1300, num=100)
dpa_range_contour = np.geomspace(1e-3, 1e03, num=100)
fpy = 3600 * 24 * 365.25

inventories = []
inventories_no_damage = []
inventories_standard_temp = []
inventories_contour = []
inventories_no_damage_contour = []
inventories_normalised_contour = []

for dpa in dpa_range:
    phi = dpa / fpy
    (
        H_retention_standard_temp,
        trap_densities_standard_temp,
        trap_filling_ratios_standard_temp,
    ) = analytical_model(phi=phi, T=700)
    inventory_per_dpa = []
    inventories_no_damage = []
    inventories_standard_temp.append(H_retention_standard_temp)
    for T in T_range:
        (H_retention, trap_densities, trap_filling_ratios) = analytical_model(
            phi=phi, T=T
        )
        (
            H_retention_no_damage,
            trap_densities_no_damage,
            trap_filling_ratios_no_damage,
        ) = analytical_model(phi=0, T=T)
        inventory_per_dpa.append(H_retention)
        inventories_no_damage.append(H_retention_no_damage)
    inventories.append(inventory_per_dpa)

for dpa in dpa_range_contour:
    phi = dpa / fpy
    inventory_contour_per_dpa = []
    inventories_no_damage_contour = []
    for T in T_range_contour:
        (H_retention, trap_densities, trap_filling_ratios) = analytical_model(
            phi=phi, T=T
        )
        (
            H_retention_no_damage,
            trap_densities_no_damage,
            trap_filling_ratios_no_damage,
        ) = analytical_model(phi=0, T=T)
        inventory_contour_per_dpa.append(H_retention)
        inventories_no_damage_contour.append(H_retention_no_damage)
    inventories_contour.append(inventory_contour_per_dpa)

for case in inventories_contour:
    inventory = case / np.array(inventories_no_damage_contour)
    inventories_normalised_contour.append(inventory)


# exporting
np.savetxt(
    "data/inventories.txt",
    inventories,
)
np.savetxt(
    "data/inventories_no_damage.txt",
    inventories_no_damage,
)
np.savetxt(
    "data/inventories_nomalised_contour.txt",
    inventories_normalised_contour,
)
