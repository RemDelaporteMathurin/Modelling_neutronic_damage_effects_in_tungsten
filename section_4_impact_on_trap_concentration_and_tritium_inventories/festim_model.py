import festim as F
import numpy as np 

# diffusion parameters
# hydrogen holtzner mulitplied by factor sqrt(3) for T
D_0 = 2.06e-7 * (3 ** 0.5)
E_D = 0.28

# damaged trap parameters
# defect type I
A_0_1 = 6.1838e-03
E_A_1 = 0.24

trap_D1_K = 9.0e26
trap_D1_n_max = 6.9e25
trap_D2_K = 4.2e26
trap_D2_n_max = 7.0e25

# defect type II
A_0_2 = 6.1838e-03
E_A_2 = 0.30

trap_D3_K = 2.5e26
trap_D3_n_max = 6.0e25
trap_D4_K = 5.0e26
trap_D4_n_max = 4.7e25

# defect type III
A_0_3 = 0
E_A_3 = 1

trap_D5_K = 1.0e26
trap_D5_n_max = 2.0e25

# general trap parameters
k_0 = 5.22e-17
p_0 = 1e13
trap_1_density = 2e22


def festim_sim(
    dpa=1,
    T=761,
    results_folder_name="Results/",
    total_time=100,
    cells=1000,
    export_retention_field=False,
):
    my_model = F.Simulation(log_level=40)

    # define materials
    tungsten = F.Material(D_0=D_0, E_D=E_D, id=1)
    my_model.materials = F.Materials([tungsten])

    # define traps
    fpy = 3600 * 24 * 365.25
    defined_absolute_tolerance = 1e07
    defined_relative_tolerance = 1e-01
    defined_maximum_iterations = 10

    trap_W_1 = F.Trap(
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=1.0,
        density=2e22,
        materials=tungsten,
    )
    trap_W_damage_1 = F.NeutronInducedTrap(
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=1.15,
        A_0=A_0_1,
        E_A=E_A_1,
        phi=dpa / fpy,
        K=trap_D1_K,
        n_max=trap_D1_n_max,
        materials=tungsten,
        absolute_tolerance=defined_absolute_tolerance,
        relative_tolerance=defined_relative_tolerance,
        maximum_iterations=defined_maximum_iterations,
    )
    trap_W_damage_2 = F.NeutronInducedTrap(
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=1.35,
        A_0=A_0_1,
        E_A=E_A_1,
        phi=dpa / fpy,
        K=trap_D2_K,
        n_max=trap_D2_n_max,
        materials=tungsten,
        absolute_tolerance=defined_absolute_tolerance,
        relative_tolerance=defined_relative_tolerance,
        maximum_iterations=defined_maximum_iterations,
    )
    trap_W_damage_3 = F.NeutronInducedTrap(
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=1.65,
        A_0=A_0_2,
        E_A=E_A_2,
        phi=dpa / fpy,
        K=trap_D3_K,
        n_max=trap_D3_n_max,
        materials=tungsten,
        absolute_tolerance=defined_absolute_tolerance,
        relative_tolerance=defined_relative_tolerance,
        maximum_iterations=defined_maximum_iterations,
    )
    trap_W_damage_4 = F.NeutronInducedTrap(
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=1.85,
        A_0=A_0_2,
        E_A=E_A_2,
        phi=dpa / fpy,
        K=trap_D4_K,
        n_max=trap_D4_n_max,
        materials=tungsten,
        absolute_tolerance=defined_absolute_tolerance,
        relative_tolerance=defined_relative_tolerance,
        maximum_iterations=defined_maximum_iterations,
    )
    trap_W_damage_5 = F.NeutronInducedTrap(
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=2.05,
        A_0=A_0_3,
        E_A=E_A_3,
        phi=dpa / fpy,
        K=trap_D5_K,
        n_max=trap_D5_n_max,
        materials=tungsten,
        absolute_tolerance=defined_absolute_tolerance,
        relative_tolerance=defined_relative_tolerance,
        maximum_iterations=defined_maximum_iterations,
    )

    my_model.traps = F.Traps(
        [
            trap_W_1,
            trap_W_damage_1,
            trap_W_damage_2,
            trap_W_damage_3,
            trap_W_damage_4,
            trap_W_damage_5,
        ]
    )

    vertices = np.linspace(0, 2e-03, num=cells)
    my_model.mesh = F.MeshFromVertices(vertices)

    # define temperature
    my_model.T = F.Temperature(value=T)

    # define boundary conditions
    my_model.boundary_conditions = [
        F.ImplantationDirichlet(
            surfaces=1,
            phi=1e20,
            R_p=3e-09,
            D_0=D_0,
            E_D=E_D,
        ),
    ]

    # define exports
    results_folder = results_folder_name
    my_derived_quantities = F.DerivedQuantities(
        filename=results_folder + "derived_quantities.csv"
    )
    my_derived_quantities.derived_quantities = [
        F.TotalVolume("solute", volume=1),
        F.TotalVolume("retention", volume=1),
        F.TotalVolume("1", volume=1),
        F.TotalVolume("2", volume=1),
        F.TotalVolume("3", volume=1),
        F.TotalVolume("4", volume=1),
        F.TotalVolume("5", volume=1),
        F.TotalVolume("6", volume=1),
    ]

   
    if export_retention_field:
        my_model.exports = F.Exports(
            [
                F.XDMFExport(
                    "retention",
                    label="retention",
                    folder=results_folder,
                    checkpoint=False,
                    mode=1,
                ),
                my_derived_quantities,
            ]
        )
    else:
        my_model.exports = F.Exports([my_derived_quantities])
    

    # define settings
    my_model.dt = F.Stepsize(
        initial_value=0.1,
        stepsize_change_ratio=1.01,
        dt_min=1e-1,
    )
    my_model.settings = F.Settings(
        transient=True,
        final_time=total_time,
        absolute_tolerance=1e10,
        relative_tolerance=1e-10,
        maximum_iterations=30,
    )

    # vrun simulation
    my_model.initialise()
    my_model.run()
