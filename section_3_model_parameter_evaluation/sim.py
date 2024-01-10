import sympy as sp
from compute_profile_depth import automatic_vertices
import festim as F
import numpy as np

fluence = 1.5e25
implantation_time = 72 * 3600
flux = fluence / implantation_time
resting_time = 0.5 * 24 * 3600
exposure_temp = 370
resting_temp = 295
ramp = 3 / 60
tds_time = (1000 - 300) / ramp
size = 8e-04
atom_density_W = 6.3222e28


def festim_sim(
    n1=1,
    n2=1,
    n3=1,
    n4=1,
    n5=1,
    initial_number_cells=100,
    results_foldername="Results/",
):
    """Runs a FESTIM simulation with a custom mesh generator created with the 
    automatic vertices function.

    Args:
        n1 (float): trap_density in m-3
        n2 (float): trap_density in m-3
        n3 (float): trap_density in m-3
        n4 (float): trap_density in m-3
        n4 (float): trap_density in m-3
        initial_number_of_cells (float): initial number of cells in the mesh
        results_foldername (str) results folder location
    """
    r = 0
    center = 0.7e-9
    width = 0.5e-9
    distribution = (
        1 / (width * (2 * 3.14) ** 0.5) * sp.exp(-0.5 * ((F.x - center) / width) ** 2)
    )

    my_model = F.Simulation(log_level=40)

    # define materials
    tungsten = F.Material(
        id=1,
        borders=[0, size],
        D_0=2.4e-07,
        E_D=0.39,
    )
    my_model.materials = F.Materials([tungsten])

    # define traps
    damage_dist = 1 / (1 + sp.exp((F.x - 2.3e-06) / 1e-07))
    my_model.traps = F.Traps(
        [
            F.Trap(
                k_0=4.1e-7 / (1.1e-10**2 * 6 * atom_density_W),
                E_k=0.39,
                p_0=1e13,
                E_p=1.0,
                density=2e22,
                materials=tungsten,
            ),
            F.Trap(
                k_0=4.1e-7 / (1.1e-10**2 * 6 * atom_density_W),
                E_k=0.39,
                p_0=1e13,
                E_p=1.15,
                density=n1 * damage_dist,
                materials=tungsten,
            ),
            F.Trap(
                k_0=4.1e-7 / (1.1e-10**2 * 6 * atom_density_W),
                E_k=0.39,
                p_0=1e13,
                E_p=1.35,
                density=n2 * damage_dist,
                materials=tungsten,
            ),
            F.Trap(
                k_0=4.1e-7 / (1.1e-10**2 * 6 * atom_density_W),
                E_k=0.39,
                p_0=1e13,
                E_p=1.65,
                density=n3 * damage_dist,
                materials=tungsten,
            ),
            F.Trap(
                k_0=2.4e-7 / (1.1e-10**2 * 6 * atom_density_W),
                E_k=0.39,
                p_0=1e13,
                E_p=1.85,
                density=n4 * damage_dist,
                materials=tungsten,
            ),
            F.Trap(
                k_0=2.4e-7 / (1.1e-10**2 * 6 * atom_density_W),
                E_k=0.39,
                p_0=1e13,
                E_p=2.05,
                density=n5 * damage_dist,
                materials=tungsten,
            ),
        ]
    )

    # define mesh
    vertices = automatic_vertices(
        r_p=center,
        size=size,
        mat=tungsten,
        traps=my_model.traps.traps,
        nb_cells=initial_number_cells,
        T=exposure_temp,
        implantation_time=implantation_time,
        flux=flux,
    )
    my_model.mesh = F.MeshFromVertices(vertices)

    # define temperature
    my_model.T = F.Temperature(
        value=exposure_temp * (F.t < implantation_time)
        + resting_temp
        * (F.t >= implantation_time)
        * (F.t < implantation_time + resting_time)
        + (
            (300 + ramp * (F.t - (implantation_time + resting_time)))
            * (F.t >= implantation_time + resting_time)
        )
    )

    # define boundary conditions
    my_model.boundary_conditions = [
        F.DirichletBC(surfaces=[1, 2], field="solute", value=0)
    ]

    # define sources
    my_model.sources = [
        F.Source(
            volume=1, field=0, value=flux * distribution * (F.t < implantation_time)
        )
    ]

    # define exports
    folder_results = results_foldername
    my_derived_quantities = F.DerivedQuantities(
        filename=folder_results + "last.csv",
        nb_iterations_between_exports=1,
    )

    average_T = F.AverageVolume("T", volume=1)
    H_flux_left = F.HydrogenFlux(surface=1)
    H_flux_right = F.HydrogenFlux(surface=2)
    solute = F.TotalVolume("solute", volume=1)
    retention = F.TotalVolume("retention", volume=1)
    trap_1 = F.TotalVolume("1", volume=1)
    trap_2 = F.TotalVolume("2", volume=1)
    trap_3 = F.TotalVolume("3", volume=1)
    trap_4 = F.TotalVolume("4", volume=1)
    trap_5 = F.TotalVolume("5", volume=1)
    trap_6 = F.TotalVolume("6", volume=1)
    my_derived_quantities.derived_quantities = [
        average_T,
        H_flux_left,
        H_flux_right,
        solute,
        retention,
        trap_1,
        trap_2,
        trap_3,
        trap_4,
        trap_5,
        trap_6,
    ]
    
    my_model.exports = F.Exports([my_derived_quantities])

    # define settings
    my_model.dt = F.Stepsize(
        1,
        stepsize_change_ratio=1.02,
        t_stop=implantation_time + resting_time * 0.5,
        dt_min=1e-1,
        stepsize_stop_max=50,
    )

    my_model.settings = F.Settings(
        absolute_tolerance=1e10,
        relative_tolerance=1e-10,
        final_time=implantation_time + resting_time + tds_time,
        transient=True,
        maximum_iterations=10,
        # linear_solver="mumps",
    )

    my_model.initialise()
    my_model.run()
    return my_derived_quantities.data


if __name__ == "__main__":
    # 0 dpa values
    # festim_sim(
    #     n1=0,
    #     n2=0,
    #     n3=0,
    #     n4=0,
    #     n5=0,
    #     results_foldername="data/damaged_sample_tds_fittings/dpa_0/",
    # )
    
    # 0.001 dpa values
    # festim_sim(
    #     n1=4.5e24,
    #     n2=1e24,
    #     n3=5e23,
    #     n4=1e24,
    #     n5=2e23,
    #     results_foldername="data/damaged_sample_tds_fittings/dpa_0.001/",
    # )
    
    # 0.005 dpa values
    # festim_sim(
    #     n1=7e24,
    #     n2=2.5e24,
    #     n3=1e24,
    #     n4=1.9e24,
    #     n5=1.6e24,
    #     results_foldername="data/damaged_sample_tds_fittings/dpa_0.005/",
    # )
    
    # 0.023 dpa values
    # festim_sim(
    #     n1=2.4e25,
    #     n2=1.4e25,
    #     n3=6e24,
    #     n4=2.1e25,
    #     n5=6e24,
    #     results_foldername="data/damaged_sample_tds_fittings/dpa_0.023/",
    # )
    
    # 0.1 dpa values
    # festim_sim(
    #     n1=5.4e25,
    #     n2=3.8e25,
    #     n3=2.8e25,
    #     n4=3.6e25,
    #     n5=1.1e25,
    #     results_foldername="data/damaged_sample_tds_fittings/dpa_0.1/",
    # )
    
    # 0.23 dpa values
    festim_sim(
        n1=5.8e25,
        n2=4.4e25,
        n3=3.5e25,
        n4=4.0e25,
        n5=1.4e25,
        results_foldername="data/damaged_sample_tds_fittings/dpa_0.23/",
    )
    
    # 0.5 dpa values
    festim_sim(
        n1=6.0e25,
        n2=4.8e25,
        n3=4.3e25,
        n4=4.3e25,
        n5=1.75e25,
        results_foldername="data/damaged_sample_tds_fittings/dpa_0.5/",
    )
    
    # 2.5 dpa values
    festim_sim(
        n1=6.8e25,
        n2=6.1e25,
        n3=5e25,
        n4=5e25,
        n5=2e25,
        results_foldername="data/damaged_sample_tds_fittings/dpa_2.5/",
    )
