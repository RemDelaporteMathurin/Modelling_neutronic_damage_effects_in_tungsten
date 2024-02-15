import sympy as sp
from compute_profile_depth import automatic_vertices
import festim as F
import sympy as sp

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

# diffusion properties holtzner D
D_0 = 1.6e-07
E_D = 0.28


def festim_sim(
    n1, n2, n3, n4, n5, initial_number_cells=100, results_foldername="Results/"
):
    """Runs a FESTIM simulation with a custom mesh generator created with the
    automatic_vertices function.

    Args:
        n1 (float): trap density in m-3
        n2 (float): trap density in m-3
        n3 (float): trap density in m-3
        n4 (float): trap density in m-3
        n4 (float): trap density in m-3
        initial_number_of_cells (float): initial number of cells in the mesh
        results_foldername (str): results folder location
    """
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
        D_0=D_0,
        E_D=E_D,
    )
    my_model.materials = F.Materials([tungsten])

    # define traps
    damage_dist = 1 / (1 + sp.exp((F.x - 2.3e-06) / 1e-07))
    traps = [
        F.Trap(
            k_0=tungsten.D_0 / (1.1e-10**2 * 6 * atom_density_W),
            E_k=tungsten.E_D,
            p_0=1e13,
            E_p=1.0,
            density=2e22,
            materials=tungsten,
        )
    ]
    if sum([n1, n2, n3, n4, n5]) > 0:  # if the traps densities are not zero add traps
        traps += [
            F.Trap(
                k_0=tungsten.D_0 / (1.1e-10**2 * 6 * atom_density_W),
                E_k=tungsten.E_D,
                p_0=1e13,
                E_p=1.15,
                density=n1 * damage_dist,
                materials=tungsten,
            ),
            F.Trap(
                k_0=tungsten.D_0 / (1.1e-10**2 * 6 * atom_density_W),
                E_k=tungsten.E_D,
                p_0=1e13,
                E_p=1.35,
                density=n2 * damage_dist,
                materials=tungsten,
            ),
            F.Trap(
                k_0=tungsten.D_0 / (1.1e-10**2 * 6 * atom_density_W),
                E_k=tungsten.E_D,
                p_0=1e13,
                E_p=1.65,
                density=n3 * damage_dist,
                materials=tungsten,
            ),
            F.Trap(
                k_0=tungsten.D_0 / (1.1e-10**2 * 6 * atom_density_W),
                E_k=tungsten.E_D,
                p_0=1e13,
                E_p=1.85,
                density=n4 * damage_dist,
                materials=tungsten,
            ),
            F.Trap(
                k_0=tungsten.D_0 / (1.1e-10**2 * 6 * atom_density_W),
                E_k=tungsten.E_D,
                p_0=1e13,
                E_p=2.05,
                density=n5 * damage_dist,
                materials=tungsten,
            ),
        ]

    my_model.traps = traps
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
        value=sp.Piecewise(
            (exposure_temp, F.t < implantation_time),
            (
                resting_temp,
                (F.t >= implantation_time) & (F.t < implantation_time + resting_time),
            ),
            (
                300 + ramp * (F.t - (implantation_time + resting_time)),
                F.t >= implantation_time + resting_time,
            ),
            (0, True),
        )
    )

    # define boundary conditions
    my_model.boundary_conditions = [
        F.DirichletBC(surfaces=[1, 2], field="solute", value=0)
    ]

    # define sources
    my_model.sources = [
        F.Source(
            volume=1,
            field=0,
            value=sp.Piecewise(
                (flux * distribution, F.t < implantation_time),
                (0, True),
            ),
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
    my_derived_quantities.derived_quantities = [
        average_T,
        H_flux_left,
        H_flux_right,
        solute,
        retention,
        *[
            F.TotalVolume(str(i), volume=1)
            for i, _ in enumerate(my_model.traps.traps, start=1)
        ],
    ]
    print(len(my_model.traps.traps))

    my_model.exports = F.Exports(
        [
            my_derived_quantities,
            F.XDMFExport(
                "retention",
                label="retention",
                folder=folder_results,
                checkpoint=False,
                mode=1,
            ),
        ]
    )

    # define settings
    my_model.dt = F.Stepsize(
        1,
        stepsize_change_ratio=1.1,
        t_stop=implantation_time + resting_time * 0.5,
        dt_min=1e-1,
        stepsize_stop_max=50,
    )

    my_model.settings = F.Settings(
        absolute_tolerance=1e10,
        relative_tolerance=1e-10,
        final_time=implantation_time + resting_time + tds_time,
        transient=True,
        maximum_iterations=30,
        # linear_solver="mumps",
    )

    my_model.initialise()
    my_model.run()
    return my_derived_quantities.data


if __name__ == "__main__":
    pass
