import numpy as np
from scipy.integrate import odeint
from neutron_trap_creation_models import neutron_trap_creation_numerical
from TDS_sim import festim_sim


def generate_fig_2_annealed_trap_fitting_data():
    """
    orginal data from A.Zaloznik et al, available at https://doi.org/10.1088/0031-8949/t167/1/014031
    subsequently fitted by E.Hodille et al, available at https://doi.org/10.1088/1741-4326/aa5aa5
    """

    defect_type_1_densities = [0.230, 0.230, 0.225, 0.153, 0.107]
    defect_type_2_densities = [0.290, 0.290, 0.280, 0.280, 0.189]
    defect_type_3_densities = [0.05, 0.05, 0.05, 0.05, 0.06]
    annealing_time = 7200

    # ##### standard variables ##### #
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


def generate_fig_4_TDS_fitting_data():
    # 0 dpa values
    festim_sim(
        n1=0,
        n2=0,
        n3=0,
        n4=0,
        n5=0,
        results_foldername="data/damaged_sample_tds_fittings/dpa_0/",
    )

    # 0.001 dpa values
    festim_sim(
        n1=4.5e24,
        n2=1e24,
        n3=5e23,
        n4=1e24,
        n5=2e23,
        results_foldername="data/damaged_sample_tds_fittings/dpa_0.001/",
    )

    # 0.005 dpa values
    festim_sim(
        n1=7e24,
        n2=2.5e24,
        n3=1e24,
        n4=1.9e24,
        n5=1.6e24,
        results_foldername="data/damaged_sample_tds_fittings/dpa_0.005/",
    )

    # 0.023 dpa values
    festim_sim(
        n1=2.4e25,
        n2=1.4e25,
        n3=6e24,
        n4=2.1e25,
        n5=6e24,
        results_foldername="data/damaged_sample_tds_fittings/dpa_0.023/",
    )

    # 0.1 dpa values
    festim_sim(
        n1=5.4e25,
        n2=3.8e25,
        n3=2.8e25,
        n4=3.6e25,
        n5=1.1e25,
        results_foldername="data/damaged_sample_tds_fittings/dpa_0.1/",
    )

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


def generate_fig_5_damaged_trap_fitting_data():
    """
    TDS data from T.Swartz-Selinger, currently unpublished
    """
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

    phi = 8.9e-05
    t_damage = int(3 / phi)
    t = np.linspace(0, t_damage, 1000)
    n_0 = 0
    T_damage = 800  # K

    # trap 1
    trap1_extra_args = (phi, trap_D1_K, trap_D1_n_max, A_0_1, E_A_1, T_damage)
    n_trap1_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap1_extra_args
    )
    # trap 2
    trap2_extra_args = (phi, trap_D2_K, trap_D2_n_max, A_0_1, E_A_1, T_damage)
    n_trap2_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap2_extra_args
    )
    # trap 3
    trap3_extra_args = (phi, trap_D3_K, trap_D3_n_max, A_0_2, E_A_2, T_damage)
    n_trap3_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap3_extra_args
    )
    # trap 4
    trap4_extra_args = (phi, trap_D4_K, trap_D4_n_max, A_0_2, E_A_2, T_damage)
    n_trap4_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap4_extra_args
    )
    # trap 5
    trap5_extra_args = (phi, trap_D5_K, trap_D5_n_max, A_0_3, E_A_3, T_damage)
    n_trap5_damaged = odeint(
        neutron_trap_creation_numerical, n_0, t, args=trap5_extra_args
    )

    # exporting
    np.savetxt("data/damage_trap_D1_fitting.txt", n_trap1_damaged)
    np.savetxt("data/damage_trap_D2_fitting.txt", n_trap2_damaged)
    np.savetxt("data/damage_trap_D3_fitting.txt", n_trap3_damaged)
    np.savetxt("data/damage_trap_D4_fitting.txt", n_trap4_damaged)
    np.savetxt("data/damage_trap_D5_fitting.txt", n_trap5_damaged)


if __name__ == "__main__":
    generate_fig_2_annealed_trap_fitting_data()
    generate_fig_4_TDS_fitting_data()
    generate_fig_5_damaged_trap_fitting_data()
