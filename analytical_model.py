import numpy as np
import csv


def trap_concentration(T, phi, K, n_max, A_0, E_A):
    """
    Function to evaluate the trap concentration at a given temperature, T, and
    damage, phi, level

    Args:
        T (float, int): temperature (K)
        phi (float, int): damage rate (dpa s-1),
        K (float int): trap creation factor (m-3 dpa-1),
        n_max (float, int):  maximum trap density (m-3),
        A_0 (float, int): trap_annealing_factor (s-1),
        E_A (float, int): annealing activation energy (eV).

    return:
        n_i (float): trap concentration (m-3)
    """
    if phi == 0:
        n_i = 0
    else:
        A = A_0 * np.exp(-E_A / k_B / T)
        n_i = 1 / ((A / (phi * K) + (1 / (n_max))))

    return n_i


def mobile_H_concentration(T, imp_flux, r_p, D_0, E_D):
    """
    Function to evaluate the mobile concentration of tritium

    Args:
        T (float, int): temperature (K)
        imp_flux (float, int): implantation flux (H m-2 s-1),
        r_p (float int): implantation depth (m),
        D_0 (float, int): diffusion coefficient pre-exponential factor (m2/s)
        E_D (float, int): diffusion coefficient activation energy (eV)

    return:
        c_m (float): mobile concentration of H (m-3)

    """
    D = D_0 * np.exp(-E_D / k_B / T)
    c_m = (imp_flux * r_p) / D

    return c_m


def trapped_H_concentration(T, c_m, k_0, E_k, p_0, E_p, n_i, A_0=0, E_A=1):
    """
    Function to evaluate the trapped concentration

    Args:
        T (float, int): temperature (K)
        c_m (float, int): implantation flux (H m-2 s-1),
        k_0 (float, int): trapping pre-exponential factor (m3 s-1)
        E_k (float, int): trapping activation energy (eV)
        p_0 (float, int): detrapping pre-exponential factor (s-1)
        E_p (float, int): detrapping activation energy (eV)
        A_0 (float, int): trap_annealing_factor (s-1),
        E_A (float, int): annealing activation energy (eV).

    return:
        c_t (float): trapped concentration of H (m-3)
        filling_ratio (float): filling ratio of trap
    """
    A = A_0 * np.exp(-E_A / k_B / T)
    v_t = k_0 * np.exp(-E_k / k_B / T)
    v_dt = p_0 * np.exp(-E_p / k_B / T)

    filling_ratio = 1 / (1 + ((v_dt + A) / (v_t * c_m)))

    c_t = filling_ratio * n_i

    return c_t, filling_ratio


def retention(c_m, c_t, L):
    """
    Function to evaluate the retention

     Args:
        c_m (float, int): mobile concentration of H (m-2)
        c_t (float, int): trapped concentration of H (m-2)
        L (float, int): Length of material (m)

    return:
        retention (float): total retention of H (m-2)

    """
    retention = L * (c_m + c_t)

    return retention


def analytical_model(phi=9 / (3600 * 24 * 365.25), T=761, L=0.002):
    trap_1_density = 8.22e25
    trap_2_density = 2.53e25
    trap_3_density = trap_concentration(
        T=T, phi=phi, K=1.5e28, n_max=5.2e25, A_0=A_0, E_A=E_A
    )
    trap_4_density = trap_concentration(
        T=T, phi=phi, K=4.0e27, n_max=4.5e25, A_0=A_0, E_A=E_A
    )
    trap_5_density = trap_concentration(
        T=T, phi=phi, K=3.0e27, n_max=4.0e25, A_0=A_0, E_A=E_A
    )
    trap_6_density = trap_concentration(
        T=T, phi=phi, K=9.0e27, n_max=4.2e25, A_0=A_0, E_A=E_A
    )

    c_m = mobile_H_concentration(T=T, imp_flux=imp_flux, r_p=r_p, D_0=D_0, E_D=E_D)

    c_t_1, filling_ratio_1 = trapped_H_concentration(
        c_m=c_m, k_0=k_0_1, E_k=E_k, p_0=p_0, E_p=0.87, n_i=trap_1_density, T=T
    )
    c_t_2, filling_ratio_2 = trapped_H_concentration(
        c_m=c_m, k_0=k_0_2, E_k=E_k, p_0=p_0, E_p=1.00, n_i=trap_2_density, T=T
    )
    c_t_3, filling_ratio_3 = trapped_H_concentration(
        c_m=c_m,
        k_0=k_0_1,
        E_k=E_k,
        p_0=p_0,
        E_p=1.15,
        A_0=A_0,
        E_A=E_A,
        n_i=trap_3_density,
        T=T,
    )
    c_t_4, filling_ratio_4 = trapped_H_concentration(
        c_m=c_m,
        k_0=k_0_1,
        E_k=E_k,
        p_0=p_0,
        E_p=1.35,
        A_0=A_0,
        E_A=E_A,
        n_i=trap_4_density,
        T=T,
    )
    c_t_5, filling_ratio_5 = trapped_H_concentration(
        c_m=c_m,
        k_0=k_0_1,
        E_k=E_k,
        p_0=p_0,
        E_p=1.65,
        A_0=A_0,
        E_A=E_A,
        n_i=trap_5_density,
        T=T,
    )
    c_t_6, filling_ratio_6 = trapped_H_concentration(
        c_m=c_m,
        k_0=k_0_1,
        E_k=E_k,
        p_0=p_0,
        E_p=1.85,
        A_0=A_0,
        E_A=E_A,
        n_i=trap_6_density,
        T=T,
    )

    trap_densities = [
        trap_1_density,
        trap_2_density,
        trap_3_density,
        trap_4_density,
        trap_5_density,
        trap_6_density,
    ]

    trap_filling_ratios = [
        filling_ratio_1,
        filling_ratio_2,
        filling_ratio_3,
        filling_ratio_4,
        filling_ratio_5,
        filling_ratio_6,
    ]

    total_trapped_H_concentration = c_t_1 + c_t_2 + c_t_3 + c_t_4 + c_t_5 + c_t_6

    total_retention = retention(c_m=c_m, c_t=total_trapped_H_concentration, L=L)

    return (total_retention, trap_densities, trap_filling_ratios)


def de_trapping_rate(E_p, T):
    v_dt = p_0 * np.exp(-E_p / k_B / T)

    return v_dt


def annealing_rate(T):
    A = A_0 * np.exp(-E_A / k_B / T)

    return A


# ##### standard values ##### #

L = 0.002
imp_flux = 1e20
r_p = 3e-09
D_0 = 4.1e-7
E_D = 0.39
k_B = 8.6173303e-5
A_0 = 6.18e-03
E_A = 0.28
k_0_1 = 5.22e-17
k_0_2 = 8.93e-17
E_k = 0.39
p_0 = 1e13

trap_1_conc = 8.22e25
trap_2_conc = 2.53e25

# ##### test ranges ##### #

T_range = np.linspace(400, 1300, num=1000)
dpa_range = np.geomspace(1e-03, 1e03, num=1000)

# ##### analytical test results gathering ##### #

trap_density_values_by_T = []
trap_density_values_standard_temp = []

for dpa in dpa_range:
    phi = dpa / (3600 * 24 * 365.25)
    (H_retention, trap_densities, trap_filling_ratios) = analytical_model(
        phi=phi, T=700
    )
    trap_density_values_standard_temp.append(trap_densities)

for T in T_range:
    trap_densities_by_T = []
    (H_retention, trap_densities, trap_filling_ratios) = analytical_model(
        phi=9 / (3600 * 24 * 365.25), T=T
    )
    trap_density_values_by_T.append(trap_densities)


# ##### Post processing ##### #

trap_1_densities_standard_temp, trap_2_densities_standard_temp = [], []
trap_3_densities_standard_temp, trap_4_densities_standard_temp = [], []
trap_5_densities_standard_temp, trap_6_densities_standard_temp = [], []
trap_1_densities_by_T, trap_2_densities_by_T = [], []
trap_3_densities_by_T, trap_4_densities_by_T = [], []
trap_5_densities_by_T, trap_6_densities_by_T = [], []

for density in trap_density_values_standard_temp:
    trap_1_densities_standard_temp.append(density[0])
    trap_2_densities_standard_temp.append(density[1])
    trap_3_densities_standard_temp.append(density[2])
    trap_4_densities_standard_temp.append(density[3])
    trap_5_densities_standard_temp.append(density[4])
    trap_6_densities_standard_temp.append(density[5])

for density in trap_density_values_by_T:
    trap_1_densities_by_T.append(density[0])
    trap_2_densities_by_T.append(density[1])
    trap_3_densities_by_T.append(density[2])
    trap_4_densities_by_T.append(density[3])
    trap_5_densities_by_T.append(density[4])
    trap_6_densities_by_T.append(density[5])


results_folder = "Results/analytical_model_testing/"
np.savetxt(
    results_folder + "trap_1_densities_standard_T.csv",
    trap_1_densities_standard_temp,
    delimiter=",",
)
np.savetxt(
    results_folder + "trap_2_densities_standard_T.csv",
    trap_2_densities_standard_temp,
    delimiter=",",
)
np.savetxt(
    results_folder + "trap_3_densities_standard_T.csv",
    trap_3_densities_standard_temp,
    delimiter=",",
)
np.savetxt(
    results_folder + "trap_4_densities_standard_T.csv",
    trap_4_densities_standard_temp,
    delimiter=",",
)
np.savetxt(
    results_folder + "trap_5_densities_standard_T.csv",
    trap_5_densities_standard_temp,
    delimiter=",",
)
np.savetxt(
    results_folder + "trap_6_densities_standard_T.csv",
    trap_6_densities_standard_temp,
    delimiter=",",
)

np.savetxt(
    results_folder + "trap_1_densities_varying_T.csv",
    trap_1_densities_by_T,
    delimiter=",",
)
np.savetxt(
    results_folder + "trap_2_densities_varying_T.csv",
    trap_2_densities_by_T,
    delimiter=",",
)
np.savetxt(
    results_folder + "trap_3_densities_varying_T.csv",
    trap_3_densities_by_T,
    delimiter=",",
)
np.savetxt(
    results_folder + "trap_4_densities_varying_T.csv",
    trap_4_densities_by_T,
    delimiter=",",
)
np.savetxt(
    results_folder + "trap_5_densities_varying_T.csv",
    trap_5_densities_by_T,
    delimiter=",",
)
np.savetxt(
    results_folder + "trap_6_densities_varying_T.csv",
    trap_6_densities_by_T,
    delimiter=",",
)
