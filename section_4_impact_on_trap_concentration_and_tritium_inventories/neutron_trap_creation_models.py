import numpy as np

T_range = np.linspace(400, 1300, num=50)
dpa_range = np.geomspace(1e-5, 1e03, num=10)
T_range_contour = np.linspace(400, 1300, num=100)
dpa_range_contour = np.geomspace(1e-3, 1e03, num=100)

# ##### standard values ##### #

L = 0.002
k_B = 8.6173303e-5
imp_flux = 1e20
r_p = 3e-09
fpy = 3600 * 24 * 365.25

# diffusion parameters
D_0 = 4.1e-7
E_D = 0.28

# damaged trap parameters
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
# A_0_3 = A_0_2
E_A_3 = 1
# E_A_3 = E_A_2

trap_D5_K = 1.2e26
trap_D5_n_max = 2.0e25

# general trap parameters
k_0 = 5.22e-17
p_0 = 1e13
trap_1_density = 2e22


def neutron_trap_creation_analytical_steady(T, phi, K, n_max, A_0, E_A):
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


def analytical_model(phi=9 / fpy, T=700, L=0.002):
    trap_d1_density = neutron_trap_creation_analytical_steady(
        T=T, phi=phi, K=trap_D1_K, n_max=trap_D1_n_max, A_0=A_0_1, E_A=E_A_1
    )
    trap_d2_density = neutron_trap_creation_analytical_steady(
        T=T, phi=phi, K=trap_D2_K, n_max=trap_D2_n_max, A_0=A_0_1, E_A=E_A_1
    )
    trap_d3_density = neutron_trap_creation_analytical_steady(
        T=T, phi=phi, K=trap_D3_K, n_max=trap_D3_n_max, A_0=A_0_2, E_A=E_A_2
    )
    trap_d4_density = neutron_trap_creation_analytical_steady(
        T=T, phi=phi, K=trap_D4_K, n_max=trap_D4_n_max, A_0=A_0_2, E_A=E_A_2
    )
    trap_d5_density = neutron_trap_creation_analytical_steady(
        T=T, phi=phi, K=trap_D5_K, n_max=trap_D5_n_max, A_0=A_0_3, E_A=E_A_3
    )

    c_m = mobile_H_concentration(T=T, imp_flux=imp_flux, r_p=r_p, D_0=D_0, E_D=E_D)

    c_t_1, filling_ratio_1 = trapped_H_concentration(
        c_m=c_m, k_0=k_0, E_k=E_D, p_0=p_0, E_p=0.87, n_i=trap_1_density, T=T
    )
    c_t_d1, filling_ratio_d1 = trapped_H_concentration(
        c_m=c_m,
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=1.15,
        A_0=A_0_1,
        E_A=E_A_1,
        n_i=trap_d1_density,
        T=T,
    )
    c_t_d2, filling_ratio_d2 = trapped_H_concentration(
        c_m=c_m,
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=1.35,
        A_0=A_0_1,
        E_A=E_A_1,
        n_i=trap_d2_density,
        T=T,
    )
    c_t_d3, filling_ratio_d3 = trapped_H_concentration(
        c_m=c_m,
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=1.65,
        A_0=A_0_2,
        E_A=E_A_2,
        n_i=trap_d3_density,
        T=T,
    )
    c_t_d4, filling_ratio_d4 = trapped_H_concentration(
        c_m=c_m,
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=1.85,
        A_0=A_0_2,
        E_A=E_A_2,
        n_i=trap_d4_density,
        T=T,
    )
    c_t_d5, filling_ratio_d5 = trapped_H_concentration(
        c_m=c_m,
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=2.05,
        A_0=A_0_3,
        E_A=E_A_3,
        n_i=trap_d5_density,
        T=T,
    )

    trap_densities = [
        trap_1_density,
        trap_d1_density,
        trap_d2_density,
        trap_d3_density,
        trap_d4_density,
        trap_d5_density,
    ]

    trap_filling_ratios = [
        filling_ratio_1,
        filling_ratio_d1,
        filling_ratio_d2,
        filling_ratio_d3,
        filling_ratio_d4,
        filling_ratio_d5,
    ]

    total_trapped_H_concentration = c_t_1 + c_t_d1 + c_t_d2 + c_t_d3 + c_t_d4 + c_t_d5

    total_retention = retention(c_m=c_m, c_t=total_trapped_H_concentration, L=L)

    return (total_retention, trap_densities, trap_filling_ratios)


def numerical_trap_creation_model(n, t, phi, K, n_max, A_0, E_A, T):

        dndt = phi * K * (1 - (n / n_max)) - A_0 * np.exp(-E_A / (k_B * T)) * n

        return dndt