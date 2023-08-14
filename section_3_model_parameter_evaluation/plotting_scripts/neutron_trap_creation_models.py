import numpy as np

k_B = 8.617333e-05


def neutron_trap_creation_numerical(
    n, t, phi=9.64e-7, K=3.5e28, n_max=1e40, A_0=6.1838e-03, E_A=0.2792, T=298
):
    """
    Temporal evolution of n with resepct to time, for a given value of n and t

    Args:
        n (float): number of traps (m-3)
        t (float): time of simulation (s)
        phi (float): damage per second (dpa s-1). Defaults to 9.64e-07
        K (float): trap creation factor (traps dpa-1). Defaults to 1e28 traps
            s-1
        n_max (float): maximum traps per unit damage (m-3). Defaults to 1e40
            m-3
        A_0 (float): trap annealing factor (s-1). Defaults to 1e-02 s-1
        E_A (float): Annealing activation energy (eV). Defaults to 0.1034 eV
        T (float): the annealing temperature (K). Defaults to 298 K

    Returns:
        float: dn/dt
    """
    dndt = phi * K * (1 - (n / n_max)) - A_0 * np.exp(-E_A / (k_B * T)) * n
    return dndt


def neutron_trap_creation_analytical_steady(T, phi, K, n_max, A_0, E_A):
    """
    Function to evaluate the trap concentration at a given temperature, T, and
    damage, phi

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


def trap_density(T, phi, K, n_max, A_0, E_A):
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
    trap_d1_density = trap_density(
        T=T, phi=phi, K=1.5e28, n_max=5.2e25, A_0=A_0, E_A=E_A
    )
    trap_d2_density = trap_density(
        T=T, phi=phi, K=4.0e27, n_max=4.5e25, A_0=A_0, E_A=E_A
    )
    trap_d3_density = trap_density(
        T=T, phi=phi, K=3.0e27, n_max=4.0e25, A_0=A_0, E_A=E_A
    )
    trap_d4_density = trap_density(
        T=T, phi=phi, K=9.0e27, n_max=4.2e25, A_0=A_0, E_A=E_A
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
        A_0=A_0,
        E_A=E_A,
        n_i=trap_d1_density,
        T=T,
    )
    c_t_d2, filling_ratio_d2 = trapped_H_concentration(
        c_m=c_m,
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=1.35,
        A_0=A_0,
        E_A=E_A,
        n_i=trap_d2_density,
        T=T,
    )
    c_t_d3, filling_ratio_d3 = trapped_H_concentration(
        c_m=c_m,
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=1.65,
        A_0=A_0,
        E_A=E_A,
        n_i=trap_d3_density,
        T=T,
    )
    c_t_d4, filling_ratio_d4 = trapped_H_concentration(
        c_m=c_m,
        k_0=k_0,
        E_k=E_D,
        p_0=p_0,
        E_p=1.85,
        A_0=A_0,
        E_A=E_A,
        n_i=trap_d4_density,
        T=T,
    )

    trap_densities = [
        trap_1_density,
        trap_d1_density,
        trap_d2_density,
        trap_d3_density,
        trap_d4_density,
    ]

    trap_filling_ratios = [
        filling_ratio_1,
        filling_ratio_d1,
        filling_ratio_d2,
        filling_ratio_d3,
        filling_ratio_d4,
    ]

    total_trapped_H_concentration = c_t_1 + c_t_d1 + c_t_d2 + c_t_d3 + c_t_d4

    total_retention = retention(c_m=c_m, c_t=total_trapped_H_concentration, L=L)

    return (total_retention, trap_densities, trap_filling_ratios)
