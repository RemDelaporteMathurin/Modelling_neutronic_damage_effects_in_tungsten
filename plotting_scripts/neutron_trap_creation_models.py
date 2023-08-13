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
