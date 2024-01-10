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


def annealing_sim(A_0, E_A, n_0, T, t):
    """
    Runs a numerical model of annealing effects on traps induced by neutron
    damage

    Args:
        A_0 (float): trap annealing factor (s-1).
        E_A (float): Annealing activation energy (eV).

    Returns:
        (list): A list of trap densities at various annealing tempertures
    """
    k_B = 8.617333e-05
    A = A_0 * np.exp((-E_A) / (k_B * T))
    annealed_trap_densities = n_0 * np.exp(-A * t)

    return annealed_trap_densities
