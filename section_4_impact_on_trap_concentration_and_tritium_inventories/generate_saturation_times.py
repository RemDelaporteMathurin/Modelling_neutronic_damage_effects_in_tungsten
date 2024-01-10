import numpy as np
from scipy.integrate import odeint

# ##### trap saturation evaluation ##### #

k_B = 8.617333e-05
fpy = 86400 * 365
A_0 = 6.1838e-03
E_A = 0.2792
n_0 = 0

t = np.geomspace(1e01, 1e09, int(1e04))

dpa_values_traps = np.geomspace(1e-03, 1e04, num=50)
T_values_traps = np.linspace(1300, 400, num=20)

dpa_values_inv = np.geomspace(1e-03, 1e03, num=50)
T_values_inv = np.linspace(1300, 400, num=50)
T_values_inv = T_values_inv[:34]

trap_densities = []

for T in T_values_traps:
    traps_per_T = []
    for dpa in dpa_values_traps:
        phi = dpa / fpy

        def numerical_trap_creation_model_trap_1(n, t, phi=phi, A_0=A_0, E_A=E_A, T=T):
            K = 1.5e28
            n_max = 5.2e25
            dndt = phi * K * (1 - (n / n_max)) - A_0 * np.exp(-E_A / (k_B * T)) * n
            return dndt

        trap_1 = odeint(
            numerical_trap_creation_model_trap_1,
            n_0,
            t,
        )
        traps_per_T.append(trap_1)
    trap_densities.append(traps_per_T)

normalised_traps = []
for T_case in trap_densities:
    normalised_traps_per_T = []
    for dpa_case in T_case:
        norm_values = dpa_case / dpa_case[-1]
        normalised_traps_per_T.append(norm_values)
    normalised_traps.append(normalised_traps_per_T)

saturation_time_traps = []
for T_case in normalised_traps:
    char_times_per_T = []
    for dpa_case in T_case:
        char_time = np.where(dpa_case > 0.99)
        first = char_time[0][0]
        char_times_per_T.append(t[first])
    saturation_time_traps.append(char_times_per_T)

np.savetxt("data/saturation_times_traps.txt", saturation_time_traps)

inventories = []
inventories_steady = []
ts = []

results_folder = "data/case_1e09s/"
steady_results_folder = "data/case_steady/"

for T in T_values_inv:
    invs_per_T = []
    ts_per_T = []
    steady_invs_per_T = []
    for dpa in dpa_values_inv:
        data_file = (
            results_folder + "dpa={:.2e}/T={:.0f}/derived_quantities.csv".format(dpa, T)
        )
        steady_data_file = (
            steady_results_folder
            + "dpa={:.2e}/T={:.0f}/derived_quantities.csv".format(dpa, T)
        )
        data = np.genfromtxt(data_file, delimiter=",", names=True)
        steady_data = np.genfromtxt(steady_data_file, delimiter=",", names=True)
        invs_per_T.append(data["Total_retention_volume_1"])
        steady_invs_per_T.append(steady_data["Total_retention_volume_1"])
        ts_per_T.append(data["ts"])
    inventories.append(invs_per_T)
    inventories_steady.append(steady_invs_per_T)
    ts.append(ts_per_T)

normalised_inventories = []
for T_case, T_case_steady in zip(inventories, inventories_steady):
    normalised_invs_per_T = []
    for dpa_case, dpa_case_steady in zip(T_case, T_case_steady):
        norm_values = dpa_case / dpa_case_steady
        # norm_values = dpa_case / dpa_case[-1]
        normalised_invs_per_T.append(norm_values)
    normalised_inventories.append(normalised_invs_per_T)

characteristic_times = []
for T_case, T_t in zip(normalised_inventories, ts):
    char_times_per_T = []
    for dpa_case, t in zip(T_case, T_t):
        char_time = np.where(dpa_case > 0.95)
        # char_time = np.where(dpa_case > 0.99)
        first = char_time[0][0]
        char_times_per_T.append(t[first])
    characteristic_times.append(char_times_per_T)


# export data
np.savetxt("data/saturation_times_inventories.txt", characteristic_times)
