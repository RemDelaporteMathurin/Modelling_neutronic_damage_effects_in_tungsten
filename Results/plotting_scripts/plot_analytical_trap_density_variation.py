import matplotlib.pyplot as plt
import numpy as np

T_range = np.linspace(400, 1300, num=1000)
dpa_range = np.geomspace(1e-03, 1e03, num=1000)

results_folder = "../analytical_model_testing/"
trap_1_densities_standard_temp = np.genfromtxt(
    results_folder + "trap_1_densities_standard_T.csv", delimiter=","
)
trap_2_densities_standard_temp = np.genfromtxt(
    results_folder + "trap_2_densities_standard_T.csv", delimiter=","
)
trap_3_densities_standard_temp = np.genfromtxt(
    results_folder + "trap_3_densities_standard_T.csv", delimiter=","
)
trap_4_densities_standard_temp = np.genfromtxt(
    results_folder + "trap_4_densities_standard_T.csv", delimiter=","
)
trap_5_densities_standard_temp = np.genfromtxt(
    results_folder + "trap_5_densities_standard_T.csv", delimiter=","
)
trap_6_densities_standard_temp = np.genfromtxt(
    results_folder + "trap_6_densities_standard_T.csv", delimiter=","
)

trap_1_densities_by_T = np.genfromtxt(
    results_folder + "trap_1_densities_varying_T.csv", delimiter=","
)
trap_2_densities_by_T = np.genfromtxt(
    results_folder + "trap_2_densities_varying_T.csv", delimiter=","
)
trap_3_densities_by_T = np.genfromtxt(
    results_folder + "trap_3_densities_varying_T.csv", delimiter=","
)
trap_4_densities_by_T = np.genfromtxt(
    results_folder + "trap_4_densities_varying_T.csv", delimiter=","
)
trap_5_densities_by_T = np.genfromtxt(
    results_folder + "trap_5_densities_varying_T.csv", delimiter=","
)
trap_6_densities_by_T = np.genfromtxt(
    results_folder + "trap_6_densities_varying_T.csv", delimiter=","
)


green_ryb = (117 / 255, 184 / 255, 42 / 255)
firebrick = (181 / 255, 24 / 255, 32 / 255)
pewter_blue = (113 / 255, 162 / 255, 182 / 255)
blue_jeans = (96 / 255, 178 / 255, 229 / 255)
electric_blue = (83 / 255, 244 / 255, 255 / 255)

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)


fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=([12, 4.8]))
# variation with temp
axs[0].plot(T_range, trap_1_densities_by_T, label="Trap 1 (0.87 eV)", color="black")
axs[0].plot(T_range, trap_2_densities_by_T, label="Trap 2 (1.00 eV)", color="grey")
axs[0].plot(T_range, trap_3_densities_by_T, label="Trap D1 (1.15 eV)", color=firebrick)
axs[0].plot(
    T_range, trap_4_densities_by_T, label="Trap D2 (1.35 eV)", color=pewter_blue
)
axs[0].plot(
    T_range, trap_5_densities_by_T, label="Trap D3 (1.65 eV)", color=electric_blue
)
axs[0].plot(T_range, trap_6_densities_by_T, label="Trap D4 (1.85 eV)", color=green_ryb)
axs[0].set_ylabel("Trap density (m$^{-3}$)")
axs[0].set_xlabel("Temperature (K)")
axs[0].set_xlim(400, 1300)
axs[0].set_ylim(1e23, 1e26)
axs[0].set_yscale("log")
axs[0].spines["top"].set_visible(False)
axs[0].spines["right"].set_visible(False)

# damage variation
axs[1].plot(
    dpa_range,
    trap_1_densities_standard_temp,
    label="Trap 1 (0.87 eV)",
    color="black",
)
axs[1].plot(
    dpa_range,
    trap_2_densities_standard_temp,
    label="Trap 2 (1.00 eV)",
    color="grey",
)
axs[1].plot(
    dpa_range,
    trap_3_densities_standard_temp,
    label="Trap D1 (1.15 eV)",
    color=firebrick,
)
axs[1].plot(
    dpa_range,
    trap_4_densities_standard_temp,
    label="Trap D2 (1.35 eV)",
    color=pewter_blue,
)
axs[1].plot(
    dpa_range,
    trap_5_densities_standard_temp,
    label="Trap D3 (1.65 eV)",
    color=electric_blue,
)
axs[1].plot(
    dpa_range,
    trap_6_densities_standard_temp,
    label="Trap D4 (1.85 eV)",
    color=green_ryb,
)
axs[1].set_xlabel(r"Damage rate (dpa/fpy)")
axs[1].set_xlim(dpa_range[0], dpa_range[-1])
axs[1].set_ylim(1e23, 1e26)
axs[1].set_yscale("log")
axs[1].set_xscale("log")
plt.legend()
axs[1].spines["top"].set_visible(False)
axs[1].spines["right"].set_visible(False)

# plt.subplots_adjust(wspace=0.4)
plt.tight_layout()
plt.subplots_adjust(wspace=0.15)

plt.show()
