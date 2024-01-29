import matplotlib.pyplot as plt
import numpy as np
import csv

"""
TDS data from T.Swartz-Selinger, currently unpublished
"""
implantation_time = 72 * 3600
resting_time = 0.5 * 24 * 3600
area = 12e-03 * 15e-03

# reference TDS
data_0_1 = np.genfromtxt(
    "../data/tds_data_schwartz_selinger/0.1_dpa.csv", delimiter=","
)
T_0_1 = data_0_1[:, 0]
flux_0_1 = data_0_1[:, 1] / area

# obtain tds fitting data
T_sim = []
flux1, flux2 = [], []
solute = []
retention = []
t = []
trap_1 = []
trap_D1 = []
trap_D2 = []
trap_D3 = []
trap_D4 = []
trap_D5 = []

data = np.genfromtxt(
    "../data/damaged_sample_tds_fittings/dpa_0.1/last.csv",
    delimiter=",",
    names=True,
)

t = data["ts"]

idx = np.where(
    t >= implantation_time + resting_time * 0.75
)  # only consider t > 0.75 * resting_time + implantation_time
t = t[idx]
T_sim = data["Average_T_volume_1"][idx]
flux1 = data["Flux_surface_1_solute"][idx]
flux2 = data["Flux_surface_2_solute"][idx]
solute = data["Total_solute_volume_1"][idx]
retention = data["Total_retention_volume_1"][idx]
trap_1 = data["Total_1_volume_1"][idx]
trap_D1 = data["Total_2_volume_1"][idx]
trap_D2 = data["Total_3_volume_1"][idx]
trap_D3 = data["Total_4_volume_1"][idx]
trap_D4 = data["Total_5_volume_1"][idx]
trap_D5 = data["Total_6_volume_1"][idx]


trap_1_contribution = (np.diff(trap_1) / np.diff(t)) * -1
trap_D1_contribution = (np.diff(trap_D1) / np.diff(t)) * -1
trap_D2_contribution = (np.diff(trap_D2) / np.diff(t)) * -1
trap_D3_contribution = (np.diff(trap_D3) / np.diff(t)) * -1
trap_D4_contribution = (np.diff(trap_D4) / np.diff(t)) * -1
trap_D5_contribution = (np.diff(trap_D5) / np.diff(t)) * -1
solute_contribution = (np.diff(solute) / np.diff(t)) * -1


# ##### plotting ##### #

# plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

green_ryb = (117 / 255, 184 / 255, 42 / 255)
firebrick = (181 / 255, 24 / 255, 32 / 255)
pewter_blue = (113 / 255, 162 / 255, 182 / 255)
electric_blue = (83 / 255, 244 / 255, 255 / 255)

plt.figure(figsize=(6.4, 5.5))

plt.scatter(T_0_1, flux_0_1, label=r"Exp", color="black")
plt.plot(
    T_sim,
    -np.asarray(flux1) - np.asarray(flux2),
    label="Simulation",
    color="orange",
    linewidth=2,
)
plt.plot(
    T_sim[1:],
    trap_1_contribution,
    linestyle="dashed",
    color="grey",
    label=r"Trap 1",
)
plt.plot(
    T_sim[1:],
    trap_D1_contribution,
    linestyle="dashed",
    color=firebrick,
    label=r"Trap D1",
)
plt.plot(
    T_sim[1:],
    trap_D2_contribution,
    linestyle="dashed",
    color=pewter_blue,
    label=r"Trap D2",
)
plt.plot(
    T_sim[1:],
    trap_D3_contribution,
    linestyle="dashed",
    color=electric_blue,
    label=r"Trap D3",
)
plt.plot(
    T_sim[1:],
    trap_D4_contribution,
    linestyle="dashed",
    color=green_ryb,
    label=r"Trap D4",
)
plt.plot(
    T_sim[1:],
    trap_D5_contribution,
    linestyle="dashed",
    color="black",
    label=r"Trap D5",
)
plt.fill_between(T_sim[1:], 0, trap_1_contribution, color="grey", alpha=0.1)
plt.fill_between(T_sim[1:], 0, trap_D1_contribution, color="grey", alpha=0.1)
plt.fill_between(T_sim[1:], 0, trap_D2_contribution, color="grey", alpha=0.1)
plt.fill_between(T_sim[1:], 0, trap_D3_contribution, color="grey", alpha=0.1)
plt.fill_between(T_sim[1:], 0, trap_D4_contribution, color="grey", alpha=0.1)
plt.fill_between(T_sim[1:], 0, trap_D5_contribution, color="grey", alpha=0.1)

plt.xlim(300, 1000)
plt.ylim(0, 1.0e17)
plt.xlabel("Temperature (K)")
plt.ylabel(r"Desorption flux (D m$ ^{-2}$ s$ ^{-1}$)")
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.legend()
plt.tight_layout()

plt.show()
