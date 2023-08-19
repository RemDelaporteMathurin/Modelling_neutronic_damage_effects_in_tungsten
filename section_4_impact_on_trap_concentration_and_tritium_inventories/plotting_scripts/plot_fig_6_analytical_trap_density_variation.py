import matplotlib.pyplot as plt
import numpy as np

# import data
trap_D1_density_variation = np.genfromtxt(
    "../data/damage_trap_D1_density_variation.txt"
)
trap_D2_density_variation = np.genfromtxt(
    "../data/damage_trap_D2_density_variation.txt"
)
trap_D3_density_variation = np.genfromtxt(
    "../data/damage_trap_D3_density_variation.txt"
)
trap_D4_density_variation = np.genfromtxt(
    "../data/damage_trap_D4_density_variation.txt"
)
dpa_values = np.geomspace(1e-3, 1e03, num=1000)

# ##### plotting ##### #
plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

green_ryb = (117 / 255, 184 / 255, 42 / 255)
firebrick = (181 / 255, 24 / 255, 32 / 255)
pewter_blue = (113 / 255, 162 / 255, 182 / 255)
electric_blue = (83 / 255, 244 / 255, 255 / 255)

plt.figure()
plt.plot(
    dpa_values,
    trap_D1_density_variation,
    label=r"Trap D1",
    color=firebrick,
)
plt.plot(
    dpa_values,
    trap_D2_density_variation,
    label=r"Trap D2",
    color=pewter_blue,
)
plt.plot(
    dpa_values,
    trap_D3_density_variation,
    label=r"Trap D3",
    color=electric_blue,
)
plt.plot(
    dpa_values,
    trap_D4_density_variation,
    label=r"Trap D4",
    color=green_ryb,
)
plt.xlabel(r"Damage rate (dpa/fpy)")
plt.ylabel(r"Trap density (m$^{-3}$)")
plt.xlim(dpa_values[0], dpa_values[-1])
plt.ylim(1e23, 1e26)
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.tight_layout()
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.show()
