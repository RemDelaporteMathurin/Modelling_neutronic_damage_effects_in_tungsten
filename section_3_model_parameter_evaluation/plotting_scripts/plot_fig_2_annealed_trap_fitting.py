from matplotlib import pyplot as plt
import numpy as np

"""
orginal data from A.Zaloznik et al, available at https://doi.org/10.1088/0031-8949/t167/1/014031
subsequently fitted by E.Hodille et al, available at https://doi.org/10.1088/1741-4326/aa5aa5
"""

test_temperatures = [298, 600, 800, 1000, 1200]
trap_3_densities = [0.09, 0.08, 0.06, 0.00, 0.00]  # at.fr
trap_4_densities = [0.28, 0.23, 0.19, 0.15, 0.05]  # at.fr
annealing_time = 3600

# convert trap densities to m-3
trap_3_densities = (np.array(trap_3_densities) / 100) * 6.3e28
trap_4_densities = (np.array(trap_4_densities) / 100) * 6.3e28

# read fitting data
annealed_trap_3_densities = np.genfromtxt("../data/annealed_trap_3_densities.txt")
annealed_trap_4_densities = np.genfromtxt("../data/annealed_trap_4_densities.txt")

T_values = np.linspace(1, 1400, num=1000)

# ##### Plotting ##### #

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

green_ryb = (117 / 255, 184 / 255, 42 / 255)
firebrick = (181 / 255, 24 / 255, 32 / 255)

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(5, 7.5))

# trap 3
plt.sca(axs[0])
# ref values
plt.scatter(
    test_temperatures,
    trap_3_densities,
    marker="x",
    color=green_ryb,
)
# fitting
plt.plot(
    T_values,
    annealed_trap_3_densities,
    color=green_ryb,
    label=r"Trap 3",
)

# trap 4
plt.sca(axs[1])
# ref values
plt.scatter(
    test_temperatures,
    trap_4_densities,
    marker="x",
    color=firebrick,
)
# fitting
plt.plot(
    T_values,
    annealed_trap_4_densities,
    color=firebrick,
    label=r"Trap 4",
)

for ax in [axs[0], axs[1]]:
    plt.sca(ax)
    plt.ylim(bottom=0)
    plt.xlim(0, 1400)
    plt.ylabel(r"Trap density, n$_{\mathrm{t}}$ (m$^{-3}$)")
    plt.legend(loc="lower left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
plt.xlabel(r"Annealing temperature (K)")
plt.subplots_adjust(wspace=0.112, hspace=0.2)
plt.tight_layout()

plt.show()
