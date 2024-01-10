import matplotlib.pyplot as plt
import numpy as np

"""
TDS data from T.Swartz-Selinger, currently unpublished
"""
dpa_values = [0, 0.001, 0.005, 0.023, 0.1, 0.23, 0.5, 2.5]

# trap density variations, fitted from TDS data
trap_D1_densities = [0, 4.5e24, 7.0e24, 2.4e25, 5.4e25, 5.8e25, 6.0e25, 6.8e25]
trap_D2_densities = [0, 1.0e24, 2.5e24, 1.4e25, 3.8e25, 4.4e25, 4.8e25, 6.1e25]
trap_D3_densities = [0, 5.0e23, 1.0e24, 6.0e24, 2.8e25, 3.5e25, 4.4e25, 5.1e25]
trap_D4_densities = [0, 1.0e24, 1.9e24, 2.1e25, 3.6e25, 4.0e25, 4.2e25, 4.9e25]
trap_D5_densities = [0, 2.0e23, 1.6e24, 6.0e24, 1.1e25, 1.4e25, 1.8e25, 2.0e25]

# read fitting data
trap_D1_fitting = np.genfromtxt("../data/damage_trap_D1_fitting.txt")
trap_D2_fitting = np.genfromtxt("../data/damage_trap_D2_fitting.txt")
trap_D3_fitting = np.genfromtxt("../data/damage_trap_D3_fitting.txt")
trap_D4_fitting = np.genfromtxt("../data/damage_trap_D4_fitting.txt")
trap_D5_fitting = np.genfromtxt("../data/damage_trap_D5_fitting.txt")
dpa_x_values = np.linspace(0, 3, num=100)

# ##### plotting ##### #
plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

green_ryb = (117 / 255, 184 / 255, 42 / 255)
firebrick = (181 / 255, 24 / 255, 32 / 255)
pewter_blue = (113 / 255, 162 / 255, 182 / 255)
electric_blue = (83 / 255, 244 / 255, 255 / 255)

fig = plt.figure(figsize=(4.5, 12))

# create one big plot to have a common y label
ax = fig.add_subplot(111)
ax.spines["top"].set_color("none")
ax.spines["bottom"].set_color("none")
ax.spines["left"].set_color("none")
ax.spines["right"].set_color("none")
ax.tick_params(labelcolor="w", top=False, bottom=False, left=False, right=False)
ax.set_ylabel(r"Trap density (m$^{-3}$)")

# highlight trap 1
ax1 = fig.add_subplot(511)
plt.sca(ax1)
plt.scatter(dpa_values, trap_D1_densities, marker="x", color=firebrick)
plt.plot(dpa_x_values, trap_D1_fitting, color=firebrick, label=r"Trap D1")
plt.plot(dpa_x_values, trap_D2_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D3_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D4_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D4_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D5_fitting, color="grey", alpha=0.2)
plt.legend(loc="lower right")

# highlight trap D2
ax2 = fig.add_subplot(512)
plt.sca(ax2)
plt.scatter(dpa_values, trap_D2_densities, color=pewter_blue, marker="x")
plt.plot(dpa_x_values, trap_D2_fitting, color=pewter_blue, label=r"Trap D2")
plt.plot(dpa_x_values, trap_D1_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D3_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D4_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D5_fitting, color="grey", alpha=0.2)
plt.legend(loc="lower right")

# highlight trap D3
ax3 = fig.add_subplot(513)
plt.sca(ax3)
plt.scatter(dpa_values, trap_D3_densities, color=electric_blue, marker="x")
plt.plot(dpa_x_values, trap_D3_fitting, color=electric_blue, label=r"Trap D3")
plt.plot(dpa_x_values, trap_D1_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D2_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D4_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D5_fitting, color="grey", alpha=0.2)
plt.legend(loc="lower right")

# highlight trap D4
ax4 = fig.add_subplot(514)
plt.sca(ax4)
plt.scatter(dpa_values, trap_D4_densities, color=green_ryb, marker="x")
plt.plot(dpa_x_values, trap_D4_fitting, color=green_ryb, label=r"Trap D4")
plt.plot(dpa_x_values, trap_D1_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D2_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D3_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D5_fitting, color="grey", alpha=0.2)
plt.legend(loc="lower right")

# highlight trap D5
ax5 = fig.add_subplot(515)
plt.sca(ax5)
plt.scatter(dpa_values, trap_D5_densities, color=green_ryb, marker="x")
plt.plot(dpa_x_values, trap_D5_fitting, color="black", label=r"Trap D5")
plt.plot(dpa_x_values, trap_D1_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D2_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D3_fitting, color="grey", alpha=0.2)
plt.plot(dpa_x_values, trap_D4_fitting, color="grey", alpha=0.2)
plt.legend(loc="lower right")
plt.xlabel(r"Damage (dpa)")

for ax in [ax1, ax2, ax3, ax4, ax5]:
    plt.sca(ax)
    ax.get_shared_x_axes().join(ax, ax5)
    plt.ylim(0, 8e25)
    plt.xlim(left=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# remove the xticks for top plots
ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.set_xticklabels([])
ax4.set_xticklabels([])

plt.tight_layout()

plt.show()
