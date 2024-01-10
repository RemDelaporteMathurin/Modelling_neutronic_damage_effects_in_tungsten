import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

dpa_values_traps = np.geomspace(1e-03, 1e04, num=50)
T_values_traps = np.linspace(1300, 400, num=20)

dpa_values_inv = np.geomspace(1e-03, 1e03, num=50)
T_values_inv = np.linspace(1300, 400, num=50)
T_values_inv = T_values_inv[:34]

saturation_time_traps = np.loadtxt("../data/saturation_times_traps.txt")
saturation_times_inventories = np.loadtxt("../data/saturation_times_inventories.txt")

# ##### Plotting ##### #

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)


def plot_both_traps_and_inventories():
    norm = Normalize(vmin=min(T_values_traps), vmax=max(T_values_traps))
    colorbar = cm.inferno
    sm = plt.cm.ScalarMappable(cmap=colorbar, norm=norm)

    colours_traps = [colorbar(norm(T)) for T in T_values_traps]
    colours_invs = [colorbar(norm(T)) for T in T_values_inv]

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=([6.4, 9.6]))

    day = 3600 * 24
    month = day * 31
    year = day * 365.25

    for T, sat_times, colour in zip(
        T_values_traps, saturation_time_traps, colours_traps
    ):
        axs[0].plot(dpa_values_traps, sat_times, color=colour)

    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_ylabel(r"Trap density saturation time (s)")
    axs[0].spines["top"].set_visible(False)
    axs[0].spines["right"].set_visible(False)
    axs[0].set_xlim(1e-03, 1e03)
    axs[0].set_ylim(1e03, 1e07)

    axs[0].vlines(5, 1e03, 1e09, color="black", alpha=0.5, linestyle="dashed")
    axs[0].vlines(20, 1e03, 1e09, color="black", alpha=0.5, linestyle="dashed")
    axs[0].annotate("DEMO", [6e-01, 7.5e06], color="black", alpha=0.5)
    axs[0].annotate("ARC", [2.5e1, 7.5e06], color="black", alpha=0.5)

    ##########################

    for T, sat_times, colour in zip(
        T_values_inv, saturation_times_inventories, colours_invs
    ):
        axs[1].plot(dpa_values_inv, sat_times, color=colour)

    axs[1].set_xscale("log")
    axs[1].set_yscale("log")
    axs[1].spines["top"].set_visible(False)
    axs[1].spines["right"].set_visible(False)
    axs[1].set_ylabel(r"T retention saturation time (s)")
    axs[1].set_xlabel(r"Damage rate (dpa/fpy)")
    axs[1].set_xlim(1e-03, 1e03)
    axs[1].set_ylim(1e03, 1e09)

    axs[1].hlines(year, 1e-03, 5e04, color="black", alpha=0.5, linestyle="dotted")
    axs[1].annotate("1 FPY", [1.5e-03, year * 1.1], color="black", alpha=0.5)
    axs[1].vlines(5, 1e03, 1e09, color="black", alpha=0.5, linestyle="dashed")
    axs[1].vlines(20, 1e03, 1e09, color="black", alpha=0.5, linestyle="dashed")

    plt.tight_layout()

    fig.colorbar(sm, label=r"Temperature (K)", ax=axs, aspect=40)


plot_both_traps_and_inventories()

plt.show()
