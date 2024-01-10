import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

T_range = np.linspace(400, 1300, num=50)
dpa_range = np.geomspace(1e-5, 1e03, num=10)
T_range_contour = np.linspace(400, 1300, num=100)
dpa_range_contour = np.geomspace(1e-3, 1e03, num=100)

# import data
inventories = np.genfromtxt("../data/inventories.txt")
inventories_no_damage = np.genfromtxt("../data/inventories_no_damage.txt")
inventories_normalised_contour = np.genfromtxt(
    "../data/inventories_nomalised_contour.txt"
)

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

plot_dpa_range = dpa_range
norm = LogNorm(vmin=min(plot_dpa_range), vmax=max(plot_dpa_range))
colorbar = cm.viridis
sm = plt.cm.ScalarMappable(cmap=colorbar, norm=norm)

colours = [colorbar(norm(dpa)) for dpa in plot_dpa_range]

fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=([6.4, 9.6]))

for inv, colour in zip(inventories, colours):
    axs[0].plot(T_range, inv, color=colour)

axs[0].plot(T_range, inventories_no_damage, label=r"undamaged", color="black")
axs[0].set_ylabel(r"T inventory (m$^{-2}$)")
axs[0].set_xlim(400, 1300)
axs[0].set_ylim(1e16, 1e24)
axs[0].set_yscale("log")
axs[0].spines["top"].set_visible(False)
axs[0].spines["right"].set_visible(False)
axs[0].legend()

plt.colorbar(sm, label=r"Damage rate (dpa fpy$^{-1}$)", ax=axs[0])

dpa_range_contour = np.array(dpa_range_contour)
X, Y = np.meshgrid(T_range_contour, dpa_range_contour)

# ##### normalised to 0 dpa ##### #
CS = axs[1].contourf(
    X,
    Y,
    inventories_normalised_contour,
    norm=LogNorm(),
    levels=np.geomspace(
        np.min(inventories_normalised_contour),
        np.max(inventories_normalised_contour),
        num=1000,
    ),
    cmap="plasma",
)
for c in CS.collections:
    c.set_edgecolor("face")
CS2 = axs[1].contour(
    X,
    Y,
    inventories_normalised_contour,
    levels=[1e00, 1e01, 1e02, 1e03, 1e04, 1e05],
    colors="black",
)
axs[1].clabel(CS2, inline=True, fontsize=10, fmt="%.0e")

plt.colorbar(
    CS,
    label=r"Normalised T inventory (inv/inv$_{0 \ \mathrm{dpa}}$)",
    format="%.0e",
    ax=axs[1],
)

axs[1].set_yscale("log")
axs[1].set_ylabel(r"Damage rate (dpa fpy$^{-1}$)")
axs[1].set_xlabel(r"Temperature (K)")

plt.tight_layout()

plt.show()
