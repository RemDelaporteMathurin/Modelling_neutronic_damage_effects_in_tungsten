import matplotlib.pyplot as plt
import numpy as np
import csv
from matplotlib import cm
from matplotlib.colors import LogNorm

"""
TDS data from T.Swartz-Selinger, currently unpublished
"""
implantation_time = 72 * 3600
resting_time = 0.5 * 24 * 3600
sample_area = 12e-03 * 15e-03
dpa_values = [2.5, 0.5, 0.23, 0.1, 0.023, 0.005, 0.001, 0]

# ##### get tds and fitting data #### #


def retrieve_tds_and_fitting_data():
    tds_T_and_flux = []
    fitting_data = []

    for dpa in dpa_values:
        tds_data_file = "../data/tds_data_schwartz_selinger/{}_dpa.csv".format(dpa)
        tds_data = np.genfromtxt(tds_data_file, delimiter=",", dtype=float)
        tds_data_T = tds_data[:, 0]
        tds_data_flux = tds_data[:, 1] / sample_area
        tds_T_and_flux.append([tds_data_T, tds_data_flux])

        results_folder = "../data/damaged_sample_tds_fittings/"
        fitting_file = results_folder + "dpa_{}/last.csv".format(dpa)
        T_sim = []
        flux_1 = []
        flux_2 = []
        with open(fitting_file, "r") as csvfile:
            plots = csv.reader(csvfile, delimiter=",")
            for row in plots:
                if "t(s)" not in row:
                    if float(row[0]) >= implantation_time + resting_time * 0.75:
                        T_sim.append(float(row[1]))
                        flux_1.append(float(row[2]))
                        flux_2.append(float(row[3]))
        flux = -np.asarray(flux_1) - np.asarray(flux_2)
        fitting_data.append([T_sim, flux])

    return tds_T_and_flux, fitting_data


tds_T_and_flux, fitting_data = retrieve_tds_and_fitting_data()


# ##### plotting ##### #

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

# colourbar
plot_dpa_values = dpa_values[:-1]
norm = LogNorm(vmin=min(plot_dpa_values), vmax=max(plot_dpa_values))
colorbar = cm.viridis
sm = plt.cm.ScalarMappable(cmap=colorbar, norm=norm)
colours = [colorbar(norm(dpa)) for dpa in dpa_values]

plt.figure(figsize=(6.4, 5.5))

plt.plot(
    tds_T_and_flux[-1][0],
    tds_T_and_flux[-1][1],
    label=r"undamaged",
    linewidth=3,
    color="black",
)

for case, colour in zip(tds_T_and_flux[:-1], colours):
    T_values = case[0]
    flux_values = case[1]
    plt.plot(T_values, flux_values, color=colour, linewidth=3)

for case in fitting_data[:-1]:
    T_values = case[0]
    flux_values = case[1]
    plt.plot(
        T_values,
        flux_values,
        color="grey",
        linewidth=2,
        linestyle="dashed",
        alpha=0.5,
    )
plt.plot(
    fitting_data[-1][0],
    fitting_data[-1][1],
    color="grey",
    linewidth=2,
    linestyle="dashed",
    alpha=0.5,
    label="fittings",
)

plt.xlim(300, 1000)
plt.ylim(0, 1e17)
plt.xlabel(r"Temperature (K)")
plt.ylabel(r"Desorption flux (D m$ ^{-2}$ s$ ^{-1}$)")
plt.legend()
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.subplots_adjust(wspace=0.112, hspace=0.071)
plt.colorbar(sm, label=r"Damage (dpa)")

plt.show()
