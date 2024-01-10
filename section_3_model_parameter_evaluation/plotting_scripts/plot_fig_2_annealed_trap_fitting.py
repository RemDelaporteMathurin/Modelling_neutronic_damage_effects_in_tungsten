from matplotlib import pyplot as plt
import numpy as np

"""
orginal data from A.Zaloznik et al, available at https://doi.org/10.1088/0031-8949/t167/1/014031
"""

test_temperatures = [370, 400, 500, 600, 800]
defect_type_1_densities = [0.230, 0.230, 0.225, 0.153, 0.107] # at.fr
defect_type_2_densities = [0.290, 0.290, 0.280, 0.280, 0.189]  # at.fr
defect_type_3_densities = [0.05, 0.05, 0.05, 0.05, 0.06]  # at.fr
annealing_time = 7200

# read fitting data
annealed_defect_type_1_densities = np.genfromtxt("../data/annealed_defect_1_densities.txt")
annealed_defect_type_2_densities = np.genfromtxt("../data/annealed_defect_2_densities.txt")
annealed_defect_type_3_densities = np.genfromtxt("../data/annealed_defect_3_densities.txt")

T_values = np.linspace(1, 900, num=1000)

# ##### Plotting ##### #

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

green_ryb = (117 / 255, 184 / 255, 42 / 255)
firebrick = (181 / 255, 24 / 255, 32 / 255)
electric_blue = (83 / 255, 244 / 255, 255 / 255)

plt.figure()

# defect type I
plt.scatter(
    test_temperatures,
    defect_type_1_densities,
    marker="x",
    color=green_ryb,
)
plt.plot(
    T_values,
    annealed_defect_type_1_densities,
    color=green_ryb,
    label=r"Defect type I",
)

# defect type II
plt.scatter(
    test_temperatures,
    defect_type_2_densities,
    marker="x",
    color=firebrick,
)
plt.plot(
    T_values,
    annealed_defect_type_2_densities,
    color=firebrick,
    label=r"Defect type II",
)

# defect type III
plt.scatter(
    test_temperatures,
    defect_type_3_densities,
    marker="x",
    color=electric_blue,
)
plt.plot(
    T_values,
    annealed_defect_type_3_densities,
    color=electric_blue,
    label=r"Defect type III",
)

plt.ylim(bottom=0)
plt.xlim(350, 850)
plt.ylabel(r"Trap density, n$_{\mathrm{t}}$ (m$^{-3}$)")
plt.legend()
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xlabel(r"Annealing temperature (K)")
plt.tight_layout()

plt.show()
