import matplotlib.pyplot as plt
import numpy as np

x_0 = 2.3e-06
dx_0 = 1.0e-07
x_values = np.linspace(0, 1e-5, num=1000)

# distrubution
traps = 1 / (1 + np.exp((x_values - x_0) / dx_0))

# ##### plotting ##### #

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

plt.figure()
plt.plot(x_values, traps, color="black")
plt.xlim(0, 1e-05)
plt.ylim(0, 1.01)
plt.ylabel("Trap density ratio")
plt.xlabel("x (m)")
ax = plt.gca()
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.tight_layout()

plt.show()
