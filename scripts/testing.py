import matplotlib.pyplot as plt
import numpy as np

dpa_values = np.geomspace(1e-05, 1e01, num=7)
gradients = [0.11, 0.11, 0.11, 0.11, 0.04, 0.0113, 0.006]
intercepts = [5.77, 5.49, 4.97, 4.49, 3.53, 2.93, 3.1]

test_dpa_values = np.geomspace(1e-05, 1e01, num=100)
gradient_fit = 0.016 * test_dpa_values ** (-0.434)
intercept_fit = -0.21 * np.log(test_dpa_values) + 3.6

# plt.figure()
# plt.scatter(dpa_values, gradients)
# plt.plot(test_dpa_values, gradient_fit)
# plt.xlim(5e-04, 1.5e01)
# plt.ylim(0, 0.15)
# plt.xscale("log")

plt.figure()
plt.scatter(dpa_values, intercepts)
plt.plot(test_dpa_values, intercept_fit)
plt.xscale("log")

plt.show()