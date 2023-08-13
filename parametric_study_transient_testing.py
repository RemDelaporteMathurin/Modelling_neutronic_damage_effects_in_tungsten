from neutron_induced_traps_model import festim_sim
import numpy as np
from os.path import exists


dpa_values = np.geomspace(1e-03, 1e03, num=50)
temperature_values = np.linspace(400, 1300, num=50)


def case_steady():
    for temperature in temperature_values:
        for dpa in dpa_values:
            print("Running case: dpa={:.2e}, T={:.0f}".format(dpa, temperature))
            festim_sim(
                T=temperature,
                dpa=dpa,
                results_folder_name="Results/parametric_studies/case_steady/dpa={:.2e}/T={:.0f}/".format(
                    dpa, temperature
                ),
                transient_run=False,
            )


def case_1e09s():
    for temperature in temperature_values:
        for dpa in dpa_values:
            print("Running case: dpa={:.2e}, T={:.0f}".format(dpa, temperature))
            festim_sim(
                T=temperature,
                dpa=dpa,
                results_folder_name="Results/parametric_studies/case_1e09s/dpa={:.2e}/T={:.0f}/".format(
                    dpa, temperature
                ),
                transient_run=True,
                total_time=1e09,
            )


def case_24h():
    dpa_values = np.geomspace(1e-05, 1e03, num=17)
    temperature_values = np.linspace(400, 1300, num=50)

    results_folder = "Results/parametric_studies/case_24h/"

    for temperature in temperature_values:
        for dpa in dpa_values:
            filename_test = (
                results_folder
                + "dpa={:.2e}/T={:.0f}/derived_quantities.csv".format(dpa, temperature)
            )
            if exists(filename_test):
                continue
            else:
                print("Running case: dpa={:.2e}, T={:.0f}".format(dpa, temperature))
                festim_sim(
                    T=temperature,
                    dpa=dpa,
                    results_folder_name=results_folder
                    + "dpa={:.2e}/T={:.0f}/".format(dpa, temperature),
                    transient_run=True,
                    total_time=24 * 3600,
                )


case_steady()
case_1e09s()
case_24h()