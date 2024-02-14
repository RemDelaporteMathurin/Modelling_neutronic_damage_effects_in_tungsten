from festim_model import festim_sim
import festim as F
import numpy as np
from scipy.integrate import odeint
from os.path import exists

# defect type I
A_0_1 = 6.1838e-03
E_A_1 = 0.24

trap_D1_K = 9.0e26
trap_D1_n_max = 6.9e25
trap_D2_K = 4.2e26
trap_D2_n_max = 7.0e25

# defect type II
A_0_2 = 6.1838e-03
E_A_2 = 0.30

trap_D3_K = 2.5e26
trap_D3_n_max = 6.0e25
trap_D4_K = 5.0e26
trap_D4_n_max = 4.7e25

# defect type III
A_0_3 = 0
E_A_3 = 1

trap_D5_K = 1.0e26
trap_D5_n_max = 2.0e25

# common values
k_B = F.k_B  # eV K-1
fpy = 3600 * 24 * 365
day = 3600 * 24

def generate_fig_7_inventory_transient_and_distribution():
    
    dpa_values = np.geomspace(1e-05, 1e02, 8)
    T = 700
    
    for dpa in dpa_values:
        n = 1000 # starting value
        for _ in range(15):
            try:
                print("running case T = {:.0f}, dpa = {:.1e}, n = {}".format(T, dpa, n))
                
                my_folder_name = "data/festim_model_results/dpa={:.1e}/T={:.0f}/".format(dpa, T)
                
                festim_sim(
                    dpa=dpa,
                    T=T,
                    results_folder_name=my_folder_name,
                    total_time=fpy,
                    cells=n,
                    export_retention_field=False
                )
                break
            except Exception as e:
                n *= 1.5
                n = int(n)
                print("increasing n to {}".format(n))
    # undamaged case
    festim_sim(
        dpa=0,
        T=T,
        results_folder_name="data/festim_model_results/dpa=0.0e+00/T={:.0f}/".format(T),
        total_time=fpy,
        cells=1000,
        export_retention_field=True
    )
    
    # profiles
    for dpa in dpa_values:
        n = 5000 # starting value
        for _ in range(15):
            try:
                print("running case dpa = {:.1e}, n = {}".format(dpa, n))
                
                my_folder_name = "data/profiles/data/dpa={:.1e}/".format(dpa)
                
                festim_sim(
                    dpa=dpa,
                    T=700,
                    results_folder_name=my_folder_name,
                    total_time=day,
                    cells=n,
                    export_retention_field=True
                )
                break
            except Exception as e:
                n *= 1.5
                n = int(n)
                print("increasing n to {}".format(n))
    
    # undamaged case
    festim_sim(
        dpa=0,
        T=700,
        results_folder_name="data/festim_model_results/dpa=0.0e+00/T={:.0f}/".format(T),
        total_time=day,
        cells=1000,
        export_retention_field=True
    )


def generate_fig_8_inventory_variataion():
    
    dpa_values = np.geomspace(1e-05, 1e02, 8)
    T_values = np.linspace(600, 1300, 50)

    for T in T_values:    
        for dpa in dpa_values:
            n = 5000 # starting value
            for _ in range(15):
                try:
                    print("running case T = {:.0f}, dpa = {:.1e}, n = {}".format(T, dpa, n))
                    
                    my_folder_name = "data/festim_model_results/dpa={:.1e}/T={:.0f}/".format(dpa, T)
                    
                    festim_sim(
                        dpa=dpa,
                        T=T,
                        results_folder_name=my_folder_name,
                        total_time=fpy,
                        cells=n,
                        export_retention_field=True
                    )
                    break
                except Exception as e:
                    n *= 1.5
                    n = int(n)
                    print("increasing n to {}".format(n))

        # undamaged case
        print("running case T = {:.0f}, dpa = 0.0e+00".format(T))
        my_folder_name = "data/festim_model_results/dpa=0.0e+00/T={:.0f}/".format(T)
        festim_sim(
            dpa=0,
            T=T,
            results_folder_name=my_folder_name,
            total_time=fpy,
            cells=1000,
            export_retention_field=False
        )

    
    
if __name__ == "__main__":
    generate_fig_7_inventory_transient_and_distribution()
    generate_fig_8_inventory_variataion()