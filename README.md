# Modelling_neutron_damage_effects_in_tungsten
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7863889.svg)](https://doi.org/10.5281/zenodo.7863889)

This repository contains all the scipts required to reproduce data for the paper "Modelling neutron damage effects in tungsten".

## Abstract
A damage-induced hydrogen trap creation model is proposed, and parameters for tungsten are identified using experimental data.
The methodology for obtaining these parameters using thermo-desorption analysis spectra data is outlined.
Self-damaged and optionally annealed tungsten samples have undergone TDS analysis, which has been analysed to identify the properties of extrinsic traps induced by the damage and to determine how they evolve with damage and annealing temperature.
A parametric study investigated the impact of the damage rate and temperature on tritium inventories in tungsten.
Tritium transport simulations have been performed with FESTIM considering a 1D model of a 2 mm sample of tungsten with damage rates and temperatures varying from 0-100 dpa/fpy and 600-1300 K, respectively. 
The results show that after 24 h simultaneous exposure to tritium implantation and neutron damage at 700 K , tritium inventories can increase by up to four orders of magnitude when damaged up to 100 dpa compared to an undamaged case, increasing further to five orders of magnitude after one full power year.
The time taken to reach saturation shows the need for kinetic models of trapping properties on time scales relevant to reactor operation.
The trap-creation model parameterisation procedure can be used to investigate neutron damage effects on other fusion-relevant materials such as EUROfer.

## Reproduce data
1. Clone this repository to your local machine.
2. Create a Conda environment using the provided `environment.yml` file:

    ```bash
    conda env create -f environment.yml
    ```

   This will set up a Conda environment named `festim-env` with all the required dependencies for running the FESTIM scripts and visualisation scripts.

3. Activate the Conda environment:

    ```bash
    conda activate festim-env
    ```

4. Execute the Python scripts using the activated Conda environment and ensure compatibility with FESTIM requirements.

5. Navigate to the desired folder based on the simulation you are interested in.

**Note**: To reproduce figure 7, [Paraveiw](https://www.paraview.org/) is required to obtain the retention profiles. This can be done using the [Plot over line](https://docs.paraview.org/en/latest/Tutorials/ClassroomTutorials/beginningPlotting.html) feature for each of the `retention.xdmf` files within `section_4_impact_on_trap_concentration/data/profiles/`.

## Contact

For any questions or issues, please contact james.dark@cea.fr.

