# VirCA (under development)
Virus Capsid Assembly (VirCA) is a computational method to predict the assembly of viral capsids (or other protein shells) based on properties of the protein units and solvent.

The project is under development by team members in the [Luque lab](https://www.luquelab.com). The project is available at the [GitHub repository](https://github.com/luquelab/VirCA).

## Python environment
Run the program using

    python3 VirCA.py

To recreate the environment used to test the code run

    conda env create -f environment.yml

To activate the environment, run

    conda activate env-VirCA

The environment was originally encoded on a MacBook Air M1 2020 running the commands

    conda create -n env-VirCA python=3.8.11 scipy=1.7.1 astroid babel

    conda install -c conda-forge matplotlib 

The environment was exported running the command

    conda env export > environment.yml

