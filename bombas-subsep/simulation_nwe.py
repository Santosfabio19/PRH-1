import win32com.client as winclt
import pandas as pd
from os import path, makedirs
from warnings import warn
import pickle
from connect_subsep import SimulationConnector

# UniSim File
file_path = r"C:\Users\rodri\OneDrive\Documentos\UniSim Design R492\Cases\nwe_sub_separtion_dynamic.usc"

# Folder to save date
folder = './'

continue_simulation = False
save_simulation = False
# Para continuar uma simulação anterior:
# 1 - executar primeira Simulação com continue_simulation = False e save_simulation = True
# 2 - executar segunda Simulação com continue_simulation = True e save_simulation = True

Simulation = SimulationConnector(file_path, folder, continue_simulation)

