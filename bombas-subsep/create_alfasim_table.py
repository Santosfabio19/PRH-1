# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 18:02:24 2024

@author: rodri
"""

import win32com.client as winclt

unisim = winclt.Dispatch("UnisimDesign.Application")

# Abrindo o arquivo de simulação
file_path = r"C:\Users\rodri\Downloads\hisep_only_flash_get_data.usc"
Case = unisim.SimulationCases.Open(file_path)


main_pfd = Case.Flowsheet  # Carregando o flowsheet principal
print(f'Os subflowsheets encontrados em {main_pfd.TaggedName} são:')

#%% Acessando todas as correntes do flowsheet principal
# Correntes de massa
material_streams = {}
print('As correntes de massa são:') 
for name in main_pfd.MaterialStreams.Names:
    material_streams[name] = main_pfd.MaterialStreams.Item(name)
    print(f"{name}")


