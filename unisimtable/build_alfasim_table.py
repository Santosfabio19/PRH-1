
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 18:02:24 2024

@author: rodri
"""

import win32com.client as winclt
from setup_gc_eos import gc_eos_class
from eos_database_new_resumed import *

from species_builder import R, Species, Mixture
#from transfer_phenone_parameters import viscosity
#from hydrate_module import *

unisim = winclt.Dispatch("UnisimDesign.Application")

# Abrindo o arquivo de simulação
#file_path = r"C:\Users\fabio\Documents\GitHub\PRH-1\unisimtable\nwe_sub_separtion_isolated - change_composition.usc" # PC protec
file_path = r"C:\Users\fabio\projects\PRH-1\unisimtable\nwe_sub_separtion_isolated - change_composition.usc"  # Pessoal 
Case = unisim.SimulationCases.Open(file_path)


main_pfd = Case.Flowsheet  # Carregando o flowsheet principal
print(f'Os subflowsheets encontrados em {main_pfd.TaggedName} são:')


list_names  = ["CH4",	  "C2H6",	  "C3H8",	  "iC4H10",  "nC4H10",  "iC5H12",  "nC5H12",  "nC6H14",  "nC7H16",  "nC8H18",	
               "nC9H20",  "nC10H22", "nC11H24", "nC12H26", "nC14H30", "N2", "H2O", "CO2", "C15+"]

nwe_3 = [0.238095238095238,0.0262608309264608,0.0261719617862697,0.00648744723394801,
         0.0183070428793601,0.0112419462341702,0.0105754276827372,0.0193290379915574,
         0.0210708965703968,0.0190209179995643,0.0173567422429736,0.01597710063164,
         0.0148135692205084,0.0138181989380497,0.0122023902400072] + \
        [0.000844257, 0.002962305, 0.446300822] + [0.066207509]

dict_composition_3 = {list_names[i]: nwe_3[i] for i in range(len(nwe_3))}

mixture_nwe_3 = Mixture(list_of_species, dict_composition_3)

volumn_desviation = [0]*19


#%% Acessando todas as correntes do flowsheet principal
# Correntes de massa
material_streams = {}
for name in main_pfd.MaterialStreams.Names:
    material_streams[name] = main_pfd.MaterialStreams.Item(name)
    
    
for op in main_pfd.Operations:
    print(op.Name)

#for it in main_pfd.Operations.Item("V-100-2"):
 #   print(it.)
    


# COLUMNS = (PT, TM, ROG, ROHL, DROGDP, DROHLDP, DROGDT, DROHLDT, RS, VISG, VISHL, CPG, CPHL, HG, HHL, TCG, TCHL, SIGGHL)


P = [i for i in range(120,801,5)]   # 120,801,5
T = [i for i in range(40,131,3)]    # 40,131, 3 

pvt_table = []

#feed_pressure = getattr(material_streams['feed'],'Pressure')
#feed_temperature = getattr(material_streams['feed'],'Temperature')

vapor_feed_pressure = getattr(material_streams['vapor_feed'],'Pressure')
vapor_feed_temperature = getattr(material_streams['vapor_feed'],'Temperature')

#feed_pressure.SetValue(10000)
#feed_temperature.SetValue(100)

vapor_feed_pressure.SetValue(10000)
vapor_feed_temperature.SetValue(100)

vap_enthalpy_ref = material_streams['vapor_vapor'].MassEnthalpyValue*1000
liq_enthalpy_ref = material_streams['vapor_liquid'].MassEnthalpyValue*1000

pinit = gc_eos_class(mixture_nwe_3, 85+273.15, 1.2e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'liquid')
liquid_, vapor_, both = 0,0,0

for Pin in P:
    for Tin in T:
        vapor_feed_pressure.SetValue(Pin*100)
        vapor_feed_temperature.SetValue(Tin)
        mass_flow_vap = material_streams['vapor_vapor'].MassFlow()
        mass_flow_liq = material_streams['vapor_liquid'].MassFlow()
        
        
        if mass_flow_vap > 0 and mass_flow_liq <= 0  :
            vap = pinit.copy_change_x_and_conditions(Tin+273.15,Pin*100,None,material_streams['vapor_vapor'].ComponentMolarFraction(),'gas')
            drogdp = vap.evaluate_drhodP()/1000
            print("Apenas Vapor")
            drogdt = vap.evaluate_drhodT()
            vap.evaluate_der_rho()
            viscosity_vap = material_streams['vapor_vapor'].Viscosity()*0.001
            cp_vap = material_streams['vapor_vapor'].MassHeatCapacityValue*material_streams['vapor_vapor'].CpCv()*1000
            rho_vap = material_streams['vapor_vapor'].MassDensity()
            vap_k = material_streams['vapor_vapor'].ThermalConductivity()
            vap_enthalpy = material_streams['vapor_vapor'].MassEnthalpyValue*1000-vap_enthalpy_ref
            ################################################
            drohldp = 0
            drohldt = 0
            viscosity_liq = 0
            cp_liq = 0
            rho_liq = 0
            liq_k = 0
            liq_enthalpy = 0
            surface_tension = 0
            vapor_ += 1 
            
        
        elif mass_flow_liq > 0 and mass_flow_vap <= 0  :
            liq = pinit.copy_change_x_and_conditions(Tin+273.15,Pin*100,None,material_streams['vapor_liquid'].ComponentMolarFraction(),'liquid')
            drohldp = liq.evaluate_drhodP()/1000
            print("Apenas Liquido")
            drohldt = liq.evaluate_drhodT()
            liq.evaluate_der_rho()
            viscosity_liq = material_streams['vapor_liquid'].Viscosity()*0.001
            cp_liq = material_streams['vapor_liquid'].MassHeatCapacityValue*material_streams['vapor_liquid'].CpCv()*1000
            rho_liq = material_streams['vapor_liquid'].MassDensity()
            liq_k = material_streams['vapor_liquid'].ThermalConductivity()
            liq_enthalpy = material_streams['vapor_liquid'].MassEnthalpyValue*1000-liq_enthalpy_ref
            ###############################
            surface_tension = material_streams["vapor_liquid"].SurfaceTension
            drogdp = 0
            drogdt = 0
            viscosity_vap = 0
            cp_vap = 0
            rho_vap = 0
            vap_k = 0
            vap_enthalpy = 0
            liquid_ += 1 
            
        else:
            vap = pinit.copy_change_x_and_conditions(Tin+273.15,Pin*100,None,material_streams['vapor_vapor'].ComponentMolarFraction(),'gas')
            drogdp = vap.evaluate_drhodP()/1000
            print("Vapor e Liquido")
            drogdt = vap.evaluate_drhodT()
            vap.evaluate_der_rho()
            viscosity_vap = material_streams['vapor_vapor'].Viscosity()*0.001
            cp_vap = material_streams['vapor_vapor'].MassHeatCapacityValue*material_streams['vapor_vapor'].CpCv()*1000
            rho_vap = material_streams['vapor_vapor'].MassDensity()
            vap_k = material_streams['vapor_vapor'].ThermalConductivity()
            vap_enthalpy = material_streams['vapor_vapor'].MassEnthalpyValue*1000-vap_enthalpy_ref
            liq = pinit.copy_change_x_and_conditions(Tin+273.15,Pin*100,None,material_streams['vapor_liquid'].ComponentMolarFraction(),'liquid')
            drohldp = liq.evaluate_drhodP()/1000
            drohldt = liq.evaluate_drhodT()
            liq.evaluate_der_rho()
            viscosity_liq = material_streams['vapor_liquid'].Viscosity()*0.001
            cp_liq = material_streams['vapor_liquid'].MassHeatCapacityValue*material_streams['vapor_liquid'].CpCv()*1000
            rho_liq = material_streams['vapor_liquid'].MassDensity()
            liq_k = material_streams['vapor_liquid'].ThermalConductivity()
            liq_enthalpy = material_streams['vapor_liquid'].MassEnthalpyValue*1000-liq_enthalpy_ref
            surface_tension = material_streams["vapor_liquid"].SurfaceTension          #0.001*main_pfd.Operations.Item('V-100-2')
            both += 1 
        
        
        
        rs = mass_flow_vap/(mass_flow_vap+mass_flow_liq)
        #rs = mass_flow_vap/rho_vap/(mass_flow_vap/rho_vap+mass_flow_liq/rho_liq)
        
        
        pvt_row = [Pin*100*1e3, Tin, rho_vap, rho_liq, drogdp, drohldp, drogdt, drohldt, rs, 
                   viscosity_vap, viscosity_liq, cp_vap, cp_liq, vap_enthalpy, liq_enthalpy,
                   vap_k, liq_k, surface_tension]
        pvt_table.append(pvt_row)
        print(vapor_,liquid_,both)
        
import csv 
path_saida = r'C:\Users\fabio\projects\PRH-1\unisimtable\tabelas\pvt_table.csv'
with open(path_saida, mode='w', newline='') as file:
    writer = csv.writer(file, delimiter=";")
    # Escrever o cabeçalho (se necessário)
    writer.writerow(["Pin", "Tin", "rho_vap", "rho_liq", "drogdp", "drohldp", "drogdt", "drohldt", "rs", 
                     "viscosity_vap", "viscosity_liq", "cp_vap", "cp_liq", "vap_enthalpy", "liq_enthalpy",
                     "vap_k", "liq_k", "surface_tension"])
    # Escrever as linhas de dados
    writer.writerows(pvt_table)

print(f'Tabela salva em {path_saida}')
        
        
        


