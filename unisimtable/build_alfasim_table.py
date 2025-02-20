
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 18:02:24 2024

@author: rodri
"""

import win32com.client as winclt
from setup_gc_eos import pinit

from species_builder import R, Species, Mixture
from transfer_phenone_parameters import viscosity
from hydrate_module import *

unisim = winclt.Dispatch("UnisimDesign.Application")

# Abrindo o arquivo de simulação
file_path = r"C:\Users\rodri\Downloads\hisep_only_flash_get_data2.usc"
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


# COLUMNS = (PT, TM, ROG, ROHL, DROGDP, DROHLDP, DROGDT, DROHLDT, RS, VISG, VISHL, CPG, CPHL, HG, HHL, TCG, TCHL, SIGGHL)


P = [i for i in range(120,801,5)]
T = [i for i in range(40,131,3)]

pvt_table = []

feed_pressure = getattr(material_streams['feed'],'Pressure')
feed_temperature = getattr(material_streams['feed'],'Temperature')

vapor_feed_pressure = getattr(material_streams['vapor_feed'],'Pressure')
vapor_feed_temperature = getattr(material_streams['vapor_feed'],'Temperature')

feed_pressure.SetValue(10000)
feed_temperature.SetValue(100)

vapor_feed_pressure.SetValue(10000)
vapor_feed_temperature.SetValue(100)

vap_enthalpy_ref = material_streams['top'].MassEnthalpyValue*1000
liq_enthalpy_ref = material_streams['down'].MassEnthalpyValue*1000

pinit = gc_eos_class(mixture_nwe_3, 85+273.15, 1.2e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'liquid')

for Pin in P:
    for Tin in T:
        feed_pressure.SetValue(Pin*100)
        feed_temperature.SetValue(Tin)
        mass_flow_vap = material_streams['top'].MassFlow()
        mass_flow_liq = material_streams['down'].MassFlow()
        
        
        
        vap = pinit.copy_change_x_and_conditions(Tin+273.15,Pin*100,None,material_streams['top'].ComponentMolarFraction(),'gas')
        drogdp = vap.drhodP/1000
        drogdt = vap.drhodT
        vap.evaluate_der_rho()
        viscosity_vap = material_streams['top'].Viscosity()*0.001
        cp_vap = material_streams['top'].MassHeatCapacityValue*material_streams['top'].CpCv()*1000
        rho_vap = material_streams['top'].MassDensity()
        vap_k = material_streams['top'].ThermalConductivity()
        vap_enthalpy = material_streams['top'].MassEnthalpyValue*1000-vap_enthalpy_ref
        
        
        liq = pinit.copy_change_x_and_conditions(Tin+273.15,Pin*100,None,material_streams['down'].ComponentMolarFraction(),'liquid')
        drohldp = liq.drhodP/1000
        drohldt = liq.drhodT
        liq.evaluate_der_rho()
        viscosity_liq = material_streams['down'].Viscosity()*0.001
        cp_liq = material_streams['down'].MassHeatCapacityValue*material_streams['down'].CpCv()*1000
        rho_liq = material_streams['down'].MassDensity()
        liq_k = material_streams['down'].ThermalConductivity()
        liq_enthalpy = material_streams['down'].MassEnthalpyValue*1000-liq_enthalpy_ref
        
        
        surface_tension = 0.001*main_pfd.Operations.Item('Varificacao_K').Cell('B2').CellValue
        
        
        
        
        #rs = mass_flow_vap/(mass_flow_vap+mass_flow_liq)
        rs = mass_flow_vap/rho_vap/(mass_flow_vap/rho_vap+mass_flow_liq/rho_liq)
        
        
        pvt_row = [Pin*100*1e3, Tin, rho_vap, rho_liq, drogdp, drohldp, drogdt, drohldt, rs, 
                   viscosity_vap, viscosity_liq, cp_vap, cp_liq, vap_enthalpy, liq_enthalpy,
                   vap_k, liq_k, surface_tension]
        pvt_table.append(pvt_row)
        


