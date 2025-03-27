
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 18:02:24 2024

@author: rodri
"""

import win32com.client as winclt
from setup_gc_eos import gc_eos_class
from eos_database_new_resumed import *
import csv
from species_builder import R, Species, Mixture
from solver_thermo import *
#from transfer_phenone_parameters import viscosity
#from hydrate_module import *
from Reorder import *  


unisim = winclt.Dispatch("UnisimDesign.Application")

# Abrindo o arquivo de simulação
#file_path = r"C:\Users\fabio\Documents\GitHub\PRH-1\unisimtable\nwe_sub_separtion_isolated - change_composition.usc" # PC protec
file_path = r"C:\Users\fabio\projects\PRH-1\unisimtable\nwe_sub_separtion_isolated - change_composition.usc"  # Pessoal 
Case = unisim.SimulationCases.Open(file_path)


main_pfd = Case.Flowsheet  # Carregando o flowsheet principal
print(f'Os subflowsheets encontrados em {main_pfd.TaggedName} são:')


list_names  = ["CH4",	  "C2H6",	  "C3H8",	  "iC4H10",  "nC4H10",  "iC5H12",  "nC5H12",  "nC6H14",  "nC7H16",  "nC8H18",	
               "nC9H20",  "nC10H22", "nC11H24", "nC12H26", "nC14H30", "N2", "H2O", "CO2", "C15+"]

list_names_unisim = ["C15+", "CH4", "C2H6", "C3H8", "iC4H10", "nC4H10", "iC5H12", "nC5H12", "nC6H14", "nC7H16", "nC8H18",
              "nC9H20", "nC10H22", "N2", "CO2","H2O", "nC11H24", "nC12H26", "nC14H30"]


nwe_3 = [0.238095238095238,0.0262608309264608,0.0261719617862697,0.00648744723394801,
         0.0183070428793601,0.0112419462341702,0.0105754276827372,0.0193290379915574,
         0.0210708965703968,0.0190209179995643,0.0173567422429736,0.01597710063164,
         0.0148135692205084,0.0138181989380497,0.0122023902400072] + \
        [0.000844257, 0.002962305, 0.446300822] + [0.066207509]
        
nwe_4 = [0.00028290304301490683,0.34587337798028883,0.02791144428568154,0.020996469524193925,
         0.0042042015790704395,0.010633837951991394,0.005135960232073656,0.004416221726969988,0.005740061651108351,
         0.004491845632858752,0.002897243207370615,0.0019277373790942429,0.0012719502139215899,0.0014454143220052808,
         0.5586533294833951,0.002389788697400429,0.0008630438912499545,0.0006126407468087032,0.00025252845150226154]        
        


#dict_composition_3 = {list_names[i]: nwe_3[i] for i in range(len(nwe_3))}

#mixture_nwe_3 = Mixture(list_of_species, dict_composition_3)

volumn_desviation = [0]*19

x0 = [0] + [4]*14 + [0] + [0] + [0] + [18]
y0 = [2]*5 + [0]*13 + [0]*2 + [7] + [0] + [40] + [0]
w0 = [1e-10]*15 + [0] + [1] + [1] + [0]







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


P = [i for i in range(120,805,3)]   # 120,801,5    # 22 
T = [i for i in range(40,135,3)]    # 40,131, 3  # 12 
print(len(T))
print(len(P))


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

liquid_, vapor_, both = 0,0,0

comp = material_streams['vapor_vapor'].ComponentMolarFraction()

reorder_obj = Reorder(list_names_unisim, list(comp))
new_data1, new_data2 = reorder_obj.reorder()


dict_composition_3 = {list_names[i]: new_data2[i] for i in range(len(new_data2))}
#dict_composition_4 = {list_names_unisim[i]: nwe_4[i] for i in range(len(nwe_4))}
#dict_composition_5 = {list_names_unisim[i]: comp[i] for i in range(len(comp))}


mixture_nwe_3 = Mixture(list_of_species, dict_composition_3)


pinit = gc_eos_class(mixture_nwe_3, 85+273.15, 1.2e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'liquid')





feed = material_streams['vapor_feed'].ComponentMolarFraction()
reorder_obj = Reorder(list_names_unisim, list(feed))
new_data1, new_data2 = reorder_obj.reorder()

vap_input = pinit.copy_change_x_and_conditions(40 +273.15, 801*100,None,feed,'gas')

DENSITY = []
Pbub = []
Tm = []

for PT in P:
    for TM in T:
        vapor_feed_pressure.SetValue(PT*100)
        vapor_feed_temperature.SetValue(TM)
        mass_flow_vap = material_streams['vapor_vapor'].MassFlow()
        mass_flow_liq = material_streams['vapor_liquid'].MassFlow()
        Tm.append(TM)
        
        
        if mass_flow_vap > 0 and mass_flow_liq <= 0  :
            comp = material_streams['vapor_vapor'].ComponentMolarFraction()
            reorder_obj = Reorder(list_names_unisim, list(comp))
            _, comp = reorder_obj.reorder()
            
            
            vap = pinit.copy_change_x_and_conditions(TM+273.15,PT*100,None,comp,'gas')
            #a = vap.mass_rho
            #DENSITY.append(a)
            print("Apenas Vapor")
            
            solver_vap = solver_eos(vap_input)
            solver_vap.set_estimated_conditions(x0, y0, w0)
            Tbub = TM
            Pbub_, _,_ = solver_vap.bubble_T(TM+273.15, 1000, y0)
            Pbub.append(Pbub_)
            
            DROGDP = vap.evaluate_drhodP()/1000
            
            DROGDT = vap.evaluate_drhodT()
            vap.evaluate_der_rho()
            CPG = material_streams['vapor_vapor'].MassHeatCapacityValue*material_streams['vapor_vapor'].CpCv()*1000
            ROG = material_streams['vapor_vapor'].MassDensityValue
            VISG = material_streams['vapor_vapor'].Viscosity() * 10e-3 
            TCG = material_streams['vapor_vapor'].ThermalConductivity()
            HG = material_streams['vapor_vapor'].MassEnthalpyValue/1000-vap_enthalpy_ref
            SEG = 100
            ################################################
            DROHLDP = vap.evaluate_drhodP()/1000
            DROHLDT = vap.evaluate_drhodT()
            VISHL = material_streams['vapor_vapor'].Viscosity() * 10e-3 
            CPHL = material_streams['vapor_vapor'].MassHeatCapacityValue*material_streams['vapor_vapor'].CpCv()*1000
            ROHL = material_streams['vapor_vapor'].MassDensityValue
            TCHL = material_streams['vapor_vapor'].ThermalConductivity()
            HHL = material_streams['vapor_vapor'].MassEnthalpyValue*1000-vap_enthalpy_ref
            SIGGHL = 0.
            SEHL = 100 
            vapor_ += 1 
            
        
        elif mass_flow_liq > 0 and mass_flow_vap <= 0  :
            reorder_obj = Reorder(list_names_unisim, list(comp))
            _, comp = reorder_obj.reorder()
            
            
            liq = pinit.copy_change_x_and_conditions(TM+273.15,PT*100,None,comp,'liquid')
            #DENSITY.append(liq.mass_rho)
            ComMolFrac = material_streams['vapor_liquid'].ComponentMolarFraction()
            DROHLDP = liq.evaluate_drhodP()/1000
            print("Apenas Liquido")
            
            
            solver_vap = solver_eos(vap_input)
            solver_vap.set_estimated_conditions(x0, y0, w0)
            Tbub = TM
            Pbub_, _,_ = solver_vap.bubble_T(TM +273.15, 1000, y0)
            Pbub.append(Pbub_)
            
            
            DROHLDT = liq.evaluate_drhodT()
            liq.evaluate_der_rho()
            VISHL = material_streams['vapor_liquid'].Viscosity()* 10e-3 
            CPHL = material_streams['vapor_liquid'].MassHeatCapacityValue*material_streams['vapor_liquid'].CpCv()*1000
            ROHL = material_streams['vapor_liquid'].MassDensityValue
            TCHL = material_streams['vapor_liquid'].ThermalConductivity()
            HHL = material_streams['vapor_liquid'].MassEnthalpyValue*1000-liq_enthalpy_ref
            SIGGHL = material_streams["vapor_liquid"].SurfaceTension()
            SEHL = 100
            ############
            DROGDP = vap.evaluate_drhodP()/1000
            DROGDT = vap.evaluate_drhodT()
            vap.evaluate_der_rho()
            VISG = material_streams['vapor_liquid'].Viscosity()*10e-2
            CPG = material_streams['vapor_liquid'].MassHeatCapacityValue*material_streams['vapor_vapor'].CpCv()*1000
            ROG = material_streams['vapor_liquid'].MassDensity()
            TCG = material_streams['vapor_liquid'].ThermalConductivity()
            HG = material_streams['vapor_liquid'].MassEnthalpyValue*1000-liq_enthalpy_ref
            SEG = 100.0
            liquid_ += 1 
        else:
            reorder_obj = Reorder(list_names_unisim, list(comp))
            _, comp = reorder_obj.reorder()
            vap = pinit.copy_change_x_and_conditions(TM+273.15,PT*100,None,comp,'gas')
            DENSITY.append(vap.mass_rho)
            
            
            DROGDP = vap.evaluate_drhodP()/1000
            print("Vapor e Liquido")
            DROGDT = vap.evaluate_drhodT()
            vap.evaluate_der_rho()
            VISG = material_streams['vapor_vapor'].Viscosity()*10e-2
            CPG = material_streams['vapor_vapor'].MassHeatCapacityValue*material_streams['vapor_vapor'].CpCv()*1000
            ROG = material_streams['vapor_vapor'].MassDensity()
            TCG = material_streams['vapor_vapor'].ThermalConductivity()
            HG = material_streams['vapor_vapor'].MassEnthalpyValue*1000-vap_enthalpy_ref
            SEG = 100.0
            
            liq = pinit.copy_change_x_and_conditions(TM+273.15,PT*100,None,comp,'liquid')
            
            DROHLDP = liq.evaluate_drhodP()/1000
            DROHLDT = liq.evaluate_drhodT()
            liq.evaluate_der_rho()
            VISHL = material_streams['vapor_liquid'].Viscosity()* 10e-3 
            CPHL = material_streams['vapor_liquid'].MassHeatCapacityValue*material_streams['vapor_liquid'].CpCv()*1000
            ROHL = material_streams['vapor_liquid'].MassDensityValue  # Pode ser densidade molar
            TCHL = material_streams['vapor_liquid'].ThermalConductivity()
            HHL = (material_streams['vapor_liquid'].MassEnthalpyValue*1000-liq_enthalpy_ref) 
            SIGGHL = material_streams["vapor_liquid"].SurfaceTension()
            SEHL = 100.0
            
            
            #surface_tension = material_streams["vapor_liquid"].SurfaceTension          #0.001*main_pfd.Operations.Item('V-100-2')
            both += 1 
        
        
        
        RS = mass_flow_vap/(mass_flow_vap+mass_flow_liq)
        #rs = mass_flow_vap/rho_vap/(mass_flow_vap/rho_vap+mass_flow_liq/rho_liq)
        
        
        #pvt_row = [Pin*100*1e3, Tin, rho_vap, rho_liq, drogdp, drohldp, drogdt, drohldt, rs, 
        #           viscosity_vap, viscosity_liq, cp_vap, cp_liq, vap_enthalpy, liq_enthalpy,
        #           vap_k, liq_k, surface_tension]
        
        pvt_row = [PT*100*1e3,TM,ROG,ROHL,DROGDP,DROHLDP,DROGDT,DROHLDT,
                   RS,VISG,VISHL,CPG,CPHL,HG,HHL,TCG,
                   TCHL,SIGGHL,SEG,SEHL]
            
        
        pvt_table.append(pvt_row)
        print(vapor_,liquid_,both)
        


# %%


path_saida = r'C:\Users\fabio\projects\PRH-1\unisimtable\tabelas\pvt_table_32x228.tab'

with open(path_saida, mode='w', newline='') as file:
    # Escrever o cabeçalho
     # Escrever o cabeçalho
    columns = ["PT", "TM", "ROG", "ROHL", "DROGDP", "DROHLDP", "DROGDT", "DROHLDT", "RS", "VISG", "VISHL", "CPG", "CPHL", "HG", "HHL", "TCG", "TCHL", "SIGGHL", "SEG", "SEHL"]
    COLUMNS_row = "COLUMNS = (" + ",".join(columns) + ")\n"

    
    # Escrever as linhas de dados no formato desejado
    COMPONENT_row = "COMPONENTS = (" + ",".join(f'"{component}"' for component in list_names) + "),\\"    
    MOLES_row = "MOLES = (" + ",".join(map(str, new_data2)) + "),\\"
    MOL = critical_table[12,:]
    MOLWEIGHT_row = "MOLWEIGHT = (" + ",".join(map(str, MOL)) + ") g/mol,\\"
    DENSITY_row = "DENSITY = (" + ",".join(map(str, DENSITY)) + ") g/cm3,\\"
    
    
    
    STDTEMPERATURE = 0.288710e+03  # K
    STDTEMPERATURE_row = f"STDTEMPERATURE = {STDTEMPERATURE} K,\\"
    GOR = .382991E+03  # Sm3/Sm3
    GOR_row = f"GOR = {GOR} Sm3/Sm3,\\"
    GLR = .382991E+03  # Sm3/Sm3
    GLR_row = f"GLR = {GLR} Sm3/Sm3,\\"
    
    STDGASDENSITY = vap.rho  # kg/m3
    
    STDGASDENSITY_row = f"STDGASDENSITY = {STDGASDENSITY} kg/m3,\\"
 
    STDOILDENSITY = liq.rho  # kg/m3
    
    STDOILDENSITY_row = f"STDOILDENSITY = {STDOILDENSITY} kg/m3,\\"
    
    CRITICALPRESSURE = .100000E+01  # ATM
    CRITICALPRESSURE_row = f"CRITICALPRESSURE = {CRITICALPRESSURE} ATM,\\"
    
    CRITICALTEMPERATURE = .288710E+03  # 
    CRITICALTEMPERATURE_row = f"CRITICALTEMPERATURE = {CRITICALTEMPERATURE} K,\\"
    
    MESHTYPE = "STANDARD"
    MESHTYPE_row = f"MESHTYPE = {MESHTYPE},\\"
    
    STDPRESSURE = 1.4444444444e+01  # ATM
    STDPRESSURE_row = f"STDPRESSURE = {STDPRESSURE} ATM,\\"
    PRESSURE_row = "PRESSURE = (" + ",".join(map(str, P)) + ") Bar,\\"
    TEMPERATURE_row = "TEMPERATURE = (" + ",".join(map(str, T)) + ") C,\\"
    #BUBBLEPRESSURES = T
    BUBBLEPRESSURES_row = "BUBBLEPRESSURES = (" + ",".join(map(str, Pbub)) + ") Pa,\\"
    #BUBBLETEMPERATURES = 
    BUBBLETEMPERATURES_row = "BUBBLETEMPERATURES= (" + ",".join(map(str, Tm)) + ") C,\\"
        
    
    
    
    
    
    PRESSURE_row = "PRESSURE = (" + ",".join(map(str, P)) + ") Pa,\\"
    TEMPERATURE_row = "TEMPERATURE = (" + ",".join(map(str, T)) + ") C,\\"
    
    
    
    
    
    
    
    file.write('PVTTABLE LABEL = "NONEQ TABSAT5", PHASE = TWO,\\ \n')    
    file.write(COMPONENT_row + "\n")
    #file.write(MOLES_row + '\n')
    #file.write(MOLWEIGHT_row + '\n')
    #file.write(DENSITY_row + '\n')
    file.write(STDPRESSURE_row + '\n')
    file.write(STDTEMPERATURE_row + '\n')
    #file.write(GOR_row + '\n')
    #file.write(GLR_row + '\n')
    #file.write(STDGASDENSITY_row + '\n')
    #file.write(STDOILDENSITY_row + '\n')
    #file.write(CRITICALPRESSURE_row + '\n')
    #file.write(CRITICALTEMPERATURE_row + '\n')
    file.write(MESHTYPE_row + '\n')
    #file.write(STDPRESSURE_row + '\n')
    file.write(PRESSURE_row + '\n')
    file.write(TEMPERATURE_row + '\n')
    file.write(BUBBLEPRESSURES_row + "\n")
    file.write(BUBBLETEMPERATURES_row + "\n")
    
    file.write(COLUMNS_row)
    
    for row in pvt_table:
        try:
            # Converte todos os valores para float (sem formatação científica)
            formatted_values = [str(float(value)) for value in row]
            formatted_row = "PVTTABLE POINT = (" + ",".join(formatted_values) + ")\n"
            file.write(formatted_row)
        except (ValueError, TypeError) as e:
            print(f"Erro ao processar a linha: {row}. Erro: {e}")

print(f'Tabela salva em {path_saida}')
        

        
        


