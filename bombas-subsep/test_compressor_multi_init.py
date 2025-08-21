# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 13:59:02 2025

@author: fabio
"""

from eos_database_new_resumed import *
from gc_eos import gc_eos_class
from solver_thermo import solver_eos
from transfer_phenone_parameters import viscosity
from compressor_class import CompressorClass
from compression import compression
import time 
from compression import multi_stage_compression
from setup_simulation import *

file_path_unisim = r"C:\Users\fabio\Downloads\bombas_novo_tema\model\nwe_sub_separtion_dynamic(1).usc"

composition = [1]*30

list_names_unisim = ["C15+", "CH4", "C2H6", "C3H8", "iC4H10", "nC4H10", "iC5H12", "nC5H12", "nC6H14", "nC7H16", "nC8H18",
              "nC9H20", "nC10H22", "N2", "CO2","H2O", "nC11H24", "nC12H26", "nC14H30"]

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

stream_nwe_3 = gc_eos_class(mixture_nwe_3, 85+273.15, 1.2e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'liquid')

solver_stream_nwe_3 = solver_eos(stream_nwe_3)

x0 = [0] + [4]*14 + [0] + [0] + [0] + [18]
y0 = [2]*5 + [0]*13 + [0]*2 + [7] + [0] + [40] + [0]
w0 = [1e-10]*15 + [0] + [1] + [0]*2

stream_nwe_33 = gc_eos_class(mixture_nwe_3, 80+273.15, 1.55e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'gas')

phi3, stream_nwe_3_vap3, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_33, [], [], 0.6)

compressor = CompressorClass()
compressor.change_parameters()

Gi1 = stream_nwe_33.copy_change_conditions(323, 'icog', stream_nwe_33.V*0.8, 'gas')







arm_pump_1 = []
arm_pump_2 = []
    
    
 
def organize_separately(list1, list2):
    # Pair etas from list1 and list2
    eta_pairs = [[float(p1[
        0]*100), float(p2[0]*100)] for p1, p2 in zip(list1, list2)]
    
    # Pair heads from list1 and list2
    PressureIncrease_pairs = [[float(p1[1]), float(p2[1])] for p1, p2 in zip(list1, list2)]
    
    return eta_pairs, PressureIncrease_pairs



    ########################################################################
name_streams = ['feed_sep','1st_heat_out','scrubber_in','pump_1_in',
                'pump_1_out','pump_2_in','pump_2_out']

name_valves = ['Valve_bypass_1', 'Valve_heat_1', 'Valve_bypass_2', 'Valve_heat_2']

name_pumps = ['pump_1', 'pump_2']


ns = [[ 90,50],
      [ 90,60],
      [ 90,70],
      [100,70],
      [100,80],
      [120,80]]

nsim = 80*6

unit_time = 'minutes'
sample_time = 1


########################################################################
unisim = winclt.Dispatch("UnisimDesign.Application")

simu = Simulator(file_path_unisim, unit_time, sample_time)    
Case = unisim.SimulationCases.Open(file_path_unisim)
#####################################################################################################
main_pfd = Case.Flowsheet
simu.get_pumps(name_pumps) 
simu.get_streams(name_streams)
simu.get_valves(name_valves)
k = 0
#a = main_pfd.Operations['pump_1'].AdiabaticHead.SetValue
simu.set_visible()

simu.start_streams()


curves= []
start = time.time()



#corrente =  gc_eos_class(mixture, T+273.15, P, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'gas')
#corrente_pump_2 =  gc_eos_class(mixture2, T_pump_2+273.15, P_pump_2, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'gas')
    
    
#phi1, corrente_vap, liqui = solver1.evaluate_flash(corrente, [],[],0.6)
#phi2, corrente_vap_pump_2, liqui2 = solver2.evaluate_flash(corrente_pump_2, [],[],0.6)
while k <= nsim:
    
    if k < 80:
        i = 0
    elif k < (80*2):
        i = 1
    elif k < (80*3):
        i = 2
    elif k < (80*4):
        i = 3
    elif k < (80*5):
        i = 4
    else:
        i = 5
    
    comp = simu.stream_states['pump_1_in']['x'][-1]
    comp_pump_2 = simu.stream_states['pump_2_in']['x'][-1]
    dict_comp =  {list_names[i]: comp[i] for i in range(len(comp))}
    dict_comp2 =  {list_names[i]: comp_pump_2[i] for i in range(len(comp_pump_2))}
    
    mixture = Mixture(list_of_species, dict_comp)
    mixture2 = Mixture(list_of_species, dict_comp2)
    
    #################################################################
    T = simu.stream_states['pump_1_in']['T'][-1]
    P = simu.stream_states['pump_1_in']['P'][-1]
    T_pump_2 = simu.stream_states['pump_2_in']['T'][-1]
    P_pump_2 = simu.stream_states['pump_2_in']['P'][-1]
    
    ###########################################3
    # Talvez criar o metodo copy and change composition
    corrente_1 =  gc_eos_class(mixture, T+273.15, P, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'gas')
    corrente_2 =  gc_eos_class(mixture2, T_pump_2+273.15, P_pump_2, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'gas')
    solver1 = solver_eos(corrente_1)
    solver2 = solver_eos(corrente_2)
############################################################################################3
    corrente =  gc_eos_class(mixture, T+273.15, P, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'gas')
    corrente_pump_2 =  gc_eos_class(mixture2, T_pump_2+273.15, P_pump_2, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'gas')
    
    ##############################################################################
    
    Gi2 = corrente.copy_change_conditions(corrente.T, "icog", corrente.V*0.8, 'gas')
    Gi3 = corrente_pump_2.copy_change_conditions(corrente_pump_2.T, "icog", corrente_pump_2.V*0.8, 'gas')



    vis_case_corrente = viscosity(corrente.mixture, [0]*19)          
    vis_case_corrente_pump_2 = viscosity(corrente.mixture, [0]*19)          
    c = multi_stage_compression(10,corrente, compressor, vis_case_corrente)
    c2 = multi_stage_compression(10,corrente_pump_2, compressor, vis_case_corrente_pump_2)
    
    res_pump_1 =c.character_online(simu.stream_states['pump_1_in']['MassFlow'][-1],ns[i][0] , corrente, Gi2) 
    res_pump_2 =c2.character_online(simu.stream_states['pump_2_in']['MassFlow'][-1],ns[i][1],corrente_pump_2, Gi3)
    
    
    #####################################'
    arm_pump_1.append(res_pump_1)

    arm_pump_2.append(res_pump_2)
    #arm_mass_flow_pump_1.append(simu.Streams['pump_1_in'].MassFlow())
    #arm_mass_flow_pump_2.append(simu.Streams['pump_2_in'].MassFlow())
    eta_values, pressure_values = organize_separately(arm_pump_1, arm_pump_2)
    simu.new_eta_head(eta_values[i],pressure_values[i])
    simu.simulate_n_save_streams()    
    
    k += 1

end = time.time()
print(end-start)

simu.plot_stream_state("scrubber_in","P")
simu.plot_stream_state("scrubber_in","T")
simu.plot_stream_state("pump_2_out","P")
simu.plot_stream_state("pump_2_out","T")
simu.plot_stream_state("pump_1_out","T")
simu.plot_stream_state("1st_heat_out","T")




