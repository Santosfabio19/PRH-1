# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 10:45:49 2025

@author: fabio
Simulação unisim setando apenas os valores da pump_2
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


compressor = CompressorClass()
compressor.change_parameters2()


def organize_separately(list1, list2):
    # Pair etas from list1 and list2
    eta_pairs = [[float(p1[0]*100), float(p2[0]*100)] for p1, p2 in zip(list1, list2)]
    
    # Pair heads from list1 and list2
    PressureIncrease_pairs = [[float(p1[1]), float(p2[1])] for p1, p2 in zip(list1, list2)]
    
    return eta_pairs, PressureIncrease_pairs



file_path_unisim =  r"C:\Users\fabio\folder\PRH-1\PRH-1\bombas-subsep\arquivos de simulação\nwe_sub_separtion_dynamic_set_pump2.usc"

name_streams = ['feed_sep','1st_heat_out','scrubber_in','pump_1_in',
                'pump_1_out','pump_2_in','pump_2_out']

name_valves = ['Valve_bypass_1', 'Valve_heat_1', 'Valve_bypass_2', 'Valve_heat_2']

name_pumps = ['pump_1', 'pump_2']


ns = [[ 70,60],
      [ 70,70],
      [80,70],
      [90,70],
      [100,80],
      [100,80]] # n1 C [70,100], n2 C [60,80]
    #ns = 
    #[[ 90,50],
    #[ 90,60],
    #[ 90,70],
    #[100,70],
    #[100,80],
    #[120,80]] # n1 C [70,100], n2 C [60,80]

nsim = (40*60)*6

unit_time = 'seconds'
sample_time = 10


########################################################################
unisim = winclt.Dispatch("UnisimDesign.Application")

simu2 = Simulator(file_path_unisim, unit_time, sample_time)    
Case = unisim.SimulationCases.Open(file_path_unisim)
#####################################################################################################
main_pfd = Case.Flowsheet
simu2.get_pumps(name_pumps) 
simu2.get_streams(name_streams)
simu2.get_valves(name_valves)
k = 0
#a = main_pfd.Operations['pump_1'].AdiabaticHead.SetValue
simu2.set_visible()

simu2.start_streams()


curves= []
arm_pump_1 = []
start = time.time()


while k <= nsim: # Alinhar k 
    
    if k < 10*60:
        i = 0
    elif k < ((10*60)*2):
        i = 1
    elif k < ((10*60)*3):
        i = 2
    elif k < ((10*60)*4):
        i = 3
    elif k < ((10*60)*5):
        i = 4
    else:
        i = 5
    
    comp = simu2.stream_states['pump_2_in']['x'][-1]
    dict_comp =  {list_names[i]: comp[i] for i in range(len(comp))}
    
    mixture = Mixture(list_of_species, dict_comp)
    
    #################################################################
    T = simu2.stream_states['pump_2_in']['T'][-1]
    P = simu2.stream_states['pump_2_in']['P'][-1]

############################################################################################3
    corrente =  gc_eos_class(mixture, T+273.15, P, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'gas')
    
    ##############################################################################
    
    Gi2 = corrente.copy_change_conditions(corrente.T, "icog", corrente.V*0.8, 'gas')



    vis_case_corrente = viscosity(corrente.mixture, [0]*19)          
    c = multi_stage_compression(10,corrente, compressor, vis_case_corrente)
    
    res_pump_1 =c.character_online(simu2.stream_states['pump_2_in']['MassFlow'][-1],ns[i][1] , corrente, Gi2) 
    
    
    arm_pump_1.append(res_pump_1)


    eta_values, pressure_values, head = res_pump_1
    simu2.new_eta_head((eta_values*100),pressure_values,"pump_2")
    simu2.pumps_new_value(ns[i],"pump_1")
    simu2.simulate_n_save_streams()    
    print(k)
    
    k += 1
    

end = time.time()
print(end-start)
simu2.plot_stream_state("scrubber_in","P")
simu2.plot_stream_state("scrubber_in","T")
simu2.plot_stream_state("pump_2_out","P")
simu2.plot_stream_state("pump_2_out","T")
simu2.plot_stream_state("pump_1_out","T")
simu2.plot_stream_state("1st_heat_out","T")
