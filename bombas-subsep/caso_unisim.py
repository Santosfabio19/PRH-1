# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:57:15 2025

@author: rodri
"""

from setup_simulation import *

file_path_unisim = r"C:\Users\fabio\Downloads\bombas_novo_tema\model\nwe_sub_separtion_dynamic(1).usc"

name_streams = ['feed_sep','1st_heat_out','scrubber_in','pump_1_in',
                'pump_1_out','pump_2_in','pump_2_out']

name_valves = ['Valve_bypass_1', 'Valve_heat_1', 'Valve_bypass_2', 'Valve_heat_2']

name_pumps = ['pump_1', 'pump_2']

nsim = 80*6

unit_time = 'minutes'
sample_time = 1

simu = Simulator(file_path_unisim, unit_time, sample_time)

simu.get_pumps(name_pumps) 
simu.get_streams(name_streams)
simu.get_valves(name_valves)

k = 0

simu.set_visible()

simu.start_streams()

ns = [[ 90,50],
      [ 90,60],
      [ 90,70],
      [100,70],
      [100,80],
      [120,80]] # n1 C [70,100], n2 C [60,80]

nv = [[50.0, 50.0, 50.0, 50.0],
      [50.0, 50.0, 50.0, 50.0],
      [50.0, 50.0, 50.0, 50.0],
      [50.0, 50.0, 50.0, 50.0],
      [50.0, 50.0, 50.0, 50.0],
      [50.0, 50.0, 50.0, 50.0]]

curves= []

while k <= nsim:
    
    if k < 80:
        i = 0
    elif k < 80*2:
        i = 1
    elif k < (80*3):
        i = 2
    elif k <80*4:
        i = 3
    elif k < 80*5:
        i = 4
    else:
        i = 5
    
    #if k in [99, 199, 299, 399, 499, 599]:
        
        #T = simu.stream_states['pump_1_in']['T'][-1]
        #P = simu.stream_states['pump_1_in']['P'][-1]
        #x = simu.stream_states['pump_1_in']['x'][-1]
        #fluid = stream_nwe.copy_change_x_and_conditions(T+273.15, P, None, x, 'gas')
        #curve = evaluate_curves(fluid, x0, y0, w0)
        #curve.build_curves_resumed([15,120], 120, 5, 1000)
        #curve.plot_curve(simu.get_streams_PT_points_in_k(k))
        #curves.append(curve)

        
    simu.pumps_new_value(ns[i])
    #simu.valves_new_value(nv[i])
    
    simu.simulate_n_save_streams()
    
    
    
    k += 1
    
simu.set_unvisible()
#simu.save_data('data_new_cobeq_final', {'curves':curves})


simu.plot_stream_state("scrubber_in","P")
simu.plot_stream_state("scrubber_in","T")
simu.plot_stream_state("pump_2_out","P")
simu.plot_stream_state("pump_2_out","T")
simu.plot_stream_state("pump_1_out","T")
simu.plot_stream_state("1st_heat_out","T")









