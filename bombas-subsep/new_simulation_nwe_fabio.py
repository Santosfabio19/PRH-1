# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:33:24 2024

@author: rodri
"""

from eos_database_new_resumed import *
from gc_eos import gc_eos_class
from solver_thermo import solver_eos

composition = [1]*30

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

stream_nwe_3 = gc_eos_class(mixture_nwe_3, 85+273.15, 1.2e4, None, 2, -1, Aij, volumn_desviation, 'liquid')

solver_stream_nwe_3 = solver_eos(stream_nwe_3)

x0 = [0] + [4]*14 + [0] + [0] + [0] + [18]
y0 = [2]*5 + [0]*13 + [0]*2 + [7] + [0] + [40] + [0]
w0 = [1e-10]*15 + [0] + [1] + [0]*2

stream_nwe_31 = gc_eos_class(mixture_nwe_3,  80+273.15, 1.3e4, None, 2, -1, Aij, volumn_desviation, 'liquid')

phi1, stream_nwe_3_vap1, liq_1 = solver_stream_nwe_3.evaluate_flash(stream_nwe_31, y0, x0)

def RGO_evaluation(liq,solver):
    fluid_std = liq.copy_change_conditions(15.57+273.15,101.325,None,'gas')
    phi_std, vap_std, liq_std = solver.evaluate_flash(fluid_std, y0, x0)
    RGO = vap_std.V*phi_std/(1-phi_std)/liq_std.V
    return RGO

liq_list = [liq_1, liq_2, liq_3, liq_4]

print([RGO_evaluation(liq, solver_stream_nwe_3) for liq in liq_list])


vap_refresh_1 = stream_nwe_3_vap1.copy_change_conditions(40+273.15, 1.3e4, None, 'gas')
vap_refresh_2 = stream_nwe_3_vap2.copy_change_conditions(40+273.15, 1.3e4, None, 'gas')
vap_refresh_3 = stream_nwe_3_vap3.copy_change_conditions(40+273.15,  1.8e4, None, 'gas')
vap_refresh_4 = stream_nwe_3_vap4.copy_change_conditions(40+273.15, 1.8e4, None, 'gas')

list_vap = [vap_refresh_1, vap_refresh_2, vap_refresh_3, vap_refresh_4]

list_solver = []

Pcrit = [1.15e4, 1.25e4, 1.60e4, 1.6e4]
delta_P = [500, 500, 400, 400]
Tcrit = [110, 110, 130, 130]

for i,vap in enumerate(list_vap):
    solver = solver_eos(vap)
    solver.set_estimated_conditions(x0, y0, w0)
    solver.isochoric_evaluation = True
    solver.critical_evaluation = True
    solver.water_evaluation = True
    solver.step_flash = 10
    solver.critical_init_value = [800, 0.08]
    solver.P_max_isochoric = 25000
    solver.P_max_water = 25000
    solver.build_phase_envelope([15,100],Tcrit[i],[1e3,1.2e4],Pcrit[i],5,delta_P[i])
    list_solver.append(solver)
    solver.plot_envelope()






