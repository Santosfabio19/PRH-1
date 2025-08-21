# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:33:24 2024

@author: rodri
"""

from eos_database_new_resumed import *
from gc_eos import gc_eos_class
from solver_thermo import solver_eos
from transfer_phenone_parameters import viscosity
from compressor_class import CompressorClass
from compression import compression
from casadi import *

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

stream_nwe_3 = gc_eos_class(mixture_nwe_3, 97.29+273.15, 2.223e4, None, 2, -1, Aij, volumn_desviation, 'liquid')

solver_stream_nwe_3 = solver_eos(stream_nwe_3)


x0 = [0] + [4]*14 + [0] + [0] + [0] + [18]
y0 = [40]*5 + [0]*13 + [0]*2 + [7] + [0] + [40] + [0]
w0 = [1e-10]*15 + [0] + [1] + [0]*2


solver_stream_nwe_3.set_estimated_conditions(x0, y0, w0)
solver_stream_nwe_3.isochoric_evaluation = False
solver_stream_nwe_3.critical_evaluation = True
solver_stream_nwe_3.water_evaluation = False
solver_stream_nwe_3.step_flash = 10
solver_stream_nwe_3.critical_init_value = [800, 0.28]

T_sep = 57.5
T_heat = 40.

PTfluid_nwe_1 = []
PTfluid_nwe_2 = []
PTfluid_nwe_3 = []

def evaluate_condition(list_of_conditions,fluid,solver):
    x0 = [0] + [34]*14 + [0] + [0] + [0] + [67]    
    phi_init, vap, liq = solver.evaluate_flash(fluid, solver.y0, solver.x0)
    Pbub,_,_ = solver.evaluate_bubble_T(vap,40+273.15,10000,x0)
    vap2 = vap.copy_change_conditions(40+273.15,fluid.P,None,'gas')
    print(Pbub,fluid.P)
    if Pbub > fluid.P:    
        phi, vap_final, liq3 = solver.evaluate_flash(vap2, solver.y0, solver.x0, 1)
    else:
        vap_final = vap2
        phi = 1
        
    print(phi)
    rho = vap_final.rho*vap_final.mixture.MM_m
    vis_case = viscosity(vap_final.mixture, [0]*19)
    mu = vis_case.evaluate_viscosity(vap_final.T, vap_final.P)
    list_of_conditions.append([fluid.T-273.15,fluid.P/100,phi_init,phi,rho,mu,vap_final])
    return list_of_conditions    

Pbub = [11528.761010734019, 12739.697047111044, 13898.309721201056, 15003.525574430816, 16053.593837765966,
        17045.961751143637, 17977.160396642863, 18842.73771954823, 19637.302370453184, 20266.13645894297,
        20770, 21243.86897044951, 21690, 22105.747054394185, 22490, 22850, 23190, 23490, 23760, 24000,
        24210, 24400, 24560, 24680, 24780, 24850, 24890, 24900, 24890, 24810, 24740, 24640, 24500, 24330,
        24120, 23880, 23600, 23280, 22910]

Tbub = [50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145,
        150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240]

lut_T_ref = interpolant('LUT','bspline',[Tbub],Pbub)

for T in range(80,95,5):
    for P in range(120,160,10):
        Pref = lut_T_ref(T).__float__()/100
        print(P<Pref)
        if P < Pref:
            fluid3 = stream_nwe_3.copy_change_conditions(T+273.15, P*100, None, 'liquid')
            evaluate_condition(PTfluid_nwe_3, fluid3, solver_stream_nwe_3)
        
list_solvers = []

for ref in PTfluid_nwe_3:
    solver_stream_nwe_ref = solver_eos(ref[-1])
    solver_stream_nwe_ref.set_estimated_conditions(x0, y0, w0)
    solver_stream_nwe_ref.isochoric_evaluation = True
    solver_stream_nwe_ref.critical_evaluation = True
    solver_stream_nwe_ref.water_evaluation = True
    solver_stream_nwe_ref.step_flash = 10
    solver_stream_nwe_ref.critical_init_value = [800, 0.15]
    solver_stream_nwe_ref.build_phase_envelope([10,100],100,[1e3,1.2e4],1.2e4,5,500)
    
    















