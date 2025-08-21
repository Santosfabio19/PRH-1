# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:33:24 2024

@author: rodri
"""

from eos_database_new import *
from gc_eos import gc_eos_class
from solver_thermo import solver_eos
from transfer_phenone_parameters import viscosity
from compressor_class import CompressorClass
from compression import compression
from casadi import *

composition = [1]*30

list_names  = ["CH4",	  "C2H6",	  "C3H8",	  "iC4H10",  "nC4H10",  "iC5H12",  "nC5H12",  "nC6H14",  "nC7H16",  "nC8H18",	
               "nC9H20",  "nC10H22", "nC11H24", "nC12H26", "nC13H28", "nC14H30", "nC15H32", "nC16H34", "nC17H36",	
               "nC18H38", "nC19H40", "nC20H42", "nC21H44", "nC22H46", "nC23H48", "nC24H50", "nC25H52", "N2", "H2O", "CO2"]

nwe_1 = [0.306200076074553,0.0449600608596425,0.0448079117535184,0.0111068847470521,
         0.0313427158615443,0.0192468619246862,0.0181057436287562,0.0330924305819703,
         0.0360745931773892,0.0325649113412092,0.0297157462288414,0.027353720023968,
         0.0253616869766103,0.0236575555040555,0.0221820345720906,0.020891197592344] +[0.113351084/11]*11 + \
        [0.001445417, 0.003042982, 0.155496386]

nwe_2 = [0.27169018312819,0.0354848393875713,0.0353647553287301,0.00876613629540678,
         0.0247373161212849,0.0151906334434104,0.0142900030021015,0.0261182827979586,
         0.0284719620124155,0.0257019369306632,0.0234532263091036,0.0215889912767973,
         0.0200167742604349,0.0186717842750411,0.0175072257250154,0.0164884294416909] +[0.089462624/11]*11 + \
        [0.001140799, 0.003002101, 0.302851996]

nwe_3 = [0.238095238095238,0.0262608309264608,0.0261719617862697,0.00648744723394801,
         0.0183070428793601,0.0112419462341702,0.0105754276827372,0.0193290379915574,
         0.0210708965703968,0.0190209179995643,0.0173567422429736,0.01597710063164,
         0.0148135692205084,0.0138181989380497,0.0129563583403746,0.0122023902400072] +[0.066207509/11]*11 + \
        [0.000844257, 0.002962305, 0.446300822]

dict_composition_1 = {list_names[i]: nwe_1[i] for i in range(len(nwe_1))}
mixture_nwe_1 = Mixture(list_of_species, dict_composition_1)

dict_composition_2 = {list_names[i]: nwe_2[i] for i in range(len(nwe_2))}
mixture_nwe_2 = Mixture(list_of_species, dict_composition_2)

dict_composition_3 = {list_names[i]: nwe_3[i] for i in range(len(nwe_3))}
mixture_nwe_3 = Mixture(list_of_species, dict_composition_3)

volumn_desviation = [0]*30

stream_nwe_1 = gc_eos_class(mixture_nwe_1, 97.29+273.15, 2.223e4, None, 2, -1, Aij, volumn_desviation, 'liquid')

stream_nwe_2 = gc_eos_class(mixture_nwe_2, 97.29+273.15, 2.223e4, None, 2, -1, Aij, volumn_desviation, 'liquid')

stream_nwe_3 = gc_eos_class(mixture_nwe_3, 97.29+273.15, 2.223e4, None, 2, -1, Aij, volumn_desviation, 'liquid')

solver_stream_nwe_1 = solver_eos(stream_nwe_1)
solver_stream_nwe_2 = solver_eos(stream_nwe_2)
solver_stream_nwe_3 = solver_eos(stream_nwe_3)


x0 = [0] + [4]*24 + [8]*2 + [0] + [1] + [0]
y0 = [1]*2 + [0]*23 + [0]*3 + [1] + [1]
w0 = [1]*2 + [1e-10]*29

solver_stream_nwe_1.set_estimated_conditions(x0, y0, w0)
solver_stream_nwe_1.isochoric_evaluation = False
solver_stream_nwe_1.critical_evaluation = True
solver_stream_nwe_1.water_evaluation = False
solver_stream_nwe_1.step_flash = 10


solver_stream_nwe_2.set_estimated_conditions(x0, y0, w0)
solver_stream_nwe_2.isochoric_evaluation = False
solver_stream_nwe_2.critical_evaluation = True
solver_stream_nwe_2.water_evaluation = False
solver_stream_nwe_2.step_flash = 10

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
    phi_init, vap, liq = solver.evaluate_flash(fluid, solver.y0, solver.x0)
    vap2 = vap.copy_change_conditions(40+273.15,fluid.P,None,'gas')
    print(phi_init)
    phi, vap3, liq3 = solver.evaluate_flash(vap2, solver.y0, solver.x0)
    rho = vap3.rho*vap3.mixture.MM_m
    vis_case = viscosity(vap3.mixture, [0]*30)
    mu = vis_case.evaluate_viscosity(vap3.T, vap3.P)
    if phi < 1 and phi > 0.9:
        list_of_conditions.append([fluid.T-273.15,fluid.P/100,phi_init,phi,vap,liq,liq3,rho,mu,vap2,vap3])
    return list_of_conditions    

Pbub = [13470,  14607.53923126884,  15688.318896336355,  16701.352481511713,
        17640,  18503.35586145648,  19200,  19710,  20150,  20520,  20840, 
        21090,  21300,  21450,  21560,  21640,  21690,  21740,  21830,  22200,
        21480,  20820, 20030]

Tbub = [50, 60,
        70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170,
        180, 190, 200, 210, 220, 230, 240, 250, 260, 270]

lut_T_ref = interpolant('LUT','bspline',[Tbub],Pbub)

for T in range(75,95,5):
    for P in range(120,180,10):
        Pref = lut_T_ref(T)
        if P < Pref:
            fluid3 = stream_nwe_3.copy_change_conditions(T+273.15, P*100, None, 'liquid')
            evaluate_condition(PTfluid_nwe_3, fluid3, solver_stream_nwe_3)
        
compressor = CompressorClass()

vis_case = viscosity(PTfluid_nwe_3[0][-1].mixture, [0]*30)
    
c = compression(PTfluid_nwe_3[0][-1], compressor, vis_case)

N = [200,300,400,500,600,700]

# c.characterization(N,0.05)

















