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
from compression import multi_stage_compression

composition = [1]*30

list_names  = ["CH4",	  "C2H6",	  "C3H8",	  "iC4H10",  "nC4H10",  "iC5H12",  "nC5H12",  "nC6H14",  "nC7H16",  "nC8H18",	
               "nC9H20",  "nC10H22", "nC11H24", "nC12H26", "nC14H30", "N2", "H2O", "CO2", "C15+"]

nwe_1 = [0.145939192003332,0.0824656393169513,0.104456476468138,0.0146605581007913,
         0.051228654727197,0.0177426072469804,0.0244064972927947,0.0313202832153269,
         0.0353419812428145,0.0312122431043531,0.027928706900501,0.0252569872780137,
         0.0230417050869188,0.0211758373070616,0.018208515913102] + \
        [0.004164931, 0.2181591, 0.00241566]+ [0.101291129]

nwe_2 = [0.27169018312819,0.0354848393875713,0.0353647553287301,0.00876613629540678,
         0.0247373161212849,0.0151906334434104,0.0142900030021015,0.0261182827979586,
         0.0284719620124155,0.0257019369306632,0.0234532263091036,0.0215889912767973,
         0.0200167742604349,0.0186717842750411,0.0164884294416909] + \
        [0.001140799, 0.003002101, 0.302851996] + [0.089462624]

nwe_3 = [0.238095238095238,0.0262608309264608,0.0261719617862697,0.00648744723394801,
         0.0183070428793601,0.0112419462341702,0.0105754276827372,0.0193290379915574,
         0.0210708965703968,0.0190209179995643,0.0173567422429736,0.01597710063164,
         0.0148135692205084,0.0138181989380497,0.0122023902400072] + \
        [0.000844257, 0.002962305, 0.446300822] + [0.066207509]

dict_composition_1 = {list_names[i]: nwe_1[i] for i in range(len(nwe_1))}
mixture_nwe_1 = Mixture(list_of_species, dict_composition_1)

dict_composition_2 = {list_names[i]: nwe_2[i] for i in range(len(nwe_2))}
mixture_nwe_2 = Mixture(list_of_species, dict_composition_2)

dict_composition_3 = {list_names[i]: nwe_3[i] for i in range(len(nwe_3))}
mixture_nwe_3 = Mixture(list_of_species, dict_composition_3)

volumn_desviation = [0]*19

stream_nwe_3 = gc_eos_class(mixture_nwe_3, 85+273.15, 1.2e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'liquid')

solver_stream_nwe_3 = solver_eos(stream_nwe_3)

x0 = [0] + [4]*14 + [0] + [0] + [0] + [18]
y0 = [2]*5 + [0]*13 + [0]*2 + [7] + [0] + [40] + [0]
w0 = [1e-10]*15 + [0] + [1] + [0]*2

stream_nwe_31 = gc_eos_class(mixture_nwe_3, 80+273.15, 1.3e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'liquid')
stream_nwe_32 = gc_eos_class(mixture_nwe_3, 100+273.15, 1.3e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'liquid')
stream_nwe_33 = gc_eos_class(mixture_nwe_3, 80+273.15, 1.55e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'liquid')
stream_nwe_34 = gc_eos_class(mixture_nwe_3, 100+273.15, 1.55e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'liquid')
stream_nwe_35 = gc_eos_class(mixture_nwe_3, 80+273.15, 1.8e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'liquid')
stream_nwe_36 = gc_eos_class(mixture_nwe_3, 100+273.15, 1.8e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'liquid')


# phi1, stream_nwe_3_vap1, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_31, y0, x0)
# phi2, stream_nwe_3_vap2, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_32, y0, x0)
# phi3, stream_nwe_3_vap3, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_33, y0, x0)
# phi4, stream_nwe_3_vap4, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_34, y0, x0)
# phi3, stream_nwe_3_vap5, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_35, y0, x0)
# phi4, stream_nwe_3_vap6, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_36, y0, x0)

phi1, stream_nwe_3_vap1, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_31, [], [], 0.6)
phi2, stream_nwe_3_vap2, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_32, [], [], 0.6)
phi3, stream_nwe_3_vap3, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_33, [], [], 0.6)
phi4, stream_nwe_3_vap4, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_34, [], [], 0.6)
phi3, stream_nwe_3_vap5, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_35, [], [], 0.6)
phi4, stream_nwe_3_vap6, liq = solver_stream_nwe_3.evaluate_flash(stream_nwe_36, [], [], 0.6)


vap_refresh_1 = stream_nwe_3_vap1.copy_change_conditions(40+273.15, 4e4, None, 'gas')
vap_refresh_2 = stream_nwe_3_vap2.copy_change_conditions(40+273.15, 4e4, None, 'gas')
vap_refresh_3 = stream_nwe_3_vap3.copy_change_conditions(40+273.15, 4e4, None, 'gas')
vap_refresh_4 = stream_nwe_3_vap3.copy_change_conditions(40+273.15, 4e4, None, 'gas')
vap_refresh_5 = stream_nwe_3_vap5.copy_change_conditions(40+273.15, 4e4, None, 'gas')
vap_refresh_6 = stream_nwe_3_vap6.copy_change_conditions(40+273.15, 4e4, None, 'gas')

list_vap = [vap_refresh_1, vap_refresh_2, vap_refresh_3, vap_refresh_4, vap_refresh_5, vap_refresh_6]

list_solver = []

Pcrit = [1.05e4, 1.05e4, 1.45e4, 1.45e4]
delta_P = [500, 500, 400, 400]
Tcrit = [100, 100, 120, 120]

N = [1000/60,2000/60,3000/60,4000/60,5000/60,6000/60,7000/60,8000/60]
list_map = []
dm = 2
m_max = 100

for i,vap in enumerate(list_vap):
    
    compressor = CompressorClass()
    compressor.change_parameters()
    
    m_max = 80
    dm = 1 
    N = [1000/60,1500/60,2000/60,2500/60,3000/60,3500/60,4000/60,4500/60,5000/60]
    vis_case = viscosity(vap.mixture, [0]*19)        
    c = multi_stage_compression(10,vap, compressor, vis_case)
    c.characterization(N,1,m_max=m_max)
    list_map.append(c)

list_map_ic = []

N = [1000/60,2000/60,3000/60,4000/60,5000/60,6000/60,7000/60,8000/60]

vap_refresh_1 = stream_nwe_3_vap1.copy_change_conditions(50+273.15, 4e4, None, 'gas')
vap_refresh_2 = stream_nwe_3_vap2.copy_change_conditions(50+273.15, 4e4, None, 'gas')
vap_refresh_3 = stream_nwe_3_vap3.copy_change_conditions(50+273.15, 4e4, None, 'gas')
vap_refresh_4 = stream_nwe_3_vap3.copy_change_conditions(50+273.15, 4e4, None, 'gas')
vap_refresh_5 = stream_nwe_3_vap5.copy_change_conditions(50+273.15, 4e4, None, 'gas')
vap_refresh_6 = stream_nwe_3_vap6.copy_change_conditions(50+273.15, 4e4, None, 'gas')

list_vap = [vap_refresh_1, vap_refresh_2, vap_refresh_3, vap_refresh_4, vap_refresh_5, vap_refresh_6]


dm = 2
m_max = 80

for i,vap in enumerate(list_vap):
    compressor = CompressorClass()

    compressor.change_parameters2()
    m_max = 80
    dm = 1 
    N = [1000/60,1500/60,2000/60,2500/60,3000/60,3500/60,4000/60,4500/60,5000/60]
    vis_case = viscosity(vap.mixture, [0]*19)        
    c = multi_stage_compression(10,vap, compressor, vis_case, compressible=False)
    #c.characterization(N,1,m_max=80)
    list_map_ic.append(c)




