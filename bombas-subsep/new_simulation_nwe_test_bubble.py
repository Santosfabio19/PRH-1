# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:33:24 2024

@author: rodri
"""

from eos_database_new import *
from gc_eos import gc_eos_class
from solver_thermo import solver_eos

composition = [1]*30

list_names  = ["CH4",	  "C2H6",	  "C3H8",	  "iC4H10",  "nC4H10",  "iC5H12",  "nC5H12",  "nC6H14",  "nC7H16",  "nC8H18",	
               "nC9H20",  "nC10H22", "nC11H24", "nC12H26", "nC13H28", "nC14H30", "nC15H32", "nC16H34", "nC17H36",	
               "nC18H38", "nC19H40", "nC20H42", "nC21H44", "nC22H46", "nC23H48", "nC24H50", "nC25H52", "N2", "H2O", "CO2"]

nwe_1 = [0.145939192003332,0.0824656393169513,0.104456476468138,0.0146605581007913,
         0.051228654727197,0.0177426072469804,0.0244064972927947,0.0313202832153269,
         0.0353419812428145,0.0312122431043531,0.027928706900501,0.0252569872780137,
         0.0230417050869188,0.0211758373070616,0.0195832943042615,0.018208515913102] +[0.101291129]*11 + \
        [0.004164931, 0.2181591, 0.00241566]

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
y0 = [4]*2 + [0]*23 + [0]*2 + [0.2] + [.120] + [4]
w0 = [1]*2 + [1e-10]*29

# solver_stream_nwe_1.set_estimated_conditions(x0, y0, w0)
# solver_stream_nwe_1.isochoric_evaluation = False
# solver_stream_nwe_1.critical_evaluation = True
# solver_stream_nwe_1.water_evaluation = False
# solver_stream_nwe_1.step_flash = 10

# solver_stream_nwe_1.build_phase_envelope([-10,100],350,[1e0,1.2e4],1.1e4,10,500)

# solver_stream_nwe_1.plot_envelope()

# solver_stream_nwe_2.set_estimated_conditions(x0, y0, w0)
# solver_stream_nwe_2.isochoric_evaluation = False
# solver_stream_nwe_2.critical_evaluation = True
# solver_stream_nwe_2.water_evaluation = False
# solver_stream_nwe_2.step_flash = 10

# solver_stream_nwe_2.build_phase_envelope([-10,280],320,[1e3,1.4e4],1.4e4,10,500)

# solver_stream_nwe_2.plot_envelope()


solver_stream_nwe_3.set_estimated_conditions(x0, y0, w0)
solver_stream_nwe_3.isochoric_evaluation = False
solver_stream_nwe_3.critical_evaluation = True
solver_stream_nwe_3.water_evaluation = False
solver_stream_nwe_3.step_flash = 10
solver_stream_nwe_3.critical_init_value = [800, 0.28]

for P in range(10000,23000,100):
    res = solver_stream_nwe_3.test_bubble(270+273.15, P, y0)
    print(res,P)






