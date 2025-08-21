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

nwe_1 = [0.145939192003332,0.0824656393169513,0.104456476468138,0.0146605581007913,
         0.051228654727197,0.0177426072469804,0.0244064972927947,0.0313202832153269,
         0.0353419812428145,0.0312122431043531,0.027928706900501,0.0252569872780137,
         0.0230417050869188,0.0211758373070616,0.018208515913102] + \
        [0.004164931, 0.2181591, 0.00241566]+ [0.101291129]

dict_composition_1 = {list_names[i]: nwe_1[i] for i in range(len(nwe_1))}
mixture_nwe_1 = Mixture(list_of_species, dict_composition_1)

volumn_desviation = [0]*19

stream_nwe_1 = gc_eos_class(mixture_nwe_1, 97.29+273.15, 1.4e4, None, 2, -1, Aij, volumn_desviation, 'liquid')


