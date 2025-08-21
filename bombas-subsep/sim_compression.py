# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 21:38:00 2025

@author: fabio
"""

from eos_database_new_resumed import *
from gc_eos import gc_eos_class
from solver_thermo import solver_eos
import matplotlib.pyplot as plt
from caso_unisim import *
from compression import *
from compressor_class import *
from transfer_phenone_parameters import viscosity

#from hydrate_module import *

composition = [1]*30

list_names  = ["CH4",	  "C2H6",	  "C3H8",	  "iC4H10",  "nC4H10",  "iC5H12",  "nC5H12",  "nC6H14",  "nC7H16",  "nC8H18",	
               "nC9H20",  "nC10H22", "nC11H24", "nC12H26", "nC14H30", "N2", "H2O", "CO2", "C15+"]



nwe_3 = [0.238095238095238,0.0262608309264608,0.0261719617862697,0.00648744723394801,
         0.0183070428793601,0.0112419462341702,0.0105754276827372,0.0193290379915574,
         0.0210708965703968,0.0190209179995643,0.0173567422429736,0.01597710063164,
         0.0148135692205084,0.0138181989380497,0.0122023902400072] + \
        [0.000844257, 0.002962305, 0.446300822] + [0.066207509]

volumn_desviation = [0] * 19
dict_composition_3 = {list_names[i]: nwe_3[i] for i in range(len(nwe_3))}
mixture_nwe_3 = Mixture(list_of_species, dict_composition_3)

# 2. Criação dos fluidos de entrada (similar ao código original)
vap1 = gc_eos_class(mixture_nwe_3, 46+273.15, 16323.3041307910, None, 2, -1, 
                   Aij, Bij, Cij, volumn_desviation, 'liquid')
vap2 = gc_eos_class(mixture_nwe_3, 48+273.15, 17605, None, 2, -1,
                   Aij, Bij, Cij, volumn_desviation, 'liquid')

pump_1_in = vap1.copy_change_conditions(46+273.15, 16323.3041307910, None, "gas")
pump_2_in = vap2.copy_change_conditions(48+273.15, 17605, None, "gas")


list_vap = [pump_1_in,pump_2_in]

for i,vap in enumerate(list_vap):

    compressor = CompressorClass()
    compressor.change_parameters2()  # Configura parâmetros específicos do compressor
    vis_case = viscosity(pump_2_in.mixture, [0]*19)        
    
    c = multi_stage_compression(10 , vap, compressor, vis_case, compressible=True)
    
    m = 8005/3600  # Conversão de kg/h para kg/s
    N = [3000/60,3500/60,4000/60]#,2500/60,3000/60,3500/60,4000/60,4500/60,5000/60]
    dm = 1 
    m_max = 100  
    c.characterization(N, dm, m_max=m_max)



    
    