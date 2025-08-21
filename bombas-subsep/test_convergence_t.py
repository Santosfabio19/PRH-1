# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 01:30:07 2024

@author: Rodrigo Meira
"""


from eos_database import *
from casadi import *
from numpy import exp, log, array, roots
from scipy.optimize import fsolve
from species_builder import R, Species, Mixture

from matplotlib.pyplot import plot, figure
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter
from gc_eos import gc_eos_class

def cp_ig_fun_pure(par,T):
    return par[0] + par[1]*T + par[2]*T**2 + par[3]*T**3 + par[4]*T**4 + par[5]*T**5

def enthapy_ig_fun_pure(par,T):
    return par[0]*T + par[1]*T**2 + par[2]*T**3 + par[3]*T**4 + par[4]*T**5 + par[5]*T**6

def pvap_fun_pure(par,T):
    return exp(par[0] + par[1]/(T+par[2]) + par[3]*log(T) + par[4]*T**par[5])

def kappa_fun(par,w):
    return par[0] + par[1]*w  + par[2]*w**2 + par[3]*w**3

list_of_names = ['H2O', 'N2_2*', 'CO2_2*', 'C1_2*', 'C2_2*', 'C3_2*', 
                'iC4_2*', 'nC4_2*', 'iC5_2*', 'nC5_2*', 'C6_2*', 'C7_2*', 'C8_2*', 'C9_2*', 'C10_2*', 'C11_2*', 
                'C12_2*', 'C13_2*', 'C14_2*', 'C15_2*', 'C16_2*', 'C17_2*', 'C18_2*', 'C19_2*', 'C20+_2*']

list_of_names.pop(0)

list_of_species = [Species(row[0], row[1], row[4]+273.15, row[5], row[6], row[7]) 
                   for row in critical_table if row[0] in list_of_names]

list_of_index = [i for i in range(len(critical_table)) if critical_table[i][0] in list_of_names]


Tcp = MX.sym('Tcp')
wcp = MX.sym('acentricity')

[list_of_species[i].set_cp_ig_function(
    Function('CP_i_IG',[Tcp],[cp_ig_fun_pure(cp_ig_p[i],Tcp)/list_of_species[i].MM]),
    Function('H_i_IG',[Tcp],[enthapy_ig_fun_pure(cp_ig_p[i],Tcp)/list_of_species[i].MM]))
                   for i in range(len(list_of_species))]

[list_of_species[i].set_p_vap_function(
    Function('Pvap',[Tcp],[pvap_fun_pure(p_vap_pol[i],Tcp)]))
                   for i in range(len(list_of_species))]

[list_of_species[i].set_kappa(
    Function('kappa',[wcp],[kappa_fun(kp[i],wcp)]))
                   for i in range(len(list_of_species))]

[species.evaluate_kappa() for species in list_of_species]


composition = [3.44945577070149e-003, 2.48960745236753e-003, 0.369253301096065, 0.350151113126578, 4.69925780689622e-002,
               3.16949914724554e-002, 5.27916552903917e-003, 1.38978027217423e-002, 4.97921253693917e-003, 7.57880126909168e-003,
               9.96842256582064e-003, 1.39977840819956e-002, 1.33978786580008e-002, 9.96842134103575e-003, 8.86859527510853e-003,
               7.47881522052249e-003, 6.77892590685044e-003, 7.07887819148949e-003, 6.27900478629509e-003, 5.67909973211496e-003,
               4.97921057183309e-003, 3.88938329490607e-003, 4.28931982002196e-003, 3.68941490974787e-003, 5.78908166003149e-002]

composition.pop(0)



Aij = [[Aij[i][j] for i in range(len(Aij)) if i > 0] for j in range(len(Aij)) if j > 0]

dict_composition = {list_of_names[i]: composition[i] for i in range(len(composition))}

mixture = Mixture(list_of_species, dict_composition)

volumn_desviation.pop(0)

pinit = gc_eos_class(mixture, 97.29+273.15, 2.223e4, None, 2, -1, Aij, volumn_desviation, 'liquid')

Tt = []
xco2 = []
sumy = []

lin = []

y0 = [0] + [0.3]*2 +[0.3/20]*20 + [0.1]

for T in range(0,750,10):
    
    soma, y1, vap, liq = pinit.test_dew_P(T+273.15, 40000, y0)
    xco2.append(liq.mixture.x[1])
    Tt.append(T)
    sumy.append(soma)
    lin.append(1)
    
    
    
fig = figure()
p1, = plot(Tt, sumy)
p2, = plot(Tt, lin)
fig.tight_layout()
fig.show()

fig1 = figure()
p1, = plot(Tt, xco2)
fig1.tight_layout()
fig1.show()
