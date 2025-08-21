# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 08:21:05 2024

@author: rodri
"""

from eos_database import *
from casadi import *
from numpy import exp, log, array, roots
from scipy.optimize import fsolve
from species_builder import R, Species, Mixture, pr_class

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
    return 0.378893 + 1.4897153*w-0.1713184*w**2+0.0196554*w**3


list_of_names = ['H2O', 'N2_2*', 'CO2_2*', 'C1_2*', 'C2_2*', 'C3_2*', 
                'iC4_2*', 'nC4_2*', 'iC5_2*', 'nC5_2*', 'C6_2*', 'C7_2*', 'C8_2*', 'C9_2*', 'C10_2*', 'C11_2*', 
                'C12_2*', 'C13_2*', 'C14_2*', 'C15_2*', 'C16_2*', 'C17_2*', 'C18_2*', 'C19_2*', 'C20+_2*']

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
    Function('kappa',[wcp],[kappa_fun(kappa_PR[i],wcp)]))
                   for i in range(len(list_of_species))]

[species.evaluate_kappa() for species in list_of_species]


composition = [3.44945577070149e-003, 2.48960745236753e-003, 0.369253301096065, 0.350151113126578, 4.69925780689622e-002,
               3.16949914724554e-002, 5.27916552903917e-003, 1.38978027217423e-002, 4.97921253693917e-003, 7.57880126909168e-003,
               9.96842256582064e-003, 1.39977840819956e-002, 1.33978786580008e-002, 9.96842134103575e-003, 8.86859527510853e-003,
               7.47881522052249e-003, 6.77892590685044e-003, 7.07887819148949e-003, 6.27900478629509e-003, 5.67909973211496e-003,
               4.97921057183309e-003, 3.88938329490607e-003, 4.28931982002196e-003, 3.68941490974787e-003, 5.78908166003149e-002]


composition = [0, 0.002498, 0.370531, 0.351363, 0.047155, 0.031805, 0.005297, 0.013946, 
               0.004996, 0.007605, 0.010003, 0.014046, 0.013444, 0.010003, 0.008899, 0.007505, 
               0.006802, 0.007103, 0.006301, 0.005699, 0.004996, 0.003903, 0.004304, 0.003702, 0.058091]

vap_c = [0, 0.00249822496890302, 0.370531432885458, 0.351363124684633, 0.0471552379767198, 0.0318047003797157, 0.00529743880991196, 0.0139459085163517, 0.00499644756181428, 0.00760503449923144, 0.0100029272208465, 0.0140462359516557, 0.0134442540176046, 0.0100029259918221, 0.00889929299267729, 0.00750470235938346, 0.00680239045184914, 0.00710338099003707, 0.00630073890647572, 0.00569875734351939, 0.00499644558990619, 0.0039028459895268, 0.00430416685321184, 0.00370218543465966, 0.0580911996240852]
liq_c = [0, 0.00217086901266447, 0.335769285172289, 0.314678368881983, 0.0435908160679091, 0.0301022470537834, 0.00511003837679921, 0.0135353530952131, 0.00493640145197559, 0.00752730858933425, 0.0103276219940831, 0.0149375793554574, 0.014492936869424, 0.0109830553930644, 0.00994154808215175, 0.00853714746019837, 0.00788937386269383, 0.00840524932473858, 0.00762596321167985, 0.00707416715087872, 0.00636856536599538, 0.00510470055118467, 0.00577221210983339, 0.00507623368904582, 0.12004295787762]

vap_c.pop(0)
liq_c.pop(0)

composition.pop(0)
list_of_species.pop(0)
list_of_names.pop(0)

Aij = [[Aij_PR[i][j] for i in range(len(Aij)) if i > 0] for j in range(len(Aij_PR)) if j > 0]
volumn_desviation.pop(0)

delta_v = [0. for i in volumn_desviation]

dict_composition = {list_of_names[i]: composition[i] for i in range(len(composition))}
mixture = Mixture(list_of_species, dict_composition)

dict_composition = {list_of_names[i]: vap_c[i] for i in range(len(composition))}
mixturev = Mixture(list_of_species, dict_composition)

dict_composition = {list_of_names[i]: liq_c[i] for i in range(len(composition))}
mixturel = Mixture(list_of_species, dict_composition)



#mixture.list_of_species[0].kappa = 0.710588

pinit = gc_eos_class(mixture, 600+273.15, 4e4, None, 2, -1, Aij, delta_v, 'gas')

pinit2 = gc_eos_class(mixturev, 600+273.15, 4e4, None, 2, -1, Aij, delta_v, 'gas')

pinit3 = gc_eos_class(mixturel, 600+273.15, 4e4, None, 2, -1, Aij, delta_v, 'liquid')


Tt = []
xco2 = []
sumy = []

lin = []

y0 = [0] + [0.13]*2 +[0.3/20]*20 + [0.55]

for T in range(400,650,1):
    
    soma, y1, vap, liq = pinit.test_dew_P(T+273.15, 429, y0)
    xco2.append(liq.mixture.x[1])
    Tt.append(T)
    sumy.append(soma)
    lin.append(1)
    
    
    
fig = figure()
p1, = plot(Tt, sumy)
p2, = plot(Tt, lin)
fig.tight_layout()
fig.show()









