# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:53:40 2024

@author: rodri
"""

from setup_gc_eos import *

test_crit = [0.0109, 0.0884, 0.8286, 0.0401, 0.0174, 0.003, 0.0055, 0.0019, 0.0012, 0.0014, 0.0006]
# . CO2 N2 CH4 C2H6 C3H8 iC4H10 nC4H10 iC5H12 nC5H12 nC6H14 nC7H16

list_of_names = ['H2O', 'N2_2*', 'CO2_2*', 'C1_2*', 'C2_2*', 'C3_2*', 
                'iC4_2*', 'nC4_2*', 'iC5_2*', 'nC5_2*', 'C6_2*', 'C7_2*', 'C8_2*', 'C9_2*', 'C10_2*', 'C11_2*', 
                'C12_2*', 'C13_2*', 'C14_2*', 'C15_2*', 'C16_2*', 'C17_2*', 'C18_2*', 'C19_2*', 'C20+_2*']

index_numbers = [2, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11]

list_dimi_names = [list_of_names[i] for i in index_numbers]

list_of_species = [Species(row[0], row[1], row[4]+273.15, row[5], row[6], row[7]) 
                   for row in critical_table if row[0] in list_dimi_names]

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

dict_composition = {list_dimi_names[i]: test_crit[i] for i in range(len(test_crit))}

mixture = Mixture(list_of_species, dict_composition)

Aij = [[Aij[i][j] for i in index_numbers] for j in index_numbers]

# Aij = np.array([
#     [0, -0.0199970, 0.1000000, 0.1298000, 0.1350000, 0.1298000, 0.1298000, 0.1250000, 0.1250000, 0.1250000, 0.1199000],
#     [-0.0199970, 0, 0.0359990, 0.0500000, 0.0799980, 0.0949990, 0.0900000, 0.0949990, 0.0990000, 0.1490000, 0.1439000],
#     [0.1000000, 0.0359990, 0, 0.0022413, 0.0068288, 0.0131134, 0.0123047, 0.0176275, 0.0179254, 0.0234741, 0.0288643],
#     [0.1298000, 0.0500000, 0.0022413, 0, 0.0012579, 0.0045736, 0.0040964, 0.0074333, 0.0076095, 0.0114138, 0.0153243],
#     [0.1350000, 0.0799980, 0.0068288, 0.0012579, 0, 0.0010406, 0.0008189, 0.0025834, 0.0027005, 0.0051420, 0.0078874],
#     [0.1298000, 0.0949990, 0.0131134, 0.0045736, 0.0010406, 0, 0.0000133, 0.0003462, 0.0003900, 0.0015653, 0.0032212],
#     [0.1298000, 0.0900000, 0.0123047, 0.0040964, 0.0008189, 0.0000133, 0, 0.0004951, 0.0005472, 0.0018663, 0.0036464],
#     [0.1250000, 0.0949990, 0.0176275, 0.0074333, 0.0025834, 0.0003462, 0.0004951, 0, 0.0000013, 0.0004400, 0.0014592],
#     [0.1250000, 0.0990000, 0.0179254, 0.0076095, 0.0027005, 0.0003900, 0.0005472, 0.0000013, 0, 0.0003934, 0.0003733],
#     [0.1250000, 0.1490000, 0.0234741, 0.0114138, 0.0051420, 0.0015653, 0.0018663, 0.0004400, 0.0003934, 0, 0.0002972],
#     [0.1199000, 0.1439000, 0.0288643, 0.0153243, 0.0078874, 0.0032212, 0.0036464, 0.0014592, 0.0003733, 0.0002972, 0]
# ])

volumn_desviation = [volumn_desviation[i] for i in index_numbers]

pinit = gc_eos_class(mixture, 204, 6584, None, 2, -1, Aij, volumn_desviation, 'liquid')

x0 = [0]*2 + [0.5] + [0]*6  +[1]*2

y0 = [0.2]*2 + [1] + [0]*8

T0 = sum(array(pinit.mixture.list_Tc)*array(pinit.mixture.x))*3
V0 = sum(array([sp.Vc for sp in pinit.mixture.list_of_species])*array(pinit.mixture.x))*1.5

crit_solver = lambda y: array([float(i) for i in pinit.critical_case_fehbet(y[0], y[1], y[2:],False)])
y0 = array(mixture.x)**(2/3)
y_ss = fsolve(crit_solver, [T0,V0]+[i for i in y0])

#case = pinit.copy_change_conditions(y_ss[0], None, y_ss[1], 'gas')

pinit.build_phase_envelope([-80,-20],0,[1e3,5.5e3],10e3,1,100,x0,x0)

pinit.plot_envelope()


