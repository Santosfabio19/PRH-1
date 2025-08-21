# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 15:57:35 2024

@author: rodri
"""

from eos_database import *
from casadi import *
from solver_thermo import solver_eos
from numpy import exp, log, array, roots, isnan
from scipy.optimize import fsolve
from species_builder import R, Species, Mixture
from transfer_phenone_parameters import viscosity
from hydrate_module import *

from matplotlib.pyplot import plot, figure, scatter
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter
from gc_eos import gc_eos_class

def cp_ig_fun_pure(par,T):
    return par[1] + 2*par[2]*T + 3*par[3]*T**2 + 4*par[4]*T**3 + 5*par[5]*T**4

def enthalpy_ig_fun_pure(par,T):
    return par[0] + par[1]*T + par[2]*T**2 + par[3]*T**3 + par[4]*T**4 + par[5]*T**5

def int_h_ig_RT2_dT_fun_pure(par,T):
    return -par[0]/R/T + par[1]/R*log(T) + par[2]/R*T + par[3]*T**2/2/R + par[4]*T**3/3/R + par[5]*T**4/4/R


def entropy_ig_fun_pure(par,T):
    return  par[1]*log(T) + 2*par[2]*T + 3/2*par[3]*T**2 + 4/3*par[4]*T**3 + 5/4*par[5]*T**4 + 1

def pvap_fun_pure(par,T):
    return exp(par[0] + par[1]/(T+par[2]) + par[3]*log(T) + par[4]*T**par[5])

def kappa_fun(par,w):
    return par[0] + par[1]*w  + par[2]*w**2 + par[3]*w**3

list_of_names = ['H2O', 'N2_2*', 'CO2_2*', 'C1_2*', 'C2_2*', 'C3_2*', 
                'iC4_2*', 'nC4_2*', 'iC5_2*', 'nC5_2*', 'C6_2*', 'C7_2*', 'C8_2*', 'C9_2*', 'C10_2*', 'C11_2*', 
                'C12_2*', 'C13_2*', 'C14_2*', 'C15_2*', 'C16_2*', 'C17_2*', 'C18_2*', 'C19_2*', 'C20+_2*']

dipoles = [1.8] + [0]*4 + [0.4] + [0, 0.1]*2 + [0]*15


data = {
    'yCO2': [10.0, 9.0, 8.0, 8.0, 8.0, 8.0, 9.0, 14.0, 13.0, 13.0, 13.0, 13.0, 12.0, 13.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 44.0, 42.0, 40.0, 39.0, 40.0, 39.0, 39.0, 44.0, 50.0, 47.0, 44.0, 40.0, 44.0, 39.0, 73.0, 70.0, 68.0, 68.0, 67.0, 79.0, 76.0, 76.0, 75.0, 74.0, 85.0],
    'T': [273.7, 275.8, 277.8, 280.2, 283.2, 285.1, 287.2, 274.6, 276.9, 278.1, 280.9, 284.0, 286.1, 287.4, 273.8, 276.8, 278.4, 283.4, 285.1, 287.4, 273.7, 276.9, 278.7, 280.7, 283.1, 285.1, 287.4, 273.7, 275.6, 278.5, 280.9, 281.4, 285.1, 287.4, 274.6, 276.4, 278.2, 280.2, 282.0, 273.7, 275.9, 277.8, 280.2, 281.6, 282.7],
    'Pexp': [2.52, 3.10, 3.83, 4.91, 6.80, 8.40, 10.78, 2.59, 3.24, 4.18, 6.38, 7.17, 9.24, 10.95, 2.12, 3.96, 4.83, 6.94, 8.82, 9.78, 1.81, 2.63, 4.03, 5.34, 6.98, 6.94, 9.78, 1.81, 1.99, 2.98, 4.17, 4.41, 6.84, 9.78, 1.66, 2.08, 2.58, 3.28, 4.12, 1.45, 1.85, 2.37, 2.97, 3.79, 4.37],
    'Pcalc': [2.51, 3.07, 3.78, 4.85, 6.80, 8.53, 10.99, 2.61, 3.29, 4.13, 6.43, 7.19, 9.37, 10.98, 2.17, 3.95, 4.78, 7.02, 8.58, 9.96, 1.82, 2.58, 4.03, 5.34, 7.01, 7.02, 9.96, 1.82, 2.12, 2.98, 4.13, 4.57, 6.85, 9.96, 1.56, 2.05, 2.56, 3.21, 4.14, 1.44, 1.86, 2.35, 2.95, 3.84, 4.32],
    'err': [-0.4, -0.9, -1.4, -1.1, 0.0, 1.5, 2.1, 0.8, 1.5, -1.3, 1.0, 0.3, 1.4, 0.3, 2.6, -0.3, -1.0, 1.1, -2.5, 1.8, 0.7, -1.9, 0.0, 0.0, 0.4, 1.1, 1.8, 0.7, 2.0, 0.0, -0.3, -2.3, 0.2, 1.8, -0.9, -1.5, -0.8, -0.5, 0.5, -0.5, 0.2, -0.7, -0.8, 1.2, -1.2]
}

list_of_species = [Species(row[0], row[1], row[4]+273.15, row[5], row[6], row[7]) 
                   for row in critical_table if row[0] in list_of_names]

list_of_index = [i for i in range(len(critical_table)) if critical_table[i][0] in list_of_names]


Tcp = MX.sym('Tcp')
wcp = MX.sym('acentricity')

[list_of_species[i].set_cp_ig_function(
    Function('CP_i_IG',[Tcp],[cp_ig_fun_pure(cp_ig_p[i],Tcp)*list_of_species[i].MM]),
    Function('H_i_IG',[Tcp],[enthalpy_ig_fun_pure(cp_ig_p[i],Tcp)*list_of_species[i].MM]))
                   for i in range(len(list_of_species))]

list_int_h_ig_dt = [Function('H_i_IG_RT2_dT',[Tcp],[int_h_ig_RT2_dT_fun_pure(cp_ig_p[i],Tcp)*list_of_species[i].MM]) \
                    for i in range(len(list_of_species))]

[list_of_species[i].set_entropy_ig_function(
    Function('S_i_IG',[Tcp],[entropy_ig_fun_pure(cp_ig_p[i],Tcp)*list_of_species[i].MM]))
                   for i in range(len(list_of_species))]

[list_of_species[i].set_p_vap_function(
    Function('Pvap',[Tcp],[pvap_fun_pure(p_vap_pol[i],Tcp)]))
                   for i in range(len(list_of_species))]

[list_of_species[i].set_kappa(
    Function('kappa',[wcp],[kappa_fun(kp[i],wcp)]))
                   for i in range(len(list_of_species))]

[species.evaluate_kappa() for species in list_of_species]

#%%  

list_nu_I  = [2/46, 6/46]
list_nu_II = [16/136, 8/136]
list_nu_H  = [3/34, 2/34, 1/34]

list_Aki_I  = [[],[]]
list_Aki_II = [[],[]]
list_Aki_H  = [[],[],[]]

list_Bki_I  = [[],[]]
list_Bki_II = [[],[]]
list_Bki_H  = [[],[],[]]


list_Aki_I[0] = [0]+[6.915e-2,2.614e1,8.287e2] + [0]*21
list_Bki_I[0] = [0]+[1740,38.60,-881.1] + [0]*21

list_Aki_I[1] = [0]+[1.530,1.113e-3,2.019e-3,8.547e-3]\
                + [0]*20
list_Bki_I[1] = [0]+[2028,3856,3405,3583]\
                + [0]*20

list_Aki_II[0] = [0]+[6.558e-2, 3.071e-3, 6.954e-3] \
                + [0]*21
list_Bki_II[0] = [0]+[1444,2652,1865]\
                + [0]*21
                         
list_Aki_II[1] = [0]+[1.530,4.824e-3,6.354e-3,9.765e-3,2.970e-5,2.372e-3,2.146e-6]\
                + [0]*18        
list_Bki_II[1] = [0]+[229,3183,2785,3770,6081,4988,6305]\
                + [0]*18



list_Aki_II[0] = [0]+[6.558e-2,3.071e-3,6.954e-3] + [0]*21


def int_fun_cp(delta_cp,T,T0):
    return delta_cp/R*log(T/T0) +delta_cp*T0/R*(1/T-1/T0)

def int_fun_cp_2(delta_cp,b_cp,T,T0):
    
    int_b_cp = b_cp*( 1/2*(T-T0) - T0*log(T/T0) - T0**2/2*(1/T-1/T0) )/R
    
    int_delta_cp = delta_cp/R*log(T/T0) +delta_cp*T0/R*(1/T-1/T0)
    return  int_delta_cp + int_b_cp

hydrate_model_I = parrish_n_prusnitz(list_nu_I,list_Aki_I,list_Bki_I)

fun_int_cp_I  = Function('int_CP',[Tcp],[int_fun_cp_2(-34.583,0.189,Tcp,273.15)])

fun_int_cp_II = Function('int_CP',[Tcp],[int_fun_cp_2(-36.8607,0.18,Tcp,273.15)])

water_ref_I = water_reference(273.15,101,-4297,1120,4.5959e-3,
                              fun_int_cp_I,fun_int_cp_I)
water_ref_II = water_reference(273.15,101,-4611,931,4.99644e-3,
                              fun_int_cp_II,fun_int_cp_II)
hydrate_model_II = parrish_n_prusnitz(list_nu_II,list_Aki_II,list_Bki_II)


def case_hidrate(T,P,hydrate_model,water_ref,fluid):
    vap = fluid.copy_change_conditions(T,P,None,'gas')
    delta_mu_H_beta = hydrate_model.evaluate_delta_mu_H_beta(vap)
    Pvap = list_of_species[0].evaluate_pvap(300)
    delta_mu_beta_0 = water_ref.evaluate_delta_mu_beta_0(P,T)
    return delta_mu_H_beta + delta_mu_beta_0




for i in range(len(data['yCO2'])):
    composition = [0]*2+[data['yCO2'][i]/100,1-data['yCO2'][i]/100] + [0]*21
    dict_composition = {list_of_names[i]: composition[i] for i in range(len(composition))}
    mixture = Mixture(list_of_species, dict_composition)
    vap = gc_eos_class(mixture, data['T'][i], data['Pexp'][i]*1000, None, 2, -1, Aij, volumn_desviation, 'liquid')



















