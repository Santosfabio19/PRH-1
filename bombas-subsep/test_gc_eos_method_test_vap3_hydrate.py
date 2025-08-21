
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


list_of_species = [Species(row[0], row[1], row[4]+273.15, row[5], row[6], row[7]) 
                   for row in critical_table if row[0] in list_of_names]

list_of_index = [i for i in range(len(critical_table)) if critical_table[i][0] in list_of_names]


Tcp = MX.sym('Tcp')
wcp = MX.sym('acentricity')

[list_of_species[i].set_cp_ig_function(
    Function('CP_i_IG',[Tcp],[cp_ig_fun_pure(cp_ig_p[i],Tcp)*list_of_species[i].MM]),
    Function('H_i_IG',[Tcp],[enthalpy_ig_fun_pure(cp_ig_p[i],Tcp)*list_of_species[i].MM]))
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

# -----------------------------------------FIX---------------------------------

composition = [3.44945577070149e-003, 2.48960745236753e-003, 0.369253301096065, 0.350151113126578, 4.69925780689622e-002,
               3.16949914724554e-002, 5.27916552903917e-003, 1.38978027217423e-002, 4.97921253693917e-003, 7.57880126909168e-003,
               9.96842256582064e-003, 1.39977840819956e-002, 1.33978786580008e-002, 9.96842134103575e-003, 8.86859527510853e-003,
               7.47881522052249e-003, 6.77892590685044e-003, 7.07887819148949e-003, 6.27900478629509e-003, 5.67909973211496e-003,
               4.97921057183309e-003, 3.88938329490607e-003, 4.28931982002196e-003, 3.68941490974787e-003, 5.78908166003149e-002]

composition = [.00552688,.00178958,.54945208,.28497686,.03133262,.02168416,
                .00387928,.00875324,.00268565,.00427715,.00557024,.0050729,
                .00775855,.00686334,.00576918,.00497343,.00457556,.00437662,
                .00358087,.00338193,.00268565,.00218831,.00238725,.00218831,.02427035]

composition = [.00773653,.00198397,.59036782,.31137302,.02411984,.01955394,
                .00317627,.0081392,.00258072,.0040696,.00148888,.00208443,
                .00248146,.00208443,.00178665,.00148888,.00138962,.00138962,
                .0011911,.0011911,.00089333,.00079407,.00079407,.00069481,.00714662]

dict_composition = {list_of_names[i]: composition[i] for i in range(len(composition))}

mixture = Mixture(list_of_species, dict_composition)

stream_S01 = gc_eos_class(mixture, 81.0324091443901+273.15, 25011.8375734781, None, 2, -1, Aij, volumn_desviation, 'liquid')

solver_stream_S01 = solver_eos(stream_S01)

y0 = [1e-10]*2 + [0.5]*2 +[1e-10]*20 + [1]
# y0 = array(mixture.x)**(2/3)
x0 = [0.1]+ [1e-10] + [0.0]*2 +[0.3/20]*20 + [1]

w0 = [1]*2 + [0.]*2 +[1e-10]*20 + [0]

phi, stream_vap, stream_liq = solver_stream_S01.evaluate_flash(stream_S01, x0, y0)

solver_stream_vap = solver_eos(stream_vap)
solver_stream_vap.set_estimated_conditions(x0, y0, w0)
solver_stream_vap.isochoric_evaluation = True
solver_stream_vap.critical_evaluation = True
solver_stream_vap.water_evaluation = True
solver_stream_vap.critical_init_value = [800, 0.078]

#solver_stream_vap.build_phase_envelope([-10,150],130,[1e3,2.2e4],2.2e4,10,500)
solver_stream_vap.build_phase_envelope([10,150],160,[1e3,2.3e4],2.2e4,10,500)

points = [[31.977475049757,63.1572826966245,32.0024931112913,57.4602861541931],
          [248.487203106415,468.904300995788,467.270640035532,731.328303693497]]


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


list_of_names = ['H2O', 'N2_2*', 'CO2_2*', 'C1_2*', 'C2_2*', 'C3_2*', 
                'iC4_2*', 'nC4_2*', 'iC5_2*', 'nC5_2*', 'C6_2*', 'C7_2*', 'C8_2*', 'C9_2*', 'C10_2*', 'C11_2*', 
                'C12_2*', 'C13_2*', 'C14_2*', 'C15_2*', 'C16_2*', 'C17_2*', 'C18_2*', 'C19_2*', 'C20+_2*']

# list_Aki_I[0] = [0]+[6.915e-2,2.614e1,8.287e2] + [0]*21
# list_Bki_I[0] = [0]+[1740,38.60,-881.1] + [0]*21

# list_Aki_I[1] = [0]+[1.530,1.113e-3,2.019e-3,8.547e-3]\
#                 + [0]*20
# list_Bki_I[1] = [0]+[2028,3856,3405,3583]\
#                 + [0]*20

# list_Aki_II[0] = [0]+[6.558e-2, 3.071e-3, 6.954e-3] \
#                 + [0]*21
# list_Bki_II[0] = [0]+[1444,2652,1865]\
#                 + [0]*21
                         
# list_Aki_II[1] = [0]+[1.530,4.824e-3,6.354e-3,9.765e-3,2.970e-5,2.372e-3,2.146e-6]\
#                 + [0]*18        
# list_Bki_II[1] = [0]+[229,3183,2785,3770,6081,4988,6305]\
#                 + [0]*18



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

#fun_int_cp = Function('int_CP',[Tcp],[int_fun_cp(-39.16,Tcp,273.15)])

fun_int_cp_I  = Function('int_CP',[Tcp],[int_fun_cp_2(-34.583,0.189,Tcp,273.15)])

fun_int_cp_II = Function('int_CP',[Tcp],[int_fun_cp_2(-36.8607,0.18,Tcp,273.15)])

water_ref_I = water_reference(273.15,101,-4297,1120,4.5959e-3,
                              fun_int_cp_I,fun_int_cp_I)

solver_hydrate_I = hydrate_solver(water_ref_I,hydrate_model_I,'H2O',solver_stream_vap,solver_stream_vap.lut_T_sat)

solver_hydrate_I.setup_duple_int_dt_fun(list_int_h_ig_dt)

water_ref_II = water_reference(273.15,101,-4611,931,4.99644e-3,
                              fun_int_cp_II,fun_int_cp_II)

hydrate_model_II = parrish_n_prusnitz(list_nu_II,list_Aki_II,list_Bki_II)

solver_hydrate_II = hydrate_solver(water_ref_II,hydrate_model_II,'H2O',solver_stream_vap,solver_stream_vap.lut_T_sat)

solver_hydrate_II.setup_duple_int_dt_fun(list_int_h_ig_dt)


points = [[34.9965888438368,64.4835359054855,59.9984891199637,86.34665498311],
          [223.905304503072,378.165085118064,376.998728116232,547.081373932586]]

solver_hydrate_I.evaluate_hydrate_curve([2.3e4,3e4,3.5e4,4e4,4.5e4,5e4,5.5e4,6e4],solver_stream_vap.fluid)

solver_stream_vap.plot_envelope(hydrate_points = [solver_hydrate_I.T_hydrade,solver_hydrate_I.P_hydrate],pt_points=points)

