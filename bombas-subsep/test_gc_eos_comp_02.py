from eos_database import *
from casadi import *
from numpy import exp, log, array, roots, isnan
from scipy.optimize import fsolve
from species_builder import R, Species, Mixture

from matplotlib.pyplot import plot, figure, scatter
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter
from gc_eos import gc_eos_class

def config_plot(axes):
    """
    Configures the axes of Figures
    """
    formatter = ScalarFormatter(useOffset=False, useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    axes.yaxis.set_major_formatter(formatter)

    axes.xaxis.set_minor_locator(AutoMinorLocator())
    axes.yaxis.set_minor_locator(AutoMinorLocator())
    axes.tick_params(which='both', direction='out', bottom=True, left=True)
    axes.tick_params(which='major', width=2)
    axes.tick_params(which='minor', width=1)
    #axes.xaxis.set_tick_params(which='both', right='off', top='off', direction='out', width=1)
    #axes.yaxis.set_tick_params(which='both', right='off', top='off', direction='out', width=1)

    #axes.spines['top'].set_visible(False)

        
def cp_ig_fun_pure(par,T):
    return par[0] + par[1]*T + par[2]*T**2 + par[3]*T**3 + par[4]*T**4 + par[5]*T**5

def enthapy_ig_fun_pure(par,T):
    return par[0]*T + par[1]*T**2 + par[2]*T**3 + par[3]*T**4 + par[4]*T**5 + par[5]*T**6

def pvap_fun_pure(par,T):
    return exp(par[0] + par[1]/(T+par[2]) + par[3]*log(T) + par[4]*T**par[5])

def kappa_fun(par,w):
    return par[0] + par[1]*w  + par[2]*w**2 + par[3]*w**3


def evaluate_bubble_P(T,Pr,y,liq_base):
    
    vapor_composition = {liq_base.mixture.list_of_names[i]: y[i] for i in range(len(y))}
    
    vapor_mixture = Mixture(liq_base.mixture.list_of_species, vapor_composition)
    
    P = Pr*liq_base.Pc
    
    print(P)
    
    liquid = liq_base.copy_change_conditions(T, P, None, 'liquid')
    
    vapor = liq_base.copy_change_x_and_conditions(T, P, None, y, 'gas')
    
    gas_fugacity = vapor.evaluate_fugacity_coef()
    
    liquid_fugacity = liquid.evaluate_fugacity_coef()
    
    equilibrium_equations = 1 - array(y)*gas_fugacity/array(liquid.mixture.x)/liquid_fugacity
    
    equilibrium_equations_eff = [float(equilibrium_equations[i]) 
                                 for i in range(equilibrium_equations.shape[0])]
    
    return equilibrium_equations_eff + [1-sum(y)]

def evaluate_dew_T(Tr,P,x,gas_base,Aij):
    
    liquid_composition = {gas_base.mixture.list_of_names[i]: x[i] for i in range(len(x))}
    
    liquid_mixture = Mixture(gas_base.mixture.list_of_species, liquid_composition)
    
    T = Tr*gas_base.Tc
    
    liquid = gc_eos_class(liquid_mixture, T, P, None, 2, -1,Aij, volumn_desviation, 'liquid')
    
    vapor = gc_eos_class(gas_base.mixture, T, P, None, 2, -1,Aij, volumn_desviation, 'gas')
    
    gas_fugacity = vapor.evaluate_fugacity_coef()
    
    liquid_fugacity = liquid.evaluate_fugacity_coef()
    
    equilibrium_equations = 1 - array(x)*liquid_fugacity/array(vapor.mixture.x)/gas_fugacity
    
    equilibrium_equations_eff = [float(equilibrium_equations[i]) 
                                 for i in range(equilibrium_equations.shape[0])]
    
    return equilibrium_equations_eff + [1-sum(x)]

def evaluate_dew(T,Pr,x,gas_base,Aij):
    
    liquid_composition = {gas_base.mixture.list_of_names[i]: x[i] for i in range(len(x))}
    
    liquid_mixture = Mixture(gas_base.mixture.list_of_species, liquid_composition)
    
    P = Pr*gas_base.Pc
    
    liquid = gc_eos_class(liquid_mixture, T, P, None, 2, -1,Aij, volumn_desviation, 'liquid')
    
    vapor = gc_eos_class(gas_base.mixture, T, P, None, 2, -1,Aij, volumn_desviation, 'gas')
    
    gas_fugacity = vapor.evaluate_fugacity_coef()
    
    liquid_fugacity = liquid.evaluate_fugacity_coef()
    
    equilibrium_equations = 1 - array(x)*liquid_fugacity/array(vapor.mixture.x)/gas_fugacity
    
    equilibrium_equations_eff = [float(equilibrium_equations[i]) 
                                 for i in range(equilibrium_equations.shape[0])]
    
    return equilibrium_equations_eff + [1-sum(x)],liquid,vapor


def get_y_from_flash(x,z,phi):
    return [(z[i]/phi - (1-phi)/phi*x[i]) for i in range(len(x))]

list_of_names = ['H2O', 'N2_2*', 'CO2_2*', 'C1_2*', 'C2_2*', 'C3_2*', 
                'iC4_2*', 'nC4_2*', 'iC5_2*', 'nC5_2*', 'C6_2*', 'C7_2*', 'C8_2*', 'C9_2*', 'C10_2*', 'C11_2*', 
                'C12_2*', 'C13_2*', 'C14_2*', 'C15_2*', 'C16_2*', 'C17_2*', 'C18_2*', 'C19_2*', 'C20+_2*']

#list_of_names.pop(0)

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


composition = [.00552688,.00178958,.54945208,.28497686,.03133262,.02168416,
               .00387928,.00875324,.00268565,.00427715,.00557024,.0050729,
               .00775855,.00686334,.00576918,.00497343,.00457556,.00437662,
               .00358087,.00338193,.00268565,.00218831,.00238725,.00218831,.02427035]



#Aij = [[Aij[i][j] for i in range(len(Aij)) if i > 0] for j in range(len(Aij)) if j > 0]



dict_composition = {list_of_names[i]: composition[i] for i in range(len(composition))}

mixture = Mixture(list_of_species, dict_composition)

#volumn_desviation.pop(0)


pinit = gc_eos_class(mixture, 97.29+273.15, 2.223e4, None, 2, -1, Aij, volumn_desviation, 'liquid')
Pt = []

sumy = []

lin = []

y0 = [0]*2 + [0.5]*2 +[0]*21

#for P in range(10000,40000,50):
    
#    soma, y0 = pinit.test_bubble_T(0+273.15, P, y0)
    
#    Pt.append(P)
#    sumy.append(soma)
#    lin.append(1)
    
#fig = figure()
#p1, = plot(Pt, sumy)
#p2, = plot(Pt, lin)
#fig.tight_layout()
#fig.show()

dt = 5
dp = 100

Porv = []
Pbub = []
Torv = []
Tbub = []
Vorvv = []
Vbubv = []
Vorvl = []
Vbubl = []

delta_c_bub = []
delta_c_orv = []

Zorvl = []
Zorvv = []
Zbubl = []
Zbubv = []

Tb0 = 50
Po0 = 1000

it = 0

Pb = 10000
To = 600

y0 = [0]*2 + [1] + [0] +[1]*21

diff_liq_vap = 2
while Tb0 < 360 :
    Pb, vapb, liqb = pinit.bubble_T(Tb0+273.15, 10000,y0)
    Pbub.append(Pb/100)
    Tbub.append(Tb0)
    Tb0 += dt
    Zbubl.append(liqb.Z)
    Zbubv.append(vapb.Z)
    Vbubv.append(vapb.V)
    Vbubl.append(liqb.V)    
    diff_liq_vap = sum(abs(array(liqb.mixture.x)-array(vapb.mixture.x)))
    

y0 = [0]*2 + [0.05]*2 +[0.3/20]*20 + [1]

diff_liq_vap = 2
while Po0<45000:
    To, vapo, liqo = pinit.dew_P(300+273.15, Po0,y0)
    Porv.append(Po0/100)
    Torv.append(To-273.15)
    Po0 += dp
    Zorvl.append(liqo.Z)
    Zorvv.append(vapo.Z)
    Vorvv.append(vapo.V)
    Vorvl.append(liqo.V) 
    diff_liq_vap = sum(abs(array(liqo.mixture.x)-array(vapo.mixture.x)))

if Torv[-1]-Tbub[-1]>=10:
    Tt = Tbub[-1]+dt
    while Tt <= Torv[-1]:
        if Tt == Torv[-1]:
            Tbub.append(Torv[-1])
            Pbub.append(Porv[-1])
            Zbubl.append(Zbubv[-1])
            Zbubv.append(Zbubl[-1])
            Vbubv.append(Vbubv[-1])
            Vbubl.append(Vbubl[-1])
        else:
            phi = 0.5
            P0 = Porv[-1]*100
            while 1-phi>0.005:
                ptest = pinit.copy_change_conditions(Tt+273.15,P0,None,'liquid')
                phi, phase1, phase2 = ptest.evaluate_flash(liqb.mixture.x,vapb.mixture.x)
                if isnan(phi):
                    phi = 0.5
                print(phi)
                if phase1.V>phase2.V:
                    vapc = phase1
                    liqc = phase2
                else:
                    vapc = phase2
                    liqc = phase1    
                P0 += 50
            Pbub.append(P0/100)
            Tbub.append(Tt)
            Zbubl.append(liqc.Z)
            Zbubv.append(vapc.Z)
            Vbubv.append(vapc.V)
            Vbubl.append(liqc.V)
            delta_c_bub.append(array(vapc.mixture.x)-array(liqc.mixture.x))
        Tt += dt



Pdexp = [202.65,429.01005337,900.06765995,1879.0524575,3957.4123494,8766.70171473,11883.68562813,19420.53346215,43021.55271472,55523.21181754,61265.34049775,66166.10702195,69570.03427229,70247.59751654,70249.84115388,70029.68145147,69413.71188891,68670.49992807,68207.70066475,68200.24671335,68208.29684438,68235.54560042,68240.20005301,68242.52739836]
Tdexp = [427.40089222,455.63452791,485.46086011,515.96826989,544.98756761,566.46826755,568.95265655,559.72985374,474.57621847,403.11450621,359.34639499,308.72660916,246.85405483,205.95752049,203.39897789,177.43202288,150.60489981,125.37503206,101.64732093,98.94277553,95.93950475,93.11827666,92.77934413,92.61810434]
Vdexp = [28.71302611,14.10248991,6.99366689,3.48527848,1.72027157,.80612374,.60226167,.37486885,.17075079,.12963328,.11484851,.10245311,.09136,.08557011,.08523764,.08202305,.07895434,.07623128,.07374252,.07345929,.07314441,.07284777,.07281207,.07279508]


fig2 = figure()
ax = fig2.add_subplot(1, 1, 1)
p1, = plot(Torv, Porv, 'r')
p3, = plot(Tbub, Pbub, 'b')
p2, = plot(pinit.T-273.15, pinit.P/100, 'x')
p4, = plot(Tdexp, [P/100 for P in Pdexp], 'ko', fillstyle = 'none')
config_plot(ax)
ax.set_ylabel('Pressure /bar', fontsize=12)
ax.set_xlabel('Temperature / °C', fontsize=12)
fig2.tight_layout()
fig2.show()

Vmax = 0.6
fig3 = figure()
ax = fig3.add_subplot(1, 1, 1)
Vgorv_l = [V for V in Vorvv if V < Vmax]
Pgorv_l = [Porv[i] for i in range(len(Porv)) if  Vorvv[i] < Vmax]
p1, = plot(Vgorv_l, Pgorv_l,'r')
p2, = plot(pinit.V, pinit.P/100, 'x')
p3, = plot(Vbubl, Pbub, 'b')
Vdexp_l = [V for V in Vdexp if V < Vmax]
Pdexp_l = [Pdexp[i]/100 for i in range(len(Pdexp)) if  Vdexp[i] < Vmax]
p4, = plot(Vdexp_l, Pdexp_l,'ko', fillstyle = 'none')
fig3.tight_layout()
fig3.show()
config_plot(ax)
ax.set_xlabel('Molar Volume /(m³/kmol)', fontsize=12)
ax.set_ylabel('Pressure /bar', fontsize=12)
fig3.tight_layout()
fig3.show()

Vbubt = [[Tbub[i], Pbub[i], Vbubv[i]-Vbubl[i]] for i in range(len(Vbubl))]
Vorvt = [[Torv[i], Porv[i], Vorvv[i]-Vorvl[i]] for i in range(len(Vorvl))]
Vt = Vbubt+Vorvt
for Vs in Vt:
    print(Vs)
