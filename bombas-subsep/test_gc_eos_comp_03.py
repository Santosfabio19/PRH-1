''

from eos_database import *
from casadi import *
from numpy import exp, log, array, roots
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
    formatter.set_scientific(False)
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

# list_of_names.pop(0)

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


composition = [.00773653,.00198397,.59036782,.31137302,.02411984,.01955394,
               .00317627,.0081392,.00258072,.0040696,.00148888,.00208443,
               .00248146,.00208443,.00178665,.00148888,.00138962,.00138962,
               .0011911,.0011911,.00089333,.00079407,.00079407,.00069481,.00714662]



# composition.pop(0)



# Aij = [[Aij[i][j] for i in range(len(Aij)) if i > 0] for j in range(len(Aij)) if j > 0]

dict_composition = {list_of_names[i]: composition[i] for i in range(len(composition))}

mixture = Mixture(list_of_species, dict_composition)

# volumn_desviation.pop(0)


pinit = gc_eos_class(mixture, 97.29+273.15, 2.223e4, None, 2, -1, Aij, volumn_desviation, 'liquid')
Pt = []

sumy = []

lin = []

y0 = [0] + [0.5]*2 +[0]*21

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

dt = 10
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

Tb0 = 60
Po0 = 1000

it = 0

Pb = 10000
To = 600


y0 = [0]*2 + [0.2]+[0] +[0.3/20]*20 + [1]

while Po0 < 48000 :
    To, vapo, liqo = pinit.dew_P(400+273.15, Po0,y0)
    y0 = liqo.mixture.x
    Porv.append(Po0/100)
    Torv.append(To-273.15)
    Po0 += dp
    Zorvl.append(liqo.Z)
    Zorvv.append(vapo.Z)
    Vorvv.append(vapo.V)
    Vorvl.append(liqo.V)
    delta_c_orv.append(array(vapo.mixture.x)-array(liqo.mixture.x))


y0 = [0]*2 + [0.5]+[0.05] +[0]*20 + [1]

while Tb0 < 360:
    Pb, vapb, liqb = pinit.bubble_T(Tb0+273.15, 10000,y0)
    Pbub.append(Pb/100)
    Tbub.append(Tb0)
    Tb0 += dt
    Zbubl.append(liqb.Z)
    Zbubv.append(vapb.Z)
    Vbubv.append(vapb.V)
    Vbubl.append(liqb.V)    
    delta_c_bub.append(array(vapb.mixture.x)-array(liqb.mixture.x))

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
            while 1-phi>0.001:
                ptest = pinit.copy_change_conditions(Tt+273.15,P0,None,'liquid')
                phi, phase1, phase2 = ptest.evaluate_flash(liqb.mixture.x,vapb.mixture.x)
                if phase1.V>phase2.V:
                    vapc = phase1
                    liqc = phase2
                else:
                    vapc = phase2
                    liqc = phase1    
                P0 += 100
            Pbub.append(P0/100)
            Tbub.append(Tt)
            Zbubl.append(liqc.Z)
            Zbubv.append(vapc.Z)
            Vbubv.append(vapc.V)
            Vbubl.append(liqc.V)
            delta_c_bub.append(array(vapc.mixture.x)-array(liqc.mixture.x))
        Tt += dt

Tenv = Tbub + Torv
Penv = Pbub + Porv

fun_crit = []

for i in range(len(Penv)):
    ptest = pinit.copy_change_conditions(Tenv[i]+273.15,Penv[i]*100,None,'liquid')
    fun_crit.append(ptest.critical_case_V(800,[ptest.V],False))    


Pdexp = [202.65,334.113365507,549.666572771,903.032921413,1904.377242326,4077.538156475,9225.567721217,11338.340506564,20873.158389363,47226.22545444,61331.664930986,73071.666192785,77951.398766377,87196.033166665,92221.702521948,97537.033591776,103158.721447581]
Tdexp = [385.457230425,401.13950209,417.417879897,434.15492308,459.49350801,483.20256755,499.289784593,500.214614366,488.012690062,397.582619495,319.995099722,203.483567006,101.491032243,64.868660593,57.733921297,52.810417487,49.226390269]
Vdexp = [27.014076307,16.774780204,10.443384648,6.512503772,3.20266893,1.550078937,.70803459,.580011515,.319005418,.137201365,.100035063,.073304897,.059021471,.053981249,.052808795,.051897652,.051143462]



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
