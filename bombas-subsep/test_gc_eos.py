

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


composition = [3.44945577070149e-003, 2.48960745236753e-003, 0.369253301096065, 0.350151113126578, 4.69925780689622e-002,
               3.16949914724554e-002, 5.27916552903917e-003, 1.38978027217423e-002, 4.97921253693917e-003, 7.57880126909168e-003,
               9.96842256582064e-003, 1.39977840819956e-002, 1.33978786580008e-002, 9.96842134103575e-003, 8.86859527510853e-003,
               7.47881522052249e-003, 6.77892590685044e-003, 7.07887819148949e-003, 6.27900478629509e-003, 5.67909973211496e-003,
               4.97921057183309e-003, 3.88938329490607e-003, 4.28931982002196e-003, 3.68941490974787e-003, 5.78908166003149e-002]

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

dt = 20
dp = 500

Porv = []
Pbub = []
Torv = []
Tbub = []
Vorvv = []
Vbubv = []
Vorvl = []
Vbubl = []

Zorvl = []
Zorvg = []
Zbubl = []
Zbubg = []

Tb0 = -30
Po0 = 1000

it = 0

Pb = 10000
To = 600

y0 = [0]*2 + [0.5]*2 +[0]*21

while Tb0 < 450 :
    Pb, vapb, liqb = pinit.bubble_T(Tb0+273.15, 10000,y0)
    Pbub.append(Pb/100)
    Tbub.append(Tb0)
    Tb0 += dt
    Zbubl.append(liqb.Z)
    Zbubg.append(vapb.Z)
    Vbubv.append(vapb.V)
    Vbubl.append(liqb.V)    

    

y0 = [0]*2 + [0.1]*2 +[0.3/20]*20 + [0.8]

while Po0 < 45000 :
    To, vapo, liqo = pinit.dew_P(400+273.15, Po0,y0)
    Porv.append(Po0/100)
    Torv.append(To-273.15)
    Po0 += dp
    Zorvl.append(liqo.Z)
    Zorvg.append(vapo.Z)
    Vorvv.append(vapo.V)
    Vorvl.append(liqo.V) 


Pbexp = [87882.75445381,87882.77855322,85770.82023567,83709.63847936,81697.98953877,77818.55736412,74123.34017546,70603.59051708,67250.97630658,64057.56111077,61015.78535832,55358.69168532,45569.371828,43480.57290125,43323.95119257,43516.20260796,44628.25455112,46597.59565864,50070.73897798,53524.7853482,58642.69970031,59730.65896627,59831.33554025,59814.21996524,59532.78514188,58948.72052478,58412.47119553,57770.44760533,57287.96202634,56764.6534188,56202.73640133,55604.32522233,54971.45242861,54306.08265371,53610.12270065,52885.4287778,52133.81163294,51357.03975628,49734.90779139,48937.79927667]
Tbexp = [-33.00351598,-33.00350937,-32.59924267,-32.17204044,-31.72074288,-30.74056593,-29.64657488,-28.42355407,-27.05203122,-25.50647115,-23.75238122,-19.40247491,-2.34790047,12.07139002,18.38978231,26.29959006,41.57970683,59.86264147,88.90283146,120.47545575,192.12062375,231.96045301,248.55191863,255.46360679,277.65529314,298.92340012,312.69316528,326.17875338,335.02287487,343.75552236,352.3799302,360.89827,369.31185018,377.62127529,385.82657501,393.92730961,401.92265615,409.81148258,425.263836,432.37312492]
Vbexp = [.08547519,.08547519,.08558411,.08569473,.08580716,.08603796,.08627763,.08652763,.0867898,.08706656,.0873612,.08802463,.09011804,.09158096,.0921766,.09289643,.09423612,.0957972,.09828457,.10111179,.10847498,.11341126,.11570367,.11670613,.12013432,.12375935,.12631061,.12898565,.13084442,.1327682,.13476108,.13682723,.13897097,.14119677,.14350929,.14591341,.14841424,.15101716,.15655233,.15932638]

Pdexp = [202.65,429.01005337,894.71043915,1843.41049031,3801.17134633,8188.14997806,10552.66732812,12701.9133587,19703.9139982,23359.42628889,27116.53426882,30956.11899415,34844.97211223,38735.6435436,42566.49298634,46262.21554919,48937.79927667]
Tdexp = [461.76508695,493.88620636,527.69289194,561.91644594,593.83560815,616.48534474,618.51739694,617.22907159,600.88473526,587.55449991,571.40886099,552.67875658,531.52481207,508.07918209,482.46747344,454.81845592,432.37312492]
Vdexp = [30.06643757,14.79593615,7.38869835,3.72913938,1.87561664,.9007143,.7060314,.59081571,.38747531,.32910404,.28524261,.25117989,.2240834,.20213468,.18410366,.1691217,.15932638]


fig2 = figure()
ax = fig2.add_subplot(1, 1, 1)
p1, = plot(Torv, Porv, 'r')
p2, = plot(Tbexp, [P/100 for P in Pbexp], 'ko', fillstyle = 'none')
p3, = plot(Tbub, Pbub, 'b')
p4, = plot(Tdexp, [P/100 for P in Pdexp], 'ko', fillstyle = 'none')
p5, = plot(pinit.T-273.15, pinit.P/100, 'x')
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
p2, = plot(Vbexp, [P/100 for P in Pbexp],'ko', fillstyle = 'none')
p3, = plot(Vbubl, Pbub, 'b')
p4, = plot(pinit.V, pinit.P/100, 'x')
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
