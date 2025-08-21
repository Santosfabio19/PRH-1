# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 13:11:08 2024

@author: rodri
"""

class pr_class:
    
    def __init__(self,mixture,T,P,V,Aij,phase):
        self.mixture = mixture
        self.T = T
        self.P = P
        self.V = V
        self.Aij = Aij
        self.phase = phase
        
        if [T,P,V].count(None) > 1:
            print("It is necessary 2 values in T, P or V")
        
        self.pr_parameters()
        
        if self.T==None:
            if self.phase == 'liquid':
                t_solver = lambda Z: array(self.evaluate_eos_T_liq(Z))
                res = fsolve(t_solver,0.25)
            elif self.phase == 'gas':
                t_solver = lambda Z: array(self.evaluate_eos_T_gas(Z))
                res = fsolve(t_solver,0.85)
            self.Z = res[0]
            self.T = self.P*self.V/R/self.Z
        elif self.P==None:
            
            self.P = self.evaluate_eos_P()
            self.Z = self.P*self.V/self.T/R
        elif self.V == None:
            if self.phase == 'liquid':
                #v_solver = lambda Z: array(self.evaluate_eos_V_liq(Z))
                #res = fsolve(v_solver,0.1)
                self.Z = min([y.real for y in self.evaluate_z_roots() if y.imag==0])
            elif self.phase == 'gas':
                # v_solver = lambda Z: array(self.evaluate_eos_V_gas(Z))
                # res = fsolve(v_solver,1)
                self.Z = max([y.real for y in self.evaluate_z_roots() if y.imag==0])
            # self.Z = res[0]
            self.V = self.Z*self.T*R/self.P
        
        self.rho = 1/self.V
        
        self.Zc = 0.307401
        self.Vc = self.b_m/0.253077
        b_l = self.b_m/self.Vc
        self.Pc = 0.816619*self.a_m*b_l/self.Vc**2
        self.Tc = self.Pc*self.Vc/self.Zc/R
        
    def evaluate_eos_T_gas(self,Z):
        
        self.T = self.P*self.V/Z/R        
        self.evaluate_par_a()
        
        beta = self.b_m*self.P/self.T/R
        
        q = self.a_m/self.b_m/self.T/R
        
        return Z - 1 - beta + q*beta*(Z-beta)/( Z**2 + 2*Z*beta - beta**2 )
    
    def evaluate_eos_T_liq(self,Z):
        
        self.T = self.P*self.V/Z/R
        self.evaluate_par_a()
        
        beta = self.b_m*self.P/self.T/R
        
        q = self.a_m/self.b_m/self.T/R
        
        return Z - beta + ( Z**2 + 2*Z*beta - beta**2 )*(1 + beta - Z)/q/beta
    
    
    def evaluate_eos_P(self):
        
        self.evaluate_par_a()
        
        return self.T*R/(self.V-self.b_m) - self.a_m/( self.V**2 + 2*self.V*self.b_m - self.b_m**2 )

    def evaluate_eos_V_gas(self,Z):
        
        self.V = self.T*Z*R/self.P        
        self.evaluate_par_a()
        
        beta = self.b_m*self.P/self.T/R
        q = self.a_m/self.b_m/self.T/R
        
        return Z - 1 - beta + q*beta*(Z-beta)/( Z**2 + 2*Z*beta - beta**2 )
    
    def evaluate_eos_V_liq(self,Z):
        
        self.V = self.T*Z*R/self.P        
        self.evaluate_par_a()
        
        beta = self.b_m*self.P/self.T/R
        q = self.a_m/self.b_m/self.T/R
        
        return Z - beta - ( Z**2 + 2*Z*beta - beta**2 )*(1 + beta - Z)/q/beta
       
    def pr_parameters(self):
        
        self.b = [R*sp.Tc*0.07780/sp.Pc for sp in self.mixture.list_of_species]
        
        self.a_c = [R**2*sp.Tc**2*0.45724/sp.Pc for sp in self.mixture.list_of_species]
        
        self.kappa = [sp.kappa for sp in self.mixture.list_of_species]
        
    def evaluate_a(self):
        Tr = [self.T/sp.Tc for sp in self.mixture.list_of_species]
        self.a = [self.a_c[i]*(1+self.kappa[i]*(1-Tr[i]**0.5))**2 for i in range(len(Tr))]
        
    def evaluate_der_a(self):
        
        der_sqrt_a= [-self.a[i]**0.5*self.kappa[i]/(self.T*self.mixture.list_Tc[i])**0.5 for i in range(len(self.kappa))]
        der_A = [[self.a[i]**.5*der_sqrt_a[j]+self.a[j]**.5*der_sqrt_a[i] for i in len(der_sqrt_a)] 
                      for j in le(der_sqrt_a)]
        
        Q = array(der_A)*(1-array(self.Aij))
        
        x = array(self.mixture.x)
        
        self.der_a = x.dot(Q.dot(x.transpose())) 
        
    def evaluate_der_sec_a(self):
        
        der_sec_A = [[0.25*self.T**-1.5*((1+self.kappa[i])*self.kappa[j]/self.mixture.list_Tc[j]**0.5+
                                         (1+self.kappa[j])*self.kappa[i]/self.mixture.list_Tc[i]**0.5) for i in len(der_sqrt_a)] 
                     for j in le(der_sqrt_a)]
        
        Q = array(der_sec_A)*(array(self.a_c).transpose().dot(array(self.a_c)))**0.5*(1-array(self.Aij))
        
        x = array(self.mixture.x)
        
        self.der_sec_a = x.dot(Q.dot(x.transpose())) 
        
    def evaluate_par_a(self):
        
        self.evaluate_a()
        
        Q = array([[(self.a[i]*self.a[j])**0.5*(1-self.Aij[i][j]) for i in range(len(self.a))] 
                   for j in range(len(self.a))])
        x = array(self.mixture.x)
        self.a_m = x.dot(Q.dot(x.transpose()))
        self.b_m = array(self.b).dot(x)
        
        self.b_bar = array(self.b)
        self.a_bar = 2*Q.dot(x.transpose())-self.a_m
    
    def evaluate_fugacity_coef(self):
        
        beta = self.b_m*self.P/self.T/R
        q = self.a_m/self.b_m/self.T/R
        
        ln_phi_c_i = array(self.b)/self.b_m*(self.Z-1) - log(self.Z-beta) \
                     + q*(1 + self.a_bar/self.a_m  - array(self.b)/self.b_m) \
                     /2**1.5*log((self.Z+beta*(1-2**.5))/(self.Z+beta*(1+2**.5)))
        
        return exp(ln_phi_c_i)
    
    def evaluate_z_roots(self):
        
        self.evaluate_par_a()
        
        beta = self.b_m*self.P/self.T/R
        q = self.a_m/self.b_m/self.T/R
        
        return roots([1,beta-1,(q-2-3*beta)*beta,(beta**3+beta**2-q*beta**2)])
        
        
def cp_ig_fun_pure(par,T):
    return par[0] + par[1]*T + par[2]*T**2 + par[3]*T**3 + par[4]*T**4 + par[5]*T**5

def enthapy_ig_fun_pure(par,T):
    return par[0]*T + par[1]*T**2 + par[2]*T**3 + par[3]*T**4 + par[4]*T**5 + par[5]*T**6

def pvap_fun_pure(par,T):
    return exp(par[0] + par[1]/(T+par[2]) + par[3]*log(T) + par[4]*T**par[5])

def kappa_fun(par,w):
    return par[0] + par[1]*w  + par[2]*w**2 + par[3]*w**4


def evaluate_drew_P(T,Pr,x,gas_base,Aij):
    
    liquid_composition = {gas_base.mixture.list_of_names[i]: x[i] for i in range(len(x))}
    
    liquid_mixture = Mixture(gas_base.mixture.list_of_species, liquid_composition)
    
    P = Pr*gas_base.Pc
    
    liquid = pr_class(liquid_mixture, T, P, None, Aij, 'liquid')
    
    vapor = pr_class(gas_base.mixture, T, P, None, Aij, 'gas')
    
    gas_fugacity = vapor.evaluate_fugacity_coef()
    
    liquid_fugacity = liquid.evaluate_fugacity_coef()
    
    equilibrium_equations = 1 - array(liquid.mixture.x)*liquid_fugacity/array(vapor.mixture.x)/gas_fugacity
    
    equilibrium_equations_eff = [float(equilibrium_equations[i]) 
                                 for i in range(equilibrium_equations.shape[0])]
    
    return equilibrium_equations_eff + [1-sum(x)]

def evaluate_flash(phi,x,Z_g,Z_l,feed,Aij):
    
    z = feed.mixture.x
    
    y = [(z[i]/phi - (1-phi)/phi*x[i]) for i in range(len(x))]
    
    liquid_composition = {feed.mixture.list_of_names[i]: x[i] for i in range(len(x))}
    liquid_mixture = Mixture(feed.mixture.list_of_species, liquid_composition)
    
    vapor_composition = {feed.mixture.list_of_names[i]: y[i] for i in range(len(x))}
    vapor_mixture = Mixture(feed.mixture.list_of_species, vapor_composition)
    
    V_l = Z_l*feed.T*R/feed.P
    
    liquid = pr_class(liquid_mixture, feed.T, None, V_l, Aij, 'liquid')
    
    V_g = Z_g*feed.T*R/liquid.P
    
    vapor = pr_class(vapor_mixture, feed.T, None, V_g, Aij, 'gas')
    
    gas_fugacity = vapor.evaluate_fugacity_coef()
    
    liquid_fugacity = liquid.evaluate_fugacity_coef()
    
    equilibrium_equations = 1 - array(liquid.mixture.x)*liquid_fugacity/array(vapor.mixture.x)/gas_fugacity
    
    pressure_equation = float((feed.P - liquid.P)*liquid.V/liquid.T/R)
    
    equilibrium_equations_eff = [float(equilibrium_equations[i]) 
                                 for i in range(equilibrium_equations.shape[0])]
    
    return  [pressure_equation] + [Z_g-float(vapor.Z), Z_l-float(liquid.Z)] + equilibrium_equations_eff[0:-1] + [1-sum(x)]


def evaluate_flash2(phi,x,feed,Aij):
    
    z = feed.mixture.x
    
    y = [(z[i]/phi - (1-phi)/phi*x[i]) for i in range(len(x))]
    y = [e/sum(y) for e in y]
    
    
    liquid_composition = {feed.mixture.list_of_names[i]: x[i] for i in range(len(x))}
    liquid_mixture = Mixture(feed.mixture.list_of_species, liquid_composition)
    
    vapor_composition = {feed.mixture.list_of_names[i]: y[i] for i in range(len(x))}
    vapor_mixture = Mixture(feed.mixture.list_of_species, vapor_composition)
    
    liquid = pr_class(liquid_mixture, feed.T, feed.P, None, Aij, 'liquid')
    
    vapor = pr_class(vapor_mixture, feed.T, feed.P, None, Aij, 'gas')
    
    gas_fugacity = vapor.evaluate_fugacity_coef()
    
    liquid_fugacity = liquid.evaluate_fugacity_coef()
    
    Ki = liquid_fugacity/gas_fugacity
    
    equilibrium_equations = 1 - array(liquid.mixture.x)*Ki/array(vapor.mixture.x)
    
    equilibrium_equations_eff = [float(equilibrium_equations[i]) 
                                 for i in range(equilibrium_equations.shape[0])]
    
    phi_equation = sum(array(feed.mixture.x*(Ki-1))/(1+phi*(Ki-1)))        
    
    return  [phi_equation] + equilibrium_equations_eff[0:-1] + [1-sum(x)]

def get_y_from_flash(x,z,phi):
    return [(z[i]/phi - (1-phi)/phi*x[i]) for i in range(len(x))]
