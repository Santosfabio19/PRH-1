# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 00:20:08 2024

@author: Rodrigo Meira
"""
from casadi import *
from numpy import exp, log, array, roots, zeros, linalg
from scipy.optimize import fsolve

R = 8.31446261815324

class water_reference:
    
    def __init__(self,T0,P0,delta_h0,delta_mu0,delta_V0,delta_Cp_fun,int_delta_Cp_fun):
        
        self.T0 = T0
        self.P0 = P0
        self.delta_h0 = delta_h0
        self.delta_V0 = delta_V0
        self.delta_mu0 = delta_mu0
        self.int_delta_Cp_fun = int_delta_Cp_fun
        
    def evaluate_delta_mu_beta_0(self,P,T):
        
        T_average = self.T0/2+T/2
        
        int_delta_h_dT = -self.delta_h0/R*(1/T-1/self.T0) \
                        - float(self.int_delta_Cp_fun(T))
        
        return R*T*(self.delta_mu0/R/self.T0 - int_delta_h_dT + P*self.delta_V0/R/T_average)

class parrish_n_prusnitz:
    
    def __init__(self,list_nu,list_Aki,list_Bki):
        
        self.Aki = list_Aki
        self.Bki = list_Bki
        self.n_i = len(list_Aki[0])
        self.n_k = len(list_Aki)
        self.nu = list_nu
        
    def evaluate_Cki(self,T,k):
        
        return [self.Aki[k][i]/100/T*exp(self.Bki[k][i]/T) for i in range(self.n_i)]
    
    def evaluate_Yki(self,fluid,k):
        
        Cki = array(self.evaluate_Cki(fluid.T,k))
        
        f_i = fluid.evaluate_fugacity_coef()*array(fluid.mixture.x)*fluid.P
        
        return Cki*f_i/(1+sum(Cki*f_i))
    
    def evaluate_delta_mu_H_beta(self,fluid):
        
        delta_mu = 0
        for k in range(self.n_k):
            Yki = self.evaluate_Yki(fluid, k)
            delta_mu += R*fluid.T*self.nu[k]*log(1-sum(Yki))
            
        return delta_mu
        
class hydrate_solver:
    
    def __init__(self,water_reference,hydrate_model,water_name,fluid_solver,T_sat_inter_f):
        
        self.water_ref = water_reference
        self.hydrate_model = hydrate_model
        self.water_name = water_name
        self.fluid_solver = fluid_solver
        self.T_sat_inter_f = T_sat_inter_f
        
        
    def setup_duple_int_dt_fun(self,int2_cp_dt):
        self.int2_cp_dt = int2_cp_dt
    
    def evaluate_delta_mu_0_alpha(self,fluid):
        
        index_water = fluid.mixture.list_of_names.index(self.water_name)
        
        pure_composition = [1 if index_water == i else 0 \
                            for i in range(len(fluid.mixture.x))]
        
        pure_water = fluid.copy_change_x_and_conditions(fluid.T,fluid.P,None,pure_composition,'liquid')
        
        phi_i_sat = fluid.evaluate_fugacity_coef()
        
        f_water_sat = fluid.mixture.x[index_water]*phi_i_sat[index_water]*fluid.P
        
        f_pure_water = pure_water.evaluate_pure_fugacity()*pure_water.P
        
        #print(index_water,f_water_sat,f_pure_water)
        
        return R*fluid.T*(log(f_pure_water/f_water_sat))
        
            
    def select_fluid(self,fluid):
         
        T_sat = self.T_sat_inter_f(fluid.P).__float__() + 273.15
        
        P_water_min = min(self.fluid_solver.water_line_P)
        
        test = (T_sat > fluid.T) or (P_water_min > fluid.P) 
        
        if test:
            
            if self.fluid_solver == None:
                vap = fluid
            else:
                
                (phi,vap,liq) = self.fluid_solver.evaluate_flash(fluid,self.fluid_solver.y0,self.fluid_solver.w0,phi0=0.99999)
                
            delta_mu_gas_sat = 0
        
        else:
            
            index_water = fluid.mixture.list_of_names.index(self.water_name)
            vap = fluid.copy_change_conditions(T_sat,fluid.P,None,'gas')
            
            phi_sat_i = vap.evaluate_fugacity_coef()
            phi_i = fluid.evaluate_fugacity_coef()
            
            x_h_ig_sat = vap.mixture.list_of_species[index_water].evaluate_enthalpy_ig(T_sat)
            int_h_i_ig_sat = x_h_ig_sat/R*(1/T_sat-1/fluid.T)
            
            int_h_ig = array([f(fluid.T)-f(T_sat) for f in self.int2_cp_dt])
            int_h_i_ig = int_h_ig[index_water]
            
            delta_mu_gas_sat = R*fluid.T*(log(phi_i[index_water]/phi_sat_i[index_water]) + int_h_i_ig - int_h_i_ig_sat)            
            
            delta_mu_gas_sat = delta_mu_gas_sat[0][0]
            
        return vap, test, delta_mu_gas_sat
    
    def evaluate_hydrate_function(self,T,P,fluid):
        
        actual_condition = fluid.copy_change_conditions(T,P,None,'gas')
        
        vap, test, delta_mu_gas_sat = self.select_fluid(actual_condition)
        
        delta_mu_0_alpha = self.evaluate_delta_mu_0_alpha(vap)
        
        delta_mu_beta_0 = self.water_ref.evaluate_delta_mu_beta_0(P,T)
        
        if test:
            delta_mu_H_beta = self.hydrate_model.evaluate_delta_mu_H_beta(actual_condition)
        else:
            delta_mu_H_beta = self.hydrate_model.evaluate_delta_mu_H_beta(vap)
        
        return delta_mu_gas_sat+delta_mu_0_alpha+delta_mu_beta_0+delta_mu_H_beta
        
    
    def solve_hydrate(self,T0,P,fluid):
        
        f0 = self.evaluate_hydrate_function(T0, P, fluid)
        
        if f0 < 0:
            sign_dt = 1
        elif f0 > 0:
            sign_dt = -1
        
        Ta = T0
        Tb = T0
        fa = f0
        fb = f0
            
        while abs(f0) > 1e-6:
            
            f0 = self.evaluate_hydrate_function(T0, P, fluid)
            
            if f0 < 0:
                Ta = T0
                fa = f0
            elif f0 >= 0:
                Tb = T0
                fb = f0
            
            if fa*fb > 0:
                T0 += sign_dt
            else:
                T0 = Ta-(Tb-Ta)/(fb-fa)*fa
            
        return T0
    
    def evaluate_hydrate_curve(self,P_list,fluid):
        
        self.T_hydrate = []
        self.P_hydrate = []
        
        dP = P_list[1]-P_list[0]
        
        P = P_list[-1]
        
        condition = True
        
        while condition:
            
            t_solver = lambda T: array(self.evaluate_hydrate_function(T,P,fluid))
            
            t_curve = self.solve_hydrate(350,P,fluid)
            
            self.T_hydrate.append(t_curve-273.15)
            self.P_hydrate.append(P)
            
            if len(self.P_hydrate) > 2:
                Pref = self.fluid_solver.lut_T_ref(t_curve-273.15).__float__()
                if P > Pref:
                    condition = True
                else:
                    condition = False
            else:
                condition = True
            
            P -= dP
            
        self.T_hydrate_C = [T-273.15 for T in self.T_hydrate]
        
        
        
        
        
        