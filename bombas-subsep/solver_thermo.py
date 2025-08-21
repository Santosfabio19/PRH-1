# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 00:05:17 2024

@author: Rodrigo Meira
"""
from eos_database import *
from casadi import *
from numpy import exp, log, array, roots, zeros, linalg, arange
from math import isnan
from scipy.optimize import fsolve
from species_builder import R, Species, Mixture
from matplotlib.pyplot import plot, figure
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter
import scipy.linalg as la
from numpy.linalg import svd, det
import sys



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
    # axes.xaxis.set_tick_params(which='both', right='off', top='off', direction='out', width=1)
    # axes.yaxis.set_tick_params(which='both', right='off', top='off', direction='out', width=1)

    # axes.spines['top'].set_visible(False)


class solver_eos:
    
    def __init__(self,fluid):
        
        self.fluid = fluid # EOS class
        self.dew_tol = 1e-5
        self.bub_tol = 1e-5
        self.delta_fun = 1e-5
        self.flash_subs = 0.0005
        self.flash_tol = 1e-6
        self.max_run = 1e5
        self.critical_initial = 'default'
        self.step_flash = 1
        T0 = sum(array(fluid.mixture.list_Tc)*array(fluid.mixture.x))*3
        V0 = sum(array([sp.Vc for sp in fluid.mixture.list_of_species])*array(fluid.mixture.x))*1.5
        self.critical_init_value = [T0, V0]
        self.critical_evaluation = False
        self.isochoric_evaluation = False
        self.water_evaluation = False
        self.P_max_isochoric = 60000
        self.P_max_water = 60000
    
    def evaluate_bubble_T(self, fluid, Tb, P0, y0):
        P = P0
        vap = fluid.copy_change_x_and_conditions(Tb, P, None, y0, 'gas')
        phi_vap = vap.evaluate_fugacity_coef()
        
        liq = fluid.copy_change_conditions(Tb, P, None, 'liquid')
        fliq = liq.evaluate_fugacity_coef()*array(self.fluid.mixture.x)
        #phi_vap = array([1.25]*len(self.fluid.mixture.x))
        it = 0
        k = liq.evaluate_fugacity_coef()/phi_vap
        k = array([0.0001]*len(self.fluid.mixture.x))
        k[2] = 1.86
        k[3] = 1.86
        y = k*array(self.fluid.mixture.x)
        fa = sum(y)-1
        
        fb = sum(y)-1
        delta_x = 5
        P_ant = P0
        f = 8
        Pa = P0
        Pb = P0
        while abs(f) > self.bub_tol and it < self.max_run  and delta_x > self.bub_tol:
            y = k*array(fluid.mixture.x)
            liq = fluid.copy_change_conditions(Tb, P, None, 'liquid')
            fliq = liq.evaluate_fugacity_coef()*array(fluid.mixture.x)
            sum_y_0 = 4.5
            it_2 = 0
            while abs(sum_y_0-sum(y)) > self.delta_fun and it_2 < 100:
                sum_y_0 = sum(y)
                it_2 += 1
                vap = fluid.copy_change_x_and_conditions(Tb, P, None, y, 'gas')
                phi_vap = vap.evaluate_fugacity_coef()
                y = fliq/phi_vap
            
            
            f = sum(y)-1
            # if f > 1000:
            #     P = P_ant-10
            #     f = 0


            if abs(f) > self.bub_tol:
                
                if fb < 0 and fa > 0:
                    delta_x = abs(Pa-Pb)
                    if sum(y)-1 < 0:
                        fb = sum(y)-1
                        Pb = P
                    elif sum(y)-1 > 0:
                        fa = sum(y)-1
                        Pa = P

                    P = -(fa-Pa*(fa-fb)/(Pa-Pb))/((fa-fb)/(Pa-Pb))
                    P_ant = P

                else:

                    if sum(y)-1 < 0:
                        fb = sum(y)-1
                        Pb = P
                    else:
                        fa = sum(y)-1
                        fb = sum(y)-1
                        Pa = P
                        Pb = P
                        P += 10
                        P_ant = P
            
            it += 1
            
        return P, vap, liq
    
    def test_bubble(self,Tb,P,y0):

        vap = self.fluid.copy_change_x_and_conditions(Tb, P, None, y0, 'gas')
        phi_vap = vap.evaluate_fugacity_coef()
        
        liq = self.fluid.copy_change_conditions(Tb, P, None, 'liquid')
        fliq = liq.evaluate_fugacity_coef()*array(self.fluid.mixture.x)
        
        k = liq.evaluate_fugacity_coef()/phi_vap

        k = array([10.1]*len(self.fluid.mixture.x))
        k[0] = 0.26
        k[-1] = 0.26
        
        it = 0

        y = k*array(self.fluid.mixture.x)
        y_init = y/sum(y)
        print(k*array(self.fluid.mixture.x))
        liq = self.fluid.copy_change_conditions(Tb, P, None, 'liquid')
        fliq = liq.evaluate_fugacity_coef()*array(self.fluid.mixture.x)
        sum_y_0 = 4.5
        while abs(sum_y_0-sum(y)) > self.delta_fun and it < 100:
            
            sum_y_0 = sum(y)
            vap = self.fluid.copy_change_x_and_conditions(Tb, P, None, y, 'gas')
            phi_vap = vap.evaluate_fugacity_coef()
            
            y = fliq/phi_vap
            it += 1
        
        f = sum(y)-1
        
        return f,sum_y_0,y,y_init,it
    
    def test_dew_T(self,Tb,P,x0):

        liq = self.fluid.copy_change_x_and_conditions(Tb, P, None, x0, 'liquid')
        phi_liq = liq.evaluate_fugacity_coef()
        
        vap = self.fluid.copy_change_conditions(Tb, P, None, 'gas')
        fvap = vap.evaluate_fugacity_coef()*array(self.fluid.mixture.x)
        
        k = vap.evaluate_fugacity_coef()/phi_liq

        #k = array([0.023]*len(self.fluid.mixture.x))
        #k[0] = 1.26
        #k[-3] = 1.26
        #k[-2] = 1.36
        #k[-1] = 1.36
        
        it = 0

        x = 1/k*array(self.fluid.mixture.x)
        x_init = x/sum(x)
        print(phi_liq*array(x0)-fvap)
        vap = self.fluid.copy_change_conditions(Tb, P, None, 'gas')
        fvap = liq.evaluate_fugacity_coef()*array(self.fluid.mixture.x)
        sum_x_0 = 4.5
        while abs(sum_x_0-sum(x)) > self.delta_fun and it < 100:
            
            sum_x_0 = sum(x)
            liq = self.fluid.copy_change_x_and_conditions(Tb, P, None, x, 'liquid')
            phi_liq = liq.evaluate_fugacity_coef()
            
            x = fvap/phi_liq
            
            it += 1
        
        f = sum(x)-1
        
        return f,sum_x_0,x[0]/sum(x),x[-1]/sum(x),x_init[0],x_init[-1],it
        
    
    def bubble_T(self, Tb, P0, y0):
        P = P0
        vap = self.fluid.copy_change_x_and_conditions(Tb, P, None, y0, 'gas')
        phi_vap = vap.evaluate_fugacity_coef()
        
        liq = self.fluid.copy_change_conditions(Tb, P, None, 'liquid')
        fliq = liq.evaluate_fugacity_coef()*array(self.fluid.mixture.x)
        
        #phi_vap = array([1.25]*len(self.fluid.mixture.x))
        it = 0
        k = liq.evaluate_fugacity_coef()/phi_vap
        #k = array([0.0001]*len(self.fluid.mixture.x))
        #k[0] = 1.86
        #k[-1] = 1.86
        y = k*array(self.fluid.mixture.x)
        fa = sum(y)-1
        
        fb = sum(y)-1
        delta_x = 5
        P_ant = P0
        f = 8
        Pa = P0
        Pb = P0
        while abs(f) > self.bub_tol and it < self.max_run  and delta_x > self.bub_tol:
            y = k*array(self.fluid.mixture.x)
            liq = self.fluid.copy_change_conditions(Tb, P, None, 'liquid')
            fliq = liq.evaluate_fugacity_coef()*array(self.fluid.mixture.x)
            sum_y_0 = 4.5
            it_2 = 0
            while abs(sum_y_0-sum(y)) > self.delta_fun and it_2 < 100:
                sum_y_0 = sum(y)
                it_2 += 1
                vap = self.fluid.copy_change_x_and_conditions(Tb, P, None, y, 'gas')
                phi_vap = vap.evaluate_fugacity_coef()
                y = fliq/phi_vap
                
            
            f = sum(y)-1
            if f > 1000:
                P = P_ant-10
                f = 0
            
            
            if abs(f) > self.bub_tol:
                
                if fb < 0 and fa > 0:
                    delta_x = abs(Pa-Pb)
                    if sum(y)-1 < 0:
                        fb = sum(y)-1
                        Pb = P
                    elif sum(y)-1 > 0:
                        fa = sum(y)-1
                        Pa = P
            
                    P = -(fa-Pa*(fa-fb)/(Pa-Pb))/((fa-fb)/(Pa-Pb))
                    P_ant = P

                else:
        
                    if sum(y)-1 < 0:
                        fb = sum(y)-1
                        Pb = P
                    else:
                        fa = sum(y)-1
                        fb = sum(y)-1
                        Pa = P
                        Pb = P
                        P += 10
                        P_ant = P
            
            it += 1
        return P, vap, liq

    def dew_P(self, T0, Pd, x0):

        #phi_liq = array([0.05]*len(self.fluid.mixture.x))
        
        tol = 1e-5
        T = T0
        liq = self.fluid.copy_change_x_and_conditions(T,6500, None, x0, 'liquid')
        phi_liq = liq.evaluate_fugacity_coef()
        #phi_liq = array([0.25]*len(self.fluid.mixture.x))
        vap = self.fluid.copy_change_conditions(T, 6500, None, 'gas')
        fvap = vap.evaluate_fugacity_coef()*array(self.fluid.mixture.x)

        it = 0
        x = fvap/phi_liq
        sum_x_0 = sum(x)
        
        fa = sum(x)-1
        fb = sum(x)-1
        Ta = T0
        Tb = T0
        delta_x = 5
        while abs(sum(x)-1) > self.dew_tol and it < self.max_run and delta_x > self.dew_tol:

            vap = self.fluid.copy_change_conditions(T, Pd, None, 'gas')
            fvap = vap.evaluate_fugacity_coef()*array(self.fluid.mixture.x)
            sum_x_0 = 4.5
            it_2 = 0
            while abs(sum_x_0-sum(x)) > 1e-4 and it_2 < 100: #self.delta_fun:
                sum_x_0 = sum(x)
                liq = self.fluid.copy_change_x_and_conditions(
                    T, Pd, None, x, 'liquid')
                phi_liq = liq.evaluate_fugacity_coef()
                x = fvap/phi_liq
                it_2 += 1
            
            f = sum(x)-1
            if abs(sum(x)-1) > self.dew_tol:
                if fb < 0 and fa > 0:
                    delta_x = abs(Ta-Tb)
                    if sum(x)-1 < 0:
                        fb = sum(x)-1
                        Tb = T
                    elif sum(x)-1 > 0:
                        fa = sum(x)-1
                        Ta = T
                    T = -(fa-Ta*(fa-fb)/(Ta-Tb))/((fa-fb)/(Ta-Tb))
                    
                else:
                    
                    if sum(x)-1 < 0:
                        fb = sum(x)-1
                        Tb = T
                    else:
                        fa = sum(x)-1
                        fb = sum(x)-1
                        Ta = T
                        Tb = T
                        T += 1
            it += 1
        return T, vap, liq
    
    def dew_T(self, Td, P0, x0):

        #phi_liq = array([0.05]*len(self.fluid.mixture.x))
        
        tol = 1e-5
        P = P0
        liq = self.fluid.copy_change_x_and_conditions(Td,P0, None, x0, 'liquid')
        phi_liq = liq.evaluate_fugacity_coef()
        #phi_liq = array([0.25]*len(self.fluid.mixture.x))
        vap = self.fluid.copy_change_conditions(Td, 1000, None, 'gas')
        fvap = vap.evaluate_fugacity_coef()*array(self.fluid.mixture.x)

        it = 0
        x = fvap/phi_liq
        sum_x_0 = sum(x)
        
        fa = sum(x)-1
        fb = sum(x)-1
        Pa = P0
        Pb = P0
        delta_x = 5
        
        while abs(sum(x)-1) > self.dew_tol and it < self.max_run and delta_x > self.dew_tol:

            vap = self.fluid.copy_change_conditions(Td, P, None, 'gas')
            fvap = vap.evaluate_fugacity_coef()*array(self.fluid.mixture.x)
            sum_x_0 = 4.5
            it_2 = 0
            while abs(sum_x_0-sum(x)) > 1e-4 and it_2 < 100: #self.delta_fun:
                sum_x_0 = sum(x)
                liq = self.fluid.copy_change_x_and_conditions(
                    Td, P, None, x, 'liquid')
                phi_liq = liq.evaluate_fugacity_coef()
                x = fvap/phi_liq
                it_2 += 1
            
            f = sum(x)-1
            if abs(sum(x)-1) > self.dew_tol:
                if fb < 0 and fa > 0:
                    delta_x = abs(Ta-Tb)
                    if sum(x)-1 < 0:
                        fb = sum(x)-1
                        Pb = T
                    elif sum(x)-1 > 0:
                        fa = sum(x)-1
                        Pa = T
                    P = -(fa-Ta*(fa-fb)/(Ta-Tb))/((fa-fb)/(Ta-Tb))
                    
                else:
                    
                    if sum(x)-1 < 0:
                        fb = sum(x)-1
                        Tb = P
                    else:
                        fa = sum(x)-1
                        fb = sum(x)-1
                        Ta = P
                        Tb = P
                        P += 100
            it += 1
        return P, vap, liq
    
    
    def test_water_fuga(self,fluid,water_name,P,T_range):
        
        case = []
        
        T_min = T_range[0]
        T_max = T_range[1]

        step_T_init = (T_max-T_min)/30
        
        for T in arange(T_min,T_max+step_T_init,step_T_init):
        
            fluid = self.fluid.copy_change_conditions(T+273.15,P,None,'gas')
            
            index_water = fluid.mixture.list_of_names.index(water_name)
            
            pure_composition = [1 if index_water == i else 0 \
                                for i in range(len(fluid.mixture.x))]
            
            pure_water = fluid.copy_change_x_and_conditions(T+273.15,P,None,pure_composition,'liquid')
            
            phi_i_sat = fluid.evaluate_fugacity_coef()
            
            f_water_sat = fluid.mixture.x[index_water]*phi_i_sat[index_water]
            
            f_pure_water = pure_water.evaluate_pure_fugacity()
            
            case.append(f_pure_water-f_water_sat)
            
            print(f_pure_water,f_water_sat,log(f_pure_water/f_water_sat)*(T+273.15)*8.314)
            
        return case
        
    

    def test_flash(self,phi,fluid,y0,x0):
        x = x0
        y = y0
        z = array(fluid.mixture.x)
        fun_phi_0 = 4.5
        while abs(fun_phi_0-sum(y)+sum(x)) > self.delta_fun:
            fun_phi_0 = sum(y)-sum(x)
            vap = fluid.copy_change_x_and_conditions(
                fluid.T, fluid.P, None, y, 'gas')
            liq = fluid.copy_change_x_and_conditions(
                fluid.T, fluid.P, None, x, 'liquid')
            phi_vap = vap.evaluate_fugacity_coef()
            phi_liq = liq.evaluate_fugacity_coef()
            ki = phi_liq/phi_vap
            x = z/(1+phi*(ki-1))
            y = z*ki/(1+phi*(ki-1))
        f = sum(y)-sum(x)
        
        return f
    
    def solve_phi(self,phi,zi,ki):
        
        return sum(zi*(ki-1)/(1+phi*(ki-1)))
        
    
    def evaulate_PT_flash(self,fluid,T,P,y0=[],x0=[],phi0=0.5):
        new_fluid = fluid.copy_change_conditions(
            T+273.15, P, None, 'gas')
        
        return self.evaluate_flash(new_fluid,y0,x0,phi0)
        
    
    def evaluate_flash(self,fluid,y0=[],x0=[],phi0=0.5):
        
        if len(y0)>0 and len(x0)>0:
            liq = fluid.copy_change_x_and_conditions(
                fluid.T, fluid.P, None, x0, 'liquid')
            vap = fluid.copy_change_x_and_conditions(
                fluid.T, fluid.P, None, y0, 'gas')
            phi_liq = liq.evaluate_fugacity_coef()
            phi_vap = vap.evaluate_fugacity_coef()
            
            ki = phi_liq/phi_vap
        else:
            omega = array([sp.omega for sp in fluid.mixture.list_of_species])
            Pc = array(fluid.mixture.list_Pc)
            Tc = array(fluid.mixture.list_Tc)
            
            ki = Pc/fluid.P*exp(5.37*(1+omega)*(1-Tc/fluid.T))
        
        
        z = array(fluid.mixture.x)
        
        fun_phi = lambda phi: self.solve_phi(phi, z, ki)
        
        phi = fsolve(fun_phi, phi0)
        
        x = z/(1+phi*(ki-1))
        y = z*ki/(1+phi*(ki-1))
        
        while abs(phi0-phi)>1e-6:
            phi0 = phi
            
            liq = fluid.copy_change_x_and_conditions(
                fluid.T, fluid.P, None, x, 'liquid')
            vap = fluid.copy_change_x_and_conditions(
                fluid.T, fluid.P, None, y, 'gas')
            
            phi_liq = liq.evaluate_fugacity_coef()
            phi_vap = vap.evaluate_fugacity_coef()
            
            ki = phi_liq/phi_vap
            
            fun_phi = lambda phi: self.solve_phi(phi, z, ki)
            
            phi = fsolve(fun_phi, phi0)
            
            x = z/(1+phi*(ki-1))
            y = z*ki/(1+phi*(ki-1))
        
        if vap.mass_rho > liq.mass_rho:
            return 1-phi, liq, vap
        else:
            return phi, vap, liq
    
    def evaluate_flash_old(self, fluid, y0, x0, phi=.1):
        
        liq = fluid.copy_change_x_and_conditions(
            fluid.T, fluid.P, None, x0, 'liquid')
        vap = fluid.copy_change_x_and_conditions(
            fluid.T, fluid.P, None, y0, 'gas')
        phi_liq = liq.evaluate_fugacity_coef()
        phi_vap = vap.evaluate_fugacity_coef()
        #phi_vap = array([1.45]*len(fluid.mixture.x))
        # ki = [1.8]*19
        # ki[-3] = 0.002
        # ki[-1] = 0.002
        # ki = array(ki)
        # ki0 = ki
        ki = phi_liq/phi_vap
        z = array(fluid.mixture.x)
        it = 0  
        
        if phi > 0.5:
            step_phi = -0.001
        else:
            step_phi = 0.001
        
        x = z/(1+phi*(ki-1))
        y = z*ki/(1+phi*(ki-1))
        
        print(ki)
        
        fun_phi_init = 4

        fa = fun_phi_init

        fb = fun_phi_init
        phib = phi
        phia = phi
        #vap = fluid
        #liq = fluid
        fun_case = 4
        delta_phi = 4
        while fun_case > self.flash_tol and it < self.max_run and phi > 0 and delta_phi > self.flash_tol:
            x = x0
            y = y0
            liq = fluid.copy_change_x_and_conditions(
                fluid.T, fluid.P, None, x, 'liquid')
            vap = fluid.copy_change_x_and_conditions(
                fluid.T, fluid.P, None, y, 'gas')
            phi_liq = liq.evaluate_fugacity_coef()
            phi_vap = vap.evaluate_fugacity_coef()
            
            
            ki = phi_liq/phi_vap
            
            x = z/(1+phi*(ki-1))
            y = z*ki/(1+phi*(ki-1))
                    
                
            fun_phi_0 = 4.5
            while abs(fun_phi_0-sum(y)+sum(x)) > self.delta_fun:
                fun_phi_0 = sum(y)-sum(x)
                vap = fluid.copy_change_x_and_conditions(
                    fluid.T, fluid.P, None, y, 'gas')
                liq = fluid.copy_change_x_and_conditions(
                    fluid.T, fluid.P, None, x, 'liquid')
                phi_vap = vap.evaluate_fugacity_coef()
                phi_liq = liq.evaluate_fugacity_coef()
                ki = phi_liq/phi_vap
                x = z/(1+phi*(ki-1))
                y = z*ki/(1+phi*(ki-1))
            
            print(phi_liq[-3],x[-3],phi)
            print(phi_vap[-3],y[-3])
            
            f = sum(y)-sum(x)
            if sum(abs(array(y)-array(x))) > 0.005:
                if abs(sum(y)-sum(x)) > self.flash_tol:
                    if fb < 0 and fa > 0:
                        if f < 0:
                            fb = f
                            phib = phi
                        elif f > 0:
                            fa = f
                            phia = phi
                        delta_phi = phi + (fa-phia*(fa-fb)/(phia-phib))/((fa-fb)/(phia-phib))
                        phi = -(fa-phia*(fa-fb)/(phia-phib))/((fa-fb)/(phia-phib))
                    else:
                        if fun_phi_init > 0:
                            if f < 0:
                                fb = f
                                phib = phi
                            else:
                                fa = f
                                fb = f
                                phia = phi
                                phib = phi
                                phi += step_phi
                        if fun_phi_init < 0:
                            if f > 0:
                                fa = f
                                phia = phi
                            else:
                                fa = f
                                fb = f
                                phia = phi
                                phib = phi
                                phi += step_phi
                fun_case = abs(sum(y)-sum(x))
            else:
                phi += step_phi
                fun_case = 3
            it += 1
        
        if vap.mass_rho > liq.mass_rho:
            return 1-phi, liq, vap
        else:
            return phi, vap, liq
    
    def set_estimated_conditions(self,x0,y0,w0):
        self.x0 = x0 # Initial Estimative for the bubble T
        self.y0 = y0 # Initial Estimative for the dew P
        self.w0 = w0 # Initial Estimative for the dew P in the water line
    
    def build_phase_envelope(self, T_range, T_criteria, P_range, P_criteria, delta_T, delta_P):

        T_min = T_range[0]
        T_max = T_range[1]
        P_min = P_range[0]
        P_max = P_range[1]

        step_T_init = (T_max-T_min)/30
        step_P_init = (P_max-P_min)/30

        self.Porv = []
        self.Torv = []
        self.Zorvl = []
        self.Vorvl = []
        self.Zorvv = []
        self.Vorvv = []

        self.Pbub = []
        self.Tbub = []
        self.Zbubl = []
        self.Vbubl = []
        self.Zbubv = []
        self.Vbubv = []

        self.vapx = []
        self.liqx = []
        
        self.vapox = []
        self.liqox = []
        
        x0 = self.x0
        y0 = self.y0
        
        Pd = P_min
        print('|----------------------dew curve---------------------|')
        T0 = T_max
        while Pd < P_criteria:
            Td, vapo, liqo = self.dew_P(T_max+273.15, Pd, x0)
            if sum(abs(array(liqo.mixture.x)-array(self.fluid.mixture.x))) != 0:
                x0 = liqo.mixture.x
                self.Porv.append(Pd)
                print(Pd)
                self.Torv.append(Td-273.15)
                Pd += delta_P
                self.Zorvl.append(liqo.Z)
                self.Zorvv.append(vapo.Z)
                self.Vorvv.append(vapo.V)
                self.Vorvl.append(liqo.V)
                self.vapox.append(vapo.mixture.x)
                self.liqox.append(liqo.mixture.x)
            if len(self.Torv):
                T0 = T_max
            else:
                if min(self.Torv) == self.Torv[-1]:
                    T0 = min(self.Torv) - 25
                else:    
                    T0 = min(self.Torv)
            
        Tb = T_criteria
        print('|--------------------bubble curve--------------------|') 
        P0 = 10000
        while Tb >= T_min:
            print(Tb)
            Pb, vapb, liqb = self.bubble_T(Tb+273.15, P0, y0)
            #y0 = vapb.mixture.x
            self.Pbub.append(Pb)
            self.Tbub.append(Tb)
            Tb -= delta_T
            self.Zbubl.append(liqb.Z)
            self.Zbubv.append(vapb.Z)
            self.Vbubv.append(vapb.V)
            self.Vbubl.append(liqb.V)
            self.vapx.append(vapb.mixture.x)
            self.liqx.append(liqb.mixture.x)
            #P0 = Pb - 1000
        self.Zbubl.reverse()
        self.Zbubv.reverse()
        self.Pbub.reverse()
        self.Tbub.reverse()
        self.vapx.reverse()
        self.liqx.reverse()
        self.Vbubl.reverse()
        self.Vbubv.reverse()
        
        bub_array = array([self.Tbub, self.Pbub, array(self.Vbubl)-array(self.Vbubv)])
        orv_array = array([self.Torv, self.Porv, array(self.Vorvl)-array(self.Vorvv)])
        
        # print('|-------------------flash completude-----------------|')
        # if self.Torv[-1]-self.Tbub[-1]>=delta_T:
        #     Tt = self.Tbub[-1]+1
        #     while Tt <= self.Torv[-1]:
        #         if Tt == self.Torv[-1]:
        #             self.Tbub.append(self.Torv[-1])
        #             self.Pbub.append(self.Porv[-1])
        #             self.Zbubl.append(self.Zorvl[-1])
        #             self.Zbubv.append(self.Zorvv[-1])
        #             self.Vbubv.append(self.Vorvv[-1])
        #             self.Vbubl.append(self.Vorvl[-1])
        #             self.vapx.append(vapb.mixture.x)
        #             self.liqx.append(liqb.mixture.x)
        #         else:
        #             phi = 0.5
        #             P0 = self.Porv[-1]
        #             print(Tt)
        #             while 1-phi>0.005:
        #                 ptest = self.fluid.copy_change_conditions(Tt+273.15,P0,None,'liquid')
        #                 phi, phase1, phase2 = self.evaluate_flash(ptest,liqb.mixture.x,vapb.mixture.x,1)
        #                 print([phi, P0, Tt])
        #                 if isnan(phi):
        #                     phi = 0.5
        #                 if phase1.mixture.MM_m/phase1.V<phase2.mixture.MM_m/phase2.V:
        #                     vapc = phase1
        #                     liqc = phase2
        #                 else:
        #                     vapc = phase2
        #                     liqc = phase1    
        #                 P0 += 10
        #             self.Pbub.append(P0)
        #             self.Tbub.append(Tt)
        #             self.Zbubl.append(liqc.Z)
        #             self.Zbubv.append(vapc.Z)
        #             self.Vbubv.append(vapc.V)
        #             self.Vbubl.append(liqc.V)
        #             self.vapx.append(vapc.mixture.x)
        #             self.liqx.append(liqc.mixture.x)
                
        #         Tt += self.step_flash
        
        T = array(self.Tbub + self.Torv)
        P = array(self.Pbub + self.Porv)
        
        (T0, V0) = self.critical_init_value
        
        self.fluid.critical_case(T0, V0)
        
        self.critcal_isotherm_line_P = []
        self.critcal_isotherm_line_V = []
        
        if self.critical_evaluation:
            print('|------------evaluating critical isotherm------------|') 
            P_it = self.fluid.crit_gas.P
            while P_it < P_max*1.5:
                iso_therm = self.fluid.copy_change_conditions(self.fluid.crit_gas.T,P_it,None,'gas')
                self.critcal_isotherm_line_P.append(iso_therm.P)
                self.critcal_isotherm_line_V.append(iso_therm.V)
                P_it += step_P_init
        
        V_ref = self.fluid.mixture.MM_m/250
        self.isochoric_line_T = []
        self.isochoric_line_P = []
        if self.isochoric_evaluation:
            print('|------------evaluating isochoric line---------------|') 
            V = array(self.Vbubl+self.Vorvv)
            lut_T_ref = interpolant('LUT','bspline',[self.Vbubl],self.Tbub)
            lut_P_ref = interpolant('LUT','bspline',[self.Vbubl],self.Pbub)
            T_ref = lut_T_ref(V_ref).__float__()
            P_ref = lut_P_ref(V_ref).__float__()
            P_line_r = P_ref
            while P_line_r < self.P_max_isochoric:
                line_ref = self.fluid.copy_change_conditions(None,P_line_r,V_ref,'gas')
                if line_ref.T-273.15 < max(self.Torv):
                    self.isochoric_line_P.append(line_ref.P)
                    self.isochoric_line_T.append(line_ref.T-273.15)
                P_line_r += step_P_init
            
        
       
        self.water_line_T = []
        self.water_line_P = []
        self.water_vapo = []
        self.water_liqo = []
        if self.water_evaluation:
            print('|-------------evaluating water formation-------------|') 
            self.lut_T_ref = interpolant('LUT','bspline',[self.Tbub],self.Pbub)
            Pref = 100
            P_line = self.P_max_water
            while Pref <= P_line:
                Td, vapo, liqo = self.dew_P(T_min+273.15, P_line, self.w0)
                Pref = self.lut_T_ref(liqo.T-273.15).__float__()
                P_line -= step_P_init/2
                
                self.water_vapo.append(vapo)
                self.water_liqo.append(liqo)
                if Pref < P_line:
                    self.water_line_P.append(liqo.P)
                    self.water_line_T.append(liqo.T-273.15)
            self.water_line_P.reverse()
            self.water_line_T.reverse()
            self.lut_T_sat = interpolant('LUT','bspline',[self.water_line_P],self.water_line_T)
    
    
    def build_phase_bubble(self, T_range, T_criteria, delta_T, delta_P):

        T_min = T_range[0]
        T_max = T_range[1]

        step_T_init = (T_max-T_min)/30
        
        self.Porv = []
        self.Torv = []
        self.Zorvl = []
        self.Vorvl = []
        self.Zorvv = []
        self.Vorvv = []

        self.Pbub = []
        self.Tbub = []
        self.Zbubl = []
        self.Vbubl = []
        self.Zbubv = []
        self.Vbubv = []

        self.vapx = []
        self.liqx = []
        
        self.vapox = []
        self.liqox = []
        
        x0 = self.x0
        y0 = self.y0
        
        Tb = T_criteria
        print('|--------------------bubble curve--------------------|') 
        P0 = 20000
        while Tb >= T_min:
            Pb, vapb, liqb = self.bubble_T(Tb+273.15, P0, y0)
            #y0 = vapb.mixture.x
            print(Tb)
            self.Pbub.append(Pb)
            self.Tbub.append(Tb)
            Tb -= delta_T
            self.Zbubl.append(liqb.Z)
            self.Zbubv.append(vapb.Z)
            self.Vbubv.append(vapb.V)
            self.Vbubl.append(liqb.V)
            self.vapx.append(vapb.mixture.x)
            self.liqx.append(liqb.mixture.x)
            #P0 = Pb - 1000
        self.Zbubl.reverse()
        self.Zbubv.reverse()
        self.Pbub.reverse()
        self.Tbub.reverse()
        self.vapx.reverse()
        self.liqx.reverse()
        self.Vbubl.reverse()
        self.Vbubv.reverse()
        
        bub_array = array([self.Tbub, self.Pbub, array(self.Vbubl)-array(self.Vbubv)])
        
        T = array(self.Tbub)
        P = array(self.Pbub)
        
        self.critcal_isotherm_line_P = []
        self.critcal_isotherm_line_V = []
        
        if self.critical_evaluation:
            (T0, V0) = self.critical_init_value
            
            self.fluid.critical_case(T0, V0)
            print('|------------evaluating critical isotherm------------|') 
            P_it = self.fluid.crit_gas.P
            while P_it < max(self.Pbub)*1.5:
                iso_therm = self.fluid.copy_change_conditions(self.fluid.crit_gas.T,P_it,None,'gas')
                self.critcal_isotherm_line_P.append(iso_therm.P)
                self.critcal_isotherm_line_V.append(iso_therm.V)
                P_it += delta_P
        
        V_ref = self.fluid.mixture.MM_m/300
        self.isochoric_line_T = []
        self.isochoric_line_P = []
        if self.isochoric_evaluation:
            print('|------------evaluating isochoric line---------------|')
            V = array(self.Vbubl)
            lut_T_ref = interpolant('LUT','bspline',[self.Vbubl],self.Tbub)
            lut_P_ref = interpolant('LUT','bspline',[self.Vbubl],self.Pbub)
            T_ref = lut_T_ref(V_ref).__float__()
            P_ref = lut_P_ref(V_ref).__float__()
            P_line_r = P_ref
            T_ref = 0
            while P_line_r < self.P_max_isochoric and T_ref < max(self.Tbub):
                line_ref = self.fluid.copy_change_conditions(None,P_line_r,V_ref,'gas')
                if line_ref.T-273.15 < max(self.Tbub):
                    self.isochoric_line_P.append(line_ref.P)
                    self.isochoric_line_T.append(line_ref.T-273.15)
                P_line_r += delta_P
                T_ref = max(self.isochoric_line_T)
            
       
        self.water_line_T = []
        self.water_line_P = []
        if self.water_evaluation:
            print('|-------------evaluating water formation-------------|') 
            self.lut_T_ref = interpolant('LUT','bspline',[self.Tbub],self.Pbub)
            Pref = 250
            P_line = self.P_max_water
            while Pref <= P_line:
                Td, vapo, liqo = self.dew_P(T_min+273.15, P_line, self.w0)
                Pref = self.lut_T_ref(liqo.T-273.15).__float__()
                P_line -= delta_P*2
                if Pref < P_line:
                    self.water_line_P.append(liqo.P)
                    self.water_line_T.append(liqo.T-273.15)
                print(Td)
            self.water_line_P.reverse()
            self.water_line_T.reverse()
            self.lut_T_sat = interpolant('LUT','bspline',[self.water_line_P],self.water_line_T)
        
        
    
    def plot_envelope(self, Pexp=[], Texp=[], Vexp=[],pt_points=[],hydrate_points = []):

        fig1 = figure(dpi = 150)
        ax = fig1.add_subplot(1, 1, 1)
        p1, = plot(self.Torv, [P/100 for P in self.Porv], 'r')
        p2, = plot(self.Tbub, [P/100 for P in self.Pbub], 'b')
        p3, = plot(Texp, Pexp, 'ko', fillstyle='none')
        p4, = plot(self.fluid.T-273.15, self.fluid.P/100, 'kx', markersize=10, markeredgewidth=2)
        p5, = plot(self.fluid.crit_gas.T-273.15, self.fluid.crit_gas.P/100, 'ko', markersize=10)
        if len(self.isochoric_line_T)>0:
            p6, = plot(array(self.isochoric_line_T), array(self.isochoric_line_P)/100,'--g')
        if len(self.water_line_P)>0:
            p7, = plot(array(self.water_line_T), array(self.water_line_P)/100,'--b')
        if len(pt_points)>0:
            p8, = plot(pt_points[0], pt_points[1], 'ko--')
        if len(hydrate_points)>0:
            p8, = plot(hydrate_points[0], [p/100 for p in hydrate_points[1]], 'm--')
        # p3, = plot(Tdexp, Pdexp, 'ro', fillstyle = 'none')
        config_plot(ax)
        ax.set_ylabel('Pressure /bar', fontsize=12)
        ax.set_xlabel('Temperature /°C', fontsize=12)
        fig1.tight_layout()
        fig1.show()

        Vmax = 0.6
        fig2 = figure(dpi = 150)
        ax = fig2.add_subplot(1, 1, 1)
        Vgorv_l = [V for V in self.Vorvv if V < Vmax]
        Vexp_l  = [V for V in Vexp if V < Vmax]
        Pgorv_l = [self.Porv[i] /
                    100 for i in range(len(self.Porv)) if self.Vorvv[i] < Vmax]
        Pexp_l  = [Pexp[i] for i in range(len(Pexp)) if Vexp[i] < Vmax]
        p1, = plot(Vgorv_l, [P for P in Pgorv_l], 'r')
        p2, = plot(self.Vbubl, [P/100 for P in self.Pbub], 'b')
        p3, = plot(Vexp_l, Pexp_l, 'ko', fillstyle='none')
        p4, = plot(self.fluid.V, self.fluid.P/100, 'kx', markersize=10, markeredgewidth=2)
        p5, = plot(self.fluid.crit_gas.V, self.fluid.crit_gas.P/100, 'ko', markersize=10)
        if len(self.critcal_isotherm_line_V)>0:
            p6, = plot(array(self.critcal_isotherm_line_V), array(self.critcal_isotherm_line_P)/100,'k')
        fig2.tight_layout()
        fig2.show()
        config_plot(ax)
        ax.set_xlabel('Molar Volume /(m³/kmol)', fontsize=12)
        ax.set_ylabel('Pressure /bar', fontsize=12)
        fig2.tight_layout()
        fig2.show()
    
    def plot_PT(self, Pexp=[], Texp=[], Vexp=[],pt_points=[],hydrate_points = []):

        fig1 = figure(dpi = 150)
        ax = fig1.add_subplot(1, 1, 1)
        p1, = plot(self.Torv, [P/100 for P in self.Porv], 'r')
        p2, = plot(self.Tbub, [P/100 for P in self.Pbub], 'b')
        p3, = plot(Texp, Pexp, 'ko', fillstyle='none')
        p5, = plot(self.fluid.crit_gas.T-273.15, self.fluid.crit_gas.P/100, 'ko', markersize=10)
        if len(self.isochoric_line_T)>0:
            p6, = plot(array(self.isochoric_line_T), array(self.isochoric_line_P)/100,'--g')
        if len(pt_points)>0:
            p8, = plot(pt_points[0], pt_points[1], 'ko--')
        if len(hydrate_points)>0:
            p8, = plot(hydrate_points[0], [p/100 for p in hydrate_points[1]], 'r--')
        # p3, = plot(Tdexp, Pdexp, 'ro', fillstyle = 'none')
        config_plot(ax)
        ax.set_ylabel('Pressure /bar', fontsize=20)
        ax.set_xlabel('Temperature /°C', fontsize=20)
        ax.set_ylim(90,650)
        fig1.tight_layout()
        fig1.show()
    
    def plot_many_envelope(self,list_of_solvers,list_of_hidrates,list_colors):
        
        fig1 = figure(dpi = 150, figsize=(4,3))
        ax = fig1.add_subplot(1, 1, 1)
        for i, solver in enumerate(list_of_solvers):
            p1, = plot(solver.Torv, [P/100 for P in solver.Porv], list_colors[i])
            p2, = plot(solver.Tbub, [P/100 for P in solver.Pbub], list_colors[i])
            #p3, = plot(solver.fluid.T-273.15, solver.fluid.P/100, list_colors[i]+'x', markersize=10, markeredgewidth=2)
            #p4, = plot(solver.fluid.crit_gas.T-273.15, solver.fluid.crit_gas.P/100, list_colors[i]+'o', markersize=10)
            if len(solver.isochoric_line_T)>0:
                p6, = plot(array(solver.isochoric_line_T), array(solver.isochoric_line_P)/100,'--'+list_colors[i])
            
            p7, = plot(array(list_of_hidrates[i].T_hydrate), array(list_of_hidrates[i].P_hydrate)/100,'--'+list_colors[i])
        config_plot(ax)
        ax.set_ylabel('Pressure /bar', fontsize=20)
        ax.set_xlabel('Temperature /°C', fontsize=20)
        fig1.tight_layout()
        fig1.show()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
