# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 01:30:58 2024

@author: Rodrigo Meira
"""
from eos_database import *
from casadi import *
from numpy import exp, log, array, roots, zeros, linalg
from math import isnan
from scipy.optimize import fsolve
from species_builder import R, Species, Mixture
from matplotlib.pyplot import plot, figure
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter
import scipy.linalg as la
from numpy.linalg import svd, det


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


class gc_eos_class:

    def __init__(self, mixture, T, P, V, u, w, Aij, delta_v, phase):
        self.mixture = mixture
        self.T = T
        self.P = P
        self.V = V
        self.Aij = Aij
        self.phase = phase
        self.csi = min([y.real for y in roots(
            [u*(w+u)-w, 3*(w+u), 3, -1]) if y.imag == 0])
        # self.csi = 0.2774
        self.u = u
        self.w = w
        self.delta_v = delta_v
        self.delta_vm = array(delta_v).dot(array(mixture.x).transpose())

        if [T, P, V].count(None) > 1:
            print("It is necessary 2 values in T, P or V")

        self.pr_parameters()

        if self.T == None:
            if self.phase == 'liquid':
                self.Veos = self.V-self.delta_vm
                def t_solver(Z): return array(self.evaluate_eos_T_liq(Z))
                res = fsolve(t_solver, 0.25)
                Z_eos = res[0]
                self.T = self.P*(self.Veos)/R/Z_eos
                self.Z = self.P*self.V/R/self.T
            elif self.phase == 'gas':
                self.Veos = self.V-self.delta_vm
                def t_solver(Z): return array(self.evaluate_eos_T_gas(Z))
                res = fsolve(t_solver, 0.85)
                Z_eos = res[0]
                self.T = self.P*(self.Veos)/R/Z_eos
                self.Z = self.P*self.V/R/self.T
        elif self.P == None:
            self.Veos = self.V-self.delta_vm
            self.P = self.evaluate_eos_P()
            self.Z = self.P*self.V/self.T/R

        elif self.V == None:
            if self.phase == 'liquid':
                self.Zeos = min(
                    [y.real for y in self.evaluate_z_roots() if y.imag == 0 and y.real > 0])
                self.Veos = self.Zeos*self.T*R/self.P
                self.V = self.Veos+self.delta_vm
                self.Z = self.P*self.V/R/self.T
            elif self.phase == 'gas':
                self.Zeos = max(
                    [y.real for y in self.evaluate_z_roots() if y.imag == 0 and y.real > 0])
                self.Veos = self.Zeos*self.T*R/self.P
                self.V = self.Veos+self.delta_vm
                self.Z = self.P*self.V/R/self.T
            # self.Z = res[0]

        #self.critical_parameters()

        self.rho = 1/self.V
    
    def critical_case_T(self,T0,V0,F,B,eval):
        crit_gas = self.copy_change_conditions(T0[0], None, V0, 'gas')
        
        x = array(crit_gas.mixture.x)
        
        A = crit_gas.Q.dot(x)/crit_gas.a_m
        mat_Q = zeros([len(x), len(x)])
        for i in range(len(x)):
            for j in range(len(x)):
                delta = 0
                if i==j:
                    delta = 1
                par_1 = 1/sum(x)*(delta/x[i]+(B[i]+B[j])*F[1]+B[i]*B[j]*F[1]**2)
                par_2 = (B[i]*B[j]-A[i]*B[j]-A[j]*B[i])*F[6];
                par_3 = crit_gas.a_m/(self.b_m*R*T0*sum(x))*(B[i]*B[j]*F[3]-crit_gas.Q[i][j]/crit_gas.a_m*F[5]+par_2);
                mat_Q[i][j] = (par_1+par_3)*T0/100;
        if eval:
            return mat_Q
        else:
            return array(det(mat_Q))
    
    def critical_case_V(self, T0, V0,eval):
        
        V0 = V0[0]
        
        D1 = self.u/2 + (self.u**2-4*self.w)**0.5/2
        D2 = self.u/2 - (self.u**2-4*self.w)**0.5/2
        crit_gas = self.copy_change_conditions(T0, None, V0, 'gas')
        K = crit_gas.V/self.b_m
        F1 = 1/(K-1)
        F2 = 2/(D1-D2)*(D1/(K+D1)-D2/(K+D2))
        F3 = 1/(D1-D2)*(D1**2/(K+D1)**2-D2**2/(K+D2)**2)
        F4 = 1/(D1-D2)*(D1**3/(K+D1)**3-D2**3/(K+D2)**3)
        F5 = 2/(D1-D2)*log((K+D1)/(K+D2))
        F6 = F2 - F5
        
        F = [0, F1, F2, F3, F4, F5, F6]
        
        B = self.b/self.b_m
        
        t_solver = lambda T: self.critical_case_T(T, V0, F, B, False)
        res = fsolve(t_solver, x0 = T0)
        T0 = res[0]
        x = array(crit_gas.mixture.x)
        crit_gas = self.copy_change_conditions(T0, None, V0, 'gas')
        
        mat_Q = self.critical_case_T(array([T0]), V0, F, B, True)
        
        U, S, Vh = np.linalg.svd(mat_Q, full_matrices=True)
        
        delta_y0 = U[:,-1]
        if abs(min(delta_y0))>max(delta_y0):
            delta_y0 = delta_y0*-1
        #print(U)
        delta_y0 = delta_y0/sum([i**2 for i in delta_y0])**0.5
        
        N_ = sum(delta_y0)
        B_ = sum(B*delta_y0)
        A = crit_gas.Q.dot(x)/crit_gas.a_m
        A_ = sum(A*delta_y0)
        
        a_ = delta_y0.transpose().dot(crit_gas.Q.dot(delta_y0))/crit_gas.a_m
        
        t1 = -sum(delta_y0**3/x**2)+3*N_*(B_*F[1])**2+2*(B_*F[1])**3
        t2 = 3*B_**2*(2*A_-B_)*(F[3]+F[6])-2*B_**3*F[4]-3*B_*a_*F[6]
        
        Cf = t1+crit_gas.a_m/(crit_gas.b_m*R*crit_gas.T)*t2
        Cf = Cf*(self.V-self.b_m)/self.b_m**2

        if eval:
            return T0
        else:
            return Cf
    
    def critical_case_fehbet(self, T0, V0, delta_y0,eval):
        
        delta_y0 = array(delta_y0)
        
        D1 = self.u/2 + (self.u**2-4*self.w)**0.5/2
        D2 = self.u/2 - (self.u**2-4*self.w)**0.5/2
        crit_gas = self.copy_change_conditions(T0, None, V0, 'gas')
        K = crit_gas.Veos/self.b_m
        F1 = 1/(K-1)
        F2 = 2/(D1-D2)*(D1/(K+D1)-D2/(K+D2))
        F3 = 1/(D1-D2)*(D1**2/(K+D1)**2-D2**2/(K+D2)**2)
        F4 = 1/(D1-D2)*(D1**3/(K+D1)**3-D2**3/(K+D2)**3)
        F5 = 2/(D1-D2)*log((K+D1)/(K+D2))
        F6 = F2 - F5
        
        F = [0, F1, F2, F3, F4, F5, F6]
        
        B = self.b/self.b_m
        
        #t_solver = lambda T: self.critical_case_T(T, V0, F, B, False)
        #res = fsolve(t_solver, x0 = T0)
        #T0 = res[0]
        x = array(crit_gas.mixture.x)
        crit_gas = self.copy_change_conditions(T0, None, V0, 'gas')
        
        mat_Q = self.critical_case_T(array([T0]), V0, F, B, True)
        
        N_ = sum(delta_y0)
        B_ = sum(B*delta_y0)
        A = crit_gas.Q.dot(x)/crit_gas.a_m
        A_ = sum(A*delta_y0)
        
        a_ = delta_y0.transpose().dot(crit_gas.Q.dot(delta_y0))/crit_gas.a_m
        
        
        #for i in range(len(delta_y0)):
        #    F = R*crit_gas.T*(delta_y0[i]/x[i]+F[1]*self.b_bar[i])
        #    d2fdxdx.append(F)
            
        t1 = -sum(delta_y0**3/x**2)+3*N_*(B_*F[1])**2+2*(B_*F[1])**3
        t2 = 3*B_**2*(2*A_-B_)*(F[3]+F[6])-2*B_**3*F[4]-3*B_*a_*F[6]
        
        Cf = R*crit_gas.T*t1+crit_gas.a_m/(crit_gas.b_m)*t2
        Cf = [Cf]
        
        d2fdxdx = mat_Q.dot(delta_y0)
        
        Cf.append(1 - sum(delta_y0**2))
        
        for i in d2fdxdx:
            Cf.append(i)
        
        return Cf
    
    
    def critical_case(self, T0, V0, delta_y0):
        
        V_solver = lambda V: self.critical_case_V(T0, V, False)
        
        res = fsolve(V_solver, x0 = V0)
        
        Tc = self.critical_case_V(T0, res, True)
        
        Vc = res[0]
        
        self.crit_gas = self.copy_change_conditions(Tc, None, Vc, 'gas')
    
    def evaluate_helmhotz_jacobian(self):

        Qc = self.eye_x+2*self.bij/(self.Veos-self.b_m) +\
            self.dbeta_ij/(self.Veos-self.b_m) +\
            self.a_m*self.dbeta_ij*self.Veos*2*(self.Veos-self.b_m)/R/self.T/self.b_m/(self.Veos**2+self.Veos*self.b_m-self.b_m**2)**2 +\
            self.B1c/R/self.T/self.b_m**2*self.Veos/(self.Veos**2+self.Veos*self.b_m-self.b_m**2) +\
            self.B2c/R/self.T/self.b_m**2/2**1.5 * \
            log((self.Veos+self.b_m*(1+2**.5))/(self.Veos+self.b_m*(1-2**.5)))

        return Qc

    def evaluate_eos_T_gas(self, Z):

        self.T = self.P*self.V/Z[0]/R
        self.evaluate_par_a()

        beta = self.b_m*self.P/self.T/R

        q = self.a_m/self.b_m/self.T/R

        return Z - 1 - beta + q*beta*(Z-beta)/(Z**2 + 2*Z*beta - beta**2)

    def evaluate_eos_T_liq(self, Z):

        self.T = self.P*(self.Veos)/Z[0]/R
        self.evaluate_par_a()

        beta = self.b_m*self.P/self.T/R

        q = self.a_m/self.b_m/self.T/R

        return Z - beta + (Z**2 + 2*Z*beta - beta**2)*(1 + beta - Z)/q/beta

    def evaluate_eos_P(self):

        self.evaluate_par_a()

        V_eos = self.V + self.delta_vm

        return self.T*R/(V_eos-self.b_m) - self.a_m/(V_eos**2 + 2*V_eos*self.b_m - self.b_m**2)
    
    def evaluate_der_eos_P(self):
        
        self.evaluate_der_a()
        
        self.dPdT = R/(self.Veos-self.b_m) - \
                    self.der_a/(self.Veos**2 + 2*self.Veos*self.b_m - self.b_m**2)
        
        self.dPdV = -R*self.T/(self.Veos-self.b_m)**2 + \
                    self.a_m/(self.Veos**2 + 2*self.Veos*self.b_m - self.b_m**2)**2*\
                    2*(self.Veos + self.b_m)
        
    def evaluate_der_rho(self):
        
        self.evaluate_der_eos_P()
        
        dPdrho = -self.Veos**2/self.mixture.MM_m*self.dPdV
        
        self.drhodP = 1/dPdrho
        
        self.drhodT = -self.dPdT/dPdrho
        
    def evaluate_eos_V_gas(self, Z):

        self.V = self.T*Z*R/self.P + self.delta_vm
        self.evaluate_par_a()

        beta = self.b_m*self.P/self.T/R
        q = self.a_m/self.b_m/self.T/R

        return Z - 1 - beta + q*beta*(Z-beta)/(Z**2 + 2*Z*beta - beta**2)

    def evaluate_eos_V_liq(self, Z):

        self.V = self.T*Z*R/self.P
        self.evaluate_par_a()

        beta = self.b_m*self.P/self.T/R
        q = self.a_m/self.b_m/self.T/R

        return Z - beta - (Z**2 + 2*Z*beta - beta**2)*(1 + beta - Z)/q/beta

    def pr_parameters(self):

        self.a_cons = (3+(self.u-self.w)*self.csi**2) / \
            (3+(self.u-1)*self.csi)+self.u*self.csi

        self.b = [self.csi*sp.Vc for sp in self.mixture.list_of_species]

        self.a_c = [self.a_cons*R*sp.Tc *
                    sp.Vc for sp in self.mixture.list_of_species]

        self.b = [R*sp.Tc*0.07780/sp.Pc for sp in self.mixture.list_of_species]

        self.a_c = [R**2*sp.Tc**2*0.45724 /
                    sp.Pc for sp in self.mixture.list_of_species]

        self.kappa = [sp.kappa for sp in self.mixture.list_of_species]

    def evaluate_a(self):
        self.Tr = [self.T/sp.Tc for sp in self.mixture.list_of_species]
        self.a = [self.a_c[i]*(1+self.kappa[i]*(1-self.Tr[i]**0.5))**2 \
                  for i in range(len(self.Tr))]

    def evaluate_der_a(self):
        
        A_c = array([[(i*j)**0.5 for i in self.a_c] for j in self.a_c])
        
        sqrt_alfa = [(1+self.kappa[i]*(1-self.Tr[i]**0.5)) \
                     for i in range(len(self.mixture.x))]
        
        der_sqrt_a = [-0.5**self.kappa[i]/(self.T*self.mixture.list_Tc[i])**0.5 \
                      for i in range(len(self.kappa))]
        
        der_A = [[sqrt_alfa[i]*der_sqrt_a[j]+sqrt_alfa[j]*der_sqrt_a[i] \
                  for i in range(len(der_sqrt_a))]
                 for j in range(len(der_sqrt_a))]
        
        Q = A_c*array(der_A)*(1-array(self.Aij))
        
        x = array(self.mixture.x)
        
        self.der_a = x.dot(Q.dot(x.transpose()))

    def evaluate_der_sec_a(self):
        
        A_c = array([[(i*j)**0.5 for i in self.a_c] for j in self.a_c])
        
        der_sec_A = [[0.25*self.T**-1.5*((1+self.kappa[i])*self.kappa[j]/self.mixture.list_Tc[j]**0.5 +
                                         (1+self.kappa[j])*self.kappa[i]/self.mixture.list_Tc[i]**0.5) \
                      for i in range(len(self.mixture.x))]
                     for j in range(len(self.mixture.x))]

        Q = A_c*array(der_sec_A)*(1-array(self.Aij))

        x = array(self.mixture.x)

        self.der_sec_a = x.dot(Q.dot(x.transpose()))

    def evaluate_par_a(self):

        self.evaluate_a()

        self.Q = array([[(self.a[i]*self.a[j])**0.5*(1-self.Aij[i][j]) for i in range(len(self.a))]
                        for j in range(len(self.a))])
        x = array(self.mixture.x)
        self.a_m = x.dot(self.Q.dot(x.transpose()))
        self.b_m = sum(array(self.b)*x)

        self.b_bar = array(self.b)

        self.a_bar = 2*self.Q.dot(x.transpose())-self.a_m

    def evaluate_fugacity_coef(self):

        beta = self.b_m*self.P/self.T/R
        q = self.a_m/self.b_m/self.T/R

        ln_phi_c_i = self.b_bar/self.b_m*(self.Zeos-1) - log(self.Zeos-beta) \
            - q*(1 + self.a_bar/self.a_m - array(self.b_bar)/self.b_m) \
            / 2**1.5*log((self.Zeos+beta*(1+2**.5))/(self.Zeos+beta*(1-2**.5)))

        return exp(ln_phi_c_i)

    def evaluate_der_fugacity_coef(self):

        beta = self.b_m*self.P/self.T/R
        q = self.a_m/self.b_m/self.T/R

        d_ln_phi_c_i_dz = array(self.b_bar)/self.b_m - 1/(self.Z-beta) \
            + q*(1 + self.a_bar/self.a_m - array(self.b_bar)/self.b_m) \
            / (self.Z**2+2*beta*self.Z-beta**2)

        return d_ln_phi_c_i_dz
    
    def ci_ideal(self):
        
        self.cpi = self.mixture.evaluate_cp_ig(self.T)
        
        self.cvi = self.cpi-R
    
    def ci_real(self):
            
        self.evaluate_der_sec_a()
        
        self.evaluate_der_eos_P()
        
        self.ci_ideal()
        
        # Cv_real = Cv_ideal + T*\int_V^Inf(d2PdT2)dV
        
        self.Cvt = self.cvi + self.T*self.der_sec_a/(2*2**0.5*self.b_m)*log((self.Veos+self.b_m*(1+2**0.5))/(self.Veos+self.b_m*(1-2**0.5)))
        
        # Cp_real = Cv_real - T*(dPdT)^2/dPdV
        
        self.Cpt = self.Cvt - self.T*(self.dPdT**2)*(self.dPdV**-1)
    
    def h_gas(self):
        
        hi = self.mixture.evaluate_enthalpy_ig(self.T)
        
        self.evaluate_der_a()
        
        hr = self.P*self.Veos - R*self.T - (self.a_m - self.der_a*self.T)/(2*2**0.5*self.b_m)\
                                           *log((self.Veos + self.b_m*(1+2**0.5))/
                                                (self.Veos + self.b_m*(1-2**0.5)))
        
        self.h = hr + hi;
    
    def s_gas(self):
        
        si = self.mixture.evaluate_entropy_ig(self.T)
        
        self.evaluate_der_a()
        
        sr = R*log((self.Veos-self.b_m)*self.P/R/self.T) - self.der_a*(1/(2*self.b_m*2**0.5))\
                                                       *log((self.Veos + self.b_m*(1+2**0.5))\
                                                           /(self.Veos + self.b_m*(1-2**0.5)))
        
        print(sr)
        print(si)
        
        self.s = sr + si - R*log(self.P)
        
    def sound_speed(self):
        
        self.ci_real()
        
        self.evaluate_der_eos_P()
        
        rhoeos = 1/self.Veos;
        
        c2 = -self.Veos**2*self.dPdV + self.T*self.Veos**2*(self.dPdT)**2/(self.Cvt);

        return (c2*1000/self.mixture.MM_m)**0.5;
    
    def evaluate_z_roots(self):

        self.evaluate_par_a()

        beta = self.b_m*self.P/self.T/R
        q = self.a_m/self.b_m/self.T/R

        return roots([1, (self.u-1)*beta-1, (q-self.u+(self.w-self.u)*beta)*beta, -(self.w*beta**3+self.w*beta**2+q*beta**2)])

    def copy_change_x_and_conditions(self, T, P, V, x, state):

        new_mixture = self.mixture.copy_and_change_composition(x)

        return gc_eos_class(new_mixture, T, P, V, self.u, self.w, self.Aij, self.delta_v, state)

    def copy_change_conditions(self, T, P, V, state):

        return gc_eos_class(self.mixture, T, P, V, self.u, self.w, self.Aij, self.delta_v, state)

    def bubble_T(self, Tb, P0, y0):

        tol = 1e-4
        P = P0

        vap = self.copy_change_x_and_conditions(Tb, P, None, y0, 'gas')
        phi_vap = vap.evaluate_fugacity_coef()
        
        liq = self.copy_change_conditions(Tb, P, None, 'liquid')
        fliq = liq.evaluate_fugacity_coef()*array(self.mixture.x)
        #phi_vap = array([1.25]*len(self.mixture.x))
        it = 0

        y = fliq/phi_vap

        fa = sum(y)-1

        fb = sum(y)-1

        while abs(sum(y)-1) > tol and it < 100000:
            liq = self.copy_change_conditions(Tb, P, None, 'liquid')
            fliq = liq.evaluate_fugacity_coef()*array(self.mixture.x)
            sum_y_0 = 4.5
            while abs(sum_y_0-sum(y)) > 1e-5:
                sum_y_0 = sum(y)
                vap = self.copy_change_x_and_conditions(Tb, P, None, y, 'gas')
                phi_vap = vap.evaluate_fugacity_coef()
                y = fliq/phi_vap

            f = sum(y)-1

            if abs(sum(y)-1) > tol:

                if fb < 0 and fa > 0:
                    if sum(y)-1 < 0:
                        fb = sum(y)-1
                        Pb = P
                    elif sum(y)-1 > 0:
                        fa = sum(y)-1
                        Pa = P

                    P = -(fa-Pa*(fa-fb)/(Pa-Pb))/((fa-fb)/(Pa-Pb))

                else:

                    if sum(y)-1 < 0:
                        fb = sum(y)-1
                        Pb = P
                    else:
                        fa = sum(y)-1
                        fb = sum(y)-1
                        Pa = P
                        Pb = P
                        P += 100
            

            it += 1
        print('|----------------------------------------------------|')
        return P, vap, liq

    def test_bubble_T(self, Tb, P0, y0):

        tol = 1e-6
        P = P0

        vap = self.copy_change_x_and_conditions(Tb, P, None, y0, 'gas')
        phi_vap = vap.evaluate_fugacity_coef()
        liq = self.copy_change_conditions(Tb, P, None, 'liquid')
        phi_liq = liq.evaluate_fugacity_coef()
        #phi_vap = array([1.45]*len(self.mixture.x))
        it = 0

        y = phi_liq*array(liq.mixture.x)/phi_vap
        # xco2 = [y[3]]
        sum_y_0 = 4.5
        while abs(sum_y_0-sum(y)) > 1e-6:
            sum_y_0 = sum(y)
            vap = self.copy_change_x_and_conditions(Tb, P, None, y, 'gas')
            phi_vap = vap.evaluate_fugacity_coef()
            y = phi_liq*array(liq.mixture.x)/phi_vap

        viable = True
        if isnan(phi_vap[0]):
            viable = False

        return sum(y), y, vap, liq, viable

    def dew_P(self, T0, Pd, x0):

        phi_liq = array([0.75]*len(self.mixture.x))
        tol = 1e-5
        T = T0
        liq = self.copy_change_x_and_conditions(T, Pd, None, x0, 'liquid')
        phi_liq = liq.evaluate_fugacity_coef()
        vap = self.copy_change_conditions(T, Pd, None, 'gas')
        fvap = vap.evaluate_fugacity_coef()*array(self.mixture.x)

        it = 0
        x = fvap/phi_liq

        sum_x_0 = sum(x)

        fa = sum(x)-1
        fb = sum(x)-1
        while abs(sum(x)-1) > tol and it < 10000:

            vap = self.copy_change_conditions(T, Pd, None, 'gas')
            fvap = vap.evaluate_fugacity_coef()*array(self.mixture.x)
            sum_x_0 = 4.5
            while abs(sum_x_0-sum(x)) > 1e-6:
                sum_x_0 = sum(x)
                liq = self.copy_change_x_and_conditions(
                    T, Pd, None, x, 'liquid')
                phi_liq = liq.evaluate_fugacity_coef()
                x = fvap/phi_liq

            f = sum(x)-1
            if abs(sum(x)-1) > tol:

                if fb < 0 and fa > 0:
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

    def test_dew_P(self, T0, Pd, y0):

        T = T0

        vap = self.copy_change_conditions(T, Pd, None, 'gas')
        phi_vap = vap.evaluate_fugacity_coef()
        liq = self.copy_change_x_and_conditions(T, Pd, None, y0, 'liquid')
        phi_liq = liq.evaluate_fugacity_coef()
        # phi_liq = array([.85]*len(self.mixture.x))

        x = phi_vap*array(self.mixture.x)/phi_liq

        sum_x_0 = 4.5
        while abs(sum_x_0-sum(x)) > 1e-6:
            sum_x_0 = sum(x)
            liq = self.copy_change_x_and_conditions(T, Pd, None, x, 'liquid')
            phi_liq = liq.evaluate_fugacity_coef()
            x = phi_vap*array(self.mixture.x)/phi_liq

        viable = True
        if isnan(phi_liq[0]):
            viable = False

        return sum_x_0, x, vap, liq, viable


    def evaluate_flash(self, y0, x0):

        tol = 1e-6
        liq = self.copy_change_x_and_conditions(
            self.T, self.P, None, x0, 'liquid')
        phi_liq = liq.evaluate_fugacity_coef()
        phi_vap = array([1.45]*len(self.mixture.x))

        ki = phi_liq/phi_vap
        z = array(self.mixture.x)

        it = 0

        phi = 1
        x = z/(1+phi*(ki-1))
        y = z*ki/(1+phi*(ki-1))

        fun_phi_init = sum(y)-sum(x)

        fa = fun_phi_init

        fb = fun_phi_init
        phib = 0.5
        vap = self
        liq = self
        while abs(sum(y)-sum(x)) > tol and it < 100000:

            fun_phi_0 = 4.5
            while abs(fun_phi_0-sum(y)+sum(x)) > 1e-6:
                fun_phi_0 = sum(y)-sum(x)
                vap = self.copy_change_x_and_conditions(
                    self.T, self.P, None, y, 'gas')
                liq = self.copy_change_x_and_conditions(
                    self.T, self.P, None, x, 'liquid')
                phi_vap = vap.evaluate_fugacity_coef()
                phi_liq = liq.evaluate_fugacity_coef()
                ki = phi_liq/phi_vap
                x = z/(1+phi*(ki-1))
                y = z*ki/(1+phi*(ki-1))
            f = sum(y)-sum(x)
            if abs(sum(y)-sum(x)) > tol:
                if fb < 0 and fa > 0:
                    if f < 0:
                        fb = f
                        phib = phi
                    elif f > 0:
                        fa = f
                        phia = phi
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
                            phi -= 0.01
                    if fun_phi_init < 0:
                        if f > 0:
                            fa = f
                            phia = phi
                        else:
                            fa = f
                            fb = f
                            phia = phi
                            phib = phi
                            phi -= 0.01
            it += 1

        return phi, vap, liq

    def build_phase_envelope(self, T_range, T_criteria, P_range, P_criteria, delta_T, delta_P,
                              y0, x0, critical_isothem = False, isovolumetrical_case = False):

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
        
        
        Pd = P_min
        while Pd < P_criteria:
            Td, vapo, liqo = self.dew_P(T_max+273.15, Pd, x0)
            #x0 = liqo.mixture.x
            self.Porv.append(Pd)
            self.Torv.append(Td-273.15)
            Pd += delta_P
            self.Zorvl.append(liqo.Z)
            self.Zorvv.append(vapo.Z)
            self.Vorvv.append(vapo.V)
            self.Vorvl.append(liqo.V)
            self.vapox.append(vapo.mixture.x)
            self.liqox.append(liqo.mixture.x)
        
        Tb = T_criteria
        print('cheguei!')
        while Tb >= T_min:
            Pb, vapb, liqb = self.bubble_T(Tb+273.15, P_min, y0)
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

        if self.Torv[-1]-self.Tbub[-1]>=delta_T:
            Tt = self.Tbub[-1]+1
            while Tt <= self.Torv[-1]:
                if Tt == self.Torv[-1]:
                    self.Tbub.append(self.Torv[-1])
                    self.Pbub.append(self.Porv[-1])
                    self.Zbubl.append(self.Zorvl[-1])
                    self.Zbubv.append(self.Zorvv[-1])
                    self.Vbubv.append(self.Vorvv[-1])
                    self.Vbubl.append(self.Vorvl[-1])
                    self.vapx.append(vapb.mixture.x)
                    self.liqx.append(liqb.mixture.x)
                else:
                    phi = 0.5
                    P0 = self.Porv[-1]
                    while 1-phi>0.0005:
                        ptest = self.copy_change_conditions(Tt+273.15,P0,None,'liquid')
                        phi, phase1, phase2 = ptest.evaluate_flash(liqb.mixture.x,vapb.mixture.x)
                        if isnan(phi):
                            phi = 0.5
                        if phase1.mixture.MM_m/phase1.V<phase2.mixture.MM_m/phase2.V:
                            vapc = phase1
                            liqc = phase2
                        else:
                            vapc = phase2
                            liqc = phase1    
                        P0 += 10
                    print(phi)
                    self.Pbub.append(P0)
                    self.Tbub.append(Tt)
                    self.Zbubl.append(liqc.Z)
                    self.Zbubv.append(vapc.Z)
                    self.Vbubv.append(vapc.V)
                    self.Vbubl.append(liqc.V)
                    self.vapx.append(vapc.mixture.x)
                    self.liqx.append(liqc.mixture.x)
                
                Tt += 1
        # else:
        #     self.Tbub.append(self.Torv[-1])
        #     self.Pbub.append(self.Porv[-1])
        #     self.Zbubl.append(self.Zorvv[-1])
        #     self.Zbubv.append(self.Zorvl[-1])
        #     self.Vbubv.append(self.Vorvv[-1])
        #     self.Vbubl.append(self.Vorvl[-1])
        #     self.vapx.append(vapo.mixture.x)
        #     self.liqx.append(liqo.mixture.x)
        
        T = array(self.Tbub + self.Torv)
        P = array(self.Pbub + self.Porv)
        
        T0 = sum(array(self.mixture.list_Tc)*array(self.mixture.x))*3
        V0 = sum(array([sp.Vc for sp in self.mixture.list_of_species])*array(self.mixture.x))*1.5

        self.critical_case(T0, V0, [])
        
        self.P_isotherm = []
        self.V_isotherm = []
        
        if critical_isothem:
            P_it = self.crit_gas.P
            while P_it < P_max*1.5:
                iso_therm = self.copy_change_conditions(self.crit_gas.T,P_it,None,'gas')
                self.P_isotherm.append(iso_therm.P)
                self.V_isotherm.append(iso_therm.V)
                P_it += step_P_init
        
        V_ref = self.mixture.MM_m/340
        
        
        self.line_T = []
        self.line_P = []
        if isovolumetrical_case:
            V = array(self.Vbubl+self.Vorvv)
            lut_T_ref = interpolant('LUT','bspline',[self.Vbubl],self.Tbub)
            lut_P_ref = interpolant('LUT','bspline',[self.Vbubl],self.Pbub)
            T_ref = lut_T_ref(V_ref).__float__()
            P_ref = lut_P_ref(V_ref).__float__()
            P_line_r = P_ref
            while P_line_r < 60000:
                print(P_line_r)
                line_ref = self.copy_change_conditions(None,P_line_r,V_ref,'gas')
                self.line_P.append(line_ref.P)
                self.line_T.append(line_ref.T-273.15)
                P_line_r += step_P_init
                
       #if water_prep_case:
            
        



    def plot_envelope(self, Pexp=[], Texp=[], Vexp=[]):

        fig1 = figure(dpi = 150)
        ax = fig1.add_subplot(1, 1, 1)
        p1, = plot(self.Torv, [P/100 for P in self.Porv], 'r')
        p2, = plot(self.Tbub, [P/100 for P in self.Pbub], 'b')
        p3, = plot(Texp, Pexp, 'ko', fillstyle='none')
        p4, = plot(self.T-273.15, self.P/100, 'kx', markersize=10, markeredgewidth=2)
        p5, = plot(self.crit_gas.T-273.15, self.crit_gas.P/100, 'ko', markersize=10)
        if len(self.line_T)>0:
            p6, = plot(array(self.line_T), array(self.line_P)/100,'--g')
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
        p4, = plot(self.V, self.P/100, 'kx', markersize=10, markeredgewidth=2)
        p5, = plot(self.crit_gas.V, self.crit_gas.P/100, 'ko', markersize=10)
        if len(self.V_isotherm)>0:
            p6, = plot(array(self.V_isotherm), array(self.P_isotherm)/100,'k')
        fig2.tight_layout()
        fig2.show()
        config_plot(ax)
        ax.set_xlabel('Molar Volume /(m³/kmol)', fontsize=12)
        ax.set_ylabel('Pressure /bar', fontsize=12)
        fig2.tight_layout()
        fig2.show()

    def plot_envelope2(self, Pexp=[], Texp=[], Vexp=[]):

        fig1 = figure(dpi = 150)
        ax = fig1.add_subplot(1, 1, 1)
        p1, = plot(self.Torv, [P/100 for P in self.Porv], 'r')
        p2, = plot(self.Tbub, [P/100 for P in self.Pbub], 'b')
        p3, = plot(Texp, Pexp, 'ko', fillstyle='none')
        p4, = plot(self.T-273.15, self.P/100, 'kx', markersize=10, markeredgewidth=2)
        p7, = plot(self.critical_point2.T-273.15, self.critical_point2.P/100, 'ko', markersize=10)
        if len(self.line_T)>0:
            p6, = plot(array(self.line_T), array(self.line_P)/100,'--g')
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
        p4, = plot(self.V, self.P/100, 'kx', markersize=10, markeredgewidth=2)
        p7, = plot(self.critical_point2.V, self.critical_point2.P/100, 'ko', markersize=10)
        if len(self.V_isotherm2)>0:
            p6, = plot(array(self.V_isotherm2), array(self.P_isotherm2)/100,'k')
        fig2.tight_layout()
        fig2.show()
        config_plot(ax)
        ax.set_xlabel('Molar Volume /(m³/kmol)', fontsize=12)
        ax.set_ylabel('Pressure /bar', fontsize=12)
        fig2.tight_layout()
        fig2.show()

class solver_eos:
    
    def __init__(self,fluid):
        
        self.fluid = fluid
        self.dew_tol = 1e-5
        self.orv_tol = 1e-5
        self.flash_subs = 0.0005
        self.flash_tol = 1e-6
        self.max_run = 1e5
        self.critical_initial = 'default'
        T0 = sum(array(fluid.mixture.list_Tc)*array(fluid.mixture.x))*3
        V0 = sum(array([sp.Vc for sp in fluid.mixture.list_of_species])*array(fluid.mixture.x))*1.5
        self.critical_init_value = [T0, V0]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    