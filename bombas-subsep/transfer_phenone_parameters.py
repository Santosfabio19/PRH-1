# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 16:14:47 2024

@author: rodri
"""

from numpy import exp, log, array, roots, zeros, linalg, abs

class viscosity:
    
    def __init__(self,mixture,dipole_list):
        self.mixture = mixture
        self.dipole_list = dipole_list
        self.dipr = [52.46*mixture.list_Pc[i]*dipole_list[i]/mixture.list_Tc[i]**2 \
                     for i in range(len(dipole_list))]
        No = 6.023*1e26
        R = 8314.472
        self.kte = (R*(No**2))**(1/6)
        
        self.Tcm = sum(array(mixture.x)*array(mixture.list_Tc))
        self.Pcm = sum(array(mixture.x)*array(mixture.list_Pc))
        
    
    def evaluate_viscosity(self,T,P):
        
        Tr = T/array(self.mixture.list_Tc)
        
        Pr = P/array(self.mixture.list_Pc)
        
        Fp = [0 for i in Pr]
        for i in range(len(Tr)):
            if self.dipr[i]>=0 and self.dipr[i]<0.022:
                Fp[i] = 1
            elif self.dipr[i] >= 0.022 and self.dipr[i]<0.075:
                Fp[i] = 1 + 30.55*((0.292 - self.mixture.list_Zc[i])**1.72)
            elif self.dipr[i] >= 0.075:
                Fp[i] = 1 + 30.55*((0.292 - self.mixture.list_Zc[i])**1.72)*abs(0.96 + 0.1*Tr[i] - 0.7);
            else:
                Fp[i] = 0;
        Fpm = sum(array(self.mixture.x)*array(Fp))
        Trm = T/self.Tcm
        Prm = P/self.Pcm
        
        pvap = self.mixture.evaluate_pvap(T)
        
        pvapm = sum(array(self.mixture.x)*array(pvap))
        
        epm = self.kte*((self.Tcm/((self.mixture.MM_m**3)*((self.Pcm*1000)**4)))**(1/6));
        
        Z1 = Fpm*(0.807*(Trm**0.618) - 0.357*exp(-0.449**Trm) + 0.340*exp(-4.058*Trm) + 0.018)
        
        
        a1 = 1.245e-3
        a2 = 5.1726
        b1 = 1.6553            
        b2 = 1.2723
        c1 = 0.4489             
        c2 = 3.0578
        d1 = 1.7368             
        d2 = 2.2310
        e = 1.3088
        f1 = 0.9425             
        f2 = -0.1853
        gama = -0.32861          
        epsilon = -7.6351
        delta = -37.7332        
        quissi = 0.4489

        a = (a1/Trm)*exp(a2*(Trm**gama))
        b = a*(b1*Trm - b2)
        c = (c1/Trm)*exp(c2*(Trm**delta))
        d = (d1/Trm)*exp(d2*(Trm**epsilon))
        f = f1*exp(f2*(Trm**quissi))
        
        alfa = 3.262 + 14.98*(Prm**5.508)
        beta = 1.390 + 5.746*Prm
        
        
        if Trm <= 1 and Prm<(pvapm/self.Pcm): 
            Z2 = 0.600 + 0.760*Prm**alfa + (6.99*Prm**beta - 0.6)*(1 - Trm)
        elif Trm > 1 and Prm > 0:
            if Trm < 40 and Prm <= 100:
                Z2 = Z1*(1 + a*(Prm**e)/(b*(Prm**f) + (1 + c*(Prm**d))**-1))
            else:
                Z2 = 1
        else:
            Z2 = 1
        Y = Z2/Z1
        
        Fqm = 1
        
        Fpp = (1 + (Fpm - 1)*(Y**-3))/Fpm

        Fqp = (1 + (Fqm - 1)*((Y**-1) - 0.007*((log(Y))**4)))/Fqm

        return Z2*Fpp*Fqp/epm;
        