# -*- coding: utf-8 -*-
"""
Created on Sat Mar  1 20:37:27 2025

@author: rodri
"""

from eos_database_new_resumed import *
from gc_eos import gc_eos_class
from solver_thermo import solver_eos
from hydrate_module import *
import win32com.client as winclt
import pandas as pd
from os import path, makedirs
from warnings import warn
from Reorder import Reorder
import pickle
from matplotlib.pyplot import plot, figure

list_names  = ["CH4",	  "C2H6",	  "C3H8",	  "iC4H10",  "nC4H10",  "iC5H12",  "nC5H12",  "nC6H14",  "nC7H16",  "nC8H18",	
               "nC9H20",  "nC10H22", "nC11H24", "nC12H26", "nC14H30", "N2", "H2O", "CO2", "C15+"]

nwe = [0.238095238095238,0.0262608309264608,0.0261719617862697,0.00648744723394801,
       0.0183070428793601,0.0112419462341702,0.0105754276827372,0.0193290379915574,
       0.0210708965703968,0.0190209179995643,0.0173567422429736,0.01597710063164,
       0.0148135692205084,0.0138181989380497,0.0122023902400072] + \
      [0.000844257, 0.002962305, 0.446300822] + [0.066207509]
        
dict_composition = {list_names[i]: nwe[i] for i in range(len(nwe))}

mixture_nwe = Mixture(list_of_species, dict_composition)

volumn_desviation = [0]*19

stream_nwe = gc_eos_class(mixture_nwe, 85+273.15, 1.2e4, None, 2, -1, Aij, Bij, Cij, volumn_desviation, 'liquid')

x0 = [0] + [4]*14 + [0] + [0] + [0] + [18]
y0 = [2]*5 + [0]*13 + [0]*2 + [7] + [0] + [40] + [0]
w0 = [1e-10]*15 + [0] + [1] + [0] + [0]

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


class evaluate_curves:
    
    def __init__(self,fluid,x0,y0,w0):
        
        self.fluid = fluid
        self.solver = solver_eos(fluid)
        
        self.set_solve_par(x0, y0, w0)
        
        
    def copy_evaluate_curves(self,fluid):
        
        return evaluate_curves(fluid)
    
    def copy_evaluate_curves_change_fluid(self,T,P,V,x):
        
        fluid = self.fluid.copy_change_x_and_conditions(T,P,V,x,self.phase)
        
        return evaluate_curves(fluid)
    
    def set_solve_par(self,x0,y0,w0,step_flash = 10,Tc0=800,Vc0=0.08,Pmax_iso = 25e3, Pmax_h20 = 6e4):
    
        self.solver.isochoric_evaluation = True
        self.solver.critical_evaluation = True
        self.solver.water_evaluation = True    
        
        self.solver.set_estimated_conditions(x0, y0, w0)
        self.solver.step_flash = step_flash
        self.solver.critical_init_value = [Tc0, Vc0]
        self.solver.P_max_isochoric = Pmax_iso
        self.solver.P_max_water = Pmax_h20
    
    def build_curves(self,T_list,T_criteria,P_list,P_criteria,delta_T,delta_P):
        
        self.solver.build_phase_enveloped(T_list,T_criteria,P_list,P_criteria,delta_T,delta_P)
        
        self.solver_hydrate_I = hydrate_solver(water_ref_I,hydrate_model_I,'H2O',solver,solver.lut_T_sat)
        self.solver_hydrate_I.setup_duple_int_dt_fun(list_int_h_ig_dt)
        
        self.solver_hydrate_I.evaluate_hydrate_curve(solver.water_line_P,solver.fluid)
        
    def build_curves_resumed(self,T_list,T_criteria,delta_T,delta_P):
        
        self.solver.build_phase_bubble(T_list,T_criteria,delta_T,delta_P)
        
        self.solver_hydrate_I = hydrate_solver(water_ref_I,hydrate_model_I,'H2O',self.solver,self.solver.lut_T_sat)
        self.solver_hydrate_I.setup_duple_int_dt_fun(list_int_h_ig_dt)
        
        self.solver_hydrate_I.evaluate_hydrate_curve(self.solver.water_line_P,self.solver.fluid)
        
    
    def plot_curve(self,pt_points=[]):
        
        self.solver.plot_PT(pt_points=pt_points, hydrate_points = [self.solver_hydrate_I.T_hydrate,self.solver_hydrate_I.P_hydrate])
    
    
    
class Simulator:
    
    def __init__(self,file_path_unisim,unit_time,sample_time):
        
        self.unit_time = unit_time
        
        self.sample_time = sample_time
        
        unisim = winclt.Dispatch("UnisimDesign.Application")
        
        self.case = unisim.SimulationCases.Open(file_path_unisim)
        
        self.pfd = self.case.Flowsheet
        
        self.integrador = self.case.Solver.Integrator
        
        self.integrador.Reset()
        
    def set_visible(self):

        self.case.visible = 1

    def set_unvisible(self):

        self.case.visible = 0
        
    def get_streams(self,list_streams):
        
        self.names_streams = list_streams
        self.Streams = {name: self.pfd.MaterialStreams[name] for name in list_streams}
    
    def get_pumps(self,list_pumps):
        
        self.names_pumps = list_pumps
        self.Pumps = {name: {'speed': self.pfd.Operations[name].SpeedInCompressor,
                             'Efficency':self.pfd.Operations[name].CompAdiabaticEff,
                             'Pressure':self.pfd.Operations[name].PressureIncrease,
                             'obj': self.pfd.Operations[name]} for name in list_pumps}
        
    def get_valves(self,list_valves):
        
        self.names_valves = list_valves
        self.Valves = {name: self.pfd.Operations[name].PercentOpen for name in list_valves}
    
    def pumps_new_value(self,values,selected_pump):
        
        if selected_pump == "pump_1":
            self.Pumps[self.names_pumps[0]]['speed'].SetValue(values[0])
        elif selected_pump == "pump_2":
            self.Pumps[self.names_pumps[1]]['speed'].SetValue(values[1])
        else:
            for i in range(len(self.names_pumps)):
                self.Pumps[self.names_pumps[i]]['speed'].SetValue(values[i])
    
    def valves_new_value(self,values):
        
        [self.Valves[self.names_valves[i]].SetValue(values[i]) for i in range(len(self.names_valves))]
    
    def rearrange(self,x):
        
        y = [i for i in x[1:13]] + [i for i in x[16:19]] + [x[13]] + [x[15]] + [x[14]] + [x[0]]
        
        return y
        
    def get_state(self,name):
        
        dict = {'P': [self.Streams[name].Pressure()],
                'T': [self.Streams[name].Temperature()],
                'dot_n': [self.Streams[name].MolarFlow()],
                'x': [self.rearrange(self.Streams[name].ComponentMolarFraction())],
                'MassFlow':[(self.Streams[name].MassFlow())],
                }
        
        return dict
    def get_pump_state(self, name):
        dict = {'PressureIncrease':[self.Pumps[name].PressureIncrease()]
            
            }
        
    
    def new_state(self,name,dict):    
        
        dict['P'].append(self.Streams[name].Pressure())
        dict['T'].append(self.Streams[name].Temperature())
        dict['dot_n'].append(self.Streams[name].MolarFlow())
        dict['x'].append(self.rearrange(self.Streams[name].ComponentMolarFraction()))
        #dict['MassFlow'].append(self.rearrange(self.Streams[name].MassFlow()))
        
    def start_streams(self):
        self.stream_states = {name: self.get_state(name) for name in self.names_streams}
        
        # Inicializa dicionário de listas para cada bomba e tipo de valor
        self.inputs_pumps = {
            name: {
                'speed': [],
                'efficiency': [],
                'pressure': []
            }
            for name in self.names_pumps
        }
    
        # Inicializa dicionário de listas para válvulas
        self.inputs_valves = {name: [] for name in self.names_valves}
    
        # Preenche o primeiro conjunto de valores
        for name in self.names_pumps:
            self.inputs_pumps[name]['speed'].append(self.Pumps[name]['speed'].Value)
            self.inputs_pumps[name]['efficiency'].append(self.Pumps[name]['Efficency'].Value)
            self.inputs_pumps[name]['pressure'].append(self.Pumps[name]['Pressure'].Value)
    
        for name in self.names_valves:
            self.inputs_valves[name].append(self.Valves[name].Value)


        
        
    def plot_stream_state(self, stream, state):
        y_name = {'P': 'Pressure /bar', 'T': 'Temperature /°C', 'dot_n': 'Molar flow rate /(kmol/h)'}
        factor_name = {'P': 100, 'T': 1, 'dot_n': 1}
        
        y = self.stream_states[stream][state]
        
        time = [i/60 for i in range(len(y))]
        
        fig1 = figure(dpi=150)
        ax = fig1.add_subplot(1, 1, 1)
        p1, = plot(time, [i/factor_name[state] for i in y], 'r')
        ax.set_ylabel(f"{stream} {y_name[state]}", fontsize=12)  # Adicionando o nome da stream aqui
        ax.set_title(stream)  
        ax.set_xlabel('Time /h', fontsize=12)
        fig1.tight_layout()

        
                # ---- Salvamento automático ----
        os.makedirs("figuras", exist_ok=True)  # cria pasta se não existir
        filename = f"figuras/{stream}_{state}.png"
        fig1.savefig(filename, dpi=300, bbox_inches="tight")
        print(f"Gráfico salvo em: {filename}")

        fig1.show()
    def get_streams_PT_points_in_k(self,k):
        
        P = [self.stream_states[stream]['P'][k]/100 for stream in self.names_streams]
        T = [self.stream_states[stream]['T'][k] for stream in self.names_streams]
        
        return [T,P]
    def get_stream_mflowrate(self,k):
        flow = [self.stream_states[stream]['MassFlow'][k] for stream in self.names_streams]
        return flow
    
    def new_eta_head(self,eta_values,PressureIncrease_values,selected_pump):
        if selected_pump == "pump_1":
            indices = [0]
        elif selected_pump == "pump_2":
            indices = [1]
        else:  
            indices = [0, 1]
        
        for i in indices:
            self.Pumps[self.names_pumps[i]]['Efficency'].SetValue(eta_values[i], "%")
            self.Pumps[self.names_pumps[i]]['Pressure'].SetValue(PressureIncrease_values[i],"kPa")
        
    def save_data(self, file_name, external_data_dict):
        data = {
            'streams': self.stream_states,
            'pumps': self.inputs_pumps,  # works with new structure
            'valves': self.inputs_valves
        }
        for key, value in external_data_dict.items():
            data[key] = value
        with open(file_name, 'wb') as arquivo:
            pickle.dump(data, arquivo)

        
    def simulate_n_save_streams(self):
        self.integrador.RunFor(self.sample_time, self.unit_time)
    
        for name in self.names_pumps:
            self.inputs_pumps[name]['speed'].append(self.Pumps[name]['speed'].Value)
            self.inputs_pumps[name]['efficiency'].append(self.Pumps[name]['Efficency'].Value)
            self.inputs_pumps[name]['pressure'].append(self.Pumps[name]['Pressure'].Value)
    
        for name in self.names_valves:
            self.inputs_valves[name].append(self.Valves[name].Value)

                
            