# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 22:20:49 2025

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
import pickle
from matplotlib.pyplot import plot, figure

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
    
def plot_curves_by_stage(dict,pump,stage):
    
    cores = [
        '#aec7e8',  # Azul claro
        '#ff7f0e',  # Laranja
        '#2ca02c',  # Verde
        '#d62728',  # Vermelho
        '#9467bd',  # Roxo
        '#8c564b',  # Marrom
        '#e377c2',  # Rosa
        '#7f7f7f',  # Cinza
        '#bcbd22',  # Verde-oliva
        '#17becf',  # Ciano
        '#1f77b4'   # Azul
    ]
    
    fig1 = figure(dpi = 150)
    ax = fig1.add_subplot(1, 1, 1)
    for j, comp in  enumerate(dict.keys()):
        for i, table in enumerate(dict[comp][pump][stage]):
            p1, = plot(table[0], table[1],'--',color = cores[j])
    ax.set_ylabel('Head /m', fontsize=20)
    ax.set_xlabel('Volumetric Flow Rate /(m³/h)', fontsize=20)
    fig1.tight_layout()
    fig1.show()
    
    fig2 = figure(dpi = 150)
    ax = fig2.add_subplot(1, 1, 1)
    for j, comp in  enumerate(dict.keys()):
        for i, table in enumerate(dict[comp][pump][stage]):
            p2, = plot(table[0], [k*100 for k in table[2]],'--',color = cores[j])
    ax.set_ylabel('Efficiency', fontsize=20)
    ax.set_xlabel('Volumetric Flow Rate /(m³/h)', fontsize=20)
    fig1.tight_layout()
    fig1.show()
    

with open('data_new_cobeq_5', 'rb') as file:
    # Carregar o objeto do arquivo pickle
    dados = pickle.load(file)

ks = [49, 99, 149, 199, 249, 299]

name_streams = ['feed_sep','1st_heat_out','scrubber_in','pump_1_in',
                'pump_1_out','pump_2_in','pump_2_out']


for i,k in enumerate(ks):
    P = [dados['streams'][stream]['P'][k]/100 for stream in name_streams]
    T = [dados['streams'][stream]['T'][k] for stream in name_streams]
    dados['curves'][i].plot_curve([T,P])
     
