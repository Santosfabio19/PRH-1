# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 23:11:03 2025

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
    
def plot_stream_state(dict,stream,state):
    
    y_name = {'P': 'Pressure /bar', 'T': 'Temperature /°C', 'dot_n': 'Molar flow rate /(kmol/h)'}
    factor_name = {'P': 100, 'T': 1, 'dot_n': 1}
    
    y = dict[stream][state]
    
    time = [i/60 for i in range(len(y))]
    
    fig1 = figure(figsize=(10, 3), dpi = 150)
    ax = fig1.add_subplot(1, 1, 1)
    p1, = plot(time, [i/factor_name[state] for i in y], 'r')
    ax.set_ylabel(y_name[state], fontsize=12)
    ax.set_xlabel('Time /h', fontsize=12)
    fig1.tight_layout()
    fig1.show()

def plot_pumps_state(dict,pump):
    
    y_name = 'Pump speed /rpm'
    factor_name = 1/60
    
    y = dict[pump][0:-1]
    
    time = [i/60 for i in range(len(y))]
    
    fig1 = figure(figsize=(10, 3), dpi = 150)
    ax = fig1.add_subplot(1, 1, 1)
    p1, = plot(time, [i/factor_name for i in y], 'r')
    ax.set_ylabel(y_name, fontsize=12)
    ax.set_xlabel('Time /h', fontsize=12)
    fig1.tight_layout()
    fig1.show()