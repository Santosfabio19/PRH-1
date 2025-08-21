# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 09:21:22 2024

@author: Rodrigo Meira
"""

import pandas as pd
from numpy import arange
from matplotlib.pyplot import plot, figure
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter

from tqdm import tqdm
from time import process_time


yaxis_title_x = ['$P_{manifold}$ /bar', '$q_{tr}$/(m$^3$ /h)',
                 '$P_{fbhp_1}$ /bar', '$P_{choke_1}$ /bar', '$q_{average_1}$ /(m$^3$/h)',
                 '$P_{fbhp_2}$ /bar', '$P_{choke_2}$ /bar', '$q_{average_2}$ /(m$^3$/h)',
                 '$P_{fbhp_3}$ /bar', '$P_{choke_3}$ /bar', '$q_{average_3}$ /(m$^3$/h)',
                 '$P_{fbhp_4}$ /bar', '$P_{choke_4}$ /bar', '$q_{average_4}$ /(m$^3$/h)']

yaxis_title_y = ['$P_{manifold}$ /bar',
                 '$P_{choke_1}$ /bar', '$P_{choke_2}$ /bar', '$P_{choke_3}$ /bar', '$P_{choke_4}$ /bar',
                 '$P_{intake_1}$ /bar','$dP_{bcs_1}$ /bar',
                 '$P_{intake_2}$ /bar','$dP_{bcs_2}$ /bar',
                 '$P_{intake_3}$ /bar','$dP_{bcs_3}$ /bar',
                 '$P_{intake_4}$ /bar','$dP_{bcs_4}$ /bar']

yaxis_title_BCS_envelope =['$\Delta P_{bcs_1}$ /bar',
                           '$\Delta P_{bcs_2}$ /bar',
                           '$\Delta P_{bcs_3}$ /bar',
                           '$\Delta P_{bcs_4}$ /bar']

yaxis_title_u = ['$f_BP$ /Hz',
                 'f_{ESP_1} /Hz','alpha_1',
                 'f_{ESP_2} /Hz','alpha_2',
                 'f_{ESP_3} /Hz','alpha_3',
                 'f_{ESP_4} /Hz','alpha_4']


yaxis_title_slack = ['$\delta_{P_{manifold}}$ /bar', '$\delta_{q_{tr}}$/(m$^3$ /h)',
                     '$\delta_{P_{fbhp_1}}$ /bar', '$\delta_{P_{choke_1}}$ /bar', '$\delta_{q_{average_1}}$ /(m$^3$/h)',
                     '$\delta_{P_{fbhp_2}}$ /bar', '$\delta_{P_{choke_2}}$ /bar', '$\delta_{q_{average_2}}$ /(m$^3$/h)',
                     '$\delta_{P_{fbhp_3}}$ /bar', '$\delta_{P_{choke_3}}$ /bar', '$\delta_{q_{average_3}}$ /(m$^3$/h)',
                     '$\delta_{P_{fbhp_4}}$ /bar', '$\delta_{P_{choke_4}}$ /bar', '$\delta_{q_{average_4}}$ /(m$^3$/h)']

group_y = [[0],
           [1,5,6],
           [2,7,8],
           [3,9,10],
           [4,11,12]]

group_x = [[0,1],
           [2,3,4],
           [5,6,7],
           [8,9,10],
           [11,12,13]]

group_u = [[0],
           [1,2],
           [3,4],
           [5,6],
           [7,8]]

# Routines
def config_plot(axes):
    """
    Configures the axes of Figures
    """
    formatter = ScalarFormatter(useOffset=False, useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    axes.yaxis.set_major_formatter(formatter)

    axes.xaxis.set_minor_locator(AutoMinorLocator())
    axes.yaxis.set_minor_locator(AutoMinorLocator())
    axes.tick_params(which='both', direction='out', bottom=True, left=True)
    axes.tick_params(which='major', width=2)
    axes.tick_params(which='minor', width=1)

# Caminho do arquivo CSV
caminho_arquivo = 'E:\Drivers\YandexDisk\YandexDisk\descent.csv'

# Lê o arquivo CSV como um DataFrame
df = pd.read_csv(caminho_arquivo, sep=';')

# Exibe as primeiras linhas do DataFrame
print(df.head())


yk  = [df['yk1'],df['yk2'],df['yk3']]
ykn = [df['yk1c'],df['yk1c.1'],df['yk1c.2']]
t  = [t*60 for t in df['t min']]
sp = [df['sp1'],df['sp2'],df['sp3']]


fig = figure(dpi=120)
axes = fig.add_subplot(1, 1, 1)
p1, = plot(t, sp[0],'r--')
p3, = plot(t,  ykn[0],'k-')
p2, = plot(t,  yk[0],'b-')
config_plot(axes)
axes.set_ylabel(yaxis_title_y[0], fontsize=12)
axes.set_xlabel('time / s', fontsize=12)
fig.tight_layout()
fig.show()

fig = figure(dpi=120)
axes = fig.add_subplot(2, 1, 1)
p1, = plot(t, sp[1],'r--')
p3, = plot(t,  ykn[1],'k-')
p2, = plot(t,  yk[1],'b-')
config_plot(axes)
axes.set_ylabel(yaxis_title_y[1], fontsize=12)
axes.set_xlabel('time / s', fontsize=12)

axes = fig.add_subplot(2, 1, 2)
p1, = plot(t, sp[2],'r--')
p3, = plot(t,  ykn[2],'k-')
p2, = plot(t,  yk[2],'b-')
config_plot(axes)
axes.set_ylabel(yaxis_title_y[6], fontsize=12)
axes.set_xlabel('time / s', fontsize=12)


fig.tight_layout()
fig.show()

fig = figure(dpi=120)
axes = fig.add_subplot(2, 1, 1)
p1, = plot(t, sp[1],'r--')
p3, = plot(t,  ykn[1],'k-')
p2, = plot(t,  yk[1],'b-')
config_plot(axes)
axes.set_ylabel(yaxis_title_y[1], fontsize=12)
axes.set_xlabel('time / s', fontsize=12)

axes = fig.add_subplot(2, 1, 2)
p1, = plot(t, sp[2],'r--')
p3, = plot(t,  ykn[2],'k-')
p2, = plot(t,  yk[2],'b-')
config_plot(axes)
axes.set_ylabel(yaxis_title_y[6], fontsize=12)
axes.set_xlabel('time / s', fontsize=12)


fig.tight_layout()
fig.show()

#Envelope
fig = figure()
axes = fig.add_subplot(1, 1, 1)
#p1, = plot(plant.trend.xmk[contq*3+4,:], plant.trend.zmk[contq*2+1,:])
#plot(plant.trend.xmk[contq * 3 + 4, 0], plant.trend.zmk[contq * 2 + 1, 0], 'ro'
#     , fillstyle='none')
#for n_change in yt_change:
#    plot(plant.trend.xmk[contq * 3 + 4, n_change], plant.trend.zmk[contq * 2 + 1, n_change], 'rx')
p2, = plot([20.77, 28.55], [58.07, 206.6] ,'b--')
p3, = plot([72.68, 106.4], [27.91, 145.6] ,'r--')
config_plot(axes)
#axes.set_ylabel(yaxis_title_BCS_envelope[cont], fontsize=12)
#axes.set_xlabel('$\dot{q}$ /(m$^3$/h)', fontsize=12)
fig.tight_layout()
fig.show()




# for cont,ele in enumerate(group_y[1:]):
#     fig = figure()
#     axes = fig.add_subplot(2, 1, 1)
#     p1, = plot(time, plant.trend.ymk[ele[0],:] )
#     p2, = plot(time, (control.trend.yt[ele[0],:]+1) * y_ss[ele[0]] ,'r--')
#     setpoint.append((control.trend.yt[ele[0],:]+1) * y_ss[ele[0]])
#     config_plot(axes)
#     axes.set_ylabel(yaxis_title_y[ele[0]], fontsize=12)

#     axes = fig.add_subplot(2, 1, 2)
#     p1, = plot(time, plant.trend.ymk[ele[1], :] )
#     p2, = plot(time, (control.trend.yt[ele[1], :]+1) * y_ss[ele[1]] , 'r--')
#     config_plot(axes)
#     axes.set_ylabel(yaxis_title_y[ele[1]], fontsize=12)
#     setpoint.append((control.trend.yt[ele[1],:]+1) * y_ss[ele[1]])