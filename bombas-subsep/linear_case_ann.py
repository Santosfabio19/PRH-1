# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:11:37 2024

@author: rodri
"""


import torch as thor
import torch.nn as ann
import numpy as np
from matplotlib.pyplot import plot, figure
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter


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

def func(x):
    return 3*x + 7 + thor.rand(10,1)

x_train = thor.rand(10,1)*10
y_train = func(x_train)

class regressor(ann.Module):
    def __init__(self):
        super().__init__()
        self.layers = ann.Sequential(
            ann.Linear(1,5),
            ann.ReLU(),
            ann.Linear(5,4),
            ann.ReLU(),
            ann.Linear(4,1)
        )
        
    def forward(self,x):
        return self.layers(x)

model = regressor()

optim_solver = thor.optim.Adam(model.parameters(), lr = 0.001)
loss_fn = ann.MSELoss()

epochs = 10

for i in range(epochs):
    y_pred = model(x_train)
    loss = loss_fn(y_pred, y_train)
    optim_solver.zero_grad()
    loss.backward()
    optim_solver.step()


y_pred_op = model(x_train)
y_pred_op = y_pred_op.detach().numpy()

fig1 = figure(dpi=250)
ax = fig1.add_subplot(1, 1, 1)
p1, = plot(np.array(x_train),np.array(y_train),'o')
p2, = plot(np.array(x_train),np.array(y_pred_op),'-r')    
























