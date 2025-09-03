# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 10:09:21 2025

@author: fabio
"""

import pickle
import matplotlib.pyplot as plt 
cami = r"C:\Users\fabio\Downloads\bombas_novo_tema\model\simulação sem sets.pkl"
aa = r"C:\Users\fabio\Downloads\bombas_novo_tema\model\simulação3008.pkl"
with open(aa,"rb") as (g):
    #ns = [[ 90,50],[ 90,60], [ 90,70],[100,70],[100,80],[120,80]] # n1 C [70,100], n2 C [60,80]
    simu_set_pump_1 = pickle.load(g)
with open(cami, "rb") as f:
    simu = pickle.load(f)
    



import matplotlib.pyplot as plt
t = [0]*360
# Definição dos pares de variáveis a serem plotados
variaveis = [
    ("scrubber_in", "P", "Scrubber In - Pressão"),
    ("scrubber_in", "T", "Scrubber In - Temperatura"),
    ("pump_2_out", "P", "Pump 2 Out - Pressão"),
    ("pump_2_out", "T", "Pump 2 Out - Temperatura"),
    ("pump_1_out", "T", "Pump 1 Out - Temperatura"),
    ("1st_heat_out", "T", "1st Heat Exchanger Out - Temperatura")
]
def plot_comparacao(simu,simu_set_pump_1):
    for stream, var, titulo in variaveis:
        plt.figure(figsize=(7,5))
        t = [i/60 for i in range(len(simu["streams"][stream][var]))]
        plt.plot(t,simu["streams"][stream][var], label="Simulação padrão")
        plt.plot(t,simu_set_pump_1["streams"][stream][var], label="Pump 1 ajustada")
        plt.title(titulo)
        plt.xlabel("tempo / min ")
        plt.ylabel(var)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()



