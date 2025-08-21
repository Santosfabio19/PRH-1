# -*- coding: utf-8 -*-
"""
Atualizado para simular variação de uma única válvula com 12 degraus (100 min cada)
Primeira perturbação ocorre após 10 minutos.
@author: rodri
"""

from setup_simulation import *
import pickle

import numpy as np
import matplotlib.pyplot as plt
import os

# CONFIGURAÇÕES DE SIMULAÇÃO

file_path_unisim = r"C:\Users\fabio\Downloads\bombas_novo_tema\model\nwe_sub_separtion_dynamic(1).usc"

name_streams = ["pump_1_in", "pump_1_out", "pump_2_in", "pump_2_out"]
name_valves = ["Valve_bypass_1", "Valve_heat_1", "Valve_bypass_2", "Valve_heat_2"]
name_pumps = ["pump_1", "pump_2"]

unit_time = "minutes"
sample_time = 1

# ALTERADO PARA 12 DEGRAUS, 100 MINUTOS CADA + 10 MINUTOS INICIAIS
tempo_inicial = 30  # minutos sem perturbação
tempo_por_degrau = 100
num_degraus = 12
nsim = tempo_inicial + num_degraus * tempo_por_degrau  # total de minutos de simulação

# ESCOLHA O QUE SIMULAR
simular_ns = True
simular_nv = False

# nv = true: ESCOLHA QUAL VÁLVULA ALTERAR (0 a 3)
valvula_a_alterar = 2  # 0 = Valve_bypass_1, 1 = Valve_heat_1, 2 = Valve_bypass_2, 3 = Valve_heat_2

# ns = true: ESCOLHA QUAL BOMBA ALTERAR ns (0 a 1)
bomba_a_variar = 1
  # 0 para pump1, 1 para pump2

# CONFIGURAÇÃO INICIAL
simu = Simulator(file_path_unisim, unit_time, sample_time)
simu.get_pumps(name_pumps)
simu.get_streams(name_streams)
simu.get_valves(name_valves)
simu.set_visible()
simu.start_streams()

# DEFINIÇÃO DOS VALORES DE NS E NV
ns_constante = [85, 70]  # valor fixo para ambas as bombas
nv_base = [50.0, 50.0, 10.0, 50.0]
nv_limites = [10, 90]

# GERA VARIAÇÕES DE NV PARA A VÁLVULA ESCOLHIDA (com variação aleatória)
lista_nv = []
if simular_nv:
    np.random.seed(42)
    passo_max = 10
    valor = nv_base[valvula_a_alterar]
    for _ in range(num_degraus):
        passo = np.random.uniform(-passo_max, passo_max)
        valor = np.clip(valor + passo, nv_limites[0], nv_limites[1])
        nv = nv_base.copy()
        nv[valvula_a_alterar] = valor
        lista_nv.append(nv)
else:
    lista_nv = [nv_base] * num_degraus

# GERA VARIAÇÕES DE NS (quando TRUE)
lista_ns = []
if simular_ns:
    np.random.seed(42)
    for _ in range(num_degraus):
        if bomba_a_variar == 0:
            n1 = np.random.uniform(70, 100)
            n2 = ns_constante[1]
        else:
            n1 = ns_constante[0]
            n2 = np.random.uniform(60, 80)
        lista_ns.append([n1, n2])
else:
    lista_ns = [ns_constante] * num_degraus

# PRIMEIROS 10 MINUTOS - SEM PERTURBAÇÃO (valores base)
simu.pumps_new_value(ns_constante)
simu.valves_new_value(nv_base)
#for _ in range(tempo_inicial):
#    simu.simulate_n_save_streams()

# LOOP DE SIMULAÇÃO COM 12 DEGRAUS DE 100 MINUTOS
for i in range(num_degraus):
    simu.pumps_new_value(lista_ns[i])
    simu.valves_new_value(lista_nv[i])

    #for _ in range(tempo_por_degrau):
    #    simu.simulate_n_save_streams()

simu.set_unvisible()

# PLOTA P e T
simu.plot_stream_state("scrubber_in", "P")  # Pressão
simu.plot_stream_state("scrubber_in", "T")  # Temperatura

# SALVAR RESULTADOS
P = simu.stream_states["scrubber_in"]["P"]
T = simu.stream_states["scrubber_in"]["T"]

nome_valvula = name_valves[valvula_a_alterar]
nome_bomba = name_pumps[bomba_a_variar]

def salvar_com_seguranca(objeto, nome_arquivo, descricao):
    try:
        with open(nome_arquivo, "wb") as f:
            pickle.dump(objeto, f)
        print(f"{descricao} salvos como: {nome_arquivo}")
    except Exception as e:
        print(f"[ERRO] Falha ao salvar {descricao} em '{nome_arquivo}': {e}")

sufixo = f"{nome_bomba}" if simular_ns else f"{nome_valvula}"

# PLOTA E SALVA GRÁFICO DE PRESSÃO
plt.figure()
simu.plot_stream_state("scrubber_in", "P")
plt.savefig(f"pressao_scrubber_in_{sufixo}.png", dpi=300)
plt.close()

# PLOTA E SALVA GRÁFICO DE TEMPERATURA
plt.figure()
simu.plot_stream_state("scrubber_in", "T")
plt.savefig(f"temperatura_scrubber_in_{sufixo}.png", dpi=300)
plt.close()

# SALVA OS DADOS
salvar_com_seguranca(P, f"pressao_scrubber_in_{sufixo}.pkl", "Dados de Pressão")
salvar_com_seguranca(T, f"temperatura_scrubber_in_{sufixo}.pkl", "Dados de Temperatura")
salvar_com_seguranca(lista_nv, f"valores_nv_pressao_scrubber_in_{sufixo}.pkl", "Valores de NV")

if simular_ns:
    salvar_com_seguranca(lista_ns, f"valores_ns_pressao_scrubber_in_{sufixo}.pkl", "Valores de NS")
