import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from casadi import *

# Função para carregar os dados, remover os primeiros pontos e ajustar o tempo
def load_and_adjust_time(file_path, time_offset=0, decimal=',', cut_points=175):
    try:
        # Ler o arquivo Excel
        data = pd.read_excel(file_path, decimal=decimal)
        
        # Remover os primeiros 'cut_points' pontos
        data = data.iloc[cut_points:]
        
        # Extrair colunas de tempo, pressões e vazão mássica
        time = data.iloc[:, 0].values  # Tempo (primeira coluna)
        pressure_head = data.iloc[:, 1].values  # Pressão na cabeça do poço (segunda coluna)
        pressure_bottom = data.iloc[:, 2].values  # Pressão no fundo (terceira coluna)
        mass_flow_rate = data.iloc[:, 3].values  # Vazão mássica (quarta coluna)
        
        # Ajustar o tempo para começar após o último tempo do arquivo anterior
        time = time - time[0] + time_offset
        
        return time, pressure_head, pressure_bottom, mass_flow_rate, time[-1]  # Retorna o último tempo para o próximo arquivo
    
    except FileNotFoundError:
        print(f"Arquivo não encontrado: {file_path}")
        return None, None, None, None, time_offset
    except Exception as e:
        print(f"Erro ao processar {file_path}: {e}")
        return None, None, None, None, time_offset

# Função para processar e interpolar os dados
def process_and_interpolate(file_paths, num_points=9000+114999, decimal=',', cut_points=175):
    all_times = []
    all_pressure_head = []
    all_pressure_bottom = []
    all_mass_flow_rates = []
    time_offset = 0  # Inicializa o deslocamento de tempo
    
    for file_path in file_paths:
        time, pressure_head, pressure_bottom, mass_flow_rate, last_time = load_and_adjust_time(file_path, time_offset, decimal, cut_points)
        
        if time is not None:
            all_times.append(time)
            all_pressure_head.append(pressure_head)
            all_pressure_bottom.append(pressure_bottom)
            all_mass_flow_rates.append(mass_flow_rate)
            time_offset = last_time  # Atualiza o deslocamento para o próximo arquivo
    
    # Concatenar todos os dados
    total_time = np.concatenate(all_times)
    total_pressure_head = np.concatenate(all_pressure_head)
    total_pressure_bottom = np.concatenate(all_pressure_bottom)
    total_mass_flow_rate = np.concatenate(all_mass_flow_rates)
    
    # Ordenar os dados pelo tempo e remover valores duplicados
    sorted_indices = np.argsort(total_time)
    total_time = total_time[sorted_indices]
    total_pressure_head = total_pressure_head[sorted_indices]
    total_pressure_bottom = total_pressure_bottom[sorted_indices]
    total_mass_flow_rate = total_mass_flow_rate[sorted_indices]
    
    unique_times, unique_indices = np.unique(total_time, return_index=True)
    total_pressure_head = total_pressure_head[unique_indices]
    total_pressure_bottom = total_pressure_bottom[unique_indices]
    total_mass_flow_rate = total_mass_flow_rate[unique_indices]
    
    # Criar interpolação com CasADi
    interpolator_pressure_head = interpolant("LUT", "bspline", [unique_times], total_pressure_head)
    interpolator_pressure_bottom = interpolant("LUT", "bspline", [unique_times], total_pressure_bottom)
    interpolator_mass_flow = interpolant("LUT", "bspline", [unique_times], total_mass_flow_rate)
    
    # Criar novos pontos de tempo uniformemente espaçados para interpolação
    new_time = np.linspace(min(unique_times), max(unique_times), num=num_points)
    interpolated_pressure_head = interpolator_pressure_head(new_time)
    interpolated_pressure_bottom = interpolator_pressure_bottom(new_time)
    interpolated_mass_flow_rate = interpolator_mass_flow(new_time)
    
    return unique_times, total_pressure_head, total_pressure_bottom, total_mass_flow_rate, new_time, interpolated_pressure_head, interpolated_pressure_bottom, interpolated_mass_flow_rate

# Lista de arquivos de entrada
file_paths = [
    r"C:\Users\fabio\projects\PRH-1\unisimtable\ALFAsim\rodadas\sorte.xls"
]

# Processar todos os arquivos
unique_times, total_pressure_head, total_pressure_bottom, total_mass_flow_rate, new_time, pressure_head, pressure_bottom, mass_flow_rate = process_and_interpolate(file_paths)

# Plotar os dados reais e interpolados
plt.figure(figsize=(20, 6))
plt.scatter(unique_times, total_mass_flow_rate, color='red', label="Dados Reais (Vazão)", alpha=0.5, s=10)
plt.plot(new_time, mass_flow_rate, '.', label="Interpolação (Vazão)", color='blue')
plt.xlabel('Tempo (s)')
plt.ylabel('Vazão Mássica (kg/s)')
plt.title('Comparação entre dados reais e interpolados (Vazão Mássica)')
plt.legend()
plt.grid(True)
plt.show()

# Plot da pressão na cabeça do poço
plt.figure(figsize=(12, 6))
plt.scatter(unique_times, total_pressure_head, color='green', label="Dados Reais (Pressão Cabeça)", alpha=0.5, s=10)
plt.plot(new_time, pressure_head, '-', label="Interpolação (Pressão Cabeça)", color='purple')
plt.xlabel('Tempo (s)')
plt.ylabel('Pressão na Cabeça do Poço (Pa)')
plt.title('Comparação entre Pressão real e interpolada (Cabeça do Poço)')
plt.legend()
plt.grid(True)
plt.show()

# Plot da pressão no fundo do poço
plt.figure(figsize=(12, 6))
plt.scatter(unique_times, total_pressure_bottom, color='orange', label="Dados Reais (Pressão Fundo)", alpha=0.5, s=10)
plt.plot(new_time, pressure_bottom, '-', label="Interpolação (Pressão Fundo)", color='brown')
plt.xlabel('Tempo (s)')
plt.ylabel('Pressão no Fundo do Poço (Pa)')
plt.title('Comparação entre Pressão real e interpolada (Fundo do Poço)')
plt.legend()
plt.grid(True)
plt.show()


print(len(new_time))
print(len(unique_times))
print(len(np.array(mass_flow_rate)))




tempo_pos_corte = np.array([

0,
22999,
23099,
25999,
26000,
28999,
29000,
31999,
32000,
34999,
35000,
37999,
38000,
40999,
41000,
43999,
44000,
46999,
47000,
49999,
50000,
52999,
53000,
55999,
56000,
58999,
59000,
61999,
62000,
64999,
65000,
67999,
68000,
70999,
71000,
73999,
74000,
76999,
77000,
79999,
80000,
82999,
83000,
85999,
86000,
88999,
89000,
91999,
92000,
94999,
95000,
97999,
98000,
100999,
101000,
103999,
104000,
106999,
107000,
109999,
110000,
112999,
113000,
115999,
116000,
119000,
120000,



    
])  # Adapte com seus tempos reais
pressao_pos_corte = np.array([

6,5e+07,
6,5e+07,
6,55e+07,
6,55e+07,
6,9e+07,
6,9e+07,
6,45e+07,
6,45e+07,
6,75e+07,
6,75e+07,
7,1e+07,
7,1e+07,
6,65e+07,
6,65e+07,
7e+07,
7e+07,
6,5e+07,
6,5e+07,
6,95e+07,
6,95e+07,
6,4e+07,
6,4e+07,
6,85e+07,
6,85e+07,
7,2e+07,
7,2e+07,
6,7e+07,
6,7e+07,
7,05e+07,
7,05e+07,
6,6e+07,
6,6e+07,
6,95e+07,
6,95e+07,
6,55e+07,
6,55e+07,
7,15e+07,
7,15e+07,
6,75e+07,
6,75e+07,
6,35e+07,
6,35e+07,
6,8e+07,
6,8e+07,
7,25e+07,
7,25e+07,
6,65e+07,
6,65e+07,
7e+07,
7e+07,
6,45e+07,
6,45e+07,
6,9e+07,
6,9e+07,
7,3e+07,
7,3e+07,
6,5e+07,
6,5e+07,
6,85e+07,
6,85e+07,
7,2e+07,
7,2e+07,
6,75e+07,
6,75e+07,
6,8e+07,
6,8e+07,
6,85e+07,
])  
# Adapte com seus dados
temperatura_pos_corte = np.array([
80,
80,
82,
82,
83,
83,
81,
81,
84,
84,
85,
85,
83,
83,
86,
86,
87,
87,
85,
85,
88,
88,
89,
89,
87,
87,
90,
90,
89,
89,
88,
88,
86,
86,
87,
87,
85,
85,
84,
84,
82,
82,
83,
83,
81,
81,
80,
80,
82,
82,
84,
84,
85,
85,
86,
86,
87,
87,
88,
88,
89,
89,
90,
90,
88,
88,
87,
]) 

# Vetor de tempo pós-interpolação
new_time = np.linspace(0, max(tempo_pos_corte), len(new_time))

# Criar os vetores P_node1 e T_transiente sem interpolação
P_node1 = np.zeros_like(new_time)
T_transiente = np.zeros_like(new_time)

# Preencher os vetores de acordo com os intervalos
for i in range(len(tempo_pos_corte) - 1):
    mask = (new_time >= tempo_pos_corte[i]) & (new_time < tempo_pos_corte[i + 1])
    P_node1[mask] = pressao_pos_corte[i]
    T_transiente[mask] = temperatura_pos_corte[i]

# Garantir que o último valor seja preenchido corretamente
P_node1[new_time >= tempo_pos_corte[-1]] = pressao_pos_corte[-1]
T_transiente[new_time >= tempo_pos_corte[-1]] = temperatura_pos_corte[-1]

# Exibir tamanhos para verificação
print((P_node1), (T_transiente))  # Deve retornar (123999, 123999)
