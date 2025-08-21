# -*- coding: utf-8 -*-

from connect_and_run_unisim import SimulationConnector
from numpy import log, array, save, load, pi, log10
from scipy.optimize import fsolve
from time import time

#This file needs specified Duty in heat exchangers

#%% Configuração

# UniSim File
file_path = r"C:\Users\ddssa\OneDrive - Universidade Federal da Bahia\Projetos\HISEP\Simulação\Trocadores\PartidaInicial_Well6A7A_WC00_inj_MAIOR_FI_MAIOR_clean_HeatExch_spec_UA_Duty_PY_continue_sync_bypass0.usc"

# Folder to save date
folder = './data6A7A-HeatExch-UA-testing-Qi/'

continue_simulation = False
save_simulation = False
# Para continuar uma simulação anterior:
# 1 - executar primeira Simulação com continue_simulation = False e save_simulation = True
# 2 - executar segunda Simulação com continue_simulation = True e save_simulation = True

Simulation = SimulationConnector(file_path, folder, continue_simulation)

# Adicionar novos objetos
Simulation.stream_names = ['S42', 'S43-2', 'S52', 'S53']
#Simulation.valve_names_main = []
#Simulation.valve_names_W6 = []
#Simulation.valve_names_W7 = []
#Simulation.control_names = ['TIC-100', 'TIC-101']
Simulation.stream_names_phase = ['S42','S43-2','S52', 'S53']

# Adicionar novar propriedades e modificar unidades
Simulation.properties_stream = {'MassEnthalpy':'kJ/kg','Viscosity':'cP','ThermalConductivity':'W/m-K','MassHeatCapacity':'J/kg-K'}
Simulation.properties_stream_phase = {'MassEnthalpy':'kJ/kg','Viscosity':'cP','ThermalConductivity':'W/m-K','MassHeatCapacity':'J/kg-K'}
#Simulation.properties_valve = {}
#Simulation.properties_controller = {}

#%% Construindo objetos
Simulation.spreadsheets['mass_Cv_C01_mean'] = Simulation.case.Flowsheet.Operations.Item('HeatExch').Cell('B3')
Simulation.spreadsheets['mass_Cv_C02_mean'] = Simulation.case.Flowsheet.Operations.Item('HeatExch').Cell('D3')

Simulation.build_objects()
#%% Atualização da carga térmica dos trocadores
def heatExh_BE_implicit(T2, T2_k_1, T1, rho, Cv, V, m, hi, hout, UA, n, nT):
    return -Cv*(n/nT)*V*rho*(T2-T2_k_1)/(1/60.) + m*(hi-hout)-(n/nT)*UA*((T1-4) - (T2-4)) / log((T1-4) / (T2-4))

def heatExh_BE_implicit_corr(T2, T2_k_1, T1, rho, Cv, data_tube, m, hi, hout, UA, ntubo):
    # Calor total, volume total
    V = data_tube['N_bend'] * data_tube['L_to_bend']*pi*(data_tube['din']/2)**2
    return -Cv*(ntubo)*V*rho*(T2-T2_k_1)/(1/60.) + m*(hi-hout)-(ntubo)*UA*((T1-4) - (T2-4)) / log((T1-4) / (T2-4))

def heatExh_BE_implicit_corr_2(T2, T2_k_1, T1, rho, Cv, data_tube, m, hi, hout, UA, ntubo):
    # Calor para 1 tubo -> divide a vazão
    V = data_tube['N_bend'] * data_tube['L_to_bend']*pi*(data_tube['din']/2)**2
    return -Cv*V*rho*(T2-T2_k_1)/(1/60.) + m/(ntubo)*(hi-hout)-UA*((T1-4) - (T2-4)) / log((T1-4) / (T2-4))
def eval_property(stream, property):
    if (Simulation.data_streams_properties[stream][property] < 0).any():
        property_value = 0
        for phase in ['HeavyLiquidPhase', 'LightLiquidPhase', 'VapourPhase']:
            property_value = property_value + array(Simulation.data_streams_phase_properties[stream][phase]['PhaseFraction'])[-1] * \
                    array(Simulation.data_streams_phase_properties[stream][phase][property])[-1]
    else:
        property_value = array(Simulation.data_streams_properties[stream][property])[-1]

    return property_value

def fun_UA(stream_in, stream_out, range_ntubo, data_tube, Qi):
    # %%
    din = data_tube['din']
    dout = data_tube['dout']
    # Volume = 69,5 L = 0.0695 m3
    L_to_bend = data_tube['L_to_bend']
    N_bend = data_tube['N_bend']
    L_total = N_bend * L_to_bend

    rho_in = eval_property(stream_in,'MassDensity')  # kg/m3
    rho_out = eval_property(stream_out,'MassDensity')  # kg/m3

    rho = (rho_in + rho_out) / 2

    mu_in = eval_property(stream_in, 'Viscosity') * (1 / 1000)  # Pa s
    mu_out = eval_property(stream_out, 'Viscosity') * (1 / 1000)  # Pa s

    mu = (mu_in + mu_out) / 2

    Re_ntubo = lambda ntubo: rho * ((Qi / ntubo) / (pi * (din / 2) ** 2)) * din / mu

    Re_ntubo_lista = [Re_ntubo(ntubo) for ntubo in range_ntubo]

    # Pr = cp*mu/kf
    cp_in = eval_property(stream_in, 'MassHeatCapacity')
    cp_out = eval_property(stream_out, 'MassHeatCapacity')

    cp = (cp_in + cp_out)/2

    kf_in = eval_property(stream_in, 'ThermalConductivity')
    kf_out = eval_property(stream_out, 'ThermalConductivity')

    kf = (kf_in + kf_out) / 2

    Pr = cp * mu / kf

    # Nu = 0,023*Re^(0,8)*Pr^(0,4)
    # Nu = h*d/kf

    fun_Nu_DittusBoelter = lambda Rei: 0.023 * Rei ** 0.8 * Pr ** 0.4
    fun_Nu_Colburn = lambda Rei: 0.023 * Rei ** 0.8 * Pr ** 0.33
    fun_f = lambda Rei: (1.82 * log10(Rei) - 1.64) ** (-2)
    fun_Nu_Gnielinksi = lambda Rei: (fun_f(Rei) / 8) * (Rei - 1000) * Pr / (
                1.07 + 12.7 * (fun_f(Rei) / 8) ** 0.5 * (Pr ** (2 / 3) - 1))

    Nu_DittusBoelter = [fun_Nu_DittusBoelter(Re_i) for Re_i in Re_ntubo_lista]
    Nu_Colburn = [fun_Nu_Colburn(Re_i) for Re_i in Re_ntubo_lista]

    Nu_Gnielinksi = [fun_Nu_Gnielinksi(Re_i) for Re_i in Re_ntubo_lista]

    h_DittusBoelter = [Nu_DittusBoelter_i * kf / din for Nu_DittusBoelter_i in Nu_DittusBoelter]
    h_Colburn = [Nu_Colburn_i * kf / din for Nu_Colburn_i in Nu_Colburn]
    h_Gnielinksi = [Nu_Gnielinksi_i * kf / din for Nu_Gnielinksi_i in Nu_Gnielinksi]

    UA = lambda h: 1 / (1 / (2 * pi * din * L_total * h) + log(dout / din) / (2 * pi * 14.2 * L_total))

    UA_DittusBoelter = [UA(h_DittusBoelter_i) for h_DittusBoelter_i in h_DittusBoelter]
    UA_Colburn = [UA(h_Colburn_i) for h_Colburn_i in h_Colburn]
    UA_Gnielinksi = [UA(h_Gnielinksi_i) for h_Gnielinksi_i in h_Gnielinksi]

    return UA_DittusBoelter, UA_Colburn, UA_Gnielinksi
#%%
ntrocadores_1 = 1
ntrocadores_2 = 1

range_ntubo = [4, 8]
ntubo_1 = 8
ntubo_2 = 6
id_corr = 0

data_tube_1 = {'din':0.03814,
               'dout':0.04830,
               'L_to_bend':3.5,
               'N_bend':15}

# data_tube_2 = {'din':0.0381,
#                'dout':0.04826,
#                'L_to_bend':3.5,
#                'N_bend':16}

#%% Simulação
sample_time = 1.0
sample_time_unit = 'minutes'
n_sim = 5 # Número de passos

start_time = time()
abort = False
k = 0

UA_C01 = []
UA_C02 = []

T_diff_01 = []
T_diff_02 = []

print(Simulation.heatExchangers['CO-01'].Duty.GetValue('kJ/h'),
      Simulation.heatExchangers['CO-02'].Duty.GetValue('kJ/h'))

while (k in range(n_sim)) and (not abort):
    print(Simulation.integrador.CurrentTime.GetValue('minutes'))

    m_C01 = array(Simulation.data_heatExchangers_properties['CO-01']['MassFlow'])[-1]
    h1_C01 = array(Simulation.data_streams_properties['S42']['MassEnthalpy'])[-1]
    h2_C01 = array(Simulation.data_streams_properties['S43-2']['MassEnthalpy'])[-1]

    T2_k_1_C01 = array(Simulation.data_streams_properties['S43-2']['Temperature'])[-1]
    T1_C01 = array(Simulation.data_streams_properties['S42']['Temperature'])[-1]
    Cv_C01 = Simulation.spreadsheets['mass_Cv_C01_mean'].CellValue
    rho_C01 = array(Simulation.data_streams_properties['S43-2'].MassDensity)[-1]

    Q_C01_k_i = []
    for ni in range(ntrocadores_1):

        Qi = (array(Simulation.data_streams_properties['S42']['ActualVolumeFlow'])[-1] * (1. / 3600))/ntrocadores_1

        UA_C01_corr = fun_UA('S42', 'S43-2', range_ntubo, data_tube_1, Qi)[id_corr]

        T2_C01_k_corr = fsolve(
            lambda T2: heatExh_BE_implicit_corr_2(T2, T2_k_1_C01, T1_C01, rho_C01, Cv_C01, data_tube_1, m_C01//ntrocadores_1,
                                                  h1_C01, h2_C01, UA_C01_corr[range_ntubo.index(ntubo_1)], ntubo_1), T2_k_1_C01)

        Q_C01_k_i.append(ntubo_1*UA_C01_corr[range_ntubo.index(ntubo_1)] * (((T1_C01 - 4.) - (T2_C01_k_corr - 4.)) / log((T1_C01 - 4.) / (T2_C01_k_corr - 4.))))

    m_C02 = array(Simulation.data_heatExchangers_properties['CO-02']['MassFlow'])[-1]
    h1_C02 = array(Simulation.data_streams_properties['S52']['MassEnthalpy'])[-1]
    h2_C02 = array(Simulation.data_streams_properties['S53']['MassEnthalpy'])[-1]

    T2_k_1_C02 = array(Simulation.data_streams_properties['S53']['Temperature'])[-1]
    T1_C02 = array(Simulation.data_streams_properties['S52']['Temperature'])[-1]
    Cv_C02 = Simulation.spreadsheets['mass_Cv_C02_mean'].CellValue
    rho_C02 = array(Simulation.data_streams_properties['S53'].MassDensity)[-1]

    Q_C02_k_i = []
    for ni in range(ntrocadores_2):

        Qi = (array(Simulation.data_streams_properties['S52']['ActualVolumeFlow'])[-1] * (1. / 3600)) / ntrocadores_2

        UA_C02_corr = fun_UA('S52', 'S53', range_ntubo, data_tube_2, Qi)[id_corr]

        T2_C02_k_corr = fsolve(lambda T2: heatExh_BE_implicit_corr_2(T2, T2_k_1_C02, T1_C02, rho_C02, Cv_C02, data_tube_2, m_C02/ntrocadores_2,
                                                                     h1_C02, h2_C02, UA_C02_corr[range_ntubo.index(ntubo_2)], ntubo_2), T2_k_1_C02)

        Q_C02_k_i.append(ntubo_2*UA_C02_corr[range_ntubo.index(ntubo_2)]*(((T1_C02 - 4.) - (T2_C02_k_corr - 4)) / log((T1_C02 - 4) / (T2_C02_k_corr - 4))))

    if abs(T2_C01_k_corr-T2_k_1_C01) < 100:
        Simulation.heatExchangers['CO-01'].Duty.SetValue(sum(Q_C01_k_i), 'kJ/h')
    if abs(T2_C02_k_corr-T2_k_1_C02) < 100:
        Simulation.heatExchangers['CO-02'].Duty.SetValue(sum(Q_C02_k_i), 'kJ/h')

    UA_C01.append(float(UA_C01_corr[range_ntubo.index(ntubo_1)]))
    UA_C01.append(float(UA_C02_corr[range_ntubo.index(ntubo_2)]))
    T_diff_01.append(T2_C01_k_corr-T2_k_1_C01)
    T_diff_02.append(T2_C02_k_corr-T2_k_1_C02)

    if k == 50:
        Simulation.controllers['PIC-SC'].SP.SetValue(22000)

    if k == 150:
        Simulation.controllers['PIC-SC'].SP.SetValue(23000)

    if k == 250:
        Simulation.controllers['PIC-SC'].SP.SetValue(24000)

    print(Q_C01_k_i, Q_C02_k_i)
    print(sum(Q_C01_k_i), sum(Q_C02_k_i))
    print(T2_k_1_C01, T2_k_1_C02)
    print(T2_C01_k_corr, T2_C02_k_corr)

    try:
        Simulation.time_step(sample_time, sample_time_unit)
        k += 1
    except Exception as error_message:
        print('Erro no Unisim: '+str(error_message))
        abort = True

stop_time = time()

delta_time = stop_time - start_time
print(f'Tempo de processamento: {delta_time} s')

#%% Fechando o arquivo
quit_unisim = False
# if abort:
#     quit_unisim = False
# else:
#     quit_unisim = True

Simulation.save(quit_unisim, save_simulation)

save(folder+'UA_C01', UA_C01)
save(folder+'UA_C02', UA_C02)
save(folder+'T_diff_01', T_diff_01)
save(folder+'T_diff_02', T_diff_02)