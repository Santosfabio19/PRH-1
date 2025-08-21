# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 19:10:41 2024

@author: Fapex
"""

import win32com.client as winclt
import pandas as pd
from os import path, makedirs
from warnings import warn

def concat_df(data, object_dict, properties, time, phases = False):
    data_i = {}
    for name in object_dict.keys():
        if not phases:
            setup_properties = {'time': [time]}
            for prop in properties.keys():
                try:
                    setup_properties.update({prop: [getattr(object_dict[name], prop).GetValue(properties[prop])]})
                except:
                    setup_properties.update({prop: [None]})

            data_i[name] = pd.concat([data[name], pd.DataFrame(setup_properties)])
        else:
            data_i[name] = {}
            setup_properties = {'time': [time]}
            for phase in phases.keys():

                for prop in properties.keys():
                    try:
                        obj = getattr(object_dict[name].GetFluid(), phase)
                        setup_properties.update({prop: [getattr(obj, prop).GetValue(properties[prop])]})
                    except:
                        setup_properties.update({prop: [None]})
                data_i[name][phase] = pd.concat([data[name][phase], pd.DataFrame(setup_properties)])

    return data_i

def build_df(object_dict, properties, time, phases=False):
    data = {}
    for name in object_dict.keys():
        if not phases:
            setup_properties = {'time': [time]}
            for prop in properties.keys():
                try:
                    setup_properties.update({prop: [getattr(object_dict[name], prop).GetValue(properties[prop])]})
                except:
                    setup_properties.update({prop: [None]})

            data[name] = pd.DataFrame(setup_properties)

        else:
            data[name] = {}
            for phase in phases.keys():
                setup_properties = {'time': [time]}
                for prop in properties.keys():
                    try:
                        obj = getattr(object_dict[name].GetFluid(), phase)
                        setup_properties.update({prop: [getattr(obj, prop).GetValue(properties[prop])]})
                    except:
                        setup_properties.update({prop: [None]})

                data[name][phase] = pd.DataFrame(setup_properties)

    return data

class SimulationConnector:

    @property
    def components(self):
        return ['H2O', 'N2', 'CO2', 'C1', 'C2', 'C3', 'I-C4', 'N-C4', 'I-C5', 'N-C5', 'C6', 'C7', 'C8', 'C9',
                'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20+']

    @property
    def stream_names(self):
        """
        :return: Lista com os nomes das correntes materiais
        """
        return self._stream_names

    @stream_names.setter
    def stream_names(self, stream):
        """
        Adiciona correntes na lista de correntes
        :param stream: Lista com os nomes das correntes materiais a serem adicionadas.
        """
        self._stream_names.extend(stream)

    @property
    def stream_names_Phase(self):
        """
        :return: Lista com os nomes das correntes materiais que terão os dados
        das fases salvos
        """
        return self._stream_names_phase

    @stream_names_Phase.setter
    def stream_names_Phase(self, stream):
        """
        Adiciona correntes na lista de correntes que terão os dados da fase salvos
        :param stream: Lista com os nomes das correntes materiais a serem adicionadas.
        """
        return self._stream_names_phase.extend(stream)

    @property
    def control_names(self):
        """
        :return: Lista com os nomes dos controladores
        """
        return self._control_names

    @control_names.setter
    def control_names(self, controller):
        """
        Adiciona controladores na lista de controladores
        :param controller: Lista com os nomes doss controladores a serem adicionados.
        :return:
        """
        self._control_names.extend(controller)

    @property
    def heatExchanger_names(self):
        """
        :return: Lista com os nomes dos trocadores de calor
        """
        return self._heatExchangers_names

    @heatExchanger_names.setter
    def heatExchanger_names(self, heatExchanger):
        """
        Adiciona trocadores na lista de trocadores
        :param trocador: Lista com os nomes dos trocadores a serem adicionados.
        :return:
        """
        self._heatExchangers_names.extend(heatExchanger)

    @property
    def valve_names_W6(self):
        """
        :return: Lista com os nomes dos válvulas do subfluxograma W6
        """
        return self._valve_names_W6

    @valve_names_W6.setter
    def valve_names_W6(self, valve):
        """
        Adiciona válvulas na lista de válvulas do subfluxograma W6
        :param valve: Lista com os nomes das válvulas
        """
        self._valve_names_W6.extend(valve)

    @property
    def valve_names_W7(self):
        """
        :return: Lista com os nomes dos válvulas do subfluxograma W7
        """
        return self._valve_names_W7

    @valve_names_W7.setter
    def valve_names_W7(self, valve):
        """
        Adiciona válvulas na lista de válvulas do subfluxograma W7
        :param valve: Lista com os nomes das válvulas
        """
        self._valve_names_W7.extend(valve)

    @property
    def valve_names_main(self):
        """
        :return: Lista com os nomes dos válvulas do fluxograma principal
        """
        return self._valve_names_main

    @valve_names_main.setter
    def valve_names_main(self, valve):
        """
        Adiciona válvulas na lista de válvulas do fluxograma principal
        :param valve: Lista com os nomes das válvulas
        """
        self._valve_names_main.extend(valve)

    @property
    def properties_stream(self):
        """
        :return: dicionário com as propriedades das correntes materiais (stream_names) a serem consultados no unisim,
        associado às respectivas unidades - propriedade:unidade.
        """
        return self._properties_stream

    @properties_stream.setter
    def properties_stream(self, property_with_unit):
        """
        Atualiza e/ou adiciona propriedades das correntes materiais.
        :param property_with_unit: dicionário com o par propriedade:unidade.
        """
        self._properties_stream.update(property_with_unit)

    @property
    def properties_stream_phase(self):
        """
        :return: dicionário com as propriedades das fases correntes materiais (stream_names_Phase) a serem consultados no unisim,
        associado às respectivas unidades - propriedade:unidade.
        """
        return self._properties_stream_Phase

    @properties_stream_phase.setter
    def properties_stream_phase(self, property_with_unit):
        """
        Atualiza e/ou adiciona propriedades das fases correntes materiais.
        :param property_with_unit: dicionário com o par propriedade:unidade.
        """
        self._properties_stream_Phase.update(property_with_unit)

    @property
    def properties_controller(self):
        """
        :return: dicionário com as propriedades dos controladores (control_names) a serem consultados no unisim,
        associado às respectivas unidades - propriedade:unidade.
        """
        return self._properties_controller

    @properties_controller.setter
    def properties_controller(self, property_with_unit):
        """
        Atualiza e/ou adiciona propriedades dos controladores.
        :param property_with_unit: dicionário com o par propriedade:unidade.
        """
        self._properties_controller.update(property_with_unit)

    @property
    def properties_heatExch(self):
        """
        :return: dicionário com as propriedades de todos os trocaodores a serem consultados no unisim,
        associado às respectivas unidades - propriedade:unidade.
        """
        return self._properties_heatExch

    @properties_heatExch.setter
    def properties_heatExch(self, property_with_unit):
        """
        Atualiza e/ou adiciona propriedades dos trocadores.
        :param property_with_unit: dicionário com o par propriedade:unidade.
        """
        self._properties_heatExch.update(property_with_unit)

    @property
    def properties_valve(self):
        """
        :return: dicionário com as propriedades de todas as válvulas a serem consultados no unisim,
        associado às respectivas unidades - propriedade:unidade.
        """
        return self._properties_valve

    @properties_valve.setter
    def properties_valve(self, property_with_unit):
        """
        Atualiza e/ou adiciona propriedades das válvulas.
        :param property_with_unit: dicionário com o par propriedade:unidade.
        """
        self._properties_valve.update(property_with_unit)

    def __init__(self, file_path_unisim, folder_save, continue_simulation=True):
        """
        Classe para executar uma simulação do Unisim, usando o arquivo HISEP como base.
        :param file_path_unisim: caminho completo do arquivo unisim
        :param folder_save: caminho diferencial (./) para a pasta onde os resultados serão salvos
        """

        # Local para salvar os resultados

        self.folder = folder_save

        self._simulation_name = file_path_unisim.split("\\")[-1]

        # Continue simulation

        self._continue_simulation = continue_simulation

        # variável auxiliar para verificar
        self._flag_build = False

        # criação do folder
        if not path.exists(self.folder):
            makedirs(self.folder)

        # Acessando o unisim
        self.unisim = winclt.Dispatch("UnisimDesign.Application")

        # Abrindo o arquivo de simulação
        self.case = self.unisim.SimulationCases.Open(file_path_unisim)
        self.case.visible = 1

        #%% Definições do integrador
        self.integrador = self.case.Solver.Integrator  # Carregando o objeto integrador
        #step_time = self.integrador.GetStepSize('seconds') # Coletando o tempo de integração em segundos
        self.integrador.Mode = 0  # Definindo modo automático

        if not self._continue_simulation:
            self.integrador.Reset()  # Zera o tempo inicial

        # Streams
        self._stream_names = ['S01', 'S06', 'S04',
                              'S10', 'S15', 'S13',
                              'S27', 'S20',
                              'S45-2', 'S46-2',
                              'S55', 'S56',
                              'Reinjection',
                              'S32', 'S33',
                              'S76', 'Header-2']

        self._stream_names_phase = []

        # Controllers
        self._control_names = ['LIC-001', 'LIC-002', 'PIC-SC', 'HIC']

        # Valves
        self._valve_names_main = ['LV-001', 'LV-002', 'LV-003', 'Choke_Topside-2',
                                  'XV-002', 'XV-003', 'XV-004', 'XV-017', 'XV-018', 'XV-005']

        self._valve_names_W6 = ['HV-01']

        self._valve_names_W7 = ['HV-02']

        # Heat Exchangers
        self._heatExchangers_names = ['CO-01', 'CO-02']

        # Properties
        self._properties_stream = dict(Pressure='bar', StdGasFlow='STD_m3/h', MassFlow='kg/h', MolarFlow='kgmole/h',
                                       StdLiqVolFlow='m3/h',
                                       Temperature='C', MassDensity='kg/m3', VapourFraction='',
                                       ActualGasFlow='ACT_m3/h',
                                       ActualLiqFlow='m3/h', ActualVolumeFlow='m3/h')

        self._properties_stream_Phase = dict(StdGasFlow='STD_m3/h', MassFlow='kg/h', MolarFlow='kgmole/h',
                                       StdLiqVolFlow='m3/h',MassDensity='kg/m3', ActualGasFlow='ACT_m3/h',
                                       ActualLiqFlow='m3/h', ActualVolumeFlow='m3/h', PhaseFraction='')

        self._properties_controller = dict(PV='', SP='', OP='')

        self._properties_valve = dict(ActuatorPosition='%', PressureDrop='bar', FeedTemperature='C', FeedPressure='bar',
                                      MolarFlow='kgmole/h',
                                      MassFlow='kg/h')

        self._properties_heatExch = dict(Duty='kJ/h', MassFlow='kg/h', FeedTemperature='C', FeedPressure='bar',
                                      ProductTemperature='C')

        # Objects que contém os Objetos Unisim
        self.material_streams = {}
        self.material_streams_phase = {}
        self.controllers = {}
        self.valves = {}
        self.heatExchangers = {}
        self.spreadsheets = {}

        # Objetos para salvar os dados
        self.data_valves_properties = {}
        self.data_controllers = {}
        self.data_streams_properties = {}
        self.data_heatExchangers_properties = {}

        self.data_streams_phase_properties = {}

        self._phases_name = {'VapourPhase':{},
                             'LightLiquidPhase':{},
                             'HeavyLiquidPhase':{}}

        self.dataRGO = {}
        self.dataCO2 = {}
        self.data_components = {'CompMolarFraction': {},
                                'CompMassFraction': {},
                                'CompMassFlow': {},
                                'CompMolarFlow': {}}

        self.data_components_phase = {'MassFractions': {},
                                      'MolarFractions': {}}

    def build_objects(self):

        # %% Carregando o Flowsheet Principal
        main_pfd = self.case.Flowsheet  # Carregando o flowsheet principal

        sub_pfd_W6 = main_pfd.Flowsheets.Item('TPL2')  # subflowSheet Well-06
        sub_pfd_W7 = main_pfd.Flowsheets.Item('TPL3')  # subflowSheet Well-07

        # %% Definindo os objetos
        for name in self.stream_names:
            try:
                self.material_streams[name] = main_pfd.MaterialStreams.Item(name)
            except:
                warn('Material stream %s not found. Skiped' % name)

        for name in self.stream_names_phase:
            try:
                self.material_streams_phase[name] = main_pfd.MaterialStreams.Item(name)
            except:
                warn('Material stream %s not found. Skiped' % name)

        for name in self.control_names:
            try:
                self.controllers[name] = main_pfd.Operations.Item(name)
            except:
                warn('Controller %s not found. Skiped' % name)

        for name in self.heatExchanger_names:
            try:
                self.heatExchangers[name] = main_pfd.Operations.Item(name)
            except:
                warn('Heat Exchanger %s not found. Skiped' % name)

        for name in self.valve_names_W6:
            try:
                self.valves[name] = sub_pfd_W6.Operations.Item(name)
            except:
                warn('Valve %s not found in subflowsheet W6. Skiped' % name)

        for name in self.valve_names_W7:
            try:
                self.valves[name] = sub_pfd_W7.Operations.Item(name)
            except:
                warn('Valve %s not found in subflowsheet W7. Skiped' % name)

        for name in self.valve_names_main:
            try:
                self.valves[name] = main_pfd.Operations.Item(name)
            except:
                warn('Valve %s not found in main. Skiped' % name)

        # SpreadSheets
        self.spreadsheets['RGO-S10'] = main_pfd.Operations.Item('RGO-S10').Cell('B3')
        self.spreadsheets['RGO-S01'] = main_pfd.Operations.Item('RGO-S01').Cell('B3')
        self.spreadsheets['RGO-TopSide'] = main_pfd.Operations.Item('TopSide AUX').Cell('B3')
        self.spreadsheets['CO2-TopSide-S20'] = main_pfd.Operations.Item('TopSide AUX').Cell('B7')
        self.spreadsheets['CO2-TopSide-InletHISEP'] = main_pfd.Operations.Item('TopSide AUX').Cell('B10')
        self.spreadsheets['CO2-TQ01-LIQ-InletTQ01'] = main_pfd.Operations.Item('CO2 Separadores').Cell('B8')
        self.spreadsheets['CO2-TQ01-VAP-InletTQ01'] = main_pfd.Operations.Item('CO2 Separadores').Cell('C8')
        self.spreadsheets['CO2-TQ02-LIQ-InletTQ02'] = main_pfd.Operations.Item('CO2 Separadores').Cell('B9')
        self.spreadsheets['CO2-TQ02-VAP-InletTQ02'] = main_pfd.Operations.Item('CO2 Separadores').Cell('C9')

        if self._continue_simulation:
            self._load_dataframes()
        else:
            self._create_dataframes()

    def _create_dataframes(self):

        # Controllers
        self.data_controllers = build_df(self.controllers, self.properties_controller, 0.0)

        # Streams
        self.data_streams_properties = build_df(self.material_streams, self.properties_stream, 0.0)

        self.data_streams_phase_properties = build_df(self.material_streams_phase, self.properties_stream_phase,
                                                      0.0, phases=self._phases_name)

        # Heat Exchangers
        self.data_heatExchangers_properties = build_df(self.heatExchangers, self.properties_heatExch, 0.0)

        # Valves
        self.data_valves_properties = build_df(self.valves, self.properties_valve, 0.0)

        # components
        for name in self.material_streams.keys():
            setup_CompMolarFraction = {'time': [0.0]}
            setup_CompMassFraction = {'time': [0.0]}
            setup_CompMassFlow = {'time': [0.0]}
            setup_CompMolarFlow = {'time': [0.0]}

            CompMolarFraction_stream = self.material_streams[name].ComponentMolarFraction.Values
            CompMassFraction_stream = self.material_streams[name].ComponentMassFraction.Values
            CompMassFlow_stream = self.material_streams[name].ComponentMassFlow.GetValues('kg/h')
            CompMolarFlow_stream = self.material_streams[name].ComponentMolarFlow.GetValues('kgmole/h')

            for id, component in enumerate(self.components):
                setup_CompMolarFraction.update({component: [CompMolarFraction_stream[id]]})
                setup_CompMassFraction.update({component: [CompMassFraction_stream[id]]})
                setup_CompMassFlow.update({component: [CompMassFlow_stream[id]]})
                setup_CompMolarFlow.update({component: [CompMolarFlow_stream[id]]})

            self.data_components['CompMolarFraction'][name] = pd.DataFrame(setup_CompMolarFraction)
            self.data_components['CompMassFraction'][name] = pd.DataFrame(setup_CompMassFraction)
            self.data_components['CompMassFlow'][name] = pd.DataFrame(setup_CompMassFlow)
            self.data_components['CompMolarFlow'][name] = pd.DataFrame(setup_CompMolarFlow)

        # components in phase
        for name in self.material_streams_phase.keys():
            self.data_components_phase['MassFractions'][name] = {}
            self.data_components_phase['MolarFractions'][name] = {}
            for phase in self._phases_name:
                setup_CompMolarFraction = {'time': [0.0]}
                setup_CompMassFraction = {'time': [0.0]}

                CompMolarFraction_stream = getattr(self.material_streams_phase[name].GetFluid(), phase).MolarFractions()
                CompMassFraction_stream = getattr(self.material_streams_phase[name].GetFluid(), phase).MassFractions()

                for id, component in enumerate(self.components):
                    setup_CompMolarFraction.update({component: [CompMolarFraction_stream[id]]})
                    setup_CompMassFraction.update({component: [CompMassFraction_stream[id]]})

                self.data_components_phase['MassFractions'][name][phase] = pd.DataFrame(setup_CompMassFraction)
                self.data_components_phase['MolarFractions'][name][phase] = pd.DataFrame(setup_CompMolarFraction)

        # EXTRAS
        self.dataRGO = pd.DataFrame({'time': [0.0],
                                     'S10': [self.spreadsheets['RGO-S10'].CellValue],
                                     'S01': [self.spreadsheets['RGO-S01'].CellValue],
                                     'Top': [self.spreadsheets['RGO-TopSide'].CellValue]})

        self.dataCO2 = pd.DataFrame({'time': [0.0],
                                     'TopSide-S20': [self.spreadsheets['CO2-TopSide-S20'].CellValue],
                                     'TopSide-InletHISEP': [self.spreadsheets['CO2-TopSide-InletHISEP'].CellValue],
                                     'TQ01-LIQ-InletTQ01': [self.spreadsheets['CO2-TQ01-LIQ-InletTQ01'].CellValue],
                                     'TQ01-VAP-InletTQ01': [self.spreadsheets['CO2-TQ01-VAP-InletTQ01'].CellValue],
                                     'TQ02-LIQ-InletTQ02': [self.spreadsheets['CO2-TQ02-LIQ-InletTQ02'].CellValue],
                                     'TQ02-VAP-InletTQ02': [self.spreadsheets['CO2-TQ02-VAP-InletTQ02'].CellValue]})

    def _load_dataframes(self):

        for name in self.control_names:
            try:
                self.data_controllers[name] = pd.read_csv('{}{}_{}.csv'.format(self.folder, 'data_controller', name))
            except:
                warn('Could not load data for controller %s' % name)

        # Streams
        for name in self.stream_names:
            try:
                self.data_streams_properties[name] = pd.read_csv('{}{}_{}.csv'.format(self.folder, 'data_stream', name))
            except:
                warn('Could not load data for stream %s' % name)

        # HeatExch
        for name in self.heatExchanger_names:
            try:
                self.data_heatExchangers_properties[name] = pd.read_csv('{}{}_{}.csv'.format(self.folder, 'data_heatExch', name))
            except:
                warn('Could not load data for Heat Exchanger %s' % name)
        # Valves
        for name in self.valves.keys():
            try:
                self.data_valves_properties[name] = pd.read_csv('{}{}_{}.csv'.format(self.folder, 'data_valve', name))
            except:
                warn('Could not load data for Valve %s' % name)

        for title in self.data_components.keys():
            for name in self.stream_names:
                try:
                    self.data_components[title][name] = pd.read_csv('{}data_{}_stream_{}.csv'.format(self.folder, title, name))
                except:
                    warn('Could not load components data for Stream %s' % title)

        # EXTRAS
        self.dataRGO = pd.read_csv('{}{}.csv'.format(self.folder, 'RGO'))

        self.dataCO2 = pd.read_csv('{}{}.csv'.format(self.folder, 'CO2'))

    def time_step(self, sample_time, sample_time_unit):

        # Simulando um passo de integração
        self.integrador.RunFor(sample_time, sample_time_unit)

        # Salvando as variáveis
        t_sim = self.integrador.CurrentTime.GetValue(sample_time_unit)

        self.data_controllers = concat_df(self.data_controllers, self.controllers, self.properties_controller, t_sim)

        self.data_heatExchangers_properties = concat_df(self.data_heatExchangers_properties, self.heatExchangers, self.properties_heatExch, t_sim)

        self.data_valves_properties = concat_df(self.data_valves_properties, self.valves, self.properties_valve, t_sim)

        self.data_streams_properties = concat_df(self.data_streams_properties, self.material_streams,
                                                 self.properties_stream, t_sim)

        self.data_streams_phase_properties = concat_df(self.data_streams_phase_properties, self.material_streams_phase,
                                                                  self.properties_stream_phase, t_sim, phases=self._phases_name)

        for name in self.material_streams.keys():
            setup_CompMolarFraction = {'time': [t_sim]}
            setup_CompMassFraction = {'time': [t_sim]}
            setup_CompMassFlow = {'time': [t_sim]}
            setup_CompMolarFlow = {'time': [t_sim]}

            CompMolarFraction_stream = self.material_streams[name].ComponentMolarFraction.Values
            CompMassFraction_stream = self.material_streams[name].ComponentMassFraction.Values
            CompMassFlow_stream = self.material_streams[name].ComponentMassFlow.GetValues('kg/h')
            CompMolarFlow_stream = self.material_streams[name].ComponentMolarFlow.GetValues('kgmole/h')

            for id, component in enumerate(self.components):
                setup_CompMolarFraction.update({component: [CompMolarFraction_stream[id]]})
                setup_CompMassFraction.update({component: [CompMassFraction_stream[id]]})
                setup_CompMassFlow.update({component: [CompMassFlow_stream[id]]})
                setup_CompMolarFlow.update({component: [CompMolarFlow_stream[id]]})

            self.data_components['CompMolarFraction'][name] = pd.concat(
                [self.data_components['CompMolarFraction'][name],
                 pd.DataFrame(setup_CompMolarFraction)])
            self.data_components['CompMassFraction'][name] = pd.concat([self.data_components['CompMassFraction'][name],
                                                                        pd.DataFrame(setup_CompMassFraction)])
            self.data_components['CompMassFlow'][name] = pd.concat([self.data_components['CompMassFlow'][name],
                                                                    pd.DataFrame(setup_CompMassFlow)])
            self.data_components['CompMolarFlow'][name] = pd.concat([self.data_components['CompMolarFlow'][name],

                                                                     pd.DataFrame(setup_CompMolarFlow)])
        # componentes in phase
        for name in self.material_streams_phase.keys():
            for phase in self._phases_name:
                setup_CompMolarFraction = {'time': [t_sim]}
                setup_CompMassFraction = {'time': [t_sim]}

                CompMolarFraction_stream = getattr(self.material_streams_phase[name].GetFluid(), phase).MolarFractions()
                CompMassFraction_stream = getattr(self.material_streams_phase[name].GetFluid(), phase).MassFractions()

                for id, component in enumerate(self.components):
                    setup_CompMolarFraction.update({component: [CompMolarFraction_stream[id]]})
                    setup_CompMassFraction.update({component: [CompMassFraction_stream[id]]})

                self.data_components_phase['MassFractions'][name][phase] = pd.concat([self.data_components_phase['MassFractions'][name][phase],
                                                                               pd.DataFrame(setup_CompMassFraction)])
                self.data_components_phase['MolarFractions'][name][phase] = pd.concat([self.data_components_phase['MolarFractions'][name][phase],
                                                                            pd.DataFrame(setup_CompMolarFraction)])

        self.dataRGO = pd.concat([self.dataRGO, pd.DataFrame({'time': [t_sim],
                                                              'S10': [self.spreadsheets['RGO-S10'].CellValue],
                                                              'S01': [self.spreadsheets['RGO-S01'].CellValue],
                                                              'Top': [self.spreadsheets['RGO-TopSide'].CellValue]})])

        self.dataCO2 = pd.concat([self.dataCO2, pd.DataFrame({'time': [t_sim],
                                                              'TopSide-S20': [
                                                                  self.spreadsheets['CO2-TopSide-S20'].CellValue],
                                                              'TopSide-InletHISEP': [self.spreadsheets[
                                                                                         'CO2-TopSide-InletHISEP'].CellValue],
                                                              'TQ01-LIQ-InletTQ01': [self.spreadsheets[
                                                                                         'CO2-TQ01-LIQ-InletTQ01'].CellValue],
                                                              'TQ01-VAP-InletTQ01': [self.spreadsheets[
                                                                                         'CO2-TQ01-VAP-InletTQ01'].CellValue],
                                                              'TQ02-LIQ-InletTQ02': [self.spreadsheets[
                                                                                         'CO2-TQ02-LIQ-InletTQ02'].CellValue],
                                                              'TQ02-VAP-InletTQ02': [self.spreadsheets[
                                                                                         'CO2-TQ02-VAP-InletTQ02'].CellValue]})])

    def save(self, quit_unisim=True, save_simulation=None):
        """
        Salva os dados de simulação como arquivos .csv na pasta, além de logs.
        """
        if save_simulation is None:
            save_simulation = True if self._continue_simulation else False

        # Summary
        if quit_unisim:
            self.case.Close(save_simulation)
            self.unisim.Quit() # Fecha o unisim

        # Saving
        filename = '_0_summary.txt'
        with open(self.folder + filename, 'w') as f:
            f.write('file:{}\n'.format(self._simulation_name))
            f.write('material streams:' + str(list(self.material_streams.keys())) + '\n')
            f.write('controllers:' + str(list(self.controllers.keys())) + '\n')
            f.write('valves:' + str(list(self.valves.keys())) + '\n')

        # exporting units
        filename = '_stream-properties.txt'
        outfile = open(self.folder + filename, 'w')
        outfile.writelines(
            ['{}:{}\n'.format(prop, self.properties_stream[prop]) for prop in self.properties_stream.keys()])
        outfile.close()

        filename = '_stream-phase-properties.txt'
        outfile = open(self.folder + filename, 'w')
        outfile.writelines(
            ['{}:{}\n'.format(prop, self.properties_stream_phase[prop]) for prop in self.properties_stream_phase.keys()])
        outfile.close()

        filename = '_valve-properties.txt'
        outfile = open(self.folder + filename, 'w')
        outfile.writelines(
            ['{}:{}\n'.format(prop, self.properties_valve[prop]) for prop in self.properties_valve.keys()])
        outfile.close()

        filename = '_controller-properties.txt'
        outfile = open(self.folder + filename, 'w')
        outfile.writelines(
            ['{}:{}\n'.format(prop, self.properties_controller[prop]) for prop in self.properties_controller.keys()])
        outfile.close()

        filename = '_data_heatExch-properties.txt'
        outfile = open(self.folder + filename, 'w')
        outfile.writelines(
            ['{}:{}\n'.format(prop, self.properties_heatExch[prop]) for prop in self.properties_heatExch.keys()])
        outfile.close()

        # saving data
        for name in self.data_controllers.keys():
            self.data_controllers[name].to_csv('{}{}_{}.csv'.format(self.folder, 'data_controller', name), index=False)

        for name in self.heatExchangers.keys():
            self.data_heatExchangers_properties[name].to_csv('{}{}_{}.csv'.format(self.folder, 'data_heatExch', name), index=False)

        for name in self.data_valves_properties.keys():
            self.data_valves_properties[name].to_csv('{}{}_{}.csv'.format(self.folder, 'data_valve', name), index=False)

        for name in self.data_streams_properties.keys():
            self.data_streams_properties[name].to_csv('{}{}_{}.csv'.format(self.folder, 'data_stream', name),
                                                      index=False)

        for stream in self.data_streams_phase_properties.keys():
            for phase in self.data_streams_phase_properties[stream].keys():
                self.data_streams_phase_properties[stream][phase].to_csv('{}data_stream_{}_{}.csv'.format(self.folder, stream, phase),
                                                          index=False)

        for title in self.data_components.keys():
            for name in self.data_components[title].keys():
                self.data_components[title][name].to_csv('{}data_{}_stream_{}.csv'.format(self.folder, title, name),
                                                         index=False)
        for title in self.data_components_phase.keys():
            for stream in self.data_components_phase[title].keys():
                for phase in self.data_components_phase[title][stream].keys():
                    self.data_components_phase[title][stream][phase].to_csv('{}data_{}_stream_{}_{}.csv'.format(self.folder, title, stream, phase),
                                                             index=False)

        self.dataRGO.to_csv('{}{}.csv'.format(self.folder, 'RGO'), index=False)
        self.dataCO2.to_csv('{}{}.csv'.format(self.folder, 'CO2'), index=False)