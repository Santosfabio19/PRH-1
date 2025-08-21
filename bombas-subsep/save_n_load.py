# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 16:13:53 2025

@author: rodri
"""

import pickle

def load_pickle(file_name):
    with open(file_name, 'rb') as arquivo:  # 'rb' = read binary
        dados_carregados = pickle.load(arquivo)
    return dados_carregados
    