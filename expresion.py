#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 18:12:35 2024
@author: andreaross

Script para el análisis de datos tras la cuantificación de la expresión génica
mediante qPCR. Se analizan los datos de expresión génica de los genes de interés
Se calcula la expresión relativa de los genes de interés siguiendo el método 2^-ΔΔCt
propuesto por Taylor et al. (2019).

"""

# Importación de librerías

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn
from scipy import stats 

# Carga de los datos
file = '/Users/andreaross/Desktop/bioinfo_andrea/ds.csv' # Ruta del archivo con los datos de expresión génica
df = pd.read_csv(file)

# Definición de los grupos experimentales
control = 2
housekeeping = 'mir16'
target = 'mir141'

# Eliminación de los valores NaN
df_copy = df[['grupo', 'mir16' , 'mir141']].copy()
df_copy.dropna(inplace=True)
df_copy.reset_index(drop=True, inplace=True)

# Cálculo de CT mean para el gen de referencia y el gen de interés
def ct_mean(df, control, housekeeping, target):
    df_control = df[df['grupo'] == control]
    
    ct_mean_housekeeping = round(np.mean(df_control[housekeeping]), 3)
    ct_mean_target = round(np.mean(df_control[target]), 3)

    return ct_mean_housekeeping, ct_mean_target
    
# Cálculo del ΔCt y ΔΔCt
def ct(df, control, housekeeping, target):
    ct_mean_housekeeping, ct_mean_target = ct_mean(df, control, housekeeping, target)
    
    # Crear una nueva columna ΔCt
    df['ΔCt'] = ct_mean_housekeeping - df['mir16']

    # Crear una nueva columna ΔΔCt
    df['ΔΔCt'] = 2 ** df['ΔCt']

    return df

print(ct(df_copy, control, housekeeping, target))
