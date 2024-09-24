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
file = '/Users/andreaross/Desktop/bioinfo_andrea/ds.xlsx' # Ruta del archivo con los datos de expresión génica
df = pd.read_csv(file)
df.head()
