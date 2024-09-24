# Prueba estadística para conocer la sigificancia de resultados en análisis de expresión
Script para el análisis de datos tras la cuantificación de la expresión génica
mediante qPCR. Se analizan los datos de expresión génica de los genes de interés
Se calcula la expresión relativa de los genes de interés siguiendo el método 2^-ΔΔCt
propuesto por Taylor et al. (2019).

El archivo de entrada debe ser un archivo csv con los valores de CT mean tras las corridas de amplificación por tiempo real (veáse archivo muestra: ds.csv)

1. Reemplazar la ruta en la cual se encuentre el archivo csv del cual se extraerán los datos
2. Reemplazar el gen target y housekeeping a analizar (mismo nombre de la columna donde se encuentran los CTs)
3. Correr en la terminal como  $ python expresion.py
