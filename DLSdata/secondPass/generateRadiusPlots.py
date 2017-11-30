from matplotlib import pyplot as plt
import pandas as pd
"""
Script to generate radius, PDI, and amplitude plots from .csv's generated with rejection criteria
-Wells with <10 acquisitions were rejected
"""

# Experiment variables
pH = [4, 5, 6, 7, 8]
time = [1, 2, 25, 49, 73]
drugs = ['sfb', 'nil', 'ptx', 'trm', 'dab', 'pbs']

#initiliaze dataframes for each hour
df_1h = pd.read_csv('../DLSdata/secondPass/Radius/t1h_allWellsRadius.csv')
df_2h = pd.read_csv('../DLSdata/secondPass/Radius/t2h_allWellsRadius.csv')
df_25h = pd.read_csv('../DLSdata/secondPass/Radius/t25h_allWellsRadius.csv')
df_49h = pd.read_csv('../DLSdata/secondPass/Radius/t49h_allWellsRadius.csv')
df_73h = pd.read_csv('../DLSdata/secondPass/Radius/t73h_allWellsRadius.csv')

# Radius vs time plots for each drug
# plot(time, pH4_R, time, pH5_R, pH_6

print(df_1h.iloc[1,2:17])