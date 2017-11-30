from matplotlib import pyplot as plt
import pandas as pd
"""
Script to generate radius, PDI, and amplitude plots from .csv's generated with rejection criteria
-Wells with <10 acquisitions were rejected
-All 3 measurements will be plotted if available
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

def makeRadiusArray(drugStr):
    """
    :param drugStr: 3 letter drug abbrev str
    :returns: array
    """
    drugDict= {'sfb': 0, 'nil': 1, 'ptx': 2, 'trm': 3, 'dab': 4, 'pbs': 5}
    array = [df_1h.loc[drugDict[drugStr], '1':'3'],
             df_2h.loc[drugDict[drugStr], '1':'3'],
             df_25h.loc[drugDict[drugStr], '1':'3'],
             df_49h.loc[drugDict[drugStr], '1':'3'],
             df_73h.loc[drugDict[drugStr], '1':'3']]
    return array

# Radius vs time plots for each drug

print(makeRadiusArray('sfb'))
plt.plot([1, 1, 1], makeRadiusArray('sfb')[0], 'bo')
plt.show()