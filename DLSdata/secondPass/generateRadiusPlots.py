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

def makeRadiusArray(drugStr, pHInt):
    """
    :param drugStr: 3 letter drug abbrev str
    :returns: array of radius values at pH 4 (m x n, where m = 3 data points at single time point, n = time point)
    """
    drugDict = {'sfb': 0, 'nil': 1, 'ptx': 2, 'trm': 3, 'dab': 4, 'pbs': 5}
    if pHInt == 4:
        array = [df_1h.loc[drugDict[drugStr], '1':'3'],
                 df_2h.loc[drugDict[drugStr], '1':'3'],
                 df_25h.loc[drugDict[drugStr], '1':'3'],
                 df_49h.loc[drugDict[drugStr], '1':'3'],
                 df_73h.loc[drugDict[drugStr], '1':'3']]
    if pHInt == 5:
        array = [df_1h.loc[drugDict[drugStr], '4':'6'],
                 df_2h.loc[drugDict[drugStr], '4':'6'],
                 df_25h.loc[drugDict[drugStr], '4':'6'],
                 df_49h.loc[drugDict[drugStr], '4':'6'],
                 df_73h.loc[drugDict[drugStr], '4':'6']]
    if pHInt == 6:
        array = [df_1h.loc[drugDict[drugStr], '7':'9'],
                 df_2h.loc[drugDict[drugStr], '7':'9'],
                 df_25h.loc[drugDict[drugStr], '7':'9'],
                 df_49h.loc[drugDict[drugStr], '7':'9'],
                 df_73h.loc[drugDict[drugStr], '7':'9']]
    if pHInt == 7:
        array = [df_1h.loc[drugDict[drugStr], '10':'12'],
                 df_2h.loc[drugDict[drugStr], '10':'12'],
                 df_25h.loc[drugDict[drugStr], '10':'12'],
                 df_49h.loc[drugDict[drugStr], '10':'12'],
                 df_73h.loc[drugDict[drugStr], '10':'12']]
    if pHInt == 8:
        array = [df_1h.loc[drugDict[drugStr], '13':'15'],
                 df_2h.loc[drugDict[drugStr], '13':'15'],
                 df_25h.loc[drugDict[drugStr], '13':'15'],
                 df_49h.loc[drugDict[drugStr], '13':'15'],
                 df_73h.loc[drugDict[drugStr], '13':'15']]
    return array

def plotRadiusOverTime(drugAbbrev):
    # pH 4
    plt.plot([1, 1, 1], makeRadiusArray(drugAbbrev, 4)[0], 'bo',
             [2, 2, 2], makeRadiusArray(drugAbbrev, 4)[1], 'bo',
             [25, 25, 25], makeRadiusArray(drugAbbrev, 4)[2], 'bo',
             [49, 49, 49], makeRadiusArray(drugAbbrev, 4)[3], 'bo',
             [73, 73, 73], makeRadiusArray(drugAbbrev, 4)[4], 'bo')
    plt.xlabel('time (h)')
    plt.ylabel('Radius (nm)')
    plt.show()

