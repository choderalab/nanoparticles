from matplotlib import pyplot as plt
import pandas as pd
from scipy.stats import chisquare
import seaborn as sns
import sys

"""Script to generate plots to compare correlogram data of different samples

CSV Prep
--------
1) CSV files exported from Wyatt must have "μ" character removed from A1 and replaced with "u", 
otherwise a different encoding mode must be used with pd.read_csv (which may ruin column names)
2) Files must contain 91 columns (90 wells + 1 time column). 
Wyatt exported CSV's were found to contain some missing columns (due to incomplete measurements),
and had the previous timepoint's reading prepended to it. 
3) " (Incomplete)" must be removed from column headers until I find a way to use wildcards for accessing df
"""
# Experiment parameters
pH = [4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8]
timepoint = {'1':'t1h', '2':'t2h', '25':'t25h', '49':'t49h', '73':'t73h'}
drugDict = {'A': 'Sorafenib', 'B': 'Nilotinib', 'C': 'Paclitaxel',
            'D': 'Trametinib', 'E': 'Dabrafenib', 'G': 'PBS'}

def plotCorrelogram(drugRow, hour):
    sns.set_style('darkgrid')
    paletteDict = sns.color_palette("bwr", n_colors=5)
    graphPalette = [paletteDict[4], # color for pH 4 trace
                    paletteDict[3], # pH 5
                    paletteDict[2], # pH 6
                    paletteDict[1], # pH 7
                    paletteDict[0]] # pH 8

    #Reading .csv
    df = pd.read_csv('../DLSData/CorrelogramData/'+timepoint[hour]+'.csv')
    x = df['Time (us)']
    if (len(df.columns)) is not 91:
        print("Reformatting necessary: The .csv has {0} columns when it should have 91!".format(len(df.columns)))
        sys.exit()
    i = 0
    for j in range(1,16):
        if j in [4,7,10,13]: # first well with a different pH so that replicate sample traces are same color
            i += 1
        plt.plot(x, df[drugRow+str(j)], color = graphPalette[i])

    plt.semilogx()
    plt.xlabel("Time (μs)")
    plt.ylabel("ACF")
    plt.title(drugDict[drugRow]+', time = '+str(hour)+' hours')
    plt.legend(pH, title='pH', loc=1)
    plt.savefig('ACF_'+drugDict[drugRow]+'_'+hour+'h.png')
    plt.clf()

def generateCorrelograms():
    """Main function to generate ACF plots for each drug at each timepoint"""
    rows = ['A', 'B', 'C', 'D', 'E', 'G']
    times = ['1', '2', '25', '49', '73']

    for i in rows:
        for j in times:
            plotCorrelogram(i, j)

generateCorrelograms()