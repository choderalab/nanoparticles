from matplotlib import pyplot as plt
import pandas as pd
"""Script to generate plots from .csv's"""

# Experiment variables
pH = [4, 5, 6, 7, 8]
time = [1,2,25,49,73]
drugs = ['sfb', 'nil', 'ptx', 'trm', 'dab', 'pbs']
pKa = [11.55, 11.86, 10.36, 12.6, 7.16, 7.0] # in order corresponding to drugs above; values from DrugBank

df_1h = pd.read_csv('../DLSdata/Radius/t1h_avgRadius.csv')
df_2h = pd.read_csv('../DLSdata/Radius/t2h_avgRadius.csv')
df_25h = pd.read_csv('../DLSdata/Radius/t25h_avgRadius.csv')
df_49h = pd.read_csv('../DLSdata/Radius/t49h_avgRadius.csv')
df_73h = pd.read_csv('../DLSdata/Radius/t73h_avgRadius.csv')


def makeArray(i, j): #i: drugIndex, j: pH index
    array = [df_1h.iloc[i,j],
                       df_2h.iloc[i,j],
                       df_25h.iloc[i,j],
                       df_49h.iloc[i,j],
                       df_73h.iloc[i,j]]
    return array

print(df_73h)
print(df_1h)
radDiff = ((df_73h - df_1h)/df_1h)*100 #size change as percentage


print(radDiff)

plt.plot(pKa, radDiff.iloc[0:6,1], 'bo', pKa,
         radDiff.iloc[0:6,2], 'go', pKa,
         radDiff.iloc[0:6,3], 'ro', pKa,
         radDiff.iloc[0:6,4], 'ko',pKa,
         radDiff.iloc[0:6,5], 'mo')
plt.legend(pH)
plt.xlabel("Drug pKa")
plt.ylabel("Percent Size Change")
plt.title('Size Change vs. Drug pKa')
plt.show()
