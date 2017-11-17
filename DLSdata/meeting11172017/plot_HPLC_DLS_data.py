from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
import sys

df = pd.read_csv('../DLSdata/meeting11172017/DLS+HPLC_data_organized.csv')

#plot Mass Feed ratio vs. Diameter
mass_feed_ratios = df.iloc[:, 1]
actual_mass_ratio = df.iloc[:,2]
percent_mass_drug = df.iloc[:,3]
avg_Diam = df.iloc[:, 4]
avg_PDI = df.iloc[:, 5]

samples = ['Constant drug feed mass', 'Constant dye feed mass']
sns.set_style('darkgrid')

# plt.plot(mass_feed_ratios[0:4], avg_Diam[0:4], 'bo',
#          mass_feed_ratios[4:9], avg_Diam[4:9], 'go')
# plt.legend(samples)
# plt.xlabel("Mass Feed Ratio (mg drug/mg dye)")
# plt.ylabel("Average Diameter (nm)")
# plt.title('Sorafenib/IR820 NP')
# plt.show()

# plt.plot(actual_mass_ratio[0:4], avg_Diam[0:4], 'bo',
#          actual_mass_ratio[4:9], avg_Diam[4:9], 'go')
# plt.legend(samples)
# plt.xlabel("Actual Mass Ratio (mg drug/mg dye)")
# plt.ylabel("Average Diameter (nm)")
# plt.title('Sorafenib/IR820 NP')
# plt.show()

plt.plot(mass_feed_ratios[0:4], actual_mass_ratio[0:4], 'bo',
        mass_feed_ratios[4:9], actual_mass_ratio[4:9], 'go')
plt.legend(samples)
plt.xlabel("Mass Feed Ratio (mg drug/mg dye)")
plt.ylabel("Actual Mass Ratio (mg drug/mg dye)")
plt.title('HPLC results vs initial feed')
plt.show()

