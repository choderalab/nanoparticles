import pandas as pd
import numpy as np

"""
Script to import Wyatt DLS data such as in 20170707_MI_INP_DLS_1s_10aq_exp_t1h_t2h_all_individual_acquisitions.xlsx

  Plots to generate
    -Avg radius with error bars over time at each pH
    -3D plot: R = f(time, pH) for each drug
    -try different descriptors: MW, pKa (ex. dR = f(MW, pH) )
    -Amplitude = f(time, pH)
    
Prep .xlsx file by:
-Ctrl+F and replace all "--" cells with blanks
-Make sure Radius(nm) is added to the bottom of column B when importing Radius values

This script uses accept/reject criteria agreed on by Mehtap and Zaki on 11/15/2017:
1. reject wells with <10 acquisitions
2. 
"""
#TODO make generalized functions for script that can be run by specifying filename and measurement type

# Setting up plate template that will be filled in
data = {'Row Label': ['A', 'B', 'C', 'D', 'E', 'G'],
        'Drug': ['sfb', 'nil', 'ptx', 'trm', 'dab', 'pbs']
        }
plateTemplate = pd.DataFrame(data,
                             columns=['Drug', 'Row Label', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
                                      '12', '13', '14', '15'])
print('Empty plate template:', plateTemplate)


def importExcelSheet(filename, sheetname, pandasHeader, measurementColumn):
    """
    Imports Excel column into list. Need to change the column heading and dataframe label each time.
    MAKE SURE MEASUREMENT HEADER IS AT BOTTOM OF EXCEL SHEET ALSO
    :param: Filename of excel in current directory in quotes
    """
    global columnVector
    global measurementPositions
    sheet = pd.read_excel(filename, sheetname) # import specified sheet as pandas dataframe
    columnVector = sheet[pandasHeader]  # place column E (Amplitude measurements) into list
    print(columnVector)
    measurementPositions = []
    for i, j in enumerate(columnVector):
        if j == measurementColumn: # name of measurement column # TODO make this generalizable for any measurement type
            measurementPositions.append(i)
    print("The index of each well measurement is at ", measurementPositions)
    print("CHECK: This value should be 91:",
          len(measurementPositions))  # 15 columns and 6 samples on plate + 1 to make sample computable


def averageWellAcquisitions(well_index):
    """
    This function returns the average well value for the measurement specified by columnVector, for the
    well specified by well_index.

    Modified on 11/22/2017 to only calculate averages for wells with 10 acquisitions (otherwise NaN).

    Parameters
    ----------
    well_index : int
        Index for list of wells on plate as flattened array given by the order in the excel sheet (column-major/Fortran-style)
        ex. for A1, well_index = 0
    Returns
    -------
    wellAvg : float
    """
    print('The measurement positions are: ', measurementPositions[well_index] + 1, measurementPositions[well_index + 1])
    print('The values being read are\n', columnVector[measurementPositions[well_index] + 1:measurementPositions[well_index + 1]])
    floatCount = 0
    for element in columnVector[measurementPositions[well_index] + 1:measurementPositions[well_index + 1]]:
        if np.isnan(element) == False:
            floatCount += 1
    print('The float count is: ', floatCount)
    if floatCount == 10:
        wellAvg = columnVector[measurementPositions[well_index] + 1:measurementPositions[well_index + 1]]
        wellAvg = np.nanmean(wellAvg.astype('float64'))
        print('The average measurement for this well is:', wellAvg)
        return wellAvg
    else:
        return None


def FillPlate():
    """
    Fills plate with measurement values calculated by averageWellAcquisitions
    """
    t1h_avgRadius = plateTemplate
    index = 0
    i = 0 # row index for data frame
    j = 2 # column index for data frame
    for i in range(0, 6):
        print("i is:", i)
        for j in range(2, 17):
            print("--j is:", j)
            t1h_avgRadius.iloc[i, j] = averageWellAcquisitions(index)
            index = index + 1
            print("The new well value is:", t1h_avgRadius.iloc[i, j])
            print("The value of index is:", index)  # this should be 90 at the end of the loop
    print(t1h_avgRadius) #TODO write to .csv
    return t1h_avgRadius

# run script to export data to .csv is below
# note that in the functions, lines 34, 35, 39 must be changed every time
# in the script below, lines 92, 108, and 109 must be changed every time
def firstPass():
    importExcelSheet('20170707_MI_INP_DLS_1s_10aq_exp_t73h_all_individual_acquisitions.xlsx')
    x = FillPlate()
    # Now that we have the plate filled, we need to average the three measurements at same pH
    pH_summary = [x.iloc[:, 2:5].mean(axis=1),
              x.iloc[:, 5:8].mean(axis=1),
              x.iloc[:, 8:11].mean(axis=1),
              x.iloc[:, 11:14].mean(axis=1),
              x.iloc[:, 14:17].mean(axis=1)]
    pH_errSummary = [x.iloc[:, 2:5].std(axis=1),
                 x.iloc[:, 5:8].std(axis=1),
                 x.iloc[:, 8:11].std(axis=1),
                 x.iloc[:, 11:14].std(axis=1),
                 x.iloc[:, 14:17].std(axis=1)]
    tXh_pH_averaged = pd.concat(pH_summary, axis=1) # dataframe with pH columns averaged for 3 well measurements
    tXh_pH_err = pd.concat(pH_errSummary, axis=1) # standard deviation between 3 well measurements

    tXh_pH_averaged.to_csv('t73h_avgPDI.csv')
    tXh_pH_err.to_csv('t73h_errPDI.csv')

    # Print statements
    print(tXh_pH_averaged)
    print(tXh_pH_err)

importExcelSheet('20170707_MI_INP_DLS_1s_10aq_exp_t73h_all_individual_acquisitions.xlsx', 'Sheet1', 'Unnamed: 10', 'PD Index')
x = FillPlate()
pH_summary = [x.iloc[:, 2:5].mean(axis=1),
              x.iloc[:, 5:8].mean(axis=1),
              x.iloc[:, 8:11].mean(axis=1),
              x.iloc[:, 11:14].mean(axis=1),
              x.iloc[:, 14:17].mean(axis=1)]
pH_errSummary = [x.iloc[:, 2:5].std(axis=1),
                 x.iloc[:, 5:8].std(axis=1),
                 x.iloc[:, 8:11].std(axis=1),
                 x.iloc[:, 11:14].std(axis=1),
                 x.iloc[:, 14:17].std(axis=1)]
tXh_pH_averaged = pd.concat(pH_summary, axis=1)  # dataframe with pH columns averaged for 3 well measurements
tXh_pH_err = pd.concat(pH_errSummary, axis=1)  # standard deviation between 3 well measurements

x.to_csv('secondPass/PDI/t73h_allWellsPDI.csv')
# tXh_pH_averaged.to_csv('secondPass/Radius/t1h_avgRadius.csv')
# tXh_pH_err.to_csv('secondPass/Radius/t1h_errRadius.csv')

print(pH_errSummary)
print(pH_errSummary)
# Print statements
print(tXh_pH_averaged)
print(tXh_pH_err)