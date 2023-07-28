'''Input file generator
'''
import pandas as pd
import numpy as np
import os
from functools import reduce

#print(os.getcwd())

# this code generates the input files needed for the inflow/outflow model
# change folder to working folder for your project
folder = os.getcwd()
print(folder)

# take original input data and add gas column
df = pd.read_csv('gas_fields_props.csv')
gas = pd.DataFrame(['H2', 'CH4'], columns = ['gas']) # change gases you want to investigate here
dataframes = [df, gas]
df = reduce(lambda left, right: pd.merge(left, right, how='cross'), dataframes)
#print(df)
print(df.columns)

# take each row and make it a new df with a range of pressures
z = 1
for index, row in df.iterrows():
    x = str(z) + '_timesteps.csv'
    b = df.iloc[[index]]
    p = b['Pr_psia'].iloc[0]
    pressures = np.arange(100, p + 1, 100)
    pressures = pressures.tolist()
    b.drop(columns = ['Pr_psia'], inplace=True)
    print(pressures)
    c = pd.DataFrame(pressures,columns = ['Pr_psia'])
    dataframes = [b, c]
    a = reduce(lambda left, right: pd.merge(left, right, how='cross'), dataframes)
    a = a.sort_values(by = 'Pr_psia', ascending=False)
    a.to_csv(x)
    z = z + 1

# make a list of all the timesteps files to refer to
li = []
for root, dirs, files in os.walk(folder, topdown=False):
    for name in files:
        if "_timesteps" in name:
            print(os.path.join(root, name))
            x = os.path.join(root, name)
            li.append(x)

len(li)
li2 = []

wellbore_diams_inches = np.arange(2, 10 + 1, 2)
wellbore_diams_inches = wellbore_diams_inches.tolist()
for i in li:
    # reads file to get input data
    gas_fields_props = pd.read_csv(i)
    well_diam_variation = pd.DataFrame(wellbore_diams_inches, columns = ['D_wb_inches']) # diameter of the wellbore in inches pd.DataFrame(['H2', 'CH4'], columns = ['gas'])
    dataframes = [gas_fields_props, well_diam_variation] # list the two dataframes
    df = reduce(lambda left, right: pd.merge(left, right, how='cross'), dataframes) # combine them
    print(df)
    #exit()
    for i in wellbore_diams_inches:
        df2 = df[df['D_wb_inches'] == i].reset_index()
        df2['rw_ft'] = df2['D_wb_inches'] / 12 / 2
        y = str(df2['field_name'].iloc[0]) + '_' + df2['gas'].iloc[0] + '_timesteps_' + str(
            df2['D_wb_inches'].iloc[0]) + '_inches.csv'
        df2.to_csv(y)