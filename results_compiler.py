'''
Results appender and extra calculations
'''
import pandas as pd
import os
import CoolProp.CoolProp
from CoolProp.CoolProp import PropsSI

# get wd to find folders
folder = os.getcwd()

# well power for plateau rate
HHV_H2 = 39.4  # kWh/kg
HHV_CH4 = 15.4  # kWh/kg
Tsc = 491.67  # temperature [Rankine]
Psc = 14.5038  # pressure [psia]

psi_to_pascal = 6894.76  # multiply for pascal / divide for psi
PaS_to_cP = 1000  # multiply for cP / divide for Pa S
kg_mol_to_lb_lbmol = 1000  # multiply for lb/lbmol / divide for kg/mol
kg_m3_to_lb_ft3 = 0.062428  # multiply for lb/ft3 / divide for kg/m3
secs_to_days = 86400  # multiply for days / divide for seconds
m3_to_ft3 = 35.3147  # multiply for feet / divide for metres
inches_to_metres = 0.0254
metres_to_feet = 3.28084  # multiply for feet / divide for metres
rankine_to_kelvin = 5 / 9  # multiply for kelvin / divide for rankine

def well_power_MW(gas, q):
    if gas == 'H2':
        return HHV_H2 * (q / m3_to_ft3 * 1000000) / 24 * PropsSI('D', 'T', Tsc * rankine_to_kelvin, 'P',
                                                                              Psc * psi_to_pascal, gas) / 1000
    elif gas == 'CH4':
        return HHV_CH4 * (q / m3_to_ft3 * 1000000) / 24 * PropsSI('D', 'T', Tsc * rankine_to_kelvin, 'P',
                                                                               Psc * psi_to_pascal, gas) / 1000
    else:
        return 0

li = []
for root, dirs, files in os.walk(folder + "\OUTPUT", topdown=False):
    for name in files:
        if "storage_props" in name:
            print(os.path.join(root, name))
            x = os.path.join(root, name)
            li.append(x)
print(li)

aggregated_results = pd.DataFrame()

for i in li:
    df = pd.read_csv(i, index_col=0)
    print(df)
    df['root'] = i
    #aggregated_results = aggregated_results.append(df.iloc[0])[df.columns.tolist()]
    aggregated_results = pd.concat([aggregated_results,df], ignore_index=True)
print(aggregated_results)

aggregated_results['well power plateau [MW]'] = aggregated_results.apply(lambda x:well_power_MW(x['gas'], x['plateau rate [MMSCF/D]']), axis = 1)

aggregated_results['per well power plateau [MW]'] = aggregated_results['well power plateau [MW]'] / aggregated_results['n_wells']

aggregated_results['per well plateau rate [MMSCF/D]'] = aggregated_results['plateau rate [MMSCF/D]'] / aggregated_results['n_wells']

aggregated_results.to_csv('results_appended.csv', index = False)
print(aggregated_results)
# END OF CODE
