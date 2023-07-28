# cushion_gas_calculator
Python code to calculate cushion gas requirement for closed gas storage reservoirs

Paper with explanation is here: https://doi.org/10.1144/SP528-2022-71

Cite as: Scafidi, J., Schirrer, L., Vervoort, I., &#38; Heinemann, N. (2023). An open-source tool for the calculation of field deliverability and cushion-gas requirements in volumetric gas-reservoir storage sites. Geological Society, London, Special Publications, 528(1). https://doi.org/10.1144/SP528-2022-71

Instructions:

1. Update the files gas_fields_props.csv with your gas field's properties. You can add multiple lines here
2. Update the file n_wells_single.csv with the names of the fields in gas_fields_props.csv and the number of wells you want to investigate
3. Run file_generator.py
4. Delete any generated files that do not contain the name of a gas field in their name e.g. "1_timesteps.csv"
5. Create two folders in the directory: one called "INPUT" and one called "OUTPUT"
6. Move the generated files with the word "timesteps" in their name into the folder called "INPUT"
7. Move the n_wells_single.csv file into the folder called "INPUT"
8. Run main.py - the results will appear in the folder named OUTPUT
9. Run results_compiler.py
10. If you want to redo results the OUTPUT folder must be empty so copy any results you wish to keep elsewhere
