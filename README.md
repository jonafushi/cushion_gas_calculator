# cushion_gas_calculator
Python code to calculate cushion gas requirement for closed gas storage reservoirs

Paper with explanation is here: https://doi.org/10.1144/SP528-2022-71

Cite as: Scafidi, J., Schirrer, L., Vervoort, I., &#38; Heinemann, N. (2023). An open-source tool for the calculation of field deliverability and cushion-gas requirements in volumetric gas-reservoir storage sites. Geological Society, London, Special Publications, 528(1). https://doi.org/10.1144/SP528-2022-71

Instructions:

1. Update the files gas_fields_props.csv with your gas field's properties. You can add multiple lines here
2. Update the file n_wells_single.csv with the names of the fields in gas_fields_props.csv and the number of wells you want to investigate
3. Run file_generator.py
4. Delete any generated files that do not contain the name of a gas field
5. Move the generated files into a folder called "INPUT"
6. Run main.py - the results will appear in a folder named OUTPUT
7. Run results_compiler.py
8. If you want to redo results the OUTPUT folder must be empty so copy any results you wish to keep elsewhere
