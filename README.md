# potentialFromDDEC
A script that calculates the potential drop, dipole correction and creates figures from the DDEC6 surface charges

Example command line use, when filename is "DDEC6_even_tempered_net_atomic_charges.xyz", the step of z-axis in figures is 0.1 Angstroms, the index on first electrode atom is 0 and the last is 448. Sigma for gaussian function is 2 
python potentialFromDDEC.py -i DDEC6_even_tempered_net_atomic_charges.xyz -b 0.1 -s 0 -e 448 -N 2

### Example

#### COMMAND TO RUN PROGRAM

##### Windows:
```
python potentialFromDDEC.py -i DDEC6_even_tempered_net_atomic_charges.xyz -b 0.1 -s 0 -e 448 -N 2 > output.txt
```

##### Linux:
```
python potentialFromDDEC.py -i DDEC6_even_tempered_net_atomic_charges.xyz -b 0.1 -s 0 -e 448 -N 2 > output.txt
```

#### OUTPUT OF EXAMPLE DDEC RESULTS

##### (Printed to terminal)
Potential drop [V]: 1.8434295130861662 \
Dipole correction [V]: 0.5710576676392137

##### (Created files)
* potential.svg - Potential profile in z direction (without dipole correction)
* chargeDens.svg - Charge density profile in z direction
* data.csv - .csv file containing data to 
