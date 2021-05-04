# potentialFromDDEC
A script that calculates the potential drop, dipole correction and creates figures from the DDEC6 surface charges

### Parameters

*  -i - Input DDEC file for calculating.
* -b - The bin size for z-axis creation [Angs]
* -s - The index of first electrode atom .xyz coordinates
* -e - The index of last electrode atom .xyz coordinates
* -N - Sigma value for Gaussian filter
* -area - Electrode area [Angs2] (If not given, the electrode area is estimated as rectangle unit cell x * unit cell y)

### Example

#### COMMAND TO RUN PROGRAM

```
python potentialFromDDEC.py -i DDEC6_even_tempered_net_atomic_charges.xyz -b 0.1 -s 0 -e 448 -N 2
```

#### OUTPUT OF EXAMPLE DDEC RESULTS

##### (Printed to terminal)
Potential drop [V]: 1.8434295130861662 \
Dipole correction [V]: 0.5710576676392137

##### (Created files)
* potential.svg - Potential profile in z direction (without dipole correction)
* chargeDens.svg - Charge density profile in z direction
* data.csv - .csv file containing data to 
