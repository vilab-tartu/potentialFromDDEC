# potentialFromDDEC
A script that calculates the potential drop, dipole correction and creates figures from the DDEC6 surface charges

Example command line use, when filename is "DDEC6_even_tempered_net_atomic_charges.xyz", the step of z-axis in figures is 0.1 Angstroms, the index on first electrode atom is 0 and the last is 448. Sigma for gaussian function is 2 
python potentialFromDDEC.py -i DDEC6_even_tempered_net_atomic_charges.xyz -b 0.1 -s 0 -e 448 -N 2

