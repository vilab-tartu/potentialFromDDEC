
########Import

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
from io import StringIO
from scipy.ndimage.filters import gaussian_filter

########Constants

k = 1.9872041E-3  # kcal/mol/K
rho = 0.03332819694513222
el_chg = 1.60217662e-19 #coulomb
Na = 6.02214129e23
epsilon0 =  8.854187817e-12  #F/m
eps0 = 0.0055263494 # el/V/A
kcal_eV = 2.611448e22  # kcal/(eV)
bohr2a = 0.529177 #A/bohr
D2Cm =  3.33564*10**-30 #C*m
au2D = 2.5417 #D

"""
A script that calculates the potential drop from the DDEC6 surface charges. 
You did similar calculation for the PRACE. Could you please create a script, analogues to the DDEC_XPS.py. Say DDEC_pot.py. 
"""



def process_command_line(argv):
    """Processes arguments

    Parameters
    ----------
    argv : list
        Command line arguments.

    Returns
    -------
    out : argparse.Namespace
        Namespace of command line arguments.
    """
    parser = argparse.ArgumentParser(description="""Create a potential drop and charge density figures in z direction from a
    given DDEC6 analysis input file.""")
    # DDEC filename
    parser.add_argument('-i',
                        help="""Input DDEC file for calculating.""")
    # Optional args
    parser.add_argument('-b',
                        help="""The bin size for z-axis creation [Angs]""")
    # Optional args
    parser.add_argument('-s',
                        help="""The index of first electrode atom .xyz coordinates""")

    # Optional args
    parser.add_argument('-e',
                        help="""The index of last electrode atom .xyz coordinates""")

    # Optional args
    parser.add_argument('-N',
                        help="""Sigma value for Gaussian filter""")
    return parser.parse_args(argv)

"""
    Method for testing giving filename instead of command line excecution
"""
class test:
    def __init__(self, filename, binSize, first, last, gaussianSigma):
        self.i = filename
        self.b = binSize
        self.s = first
        self.e = last
        self.N = gaussianSigma

"""
    Function for reading in DDEC6 Charge analysis results
    Arguments:
        args - parsed command line input
    Return:
        2D list [unit Cell in Angs, AtomType, Z-Coord in Angs, Atomic charge in e, dipole in z-dir Au]
"""
def readInData(args):
    file = open(args.i, 'r')
    line = file.readline() #First line
    # Remove unnecessary components, get nr of atoms
    nr_of_atoms = int(line.replace(' ', '').replace('\n', '').replace('\r', ''))

    #Get unit cell measurements in Angstroms [x,y,z]
    line = file.readline() #Read line with unit cell data
    line = line.strip("\n").replace("[", "").replace("]", "").replace("{", "").replace("}", "")
    line = line.split("unitcell")[1].split(",")
    unitCell = [0, 0, 0]
    for i in range(len(line)):
        unitCell[i] = float(line[i].split()[i])

    #Skip uninteresting lines
    line = ""
    while not line.startswith("atom number, atomic symbol, x, y, z, net_charge, dipole_x, dipole_y, dipole_z, dipole_mag"):
        line = file.readline()  # readline

    data = str()
    for i in range(nr_of_atoms):  # Read in rows of geometry
        line = file.readline()
        data += line

    data1 = StringIO(data)  # String to behave like IO object
    atomtype = np.loadtxt(data1, usecols=(1), dtype='str')  # Atom types
    data1 = StringIO(data)
    coordsChargeDipole = np.loadtxt(data1, usecols=(4,5,8))  # Coordinates and charges and z-dipoles

    return unitCell, atomtype, coordsChargeDipole[:, 0], coordsChargeDipole[:, 1], coordsChargeDipole[:, 2]



"""
    Method for calculating charge density from read in data
    
    Arguments:
        zCoords - z coordinates of atoms read in [Angs]
        charges - charges of atoms read in [e]
        unitCell - Lengths of unit cell [x,y,z] in [Angs]
        binSize - The size of bin to which charges are divided [Angs]
    Returns:
        [z_axis with step equal to binSize i [Angs],
         charge density in these bins [e/Angs^3]]
"""
def calculateChargeDensity(charges, unitCell, binSize):

    z_bin = np.arange(0, unitCell[2], binSize) #Create z axis
    charges_bin = [0] * len(z_bin)  #Create bins for charge density axis
    for l in range(len(zCoords)):
        #Separate charges into bins
        charges_bin[int(zCoords[l] // binSize)] += charges[l]
    #Divide binned charges by the volume of bin
    charges_bin = np.asarray(charges_bin) / (binSize * unitCell[0] * unitCell[1])
    return [z_bin, charges_bin]

"""
    Returns potential at point a, given a 1D-chg distribution (rho_z)
    and points at which it was evaluated (z_array) 
    
    Arguments:
        z_array - z-axis [Angs]
        rho_z - Charge density [e/Angs^3]

"""
def potential_at_a(a, z_array, rho_z):

    z_array = np.asarray(z_array)
    rho_z = np.asarray(rho_z)
    to_integrate = z_array < a
    dz = z_array[1] - z_array[0]
    integrand = 0
    for z_i, rho_zi in zip(z_array[to_integrate], rho_z[to_integrate]):
        integrand += (a - z_i)*rho_zi
    return -1/eps0*integrand*dz

"""
    Method for calculating the potential change in z-axis
    Arguments
        zCoords - z-axis with certain step [Angs]
        chargeDens - Charge density at a point given by z-axis [e/Angs^3]
        
    Returns:
        [zCoords [Angs],
         Potential change in Z axis [V]]
"""
def calculatePotential(zCoords, chargeDens):
    potZ = [potential_at_a(a, zCoords, chargeDens) for a in zCoords]
    return [zCoords, potZ]

"""
    Method for creating and saving Figure
    
    Arguments
        xAxis - x-axis to the plot
        Values - values in given points of x-axis
        xAxisName - The name on x-axis
        valueAxisName - The name of value (y) axis
        filename - The filename of created figure
"""
def createFigure(xAxis, Values, gaussianSigma, xAxisName, valueAxisName, filename):
    # Plotting and saving data
    f1, (ax) = plt.subplots(1, 1, sharey=False, sharex=True, figsize=(16.6 / 2.54, 8.3 / 2.54))
    f1.subplots_adjust(hspace=0)

    ax.plot(xAxis, gaussian_filter(Values, gaussianSigma), 'k')
    ax.tick_params(axis='both', labelsize=9)

    ax.set_xlabel(xAxisName, fontsize=12)
    ax.set_ylabel(valueAxisName, fontsize=12)

    f1.savefig(filename, format="svg" ,bbox_inches='tight')

"""
    Method for creating .csv file from data
    Arguments
        zCoords - z-axis with certain step [Angs]
        chargeDens - Charge density [e/Angs^3]
        potential - Potential change in Z axis
        filename - name of .csv file
"""
def dataToCsv(zCoords, chargeDens, potential, filename):
    df = pd.DataFrame({"z [Angs]":zCoords, "Charge Dens [e/Angs3]": chargeDens, "Potential [V]":potential})
    df.to_csv(filename, index=False)

"""
    Method to calculate dipole correction for 2D electrode
    Arguments
        args - parsed command line input
        unitCell - Lengths of unit cell [x,y,z] in [Angs]
        dipoles - dipoles in z direction [a.u.]
    :return
        Potential correction calculated from dipoles [V]
"""
def calculateDipoleCorr(args, unitCell, dipoles):
    startIndex = int(args.s)
    endIndex = int(args.e)+1
    return sum(dipoles[startIndex:endIndex])*au2D*D2Cm/(unitCell[0]*unitCell[1]*epsilon0*10**-20)

##For use without command line
#args = test("DDEC6_even_tempered_net_atomic_charges.xyz", "0.1",  "0", "447", "2")

############# Main excecution starts here ##########################


args = process_command_line(sys.argv[1:])
unitCell, atomtype, zCoords, charges, dipoles = readInData(args)
zBinned, chargeDens = calculateChargeDensity(charges, unitCell, float(args.b))
potential = calculatePotential(zBinned, chargeDens)[1]
print("Potential drop [V]: " + str(potential[0] - np.mean(potential[-int(len(potential)/100):-1])))
print("Dipole correction [V]: " + str(calculateDipoleCorr(args, unitCell, dipoles)))
createFigure(zBinned, chargeDens, float(args.N), r"z [$\AA$]", r"$\rho$ [$e/\AA^3$]", "chargeDens.svg")
createFigure(zBinned, potential, float(args.N), r"z [$\AA$]", r"$\phi$ [$V$]", "potential.svg")
dataToCsv(zBinned, chargeDens, potential, "data.csv")