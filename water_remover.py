#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog = 'water_remover', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
************************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/water_remover
DOI: 
************************************************

DESCRIPTION:
 This script removes water particles within the membrane following solvatation.

REQUIREMENTS:
 MDAnalysis, numpy
 
USAGE:

Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-b 		['PO4']	: name of beads delimiting the membrane boundaries
-d		[5]	: distance from membranre boundaries along z axis (Angstrom)
-w		['W']	: resname of water molecules
 
Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-b', nargs=1, dest='bead_name', default=['PO4'], help=argparse.SUPPRESS)
parser.add_argument('-d', nargs=1, dest='distance', default=[5], type=float, help=argparse.SUPPRESS)
parser.add_argument('-w', nargs=1, dest='w_resname', default=['W'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.grofilename = args.grofilename[0]
args.bead_name = args.bead_name[0]
args.distance = args.distance[0]
args.w_resname = args.w_resname[0]

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================
#generic science modules
try:
	import numpy as np
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.density
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=========================================================================================
# sanity check
#=========================================================================================
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#load universe
U = Universe(args.grofilename)
	
#select leaflets
bead_sele = U.selectAtoms("name " + args.bead_name)
if bead_sele.numberOfAtoms() == 0:
	print "Error: no atoms " + str(args.bead_name) + " found."
	sys.exit(1)
bead_z_avg = np.average(bead_sele.coordinates()[:,2])
bead_upper = bead_sele.selectAtoms("prop z > " + str(bead_z_avg))
bead_lower = bead_sele.selectAtoms("prop z < " + str(bead_z_avg))
bead_upper_z = np.average(bead_upper.coordinates()[:,2])
bead_lower_z = np.average(bead_lower.coordinates()[:,2])

#select waters
water_sele = U.selectAtoms("resname " + str(args.w_resname))
water_upper = water_sele.selectAtoms("prop z > " + str(bead_upper_z + args.distance))
water_lower = water_sele.selectAtoms("prop z < " + str(bead_lower_z - args.distance))
water_to_keep = water_upper + water_lower

#write file
not_water = U.selectAtoms("not resname " + str(args.w_resname))
to_keep = not_water + water_to_keep
to_keep.write(args.grofilename[:-4] + "_rm.gro")
to_keep.write(args.grofilename[:-4] + "_rm.pdb")

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check output in " + args.grofilename[:-4] + "_rm.gro file"
print ""
sys.exit(0)
