##!/usr/bin/env python3.5
# This python script is inspired by pipeMeshNek.f90 from Jacopo Canton.
# However, I am not that familiar with Fortran 90 and I wanted to 
# improve my python skills and understand the mesh generation better.
# Hence, I rewrite the code from scratch
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# Steffen Straub
# 2017-02-08
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

## Import modules
#----------------------------------------------------------------------
import nek_utils
#import numpy as np
import elementmodule 
import pdb
import sys
import re


## Generate a rea skeleton file
#----------------------------------------------------------------------
nek_utils.rea_skel()

# Create the mesh
#----------------------------------------------------------------------

# Input Variables
#----------------------------------------------------------------------
R = 0.5         # radius
nR = 28         # nel in radial direction
nPhi = 152       # nel in circumferential direction

# Calculations
#----------------------------------------------------------------------
nSq = nPhi/8    # nel in square region along one side of the square
dx = R/nSq       # length of one element
dy = R/nSq       # height of one element
n_tot = nR * nPhi   # total number of elements
spatial_dim = 2     # spatial dimension

# Do some checks
if ( (nPhi % 8) != 0 ):
    print("Number of circumferential elements needs to be a multiple of 8")
    sys.exit(1)
else:
    nSq = int(nSq)
    


# Create square region
#----------------------------------------------------------------------
number = 1  # number of the elements
elem_sq = list()    # collect all elements in a list
for i in range(nSq):    # go through each column (x direction)
    for j in range(nSq):    # go through each row (y direction)
        element = elementmodule.Element()   # instantiate an element object
        element.x = [i*dx, (i+1)*dx, (i+1)*dx, i*dx]  # set vertex coordinates
        element.y = [i*dy, i*dy, (i+1)*dy, (i+1)*dy]  # set vertex coordinates
        element.number = number 
        elem_sq.append(element)
        number = number + 1

    
# Write the mesh part into rea file
#----------------------------------------------------------------------
nek_utils.write_mesh(nSq**2, 2, elem_sq)




# Set the boundary conditions
#----------------------------------------------------------------------
for el in elem_sq:
#    print(el.number)
    if (el.number == 1):  # we are on the first element on the lower left
        el.bc = ['W  ','E  ','E  ','W  ']
        el.bc_params = [0, el.number+1, el.number+nSq, 0]
    elif (el.number == nSq):    # we are on the lower right corner
        el.bc = ['W  ','E  ','E  ','E  ']
        el.bc_params = [0, el.number+1, el.number+nSq, el.number-1] 
    elif (el.number == nSq**2-nSq+1):     # we are on the upper left corner
        el.bc = ['E  ','E  ','W  ','W  ']
        el.bc_params = [el.number-nSq, el.number+1, 0, 0]
    elif (el.number == nSq**2):   # last element on the upper right
        el.bc = ['E  ','W  ','W  ','E  ']
        el.bc_params = [el.number-nSq, 0, 0, el.number-1]
    elif (el.number < nSq):  # southern row
        el.bc = ['W  ','E  ','E  ','E  ']
        el.bc_params = [0, el.number+1, el.number+nSq, el.number-1]
    elif ((el.number%nSq+1) == 0):  #  western column
        el.bc = ['E  ','E  ','E  ','W  ']
        el.bc_params = [el.number-nSq, el.number+1, el.number+nSq, 0]
    elif (el.number > nSq**2-nSq):   # northern row
        el.bc = ['E  ','E  ','W  ','E  ']
        el.bc_params = [el.number-nSq, el.number+1, 0, el.number-1]
    elif ((el.number%nSq) == 0): # eastern column
        el.bc = ['E  ','W  ','E  ','E  ']
        el.bc_params = [el.number-nSq, 0, el.number+nSq, el.number-1]
    else: # inside
        el.bc = ['E  ','E  ','E  ','E  ']
        el.bc_params = [el.number-nSq, el.number+1, el.number+nSq, el.number-1]








## Write the boundary conditions
#f = open('base.rea','r')
#contents = f.readlines()
#f.close()
#
## find row for fluid boundary
#line_bc = 0
#while (not 'FLUID BOUNDARY' in contents[line_bc]):
#    line_bc = line_bc + 1
#
#bc = list()
#for elem_number in range(n_tot):    # loop through all elements
#    bc.append({bc_type: cbc
#
