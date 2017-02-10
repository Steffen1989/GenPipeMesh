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
#import elementmodule 
import pdb
import sys
import re

# Input Variables
#----------------------------------------------------------------------
R = 0.5         # radius
#nR = 28         # nel in radial direction
nR = 8
nPhi = 64       # nel in circumferential direction

# Define global variables here:
#----------------------------------------------------------------------
nSq = nPhi/8    # nel in square region along one side of the square
if ( (nPhi % 8) != 0 ):
    print("Number of circumferential elements needs to be a multiple of 8")
    sys.exit(1)
else:
    nSq = int(nSq)
dx = R/nSq       # length of one element
dy = R/nSq       # height of one element
n_tot = nR * nPhi   # total number of elements
spatial_dim = 2     # spatial dimension

class Element:
    """ A class for 2d elements """

    def __init__(self):
        # element number
        self.number = 0
        # vertices
        self.x = []
        self.y = []

        # boundary conditions for each face
        # convention is [south, east, north, west]
        self.bc = []

        # boundary condition parameters
        self.bc_con_f = []  # connected face: 1: south, 2:east, 3:north, 4:west
        self.bc_con_el = []  # number of the connected element

el_list = []    # list of all elements
number = 1
for i in range(nSq):
    for j in range(nSq):
        # populate list of elements
        el = Element()
        el.number = number
        el_list.append(el)
        number = number + 1

## A: Generate the mesh
#----------------------------------------------------------------------
## A.1: Generate the mesh for a quarter section
#----------------------------------------------------------------------
## A.1.1: Set vertex positions of elements
#----------------------------------------------------------------------
nek_utils.set_vertices(el_list,nSq,dx,dy)

## A.1.2: Set boundary conditions for faces 
#----------------------------------------------------------------------
nek_utils.set_bc(el_list,nSq)

## A.1.3: Set curved edges
#----------------------------------------------------------------------


## A.2: Generate the complete mesh 
#----------------------------------------------------------------------
## A.2.1: Set vertex positions
#----------------------------------------------------------------------
## A.2.2: Set curved edges
#----------------------------------------------------------------------
## A.2.3: Set boundary conditions
#----------------------------------------------------------------------

## B: Write the mesh to rea file
#----------------------------------------------------------------------
# generate a rea skeleton file
#----------------------------------------------------------------------
nek_utils.rea_skel()

## B.1: Write vertex positions
#----------------------------------------------------------------------
n_tot = nSq
nek_utils.write_mesh(el_list)

## B.2: Write curved edges
#----------------------------------------------------------------------
## B.3: Write boundary conditions
#----------------------------------------------------------------------
#pdb.set_trace()
nek_utils.write_bc(el_list)





