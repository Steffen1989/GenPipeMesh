#!/usr/bin/env python3.5
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
import numpy as np
import math as m
import my_math
import elementclass 
import pdb
import sys
import re

# Input Variables
#----------------------------------------------------------------------
R = 0.5         # radius
#nR = 28         # nel in radial direction
nR = 8 
#nSq = 19
nSq = 4         # nel in square region along one side of the square
#nSq = 8
L_z = 5.00      # Length in streamwise direction z
nz = 20          # Number of elements in streamwise direction z

# For Resolution
N = 7           # Polynomial order
Re_t = 180      # Friction Reynolds number


# Some Input for tuning the mesh
#----------------------------------------------------------------------

# Stretch nominal value of dr in square by this factor
stretch_sq = 1.3    
# Min to max element length along axis in square
dr_sq_ratio = 0.80  

# Ratio of min to max element x (resp. y) component along intersection
# Note that this is not the real length but its projection along x-, 
# respective y-axis
dr_sq_int_ratio = 0.9    
# First xx in onion region is increasing and (xx-1) is decreasing
distri_on = 0.0     

# Semi-major axis at the interface betwenn square and onion region
# Note: semi-minor axis is defined by position of element along y-axis
a_interf = 0.57

# Keep the outermost onion layer constant: 
# 0 = not const., 1 = const.
tog_r_out_const = 0
# Use exp. or sin distribution for semi-major axis in onion region
# exp gives a sharper decrease, hence more circle like shape in 
# the first onion layers.
# 0 = exp, 1 = sin
tog_a_on_dist = 0


# Define some global variables here:
#----------------------------------------------------------------------
dr_nominal = R/nR       # nominal length of one element
dr = dr_nominal
dz = L_z/nz

el_list = []    # list of all elements
# number of elements in one cross section
nel_quarter = (nSq**2+(nR-nSq)*nSq*2)
nel_cross_section = nel_quarter*4       

number = 1
# Populate list of elements: first, the square region
for i in range(nSq):
    for j in range(nSq):
        el = elementclass.Element()
        el.number = number
#        el.c = np.zeros(4)
        el_list.append(el)
        number = number + 1
# Populate list of elements: second, the curved region outside (onion region)
for i in range(nR-nSq):     # loop through each onion like layer outwards
    for j in range(nSq*2): # loop in clockwise direction through each layer
        el = elementclass.Element()
        el.number = number
#        el.c = np.zeros(4)
        el_list.append(el)
        number = number + 1
#### Populate list of elements: third, streamwise direction
########
#### Note that this leads to a different ordering of the elements
########
###for iz in range(1,nz):          # loop through streamwise cross sections
###    # loop over all elements in the previous cross section
###    for prev_el in range(nel_quarter):    
###        el = elementclass.Element()
###        el.number = prev_el+1 + iz*nel_cross_section
###        el_list.append(el)
###
###pdb.set_trace()        




## A: Generate the mesh
#----------------------------------------------------------------------
## A.1: Generate the mesh for a quarter section
#----------------------------------------------------------------------
## A.1.1: Set vertex positions of elements
# (This is the essential part of the code)
#----------------------------------------------------------------------
nek_utils.set_vertices(el_list, nR, nSq, dr, dz, dr_sq_ratio,\
        dr_sq_int_ratio, stretch_sq, distri_on, a_interf,\
        tog_r_out_const, tog_a_on_dist)

## A.1.2: Set boundary conditions for faces 
#----------------------------------------------------------------------
nek_utils.set_bc_q1(el_list,nR,nSq)


## A.2: Generate the complete mesh 
#----------------------------------------------------------------------
## A.2.1: Set vertex positions
#----------------------------------------------------------------------
nek_utils.compl_mesh(el_list,nR,nSq)
nek_utils.extrude(el_list,nR,nSq,nz,dz)


## A.2.2: Set boundary conditions
# (for each quarter separately)
#----------------------------------------------------------------------
#nek_utils.set_bc_q1(el_list,nR,nSq)
nek_utils.set_bc_q2(el_list,nR,nSq)
nek_utils.set_bc_q3(el_list,nR,nSq)
nek_utils.set_bc_q4(el_list,nR,nSq)


## B: Write the mesh to rea file
#----------------------------------------------------------------------
# generate a rea skeleton file
#----------------------------------------------------------------------
nek_utils.rea_skel()
## B.1: Write vertex positions
#----------------------------------------------------------------------
nek_utils.write_mesh(el_list)
## B.2: Write curved edges
#----------------------------------------------------------------------
nek_utils.write_curv(el_list)
## B.3: Write boundary conditions
#----------------------------------------------------------------------
nek_utils.write_fl_bc(el_list, nR, nSq)
nek_utils.write_th_bc(el_list, nR, nSq)

## C: Do some checks and write a little output
#----------------------------------------------------------------------
nek_utils.dump_input_vars(R, nR, nSq, N, Re_t, stretch_sq,\
        dr_sq_ratio, dr_sq_int_ratio, distri_on, a_interf,\
        tog_r_out_const, tog_a_on_dist)

nek_utils.check_mesh_quality(el_list, nR, nSq, R, N, Re_t)
