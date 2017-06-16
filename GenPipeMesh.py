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
# Radius
R = 0.5
# Number of elements in radial direction
nR = 9
# Number of elements along one side of the "square" region
nSq = 5
# Length of the pipe in streamwise (z) direction
L_z = 5.00
# Number of elements in streamwise (z) direction
nZ =  30

# Type of thermal BC
th_bc_type = 't  '

# In order to check the resolution requirements
N = 11           # Polynomial order
Re_t = 360      # Friction Reynolds number

# Do you want thermal boundary conditions?  
# "True" or "False"
if_therm = False



# Some Input for tuning the mesh
#----------------------------------------------------------------------

# Stretch nominal value of dr in square by this factor
stretch_sq = 1.2    
# Min to max element length along axis in square
dr_sq_ratio = 0.9  

# Ratio of min to max element of x (resp. y) component along intersection
# Note that this is not the real length but its projection along x-, 
# respective y-axis
dr_sq_int_ratio = 0.8    
# First xx in onion region is increasing and (xx-1) is decreasing
distri_on = 0.0     

# Semi-major axis at the interface betwenn square and onion region
# Note: semi-minor axis is defined by position of element along y-axis
#a_interf = 0.57
a_interf = 0.57

# Keep the outermost onion layer constant: 
# 0 = not const., 1 = const.
tog_r_out_const = 0
# Use exp. or sin distribution for semi-major axis in onion region
# exp gives a sharper decrease, hence more circle like shape in 
# the first onion layers.
# 0 = exp, 1 = sin
tog_a_on_dist = 0



## 0: Check if input variables are OK
#----------------------------------------------------------------------
nek_utils.check_input(nR, nSq, nZ, R, L_z, th_bc_type, N, Re_t,\
        if_therm)



# Define some global variables here:
#----------------------------------------------------------------------
dr_nominal = R/nR       # nominal length of one element
dr = dr_nominal
dz = L_z/nZ

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
#nek_utils.set_bc_q1(el_list,nR,nSq)


## A.2: Generate the complete mesh 
#----------------------------------------------------------------------
## A.2.1: Set vertex positions
#----------------------------------------------------------------------
nek_utils.compl_mesh(el_list,nR,nSq)
nek_utils.extrude(el_list,nR,nSq,nZ,dz)


## A.2.2: Set boundary conditions
# (for each quarter separately)
#----------------------------------------------------------------------
nek_utils.set_bc_q1(el_list,nR,nSq,th_bc_type)
nek_utils.set_bc_q2(el_list,nR,nSq,th_bc_type)
nek_utils.set_bc_q3(el_list,nR,nSq,th_bc_type)
nek_utils.set_bc_q4(el_list,nR,nSq,th_bc_type)


## B: Write the mesh to rea file
#----------------------------------------------------------------------
# generate a rea skeleton file
#----------------------------------------------------------------------
nek_utils.rea_skel(2, if_therm, 'base2d.rea')
nek_utils.rea_skel(3, if_therm, 'base3d.rea')
## B.1: Write vertex positions
#----------------------------------------------------------------------
nek_utils.write_mesh(el_list, nR, nSq, 2, 'base2d.rea')
nek_utils.write_mesh(el_list, nR, nSq, 3, 'base3d.rea')
## B.2: Write curved edges
#----------------------------------------------------------------------
nek_utils.write_curv(el_list, nR, nSq, 2, 'base2d.rea')
nek_utils.write_curv(el_list, nR, nSq, 3, 'base3d.rea')
## B.3: Write boundary conditions
#----------------------------------------------------------------------
nek_utils.write_fl_bc(el_list, nR, nSq, 2, 'base2d.rea')
nek_utils.write_fl_bc(el_list, nR, nSq, 3, 'base3d.rea')
if (if_therm):
    nek_utils.write_th_bc(el_list, nR, nSq, 2, 'base2d.rea')
    nek_utils.write_th_bc(el_list, nR, nSq, 3, 'base3d.rea')

## C: Do some checks and write a little output
#----------------------------------------------------------------------
nek_utils.dump_input_vars(R, nR, nSq, nZ, L_z, N, Re_t, stretch_sq,\
        dr_sq_ratio, dr_sq_int_ratio, distri_on, a_interf,\
        tog_r_out_const, tog_a_on_dist)

nek_utils.check_mesh_quality(el_list, nR, nSq, nZ, R, L_z, N, Re_t)
