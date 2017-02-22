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
#nR = 16
nR = 8
#nSq = 19
nSq = 4         # nel in square region along one side of the square
#nSq = 8

# Define global variables here:
#----------------------------------------------------------------------
nPhi = nSq*8    # nel in circumferential direction
dr = R/nR       # length of one element

# Set distribution of elements in a certain way
# OPEN
#----------------------------------------------------------------------
# "Square" region:
dr_sq_ratio = 0.90   # minium to maximum element length
fact_sq = dr_sq_ratio**(1/(nSq-1))
fact_sq_sum = 0
for i in range(0,nSq):
    fact_sq_sum = fact_sq_sum + fact_sq**i
dr_sq_max = nSq*dr/fact_sq_sum
dr_sq_min = dr_sq_max*dr_sq_ratio
dr_sq = np.zeros(nSq)
for i in range(0,nSq):
    dr_sq[i] = my_math.geom_prog(nSq, dr_sq_max, dr_sq_min, i)

# "Onion" region:
# Use geometric progression for a small increase and then cosine for 
# sharp decrease of element size
#----------------------------------------------------------------------
# NOTE: This needs some more work
#----------------------------------------------------------------------
n_on_1 = int(m.floor(1/2*(nR-nSq)))
n_on_2 = (nR-nSq) - n_on_1
dr_on_interface = dr_sq_min

# Determine element size at the transition between geometric progression
# and cosine distribution such that the total length of the onion region
# in radial direction is again (nR-nSq)*dr and we end up at r = R
def x_transition(x):
    """ This function is defined by the requirement that the total length
    of both onion regions 1 and two need to be (nR-nSq)*dr
    Note that the first element of the second region is not the same as 
    the last of the first region, i.e. both regions "share" dr_transition
    """

    sum_cos = 0
    for i in range(1,n_on_2+1):
        sum_cos = sum_cos + m.cos(i/(n_on_2+1)*m.pi/2)
    ret = -dr*(nR-nSq) + dr_sq_min * (1-(x/dr_sq_min)**(n_on_1/(n_on_1-1)))/\
            (1-(x/dr_sq_min)**(1/(n_on_1-1))) + x*sum_cos
    return ret

def x_transition_prime(x):
    """ First derivative of x_transition function
    """

    sum_cos = 0
    for i in range(1,n_on_2+1):
        sum_cos = sum_cos + m.cos(i/(n_on_2+1)*m.pi/2)
    ret = sum_cos + dr_sq_min * ( (-n_on_1/(n_on_1-1)*x**(n_on_1/(n_on_1-1)-1)*\
            dr_sq_min**(n_on_1/(1-n_on_1))) * (1-(x/dr_sq_min)**(1/(n_on_1-1)))**(-1)+\
            (1-(x/dr_sq_min)**(n_on_1/(n_on_1-1)))*(-1)*(1-(x/dr_sq_min)**(1/(n_on_1-1)))**(-2)*\
            (-dr_sq_min**(1/(1-n_on_1))*1/(n_on_1-1)*x**(1/(n_on_1-1)-1)))
    return ret
    


dr_on_transition = my_math.newton_raphson(dr,x_transition, x_transition_prime)
dr_on_1 = np.zeros(n_on_1)
for i in range(0,n_on_1):
    dr_on_1[i] = my_math.geom_prog(n_on_1,dr_on_interface, dr_on_transition, i)
dr_on_2 = np.zeros(n_on_2)
for i in range(0,n_on_2):
    dr_on_2[i] = dr_on_transition*m.cos((i+1)/(n_on_2+1)*m.pi/2)
dr_on = np.concatenate([dr_on_1, dr_on_2])
# END
#----------------------------------------------------------------------

n_tot = nR * nPhi   # total number of elements
spatial_dim = 2     # spatial dimension


el_list = []    # list of all elements
number = 1
# Populate list of elements: first, the square region
for i in range(nSq):
    for j in range(nSq):
        el = elementclass.Element()
        el.number = number
        el_list.append(el)
        number = number + 1
# Populate list of elements: second, the curved region outside (onion region)
for i in range(nR-nSq):     # loop through each onion like layer outwards
    for j in range(nSq*2): # loop in clockwise direction through each layer
        el = elementclass.Element()
        el.number = number
        el_list.append(el)
        number = number + 1

## A: Generate the mesh
#----------------------------------------------------------------------
## A.1: Generate the mesh for a quarter section
#----------------------------------------------------------------------
## A.1.1: Set vertex positions of elements
#----------------------------------------------------------------------
nek_utils.set_vertices(el_list,nR,nSq,dr_sq, dr_on)

## A.1.2: Set boundary conditions for faces 
#----------------------------------------------------------------------
nek_utils.set_bc_q1(el_list,nR,nSq)

## A.1.3: Set curved edges
#----------------------------------------------------------------------


## A.2: Generate the complete mesh 
#----------------------------------------------------------------------
## A.2.1: Set vertex positions
#----------------------------------------------------------------------
nek_utils.compl_mesh(el_list,nR,nSq)
## A.2.2: Set curved edges
#----------------------------------------------------------------------
## A.2.3: Set boundary conditions
#----------------------------------------------------------------------
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
## B.3: Write boundary conditions
#----------------------------------------------------------------------
nek_utils.write_bc(el_list, nR, nSq)

## C: Do some checks
#----------------------------------------------------------------------
nek_utils.check_mesh_quality(el_list, nR, nSq)





