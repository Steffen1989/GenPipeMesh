#!/usr/bin/env python3.5
# This python scrip is inspired by pipeMeshNek.f90 from Jacopo Canton.
# However, I am not that familiar with Fortran 90 and I wanted to 
# improve my python skills and understand the mesh generation better.
# Hence, I rewrite the code from scratch
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# Steffen Straub
# 2017-02-08
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

# Import modules
#----------------------------------------------------------------------
import nek_utils


# Generate a rea skeleton file
#----------------------------------------------------------------------
nek_utils.rea_skel()
