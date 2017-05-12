# Extra file for class definitions
import numpy as np


class Element:
    """ A class for 2d elements """

    def __init__(self):

        # element number
        self.number = 0

        # vertices
        self.x = np.zeros(8)
        self.y = np.zeros(8)
        self.z = np.zeros(8)
        
        # position within quarter section
        self.pos = ''

        # fluid and thermal boundary conditions for each face
        # convention is [south, east, north, west]
        self.fl_bc = [] 
        self.th_bc = []

        # boundary condition parameters
        self.bc_con_f = np.zeros(6)  # connected face: 1: south, 2:east, 3:north, 4:west
        # 5: back, 6: front
        self.bc_con_el = np.zeros(6)  # number of the connected element

        # curvature
        self.c = np.zeros(8) # curvature information 1: south, 2:east, 3:north, 4:west
        # on the next layer: 5: south, 6: east, 7: north, 8: west


