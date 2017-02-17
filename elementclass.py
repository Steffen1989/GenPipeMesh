# Extra file for class definitions
import numpy as np


class Element:
    """ A class for 2d elements """

    def __init__(self):

        # element number
        self.number = 0

        # vertices
        self.x = np.zeros(4)
        self.y = np.zeros(4)
        
        # position within quarter section
        self.pos = ''

        # boundary conditions for each face
        # convention is [south, east, north, west]
        self.bc = []

        # boundary condition parameters
        self.bc_con_f = np.zeros(4)  # connected face: 1: south, 2:east, 3:north, 4:west
        self.bc_con_el = np.zeros(4)  # number of the connected element

