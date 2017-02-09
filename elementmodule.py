#import numpy as np

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
        self.bc_params = []
