# A collection of important functions
import sys
import pdb
import numpy as np
import math as m
import my_math
import elementclass


def set_vertices(elements,nR,nSq,dr_sq, dr_on):
    """ Set vertex location for each element. 

    The vertices are set in a special way. For now the inner section, 
    called square section, is just a regular square. The outer part, 
    called onion region, is built up by ellipses and straight lines.
    The semi-major axis a is decreasing each layer outwards so that a=1
    in the outermost layer which gives a circle. The constant on the right
    hand side is kept constant and semi-minor and semi-major axis are
    varied.
    """
    
    # Variable definitions
    ntheta = nSq*2  # number of elements in one onion layer
    r_const = 1         # constant in ellipses
    a_on = np.zeros(2)    # semi-major axis
    a_row = np.zeros(2)    
    a_col = np.zeros(2)
    b_on = np.zeros(2)    # semi-minor axis
    b_row = np.zeros(2)    
    b_col = np.zeros(2)
    slope_on = np.zeros(2)     # slope of straight lines
    slope_row = np.zeros(2)
    slope_col = np.zeros(2)
    # intersection between straight lines and ellipses
    # at the last curve of square region
    x_inters_row = np.zeros(2)     
    x_inters_col = np.zeros(2)      
    y_inters_col = np.zeros(2)
    y_inters_row = np.zeros(2)
 
    for el in elements:
        if (el.number <= nSq**2):   
            # This is the inner, "square" section
            # OPEN square
            #--------------------------------------------------
            i = (el.number-1)%nSq     # column number
            j = int((el.number-1)/nSq)    # row number

            # Determine the semi-major axis for "square" section
            #----------------------------------------------------------------------
            # determine the minimum semi-major axis at the edge to the onion region
            drop = 0.75         # max drop in percent compared to all squares
            n_ellip_sq = nSq  # number of ellipses
            x_interface = np.sum(dr_sq)
#            x_interface = nSq*dr_sq

            b_interface = x_interface # semi-minor axis
            y_sq_interface = x_interface*drop

            # Option1: Drop by a certain percentage
#            a_interface = x_interface*b_interface/( (r_const**2*b_interface**2-y_sq_interface**2)**0.5 )
            # Opiton2: Set corner angle = 120°
#            a_interface = b_interface/(m.tan(m.pi/12))**0.5
            # Option3: Set a = 0.5
            a_interface = 0.5

            # Idea: define ellipses in the inner region such that they coincide with the 
            # element in the onion region.
            # OPEN def semi-major axis
            #--------------------------------------------------
            slope_row[0] = m.tan(m.pi/2*(j/ntheta))  # slope of the straight line on the bottom side
            slope_row[1] = m.tan(m.pi/2*((j+1)/ntheta)) # slope of the straight line on the top side
            slope_col[0] = m.tan(m.pi/2*((ntheta-i)/ntheta))  # slope of the straight line on the left side
            slope_col[1] = m.tan(m.pi/2*((ntheta-i-1)/ntheta)) # slope of the straight line on the right side

#            b_row[0] = j*dr    # small semi-minor axis  
#            b_row[1] = (j+1)*dr   # large semi-minor axis
#            b_col[0] = i*dr    # small semi-minor axis  
#            b_col[1] = (i+1)*dr   # large semi-minor axis

            b_row[0] = np.sum(dr_sq[:j])    # small semi-minor axis  
            b_row[1] = np.sum(dr_sq[:j+1])   # large semi-minor axis
            b_col[0] = np.sum(dr_sq[:i])    # small semi-minor axis  
            b_col[1] = np.sum(dr_sq[:i+1])   # large semi-minor axis

            x_inters_row[0] = my_math.intersec_ellip_line(b_interface,a_interface,r_const,slope_row[0],0)
            x_inters_row[1] = my_math.intersec_ellip_line(b_interface,a_interface,r_const,slope_row[1],0)
            y_inters_row[0] = my_math.line(slope_row[0],x_inters_row[0],0)
            y_inters_row[1] = my_math.line(slope_row[1],x_inters_row[1],0)

            x_inters_col[0] = my_math.intersec_ellip_line(a_interface,b_interface,r_const,slope_col[0],0)
            x_inters_col[1] = my_math.intersec_ellip_line(a_interface,b_interface,r_const,slope_col[1],0)
            y_inters_col[0] = my_math.line(slope_col[0],x_inters_col[0],0)
            y_inters_col[1] = my_math.line(slope_col[1],x_inters_col[1],0)

            if (j==0):
                a_row[0] = 0   # this is reset later
                a_row[1] = x_inters_row[1]*b_row[1]/( (r_const**2*b_row[1]**2 - y_inters_row[1]**2)**0.5 )
            else: 
                a_row[0] = x_inters_row[0]*b_row[0]/( (r_const**2*b_row[0]**2 - y_inters_row[0]**2)**0.5 )
                a_row[1] = x_inters_row[1]*b_row[1]/( (r_const**2*b_row[1]**2 - y_inters_row[1]**2)**0.5 )

            if (i==0):                # note that x and y need to be switched here
                a_col[0] = 0   # this is reset later 
                a_col[1] = y_inters_col[1]*b_col[1]/( (r_const**2*b_col[1]**2 - x_inters_col[1]**2)**0.5 )
            else:
                a_col[0] = y_inters_col[0]*b_col[0]/( (r_const**2*b_col[0]**2 - x_inters_col[0]**2)**0.5 )
                a_col[1] = y_inters_col[1]*b_col[1]/( (r_const**2*b_col[1]**2 - x_inters_col[1]**2)**0.5 )

            # CLOSE def semi-major axis
            #--------------------------------------------------

            if (j==0): # first row
                if (i==0):  # first col
                    x2 = my_math.intersec_ellip_ellip(a_row[1],b_row[1],r_const**2,\
                            b_col[1],a_col[1],r_const**2)
                    y2 = my_math.ellipse(a_row[1],b_row[1],r_const**2,x2)
                    el.x = np.array([np.sum(dr_sq[:i]), np.sum(dr_sq[:i+1]), x2, np.sum(dr_sq[:i])])
                    el.y = np.array([np.sum(dr_sq[:j]), np.sum(dr_sq[:j]), y2, np.sum(dr_sq[:j+1])])
                else:
                    x2 = my_math.intersec_ellip_ellip(a_row[1],b_row[1],r_const**2,\
                            b_col[1],a_col[1],r_const**2)
                    x3 = my_math.intersec_ellip_ellip(a_row[1],b_row[1],r_const**2,\
                            b_col[0],a_col[0],r_const**2)
                    y2 = my_math.ellipse(a_row[1],b_row[1],r_const**2,x2)
                    y3 = my_math.ellipse(a_row[1],b_row[1],r_const**2,x3)
                    el.x = np.array([np.sum(dr_sq[:i]), np.sum(dr_sq[:i+1]), x2, x3])
                    el.y = np.array([np.sum(dr_sq[:j]), np.sum(dr_sq[:j]), y2, y3])
            elif (j>0 and i==0): # first col
                x1 = my_math.intersec_ellip_ellip(a_row[0],b_row[0],r_const**2,\
                        b_col[1],a_col[1],r_const**2)
                x2 = my_math.intersec_ellip_ellip(a_row[1],b_row[1],r_const**2,\
                        b_col[1],a_col[1],r_const**2)
                y1 = my_math.ellipse(a_row[0],b_row[0],r_const**2,x1)
                y2 = my_math.ellipse(a_row[1],b_row[1],r_const**2,x2)
                el.x = np.array([np.sum(dr_sq[:i]), x1, x2, np.sum(dr_sq[:i])])
                el.y = np.array([np.sum(dr_sq[:j]), y1, y2, np.sum(dr_sq[:j+1])])
            elif (i> 0 and j>0):    # inside
                #find intersection between both ellipses
                x0 = my_math.intersec_ellip_ellip(a_row[0],b_row[0],r_const**2,\
                        b_col[0],a_col[0],r_const**2)
                x1 = my_math.intersec_ellip_ellip(a_row[0],b_row[0],r_const**2,\
                        b_col[1],a_col[1],r_const**2)
                x2 = my_math.intersec_ellip_ellip(a_row[1],b_row[1],r_const**2,\
                        b_col[1],a_col[1],r_const**2)
                x3 = my_math.intersec_ellip_ellip(a_row[1],b_row[1],r_const**2,\
                        b_col[0],a_col[0],r_const**2)
                y0 = my_math.ellipse(a_row[0],b_row[0],r_const**2,x0)
                y1 = my_math.ellipse(a_row[0],b_row[0],r_const**2,x1)
                y2 = my_math.ellipse(a_row[1],b_row[1],r_const**2,x2)
                y3 = my_math.ellipse(a_row[1],b_row[1],r_const**2,x3)
                el.x = np.array([x0, x1, x2, x3])
                el.y = np.array([y0, y1, y2, y3])
            else:
                sys.exit(1)
            #---------------------------------------------------
            # END square
        else:                       
            # This is the outer, "onion" section
            # OPEN onion
            #--------------------------------------------------
            i = ((el.number-1)-nSq**2)%(nSq*2) # position in clockwise manner through each layer
            k = abs(i-((nSq*2)-1))                  # position in anticlockwise manner
            j = int(((el.number-1)-nSq**2)/(nSq*2)) # onion like layer number, inner one is first,
            # starting from j=0
            l = (nR - nSq) - (j+1)                  # onion like layer number, outer one is last l=0

            # Determine the semi-major axis for "onion" section
            #----------------------------------------------------------------------
            a_wall = 0.5    # semi-major axis at last layer (wall)
            
            a_on[0] = my_math.geom_prog(nR-nSq, a_interface, a_wall, j)
            a_on[1] = my_math.geom_prog(nR-nSq, a_interface, a_wall, j+1)

#            if (j < (nR-nSq-1)):
#                a_on[0] = my_math.geom_prog(nR-nSq, a_interface, a_wall, j)
#                a_on[1] = my_math.geom_prog(nR-nSq, a_interface, a_wall, j+1)
#
#            else:   # last element has a=0.5 on both sides
#                a_on[0] = my_math.geom_prog(nR-nSq, a_interface, a_wall, j)
#                a_on[1] = a_on[0]


#            dr = dr_on
#            b_on[0] = (j+nSq)*dr
#            b_on[1] = (j+1+nSq)*dr

            b_on[0] = np.sum(dr_sq)+np.sum(dr_on[:j])
            b_on[1] = np.sum(dr_sq)+np.sum(dr_on[:j+1])
            slope_on[0] = m.tan(m.pi/2*(k/ntheta))  # slope of the straight line on the right side
            # of the element (upper part) or bottom side (lower part)
            slope_on[1] = m.tan(m.pi/2*((k+1)/ntheta)) # slope of the straight line on the left side
            # of the element (upper part) or top side (lower part)
            if (i <= (nSq-1)):  # upper part, including border /
                x0 = my_math.intersec_ellip_line(a_on[0],b_on[0],r_const**2,slope_on[1],0)
                x1 = my_math.intersec_ellip_line(a_on[0],b_on[0],r_const**2,slope_on[0],0)
                x2 = my_math.intersec_ellip_line(a_on[1],b_on[1],r_const**2,slope_on[0],0)
                x3 = my_math.intersec_ellip_line(a_on[1],b_on[1],r_const**2,slope_on[1],0)
                el.x = np.array([x0, x1, x2, x3])
                el.y[0:2] = my_math.ellipse(a_on[0],b_on[0],r_const**2,el.x[0:2])
                el.y[2:4] = my_math.ellipse(a_on[1],b_on[1],r_const**2,el.x[2:4])
            elif (i >= nSq):     # lower part, including border /
                x0 = my_math.intersec_ellip_line(b_on[0],a_on[0],r_const**2,slope_on[0],0)
                x1 = my_math.intersec_ellip_line(b_on[1],a_on[1],r_const**2,slope_on[0],0)
                x2 = my_math.intersec_ellip_line(b_on[1],a_on[1],r_const**2,slope_on[1],0)
                x3 = my_math.intersec_ellip_line(b_on[0],a_on[0],r_const**2,slope_on[1],0)
                y0 = my_math.line(slope_on[0],x0,0)
                y1 = my_math.line(slope_on[0],x1,0)
                y2 = my_math.line(slope_on[1],x2,0)
                y3 = my_math.line(slope_on[1],x3,0)
                el.y = np.array([y0, y1, y2, y3])
                el.x = np.array([x0, x1, x2, x3])


def compl_mesh(elements, nR, nSq):
    """ Complete the quarter mesh to a whole cross section """

    # Mirror the first quarter along y-axis to create mesh in 2nd quadrant
    nel_quarter = nSq**2 + (nR-nSq)*2*nSq       # number of elements in one quarter
    el_list_2nd = []     # list for the elements in the 2nd quadrant
    el_list_3rd = []     # list for the elements in the 3rd quadrant
    el_list_4th = []     # list for the elements in the 4th quadrant

    for el in elements:
        # First, create elements in 2nd quadrant
        # Define the mirrored element
        mirr_el = elementclass.Element()

        # Element number
        mirr_el.number = el.number+nel_quarter

        # Vertices 
        # Note that vertex numbering needs to be adjusted so that "right-handed" elements are created
        mirr_el.x = el.x*(-1) 
        # swap values to get right handed element
        mirr_el.x = np.array([mirr_el.x[1], mirr_el.x[0], mirr_el.x[3], mirr_el.x[2]])
        mirr_el.y = el.y
        mirr_el.y = np.array([mirr_el.y[1], mirr_el.y[0], mirr_el.y[3], mirr_el.y[2]])
        
        # Position
        mirr_el.pos = el.pos

        # Boundary condition
        mirr_el.bc = ['E  ','E  ','E  ','E  ']

        mirr_el.bc_con_f = np.zeros(4)  # connected face: 1: south, 2:east, 3:north, 4:west
        mirr_el.bc_con_el = np.zeros(4)  # number of the connected element

        # Add mirrored element to the list of elements
        el_list_2nd.append(mirr_el)

    for el in elements:
        # Second, create elements in third quadrant
        # Define the mirrored element
        mirr_el = elementclass.Element()

        # Element number
        mirr_el.number = el.number+nel_quarter*2

        # Vertices 
        # Note that vertex numbering needs to be adjusted so that "right-handed" elements are created
        mirr_el.x = el.x*(-1) 
        # swap values to get right handed element
        mirr_el.x = np.array([mirr_el.x[2], mirr_el.x[3], mirr_el.x[0], mirr_el.x[1]])
        mirr_el.y = el.y*(-1)
        mirr_el.y = np.array([mirr_el.y[2], mirr_el.y[3], mirr_el.y[0], mirr_el.y[1]])
        
        # Position
        mirr_el.pos = el.pos

        # Boundary condition
        mirr_el.bc = ['E  ','E  ','E  ','E  ']

        mirr_el.bc_con_f = np.zeros(4)  # connected face: 1: south, 2:east, 3:north, 4:west
        mirr_el.bc_con_el = np.zeros(4)  # number of the connected element

        
        # Add mirrored element to the list of elements
        el_list_3rd.append(mirr_el)

    for el in elements:
        # Third, create elements in fourth quadrant
        # Define the mirrored element
        mirr_el = elementclass.Element()

        # Element number
        mirr_el.number = el.number+nel_quarter*3

        # Vertices 
        # Note that vertex numbering needs to be adjusted so that "right-handed" elements are created
        mirr_el.x = el.x 
        # swap values to get right handed element
        mirr_el.x = np.array([mirr_el.x[3], mirr_el.x[2], mirr_el.x[1], mirr_el.x[0]])
        mirr_el.y = el.y*(-1)
        mirr_el.y = np.array([mirr_el.y[3], mirr_el.y[2], mirr_el.y[1], mirr_el.y[0]])
        
        # Position
        mirr_el.pos = el.pos

        # Boundary condition
        mirr_el.bc = ['E  ','E  ','E  ','E  ']

        mirr_el.bc_con_f = np.zeros(4)  # connected face: 1: south, 2:east, 3:north, 4:west
        mirr_el.bc_con_el = np.zeros(4)  # number of the connected element

        
        # Add mirrored element to the list of elements
        el_list_4th.append(mirr_el)

    # Add all new elements to the element list
    elements.extend(el_list_2nd)
    elements.extend(el_list_3rd)
    elements.extend(el_list_4th)



def set_bc_q1(elements,nR,nSq):
    """ Set boundary conditions for each face. 
    
    This is for quadrant 1.
    """

    nel_quarter = nSq**2 + (nR-nSq)*2*nSq       # number of elements in one quarter

    # only consider elements in the first quadrant
    elements = elements[0:nel_quarter]

    for el in elements:
        n = el.number
#         position = check_position(n, nSq)
        check_position(el,nR,nSq)
        position = el.pos

        if (n <= nSq**2):   # we are in the square section
            i = (n-1)%nSq     # column number
            j = int((n-1)/nSq)    # row number

            if (position == 'sq_low_left'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nel_quarter*3, n+1, n+nSq, n+nel_quarter])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_low_right'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nel_quarter*3, nSq**2+(2*nSq-j), n+nSq, n-1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_up_left'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n-nSq, n+1, n+nSq, n+nel_quarter])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_up_right'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n-nSq, nSq**2+(2*nSq-j), n+nSq, n-1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_north_row'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n-nSq, n+1, n+nSq, n-1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_east_col'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n-nSq, nSq**2+(2*nSq-j), n+nSq, n-1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_south_row'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nel_quarter*3, n+1, n+nSq, n-1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_west_col'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n-nSq, n+1, n+nSq, n+nel_quarter])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_internal'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n-nSq, n+1, n+nSq, n-1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            else:
                print('position of element not found!')
                sys.exit(1)
        else:                       # we are in the outer onion like region :-)
            i = ((n-1)-nSq**2)%(nSq*2) # position in clockwise manner through each layer
            k = abs(i-((nSq*2)-1))                  # position in anticlockwise manner
            j = int(((n-1)-nSq**2)/(nSq*2)) # onion like layer number, inner one is first

            if ('on_up' in position):   # we are in the upper onion part
                if (position == 'on_up_south_sq_west_y'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n+1, n+2*nSq, n+nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_up_south_sq_east_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n+1, n+2*nSq, n-1])
                    el.bc_con_f = np.array([3, 3, 1, 2])
                elif (position == 'on_up_east_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-2*nSq, n+1, n+2*nSq, n-1])
                    el.bc_con_f = np.array([3, 3, 1, 2])
                elif (position == 'on_up_south_sq'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n+1, n+2*nSq, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_up_west_y_north_wall'):
                    el.bc = ['E  ','E  ','W  ','E  ']
                    el.bc_con_el = np.array([n-2*nSq, n+1, 0, n+nel_quarter])
                    el.bc_con_f = np.array([3, 4, 0, 2])
                elif (position == 'on_up_west_y'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-2*nSq, n+1, n+2*nSq, n+nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_up_north_wall_east_edge'):
                    el.bc = ['E  ','E  ','W  ','E  ']
                    el.bc_con_el = np.array([n-2*nSq, n+1, 0, n-1])
                    el.bc_con_f = np.array([3, 3, 0, 2])
                elif (position == 'on_up_north_wall'):
                    el.bc = ['E  ','E  ','W  ','E  ']
                    el.bc_con_el = np.array([n-2*nSq, n+1, 0, n-1])
                    el.bc_con_f = np.array([3, 4, 0, 2])
                elif (position == 'on_up_intern'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-2*nSq, n+1, n+2*nSq, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
            elif ('on_low' in position):    # we are in the lower onion part
                if (position == 'on_low_west_sq_south_x'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nel_quarter*3, n+nSq*2, n-1, (k+1)*nSq])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_low_west_sq_north_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+1, n+nSq*2, n-1, (k+1)*nSq])
                    el.bc_con_f = np.array([3, 4, 2, 2])
                elif (position == 'on_low_west_sq'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+1, n+nSq*2, n-1, (k+1)*nSq])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_low_south_x_east_wall'):
                    el.bc = ['E  ','W  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nel_quarter*3, 0, n-1, n-nSq*2])
                    el.bc_con_f = np.array([3, 0, 1, 2])
                elif (position == 'on_low_south_x'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nel_quarter*3, n+nSq*2, n-1, n-nSq*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_low_east_wall_north_edge'):
                    el.bc = ['E  ','W  ','E  ','E  ']
                    el.bc_con_el = np.array([n+1, 0, n-1, n-nSq*2])
                    el.bc_con_f = np.array([3, 0, 2, 2])
                elif (position == 'on_low_east_wall'):
                    el.bc = ['E  ','W  ','E  ','E  ']
                    el.bc_con_el = np.array([n+1, 0, n-1, n-nSq*2])
                    el.bc_con_f = np.array([3, 0, 1, 2])
                elif (position == 'on_low_north_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+1, n+nSq*2, n-1, n-nSq*2])
                    el.bc_con_f = np.array([3, 4, 2, 2])
                elif (position == 'on_low_intern'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+1, n+nSq*2, n-1, n-nSq*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
            else:
                print('position assignment was not correct!')
                sys.exit(3)


def set_bc_q2(elements,nR,nSq):
    """ Set boundary conditions for each face. 
    
    This is for quadrant 2.
    """

    nel_quarter = nSq**2 + (nR-nSq)*2*nSq       # number of elements in one quarter
    # only consider elements in the second quadrant
    elements = elements[nel_quarter:(nel_quarter*2)]

    for el in elements:
        n = el.number
        position = el.pos

        if (n-nel_quarter <= nSq**2):   # we are in the square section
            i = (n-nel_quarter-1)%nSq     # column number
            j = int((n-nel_quarter-1)/nSq)    # row number

            if (position == 'sq_low_left'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nel_quarter, n-nel_quarter, n+nSq, n+1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_low_right'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nel_quarter, n-1, n+nSq, nSq**2+(2*nSq-j)+nel_quarter])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_up_left'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n-nSq, n-nel_quarter, n+nSq, n+1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_up_right'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n-nSq, n-1, n+nSq, nSq**2+(2*nSq-j)+nel_quarter])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_north_row'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n-nSq, n-1, n+nSq, n+1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_east_col'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n-nSq, n-1, n+nSq, nSq**2+(2*nSq-j)+nel_quarter])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_south_row'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nel_quarter, n-1, n+nSq, n+1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_west_col'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n-nSq, n-nel_quarter, n+nSq, n+1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_internal'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n-nSq, n-1, n+nSq, n+1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            else:
                print('position of element not found!')
                sys.exit(1)
        else:                       # we are in the outer onion like region :-)
            i = ((n-nel_quarter-1)-nSq**2)%(nSq*2) # position in clockwise manner through each layer
            k = abs(i-((nSq*2)-1))                  # position in anticlockwise manner
            j = int(((n-nel_quarter-1)-nSq**2)/(nSq*2)) # onion like layer number, inner one is first

            if ('on_up' in position):   # we are in the upper onion part
                if (position == 'on_up_south_sq_west_y'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n-nel_quarter, n+2*nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_up_south_sq_east_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n-1, n+2*nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 3])
                elif (position == 'on_up_east_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-2*nSq, n-1, n+2*nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 3])
                elif (position == 'on_up_south_sq'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n-1, n+2*nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_up_west_y_north_wall'):
                    el.bc = ['E  ','E  ','W  ','E  ']
                    el.bc_con_el = np.array([n-2*nSq, n-nel_quarter, 0, n+1])
                    el.bc_con_f = np.array([3, 4, 0, 2])
                elif (position == 'on_up_west_y'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-2*nSq, n-nel_quarter, n+2*nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_up_north_wall_east_edge'):
                    el.bc = ['E  ','E  ','W  ','E  ']
                    el.bc_con_el = np.array([n-2*nSq, n-1, 0, n+1])
                    el.bc_con_f = np.array([3, 4, 0, 3])
                elif (position == 'on_up_north_wall'):
                    el.bc = ['E  ','E  ','W  ','E  ']
                    el.bc_con_el = np.array([n-2*nSq, n-1, 0, n+1])
                    el.bc_con_f = np.array([3, 4, 0, 2])
                elif (position == 'on_up_intern'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-2*nSq, n-1, n+2*nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
            elif ('on_low' in position):    # we are in the lower onion part
                if (position == 'on_low_west_sq_south_x'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nel_quarter, (k+1)*nSq+nel_quarter, n-1, n+nSq*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_low_west_sq_north_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+1, (k+1)*nSq+nel_quarter, n-1, n+nSq*2])
                    el.bc_con_f = np.array([3, 4, 4, 2])
                elif (position == 'on_low_west_sq'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+1, (k+1)*nSq+nel_quarter, n-1, n+nSq*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_low_south_x_east_wall'):
                    el.bc = ['E  ','E  ','E  ','W  ']
                    el.bc_con_el = np.array([n+nel_quarter, n-nSq*2, n-1, 0])
                    el.bc_con_f = np.array([3, 4, 1, 0])
                elif (position == 'on_low_south_x'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nel_quarter, n-nSq*2, n-1, n+nSq*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_low_east_wall_north_edge'):
                    el.bc = ['E  ','E  ','E  ','W  ']
                    el.bc_con_el = np.array([n+1, n-nSq*2, n-1, 0])
                    el.bc_con_f = np.array([3, 4, 4, 0])
                elif (position == 'on_low_east_wall'):
                    el.bc = ['E  ','E  ','E  ','W  ']
                    el.bc_con_el = np.array([n+1, n-nSq*2, n-1, 0])
                    el.bc_con_f = np.array([3, 4, 1, 0])
                elif (position == 'on_low_north_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+1, n-nSq*2, n-1, n+nSq*2])
                    el.bc_con_f = np.array([3, 4, 4, 2])
                elif (position == 'on_low_intern'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+1, n-nSq*2, n-1, n+nSq*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
            else:
                print('position assignment was not correct!')
                sys.exit(3)

def set_bc_q3(elements,nR,nSq):
    """ Set boundary conditions for each face. 
    
    This is for quadrant 3.
    """

    nel_quarter = nSq**2 + (nR-nSq)*2*nSq       # number of elements in one quarter
    # only consider elements in the third quadrant
    elements = elements[nel_quarter*2:(nel_quarter*3)]

    for el in elements:
        n = el.number
        position = el.pos

        if (n-nel_quarter*2 <= nSq**2):   # we are in the square section
            i = (n-nel_quarter*2-1)%nSq     # column number
            j = int((n-nel_quarter*2-1)/nSq)    # row number

            if (position == 'sq_low_left'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n+nel_quarter, n-nel_quarter, n+1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_low_right'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n-1, n-nel_quarter, nSq**2+(2*nSq-j)+nel_quarter*2])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_up_left'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n+nel_quarter, n-nSq, n+1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_up_right'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n-1, n-nSq, nSq**2+(2*nSq-j)+nel_quarter*2])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_north_row'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n-1, n-nSq, n+1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_east_col'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n-1, n-nSq, nSq**2+(2*nSq-j)+nel_quarter*2])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_south_row'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n-1, n-nel_quarter, n+1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_west_col'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n+nel_quarter, n-nSq, n+1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_internal'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n-1, n-nSq, n+1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            else:
                print('position of element not found!')
                sys.exit(1)
        else:                       # we are in the outer onion like region :-)
            i = ((n-nel_quarter*2-1)-nSq**2)%(nSq*2) # position in clockwise manner through each layer
            k = abs(i-((nSq*2)-1))                  # position in anticlockwise manner
            j = int(((n-nel_quarter*2-1)-nSq**2)/(nSq*2)) # onion like layer number, inner one is first

            if ('on_up' in position):   # we are in the upper onion part
                if (position == 'on_up_south_sq_west_y'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+2*nSq, n+nel_quarter, n-nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_up_south_sq_east_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+2*nSq, n-1, n-nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 1])
                elif (position == 'on_up_east_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+2*nSq, n-1, n-2*nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 1])
                elif (position == 'on_up_south_sq'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+2*nSq, n-1, n-nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_up_west_y_north_wall'):
                    el.bc = ['W  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([0, n+nel_quarter, n-2*nSq, n+1])
                    el.bc_con_f = np.array([0, 4, 1, 2])
                elif (position == 'on_up_west_y'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+2*nSq, n+nel_quarter, n-2*nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_up_north_wall_east_edge'):
                    el.bc = ['W  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([0, n-1, n-2*nSq, n+1])
                    el.bc_con_f = np.array([0, 4, 1, 1])
                elif (position == 'on_up_north_wall'):
                    el.bc = ['W  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([0, n-1, n-2*nSq, n+1])
                    el.bc_con_f = np.array([0, 4, 1, 2])
                elif (position == 'on_up_intern'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+2*nSq, n-1, n-2*nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
            elif ('on_low' in position):    # we are in the lower onion part
                if (position == 'on_low_west_sq_south_x'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, (k+1)*nSq+nel_quarter*2, n-nel_quarter, n+nSq*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_low_west_sq_north_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, (k+1)*nSq+nel_quarter*2, n+1, n+nSq*2])
                    el.bc_con_f = np.array([4, 4, 1, 2])
                elif (position == 'on_low_west_sq'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, (k+1)*nSq+nel_quarter*2, n+1, n+nSq*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_low_south_x_east_wall'):
                    el.bc = ['E  ','E  ','E  ','W  ']
                    el.bc_con_el = np.array([n-1, n-nSq*2, n-nel_quarter, 0])
                    el.bc_con_f = np.array([3, 4, 1, 0])
                elif (position == 'on_low_south_x'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, n-nSq*2, n-nel_quarter, n+nSq*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_low_east_wall_north_edge'):
                    el.bc = ['E  ','E  ','E  ','W  ']
                    el.bc_con_el = np.array([n-1, n-nSq*2, n+1, 0])
                    el.bc_con_f = np.array([4, 4, 1, 0])
                elif (position == 'on_low_east_wall'):
                    el.bc = ['E  ','E  ','E  ','W  ']
                    el.bc_con_el = np.array([n-1, n-nSq*2, n+1, 0])
                    el.bc_con_f = np.array([3, 4, 1, 0])
                elif (position == 'on_low_north_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, n-nSq*2, n+1, n+nSq*2])
                    el.bc_con_f = np.array([4, 4, 1, 2])
                elif (position == 'on_low_intern'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, n-nSq*2, n+1, n+nSq*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
            else:
                print('position assignment was not correct!')
                sys.exit(3)

def set_bc_q4(elements,nR,nSq):
    """ Set boundary conditions for each face. 
    
    This is for quadrant 4.
    """

    nel_quarter = nSq**2 + (nR-nSq)*2*nSq       # number of elements in one quarter
    # only consider elements in the fourth quadrant
    elements = elements[nel_quarter*3:(nel_quarter*4)]

    for el in elements:
        n = el.number
        position = el.pos

        if (n-nel_quarter*3 <= nSq**2):   # we are in the square section
            i = (n-nel_quarter*3-1)%nSq     # column number
            j = int((n-nel_quarter*3-1)/nSq)    # row number

            if (position == 'sq_low_left'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n+1, n-nel_quarter*3, n-nel_quarter])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_low_right'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, nSq**2+(2*nSq-j)+nel_quarter*3, n-nel_quarter*3, n-1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_up_left'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n+1, n-nSq, n-nel_quarter])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_up_right'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, nSq**2+(2*nSq-j)+nel_quarter*3, n-nSq, n-1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_north_row'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n+1, n-nSq, n-1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_east_col'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, nSq**2+(2*nSq-j)+nel_quarter*3, n-nSq, n-1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_south_row'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n+1, n-nel_quarter*3, n-1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_west_col'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n+1, n-nSq, n-nel_quarter])
                el.bc_con_f = np.array([3, 4, 1, 2])
            elif (position == 'sq_internal'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = np.array([n+nSq, n+1, n-nSq, n-1])
                el.bc_con_f = np.array([3, 4, 1, 2])
            else:
                print('position of element not found!')
                sys.exit(1)
        else:                       # we are in the outer onion like region :-)
            i = ((n-nel_quarter*3-1)-nSq**2)%(nSq*2) # position in clockwise manner through each layer
            k = abs(i-((nSq*2)-1))                  # position in anticlockwise manner
            j = int(((n-nel_quarter*3-1)-nSq**2)/(nSq*2)) # onion like layer number, inner one is first

            if ('on_up' in position):   # we are in the upper onion part
                if (position == 'on_up_south_sq_west_y'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+2*nSq, n+1, n-nSq, n-nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_up_south_sq_east_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+2*nSq, n+1, n-nSq, n-1])
                    el.bc_con_f = np.array([3, 1, 1, 4])
                elif (position == 'on_up_east_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+2*nSq, n+1, n-2*nSq, n-1])
                    el.bc_con_f = np.array([3, 1, 1, 4])
                elif (position == 'on_up_south_sq'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+2*nSq, n+1, n-nSq, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_up_west_y_north_wall'):
                    el.bc = ['W  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([0, n+1, n-2*nSq, n-nel_quarter])
                    el.bc_con_f = np.array([0, 4, 1, 2])
                elif (position == 'on_up_west_y'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+2*nSq, n+1, n-2*nSq, n-nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_up_north_wall_east_edge'):
                    el.bc = ['W  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([0, n+1, n-2*nSq, n-1])
                    el.bc_con_f = np.array([0, 1, 1, 4])
                elif (position == 'on_up_north_wall'):
                    el.bc = ['W  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([0, n+1, n-2*nSq, n-1])
                    el.bc_con_f = np.array([0, 4, 1, 2])
                elif (position == 'on_up_intern'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+2*nSq, n+1, n-2*nSq, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
            elif ('on_low' in position):    # we are in the lower onion part
                if (position == 'on_low_west_sq_south_x'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, n+nSq*2, n-nel_quarter*3, (k+1)*nSq+nel_quarter*3])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_low_west_sq_north_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, n+nSq*2, n+1, (k+1)*nSq+nel_quarter*3])
                    el.bc_con_f = np.array([2, 4, 1, 4])
                elif (position == 'on_low_west_sq'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, n+nSq*2, n+1, (k+1)*nSq+nel_quarter*3])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_low_south_x_east_wall'):
                    el.bc = ['E  ','W  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, 0, n-nel_quarter*3, n-nSq*2])
                    el.bc_con_f = np.array([3, 0, 1, 4])
                elif (position == 'on_low_south_x'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, n+nSq*2, n-nel_quarter*3, n-nSq*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'on_low_east_wall_north_edge'):
                    el.bc = ['E  ','W  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, 0, n+1, n-nSq*2])
                    el.bc_con_f = np.array([2, 0, 1, 1])
                elif (position == 'on_low_east_wall'):
                    el.bc = ['E  ','W  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, 0, n+1, n-nSq*2])
                    el.bc_con_f = np.array([3, 0, 1, 4])
                elif (position == 'on_low_north_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, n+nSq*2, n+1, n-nSq*2])
                    el.bc_con_f = np.array([2, 4, 1, 2])
                elif (position == 'on_low_intern'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-1, n+nSq*2, n+1, n-nSq*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
            else:
                print('position assignment was not correct!')
                sys.exit(3)

    
def check_position(element, nR, nSq):
    """ Check position of the given element within the quarter region. """

    el = element
    n = el.number
    if (n <= nSq**2):   # we are in the square section
        if (n == 1 or n == nSq or n == nSq**2-nSq+1 or n == nSq**2):    # corners
            if (n == 1):  # we are on the first element on the lower left
                el.pos = 'sq_low_left'
            elif (n == nSq):    # we are on the lower right corner
                el.pos = 'sq_low_right'
            elif (n == nSq**2-nSq+1):     # we are on the upper left corner
                el.pos = 'sq_up_left'
            elif (n == nSq**2):   # last element on the upper right
                el.pos = 'sq_up_right'
            return
        elif (n > nSq**2-nSq or n%nSq == 0 or n < nSq or n%(nSq+1) == 0):  # edges
            if (n > nSq**2-nSq):   # northern row
                el.pos = 'sq_north_row'
            elif ((n%nSq) == 0): # eastern column
                el.pos = 'sq_east_col'
            elif (n < nSq):  # southern row
                el.pos = 'sq_south_row'
            elif ((n%(nSq+1)) == 0):  #  western column
                el.pos = 'sq_west_col'
            return
        else:   # interior
            el.pos = 'sq_internal'
            return
    else:   # we are in the onion region
        i = ((n-1)-nSq**2)%(nSq*2) # position in clockwise manner through each layer
        k = abs(i-((nSq*2)-1))                  # position in anticlockwise manner
        j = int(((n-1)-nSq**2)/(nSq*2)) # onion like layer number, inner one is first
        if (i<nSq):    ## Upper part
            ## 1-sided special treatment
            # southern faces with square section
            if (j==0):
                if (i==0):  # south square and west y=0
                    el.pos = 'on_up_south_sq_west_y'
                    return
                elif (i==nSq-1):    # south square and east edge    
                    el.pos = 'on_up_south_sq_east_edge'
                    return
                else:
                    el.pos = 'on_up_south_sq'
                    return
            # western faces with y=0
            elif (i==0):
                if (j==nR-nSq-1):   # western face y=0 and northern wall
                    el.pos = 'on_up_west_y_north_wall'
                    return
                else:
                    el.pos = 'on_up_west_y'
                    return
            # northern faces at wall
            elif (j==(nR-nSq)-1):
                if (i==nSq-1):  # north wall east edge
                    el.pos = 'on_up_north_wall_east_edge'
                    return
                else:
                    el.pos = 'on_up_north_wall'
                    return
            elif (i==nSq-1):    # east edge
                el.pos = 'on_up_east_edge'
                return
            ## internal upper part
            else:
                el.pos = 'on_up_intern'
                return
        elif (k<nSq):   ## Lower part
            # western faces with square section
            if (j==0):
                if (k==0):  # western face square and south x=0
                    el.pos = 'on_low_west_sq_south_x'
                    return
                elif (k==nSq-1):    # western face square and northern edge
                    el.pos = 'on_low_west_sq_north_edge'
                    return
                else:
                    el.pos = 'on_low_west_sq'
                    return
            # southern faces with x=0
            elif (k==0):
                if (j==nR-nSq-1):   # south x=0 east wall
                    el.pos = 'on_low_south_x_east_wall'
                    return
                else:
                    el.pos = 'on_low_south_x'
                    return
            # eastern faces at wall
            elif (j==(nR-nSq)-1):
                if (k==nSq-1):   # east wall north edge
                    el.pos = 'on_low_east_wall_north_edge'
                    return
                else: 
                    el.pos = 'on_low_east_wall'
                    return
            elif (k==nSq-1):     # north edge
                    el.pos = 'on_low_north_edge'
                    return
            ## interal lower part
            else:
                el.pos = 'on_low_intern'
                return
        else:
            print('Error in Position onion region.')
            sys.exit(2)


def write_mesh(elements):
    """ Write vertex locations to rea file. """
    
    mesh = []
    n_tot = len(elements)
    spatial_dim = 2
    mesh.append('{0:10d} {1:10d} {2:10d} NEL,NDIM,NELV\n'.format(n_tot,spatial_dim,n_tot))
    for el in elements:      # loop through all elements
        x = el.x
        y = el.y
        n = el.number

        mesh.append('{0:>19s} {1:10d} {2:6s}{3:1s}{4:12s}'.format\
            ('ELEMENT',n,'[    1','a',']  GROUP  0\n'))
        mesh.append('{0: 10.6f}{1: 14.6f}{2: 14.6f}{3: 14.6f}   {4:s}'.format\
            (x[0], x[1], x[2], x[3], '\n'))   # x coordinates
        mesh.append('{0: 10.6f}{1: 14.6f}{2: 14.6f}{3: 14.6f}   {4:s}'.format\
            (y[0], y[1], y[2], y[3], '\n'))  # y coordinates

    f = open('base.rea','r')
    contents = f.readlines()
    f.close()
    
    # find row for mesh data
    line_mesh = 0
    while (not 'MESH DATA' in contents[line_mesh]):
            line_mesh = line_mesh + 1
    
    # and write it to rea file
    contents[line_mesh+1:line_mesh+1] = mesh
    contents = "".join(contents)
    f = open('base.rea','w')
    f.write(contents)
    f.close()

def write_bc(elements, nR, nSq):
    """ Write boundary conditions to rea file. """

    bc = []
    dig_n_tot = len(str(elements[-1].number))   # size of element number
    for el in elements:
        for f in range(4):
            bc.append(' {boundary:3s}{current_el:{digits_n_tot}d} {face:2d}   \
{con_el:<07.1f}{con_f:14.5f}{zero1:14.5f}{zero2:14.5f}{zero3:14.5f}    {newline:s}'\
            .format(boundary=el.bc[f], current_el=el.number, digits_n_tot=dig_n_tot, face=(f+1),\
            con_el=el.bc_con_el[f], con_f=el.bc_con_f[f],\
            zero1=0.0,zero2=0.0,zero3=0.0,newline='\n'))

    f = open('base.rea','r')
    contents = f.readlines()
    f.close()
 
    # find row for bc data
    line_bc = 0
    while (not 'FLUID BOUNDARY CONDITIONS' in contents[line_bc]):
        line_bc = line_bc + 1
 
    # and write it to rea file
    contents[line_bc+1:line_bc+1] = bc
    contents = "".join(contents)
    f = open('base.rea','w')
    f.write(contents)
    f.close()



def rea_skel():
    """ Create a skeleton base.rea file. """
    reafile = 'base.rea'
    f = open(reafile, 'w')
    # write some default parameters
    f.write('****** PARAMETERS ******\n')
    f.write('   2.6100     NEKTON VERSION\n')
    f.write('   2 DIMENSIONAL RUN\n')
    f.write('         118 PARAMETERS FOLLOW\n')
    f.write('   1.00000     P001: DENSITY\n')
    f.write('  -5300.00     P002: VISCOSITY\n')
    f.write('   0.00000     P003: BETAG\n')
    f.write('   0.00000     P004: GTHETA\n')
    f.write('   0.00000     P005: PGRADX\n')
    f.write('   0.00000     P006: \n')
    f.write('   1.00000     P007: RHOCP\n')
    f.write('   1.00000     P008: CONDUCT\n')
    f.write('   0.00000     P009: \n')
    f.write('   0.00000     P010: FINTIME\n')
    f.write('   103.000     P011: NSTEPS\n')
    f.write('  -1.00000E-03 P012: DT\n')
    f.write('   20.0000     P013: IOCOMM\n')
    f.write('   0.00000     P014: IOTIME\n')
    f.write('   50.0000     P015: IOSTEP\n')
    f.write('   0.00000     P016: PSSOLVER: 0=default\n')
    f.write('   0.00000     P017: \n')
    f.write('   0.00000     P018: GRID <0 --> # cells on screen\n')
    f.write('   0.00000     P019: INTYPE\n')
    f.write('   0.00000     P020: NORDER\n')
    f.write('   1.00000E-08 P021: DIVERGENCE\n')
    f.write('   1.00000E-08 P022: HELMHOLTZ\n')
    f.write('   0.00000     P023: NPSCAL\n')
    f.write('   1.00000E-02 P024: TOLREL\n')
    f.write('   1.00000E-02 P025: TOLABS\n')
    f.write('   0.50000     P026: COURANT/NTAU\n')
    f.write('   3.00000     P027: TORDER\n')
    f.write('   0.00000     P028: TORDER: mesh velocity (0: p28=p27)\n')
    f.write('   0.00000     P029: = magnetic visc if > 0, = -1/Rm if < 0\n')
    f.write('   0.00000     P030: > 0 ==> properties set in uservp()\n')
    f.write('   0.00000     P031: NPERT: #perturbation modes\n')
    f.write('   0.00000     P032: #BCs in re2 file, if > 0\n')
    f.write('   0.00000     P033: \n')
    f.write('   0.00000     P034: \n')
    f.write('   0.00000     P035: \n')
    f.write('   0.00000     P036: XMAGNET\n')
    f.write('   0.00000     P037: NGRIDS\n')
    f.write('   0.00000     P038: NORDER2\n')
    f.write('   0.00000     P039: NORDER3\n')
    f.write('   0.00000     P040: \n')
    f.write('   0.00000     P041: 1-->multiplicattive SEMG\n')
    f.write('   0.00000     P042: 0=gmres/1=pcg\n')
    f.write('   0.00000     P043: 0=semg/1=schwarz\n')
    f.write('   0.00000     P044: 0=E-based/1=A-based prec.\n')
    f.write('   0.00000     P045: Relaxation factor for DTFS\n')
    f.write('   0.00000     P046: reserved\n')
    f.write('   0.00000     P047: vnu: mesh material prop.\n')
    f.write('   0.00000     P048: \n')
    f.write('   0.00000     P049: \n')
    f.write('   0.00000     P050: \n')
    f.write('   0.00000     P051: \n')
    f.write('   0.00000     P052: IOHIS\n')
    f.write('   0.00000     P053: \n')
    f.write('  -3.00000     P054: fixed flow rate dir: |p54|=1,2,3=x,y,z\n')
    f.write('   1.00000     P055: vol.flow rate (p54>0) or Ubar (p54<0)\n')
    f.write('   0.00000     P056: \n')
    f.write('   0.00000     P057: \n')
    f.write('   0.00000     P058: \n')
    f.write('   0.00000     P059: !=0 --> full Jac. eval. for each el.\n')
    f.write('   0.00000     P060: !=0 --> init. velocity to small nonzero\n')
    f.write('   0.00000     P061: \n')
    f.write('   0.00000     P062: >0 --> force byte_swap for output\n')
    f.write('   0.00000     P063: =8 --> force 8-byte output\n')
    f.write('   0.00000     P064: =1 --> perturbation restart\n')
    f.write('   0.00000     P065: #iofiles (eg, 0 or 64); <0 --> sep. dirs\n')
    f.write('   6.00000     P066: output : <0=ascii, else binary\n')
    f.write('   6.00000     P067: restart: <0=ascii, else binary\n')
    f.write('   10.0000     P068: STAT_COMP: how often you compute stats\n')
    f.write('   50.0000     P069: STAT_OUTP: how often you write stats\n')
    f.write('   50.0000     P070: CHKPTSTP: how often you write restart files (rs8)\n')
    f.write('   0.00000     P071: IFCHKPTRST: restart (1) or start from inital cond.\n')
    f.write('   0.00000     P072: \n')
    f.write('   0.00000     P073: \n')
    f.write('   0.00000     P074: \n')
    f.write('   0.00000     P075: \n')
    f.write('   0.00000     P076: \n')
    f.write('   0.00000     P077: \n')
    f.write('   0.00000     P078: \n')
    f.write('   0.00000     P079: \n')
    f.write('   0.00000     P080: \n')
    f.write('   0.00000     P081: \n')
    f.write('   0.00000     P082: \n')
    f.write('   0.00000     P083: \n')
    f.write('   0.00000     P084: != 0 --> sets initial timestep if p12>0\n')
    f.write('   0.00000     P085: dt retio of p84 !=0, for timesteps>0\n')
    f.write('   0.00000     P086: reserved\n')
    f.write('   0.00000     P087: \n')
    f.write('   0.00000     P088: \n')
    f.write('   0.00000     P089: \n')
    f.write('   0.00000     P090: \n')
    f.write('   0.00000     P091: \n')
    f.write('   0.00000     P092: \n')
    f.write('   20.0000     P093: Number of previous pressure solns saved\n')
    f.write('   9.00000     P094: start projecting velocity after p94 step\n')
    f.write('   9.00000     P095: start projecting pressure after p95 step\n')
    f.write('   0.00000     P096: \n')
    f.write('   0.00000     P097: \n')
    f.write('   0.00000     P098: \n')
    f.write('   4.00000     P099: dealiasing: <0--> off /3--> old /4-->new\n')
    f.write('   0.00000     P100: \n')
    f.write('   0.00000     P101: Number of additional modes to filter\n')
    f.write('   1.00000     P102: Dump out divergence at each time step\n')
    f.write('   0.01000     P103: weight of stabilizing filter\n')
    f.write('   0.00000     P104: \n')
    f.write('   0.00000     P105: \n')
    f.write('   0.00000     P106: \n')
    f.write('   0.00000     P107: !=0 --> add h2 array in hmholtz eqn\n')
    f.write('   0.00000     P108: \n')
    f.write('   0.00000     P109: \n')
    f.write('   0.00000     P110: \n')
    f.write('   0.00000     P111: \n')
    f.write('   0.00000     P112: \n')
    f.write('   0.00000     P113: \n')
    f.write('   0.00000     P114: \n')
    f.write('   0.00000     P115: \n')
    f.write('   0.00000     P116: =nelx for gtp solver\n')
    f.write('   0.00000     P117: =nely for gtp solver\n')
    f.write('   0.00000     P118: =nelz for gtp solver\n')
    f.write('      4  Lines of passive scalar data follows2 CONDUCT, 2RHOCP\n')
    f.write('   1.00000        1.00000        1.00000        1.00000        1.00000\n')
    f.write('   1.00000        1.00000        1.00000        1.00000\n')
    f.write('   1.00000        1.00000        1.00000        1.00000        1.00000\n')
    f.write('   1.00000        1.00000        1.00000        1.00000\n')
    f.write('         13   LOGICAL SWITCHES FOLLOW\n')
    f.write(' T      IFFLOW\n')
    f.write(' F      IFHEAT\n')
    f.write(' T      IFTRAN\n')
    f.write(' T F F F F F F F F F F  IFNAV & IFADVC (convection in P.S. fields)\n')
    f.write(' F F T T T T T T T T T T  IFTMSH (IF mesh for this field is T mesh)\n')
    f.write(' F      IFAXIS\n')
    f.write(' F      IFSTRS\n')
    f.write(' F      IFSPLIT\n')
    f.write(' F      IFMGRID\n')
    f.write(' F      IFMODEL\n')
    f.write(' F      IFKEPS\n')
    f.write(' F      IFMVBD\n')
    f.write(' F      IFCHAR\n')
    f.write('   2.00000       2.00000      -1.00000      -1.00000     XFAC,YFAC,XZERO,YZERO\n')
    f.write('  ***** MESH DATA *****  6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8\n')
    f.write(' ***** CURVED SIDE DATA *****\n')
    f.write('       0 Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n')
    f.write('  ***** BOUNDARY CONDITIONS *****\n')
    f.write('  ***** FLUID BOUNDARY CONDITIONS *****\n')
    f.write('   ***** NO THERMAL BOUNDARY CONDITIONS *****\n')
    f.write('    0 PRESOLVE/RESTART OPTIONS  *****\n')
    f.write('    7         INITIAL CONDITIONS *****\n')
    f.write(' C Default\n')
    f.write(' C Default\n')
    f.write(' C Default\n')
    f.write(' C Default\n')
    f.write(' C Default\n')
    f.write(' C Default\n')
    f.write(' C Default\n')
    f.write('   ***** DRIVE FORCE DATA ***** BODY FORCE, FLOW, Q\n')
    f.write('    4                 Lines of Drive force data follow\n')
    f.write(' C\n')
    f.write(' C\n')
    f.write(' C\n')
    f.write(' C\n')
    f.write('   ***** Variable Property Data ***** Overrrides Parameter data.\n')
    f.write('    1 Lines follow.\n')
    f.write('    0 PACKETS OF DATA FOLLOW\n')
    f.write('   ***** HISTORY AND INTEGRAL DATA *****\n')
    f.write('    0   POINTS.  Hcode, I,J,H,IEL\n')
    f.write('   ***** OUTPUT FIELD SPECIFICATION *****\n')
    f.write('    6 SPECIFICATIONS FOLLOW\n')
    f.write('    T      COORDINATES\n')
    f.write('    T      VELOCITY\n')
    f.write('    T      PRESSURE\n')
    f.write('    T      TEMPERATURE\n')
    f.write('    F      TEMPERATURE GRADIENT\n')
    f.write('    0      PASSIVE SCALARS\n')
    f.write('   ***** OBJECT SPECIFICATION *****\n')
    f.write('        0 Surface Objects\n')
    f.write('        0 Volume  Objects\n')
    f.write('        0 Edge    Objects\n')
    f.write('        0 Point   Objects\n')
    f.close()


def check_mesh_quality(elements, nR, nSq):
    """ Find minimum and maximum radial and circumferential 
    element lengths and element angles (distortion from 90°).
    """
    
    # only check first quadrant
    nel_quarter = nSq**2 + (nR-nSq)*2*nSq       # number of elements in one quarter
    elements = elements[0:nel_quarter]
    l_r_max = 0
    l_r_min = 1e5
    l_p_max = 0
    l_p_min = 1e5
    alph_max = m.pi/4
    alph_min = m.pi/2


    for el in elements:
        n = el.number
        i = ((n-1)-nSq**2)%(nSq*2) # position in clockwise manner through each layer


        if (n <= nSq**2 or i < nSq):    # either in "square" section or upper onion part
            vec1 = np.array([el.x[1]-el.x[0], el.y[1]-el.y[0]])  # corresponds to face 1
            vec1_norm = np.linalg.norm(vec1)
            vec2 = np.array([el.x[2]-el.x[1], el.y[2]-el.y[1]])  # corresponds to face 2
            vec2_norm = np.linalg.norm(vec2)
            vec3 = np.array([el.x[3]-el.x[2], el.y[3]-el.y[2]])  # corresponds to face 3
            vec3_norm = np.linalg.norm(vec3)
            vec4 = np.array([el.x[0]-el.x[3], el.y[0]-el.y[3]])  # corresponds to face 4
            vec4_norm = np.linalg.norm(vec4)

            l_rad = np.array([ vec4_norm, vec2_norm ])
            l_phi = np.array([ vec1_norm, vec3_norm ])


            alpha_12 = my_math.vec_angle(-vec1, vec2)
            alpha_23 = my_math.vec_angle(vec2, -vec3)
            alpha_34 = my_math.vec_angle(-vec3, vec4)
            alpha_41 = my_math.vec_angle(-vec4, vec1)
 

            alpha = np.array([ alpha_12, alpha_23, alpha_34, alpha_41 ])
#            print(alpha/m.pi*180, np.sum(alpha)/m.pi*180, n)

            l_rad_max = max(l_rad)
            l_rad_min = min(l_rad)
            l_phi_max = max(l_phi)
            l_phi_min = min(l_phi)

            alpha_max = max(alpha)
            alpha_min = min(alpha)

            # update previous values if necessary
            if (l_rad_max > l_r_max):
                l_r_max = l_rad_max
                el_r_max = n
            if (l_rad_min < l_r_min):
                l_r_min = l_rad_min
                el_r_min = n
            if (l_phi_max > l_p_max):
                l_p_max = l_phi_max
                el_p_max = n
            if (l_phi_min < l_p_min):
                l_p_min = l_phi_min
                el_p_min = n
            if (alpha_max > alph_max):
                alph_max = alpha_max
                el_alph_max = n
            if (alpha_min < alph_min):
                alph_min = alpha_min
                el_alph_min = n


    # Write a little output to stdout
    print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
    print('Some information about the mesh size:')
    print('Delta R max = {0:12.5f} at {2:d}\nDelta R min = {1:12.5f} at {3:d}'.format(l_r_max, l_r_min, el_r_max, el_r_min))
    print('Delta phi max = {0:10.5f} at {2:d}\nDelta phi min = {1:10.5f} at {3:d}'.format(l_p_max, l_p_min, el_p_max, el_p_min))
    print('alpha max = {0:10.5f}° at {2:d}\nalpha min = {1:10.5f}° at {3:d}'.format(alph_max/m.pi*180, alph_min/m.pi*180, el_alph_max, el_alph_min))
    print('Note that curvature is not considered here!')
    print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')





