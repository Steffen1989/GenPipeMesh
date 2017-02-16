# A collection of important functions
import sys
import pdb
import numpy as np
import math as m


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# The vertices are set in a special way. For now the inner section, 
# called square section, is just a regular square. The outer part, 
# called onion region, is built up by ellipses and straight lines.
# The semi-major axis a is decreasing each layer outwards so that a=1
# in the outermost layer which gives a circle. The semi-minor axis is 
# kept at b=1. The radius (constant) c is changed according to the 
# radius at the y axis (x=0) 
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

def set_vertices(elements,nR,nSq,dr):
    """ Set vertex location for each element. 

    The vertices are set in a special way. For now the inner section, 
    called square section, is just a regular square. The outer part, 
    called onion region, is built up by ellipses and straight lines.
    The semi-major axis a is decreasing each layer outwards so that a=1
    in the outermost layer which gives a circle. The semi-minor axis is 
    kept at b=1. The radius (constant) c is changed according to the 
    radius at the y axis (x=0) 
    """
    
    # Variable definitions
    ntheta = nSq*2  # number of elements in one onion layer
    rad_on = np.zeros(2)       # "radius": constant in ellipse equation
    rad_row = np.zeros(2)       # "radius": constant in ellipse equation
    rad_col = np.zeros(2)       # "radius": constant in ellipse equation
    semi_major_on = np.zeros(2)    # semi-major axis
    semi_major_row = np.zeros(2)    
    semi_major_col = np.zeros(2)
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
            b_sq = 1
            drop_sq_max = 0.80         # max drop in percent compared to all squares
            n_ellip_sq = nSq  # number of ellipses
            x_max_sq = nSq*dr
            r_max_sq = x_max_sq
            y_max_sq = x_max_sq*drop_sq_max
            a_min_sq = x_max_sq*b_sq/( (r_max_sq**2*b_sq**2-y_max_sq**2)**0.5 )

            # Idea: define ellipses in the inner region such that they coincide with the 
            # element in the onion region.
            # OPEN def semi-major axis
            #--------------------------------------------------
            slope_row[0] = m.tan(m.pi/2*(j/ntheta))  # slope of the straight line on the bottom side
            slope_row[1] = m.tan(m.pi/2*((j+1)/ntheta)) # slope of the straight line on the top side
            slope_col[0] = m.tan(m.pi/2*((ntheta-i)/ntheta))  # slope of the straight line on the left side
            slope_col[1] = m.tan(m.pi/2*((ntheta-i-1)/ntheta)) # slope of the straight line on the right side

            rad_row[0] = j*dr    # small "radius"  
            rad_row[1] = (j+1)*dr   # large "radius"
            rad_col[0] = i*dr    # small "radius"  
            rad_col[1] = (i+1)*dr   # large "radius"

            x_inters_row[0] = intersec_ellip_line(b_sq,a_min_sq,r_max_sq**2,slope_row[0],0)
            x_inters_row[1] = intersec_ellip_line(b_sq,a_min_sq,r_max_sq**2,slope_row[1],0)
            y_inters_row[0] = line(slope_row[0],x_inters_row[0],0)
            y_inters_row[1] = line(slope_row[1],x_inters_row[1],0)

            x_inters_col[0] = intersec_ellip_line(a_min_sq,b_sq,r_max_sq**2,slope_col[0],0)
            x_inters_col[1] = intersec_ellip_line(a_min_sq,b_sq,r_max_sq**2,slope_col[1],0)
            y_inters_col[0] = line(slope_col[0],x_inters_col[0],0)
            y_inters_col[1] = line(slope_col[1],x_inters_col[1],0)

            if (j==0):
                semi_major_row[0] = 0   # this is reset later
                semi_major_row[1] = x_inters_row[1]*b_sq/( (rad_row[1]**2*b_sq**2 - y_inters_row[1]**2)**0.5 )
            else: 
                semi_major_row[0] = x_inters_row[0]*b_sq/( (rad_row[0]**2*b_sq**2 - y_inters_row[0]**2)**0.5 )
                semi_major_row[1] = x_inters_row[1]*b_sq/( (rad_row[1]**2*b_sq**2 - y_inters_row[1]**2)**0.5 )

            if (i==0):                # note that x and y need to be switched here
                semi_major_col[0] = 0   # this is reset later 
                semi_major_col[1] = y_inters_col[1]*b_sq/( (rad_col[1]**2*b_sq**2 - x_inters_col[1]**2)**0.5 )
            else:
                semi_major_col[0] = y_inters_col[0]*b_sq/( (rad_col[0]**2*b_sq**2 - x_inters_col[0]**2)**0.5 )
                semi_major_col[1] = y_inters_col[1]*b_sq/( (rad_col[1]**2*b_sq**2 - x_inters_col[1]**2)**0.5 )

            # CLOSE def semi-major axis
            #--------------------------------------------------

            if (j==0): # first row
                if (i==0):  # first col
                    x2 = intersec_ellip_ellip(semi_major_row[1],b_sq,rad_row[1]**2,\
                            b_sq,semi_major_col[1],rad_col[1]**2)
                    y2 = ellipse(semi_major_row[1],b_sq,rad_row[1]**2,x2)
                    el.x = np.array([i*dr, (i+1)*dr, x2, i*dr])
                    el.y = np.array([j*dr, j*dr, y2, (j+1)*dr])
                else:
                    x2 = intersec_ellip_ellip(semi_major_row[1],b_sq,rad_row[1]**2,\
                            b_sq,semi_major_col[1],rad_col[1]**2)
                    x3 = intersec_ellip_ellip(semi_major_row[1],b_sq,rad_row[1]**2,\
                            b_sq,semi_major_col[0],rad_col[0]**2)
                    y2 = ellipse(semi_major_row[1],b_sq,rad_row[1]**2,x2)
                    y3 = ellipse(semi_major_row[1],b_sq,rad_row[1]**2,x3)
                    el.x = np.array([i*dr, (i+1)*dr, x2, x3])
                    el.y = np.array([j*dr, j*dr, y2, y3])
            elif (j>0 and i==0): # first col
                x1 = intersec_ellip_ellip(semi_major_row[0],b_sq,rad_row[0]**2,\
                        b_sq,semi_major_col[1],rad_col[1]**2)
                x2 = intersec_ellip_ellip(semi_major_row[1],b_sq,rad_row[1]**2,\
                        b_sq,semi_major_col[1],rad_col[1]**2)
                y1 = ellipse(semi_major_row[0],b_sq,rad_row[0]**2,x1)
                y2 = ellipse(semi_major_row[1],b_sq,rad_row[1]**2,x2)
                el.x = np.array([i*dr, x1, x2, i*dr])
                el.y = np.array([j*dr, y1, y2, (j+1)*dr])
            elif (i> 0 and j>0):    # inside
                #find intersection between both ellipses
                x0 = intersec_ellip_ellip(semi_major_row[0],b_sq,rad_row[0]**2,\
                        b_sq,semi_major_col[0],rad_col[0]**2)
                x1 = intersec_ellip_ellip(semi_major_row[0],b_sq,rad_row[0]**2,\
                        b_sq,semi_major_col[1],rad_col[1]**2)
                x2 = intersec_ellip_ellip(semi_major_row[1],b_sq,rad_row[1]**2,\
                        b_sq,semi_major_col[1],rad_col[1]**2)
                x3 = intersec_ellip_ellip(semi_major_row[1],b_sq,rad_row[1]**2,\
                        b_sq,semi_major_col[0],rad_col[0]**2)
                y0 = ellipse(semi_major_row[0],b_sq,rad_row[0]**2,x0)
                y1 = ellipse(semi_major_row[0],b_sq,rad_row[0]**2,x1)
                y2 = ellipse(semi_major_row[1],b_sq,rad_row[1]**2,x2)
                y3 = ellipse(semi_major_row[1],b_sq,rad_row[1]**2,x3)
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
            j = int(((el.number-1)-nSq**2)/(nSq*2)) # onion like layer number, inner one is first, start from j=0
            l = (nR - nSq) - (j+1)                              # onion like layer number, outer one is last l=0

            # Determine the semi-major axis for "onion" section
            #----------------------------------------------------------------------
            a_min_on = 1
            b_on = 1

            # determine the maxim semi-major axis at the edge to the square region
            r_min_on = (nSq)*dr

            # a_max is found at the lowest onion region where the elements at the border
            # need to be outside of the square region
            el_square = elements[nSq**2-1]  # last element in square region
            # in square region
            x_min_on = el_square.x[2]
            y_min_on = x_min_on
            a_max_on = x_min_on*b_on/( (r_min_on**2*b_on**2-y_min_on**2)**(0.5) )
            
            a_max_on = a_min_sq
            semi_major_on[0] = geom_prog(nR-nSq, a_min_on, a_max_on, j)
            semi_major_on[1] = geom_prog(nR-nSq, a_min_on, a_max_on, j+1)
            semi_major_on[0] = lin_dist(nR-nSq, a_min_on, a_max_on, j)
            semi_major_on[1] = lin_dist(nR-nSq, a_min_on, a_max_on, j+1)
            semi_major_on[0] = quad_dist(nR-nSq, a_min_on, a_max_on, j)
            semi_major_on[1] = quad_dist(nR-nSq, a_min_on, a_max_on, j+1)


            rad_on[0] = (j+nSq)*dr    
            rad_on[1] = (j+1+nSq)*dr
            slope_on[0] = m.tan(m.pi/2*(k/ntheta))  # slope of the straight line on the right side
            # of the element (upper part) or bottom side (lower part)
            slope_on[1] = m.tan(m.pi/2*((k+1)/ntheta)) # slope of the straight line on the left side
            # of the element (upper part) or top side (lower part)
            if (i <= (nSq-1)):  # upper part, including border /
                x0 = intersec_ellip_line(semi_major_on[0],1,rad_on[0]**2,slope_on[1],0)
                x1 = intersec_ellip_line(semi_major_on[0],1,rad_on[0]**2,slope_on[0],0)
                x2 = intersec_ellip_line(semi_major_on[1],1,rad_on[1]**2,slope_on[0],0)
                x3 = intersec_ellip_line(semi_major_on[1],1,rad_on[1]**2,slope_on[1],0)
                el.x = np.array([x0, x1, x2, x3])
                el.y[0:2] = ellipse(semi_major_on[0],1,rad_on[0]**2,el.x[0:2])
                el.y[2:4] = ellipse(semi_major_on[1],1,rad_on[1]**2,el.x[2:4])
            elif (i >= nSq):     # lower part, including border /
                x0 = intersec_ellip_line(1,semi_major_on[0],rad_on[0]**2,slope_on[0],0)
                x1 = intersec_ellip_line(1,semi_major_on[1],rad_on[1]**2,slope_on[0],0)
                x2 = intersec_ellip_line(1,semi_major_on[1],rad_on[1]**2,slope_on[1],0)
                x3 = intersec_ellip_line(1,semi_major_on[0],rad_on[0]**2,slope_on[1],0)
                y0 = line(slope_on[0],x0,0)
                y1 = line(slope_on[0],x1,0)
                y2 = line(slope_on[1],x2,0)
                y3 = line(slope_on[1],x3,0)
                el.y = np.array([y0, y1, y2, y3])
                el.x = np.array([x0, x1, x2, x3])


def lin_dist(N, a_min, a_max, j):
    """ Linear distribution 

    N     : number of steps
    a_min : starting value at j=1
    a_max : maximum value at j=N
    j     : step
    """

    ret = a_max - (a_max-a_min)/(N) * (j)
    return ret


def quad_dist(N, a_min, a_max, j):
    """ Quadratic distribution 

    N     : number of steps
    a_min : starting value at j=1
    a_max : maximum value at j=N
    j     : step
    """

    ret = (a_max-a_min)/(N**2) * (j-N)**2 + a_min
    return ret



def geom_prog(N, a_min, a_max, j):
    """ Geometric progression a(j) = a_min * r**j 

    with a(0) = a_max 
    and a(N) = a_min*r**N = a_max

    N     : number of steps
    a_min : starting value
    a_max : maximum value
    j     : step
    """
    r = m.exp(m.log(a_min/a_max)/(N))
    ret  = a_max * r**(j)

    return ret


def ellipse(a,b,c,x):
    """ Ellipse with x**2/a**2 + y**2/b**2 = c 
    
    a : semi-major axis
    b : semi-minor axis
    c : "radius" or constant on the rhs
    x : dependent variable
    """

    ret = b*(c-x**2/a**2)**(0.5)
    return ret


def line(m,x,d):
    """ Straight line with y = m*x + d 
    
    m : slope
    d : y-intercept
    """

    ret = m*x+d
    return ret


def intersec_ellip_line(a,b,c,m,d):
    """ Intersection between ellipse and straight line. 
    
    a : semi-major axis
    b : semi-minor axis
    c : "radius" or constant on the rhs 
    m : slope of straight line
    d : y-intercept of straight line
    """

    # x value from b(c-x**2/a**2)**(0.5) = m*x + d
    A = 2*m*d/(m**2+b**2/a**2)
    B = (d**2-b**2*c)/(m**2+b**2/a**2)
    x = - A/2 + ( A**2/4 - B )**(0.5)
    return x


def intersec_ellip_ellip(a1,b1,rhs1,a2,b2,rhs2):
    """ Intersection between two ellipses.

    a1   : semi-major axis of ellipse 1
    b1   : semi-minor axis of ellipse 1
    rhs1 : constant on the rhs
    a2   : semi-major axis of ellipse 2
    b2   : semi-minor axis of ellipse 2
    rhs2 : constant on the rhs

    x**2/a**2 + y**2/b**2 = R**2
    """

    # x value from b1*(R1**2-x**2/a1**2)**0.5 = b2*(R2**2-x**2/a2**2)**0.5
    x = ( ( b1**2*rhs1 - b2**2*rhs2 ) / ( (b1/a1)**2 - (b2/a2)**2 ) )**0.5
#    pdb.set_trace()
    return x


def set_bc(elements,nR,nSq):
    """ Set boundary conditions for each face. """

    for el in elements:
        n = el.number
#         position = check_position(n, nSq)
        check_position(el,nR,nSq)
        position = el.pos

        i = ((n-1)-nSq**2)%(nSq*2) # position in clockwise manner through each layer
        k = abs(i-((nSq*2)-1))                  # position in anticlockwise manner
        j = int(((n-1)-nSq**2)/(nSq*2)) # onion like layer number, inner one is first

        if (n <= nSq**2):   # we are in the square section
            if (position == 'sq_low_left'):
                el.bc = ['W  ','E  ','E  ','W  ']
                el.bc_con_el = [0, n+1, n+nSq, 0]
                el.bc_con_f = [0, 4, 1, 0]
            elif (position == 'sq_low_right'):
                el.bc = ['W  ','W  ','E  ','E  ']
                el.bc_con_el = [0, n+1, n+nSq, n-1] 
                el.bc_con_f = [0, 4, 1, 2]
            elif (position == 'sq_up_left'):
                el.bc = ['E  ','E  ','W  ','W  ']
                el.bc_con_el = [n-nSq, n+1, 0, 0]           
                el.bc_con_f = [3, 4, 0, 0]
            elif (position == 'sq_up_right'):
                el.bc = ['E  ','W  ','W  ','E  ']
                el.bc_con_el = [n-nSq, 0, 0, n-1]
                el.bc_con_f = [3, 0, 0, 2]
            elif (position == 'sq_north_row'):
                el.bc = ['E  ','E  ','W  ','E  ']
                el.bc_con_el = [n-nSq, n+1, 0, n-1]
                el.bc_con_f = [3, 4, 0, 2]
            elif (position == 'sq_east_col'):
                el.bc = ['E  ','W  ','E  ','E  ']
                el.bc_con_el = [n-nSq, 0, n+nSq, n-1]
                el.bc_con_f = [3, 0, 1, 2]
            elif (position == 'sq_south_row'):
                el.bc = ['W  ','E  ','E  ','E  ']
                el.bc_con_el = [0, n+1, n+nSq, n-1]
                el.bc_con_f = [0, 4, 1, 2]
            elif (position == 'sq_west_col'):
                el.bc = ['E  ','E  ','E  ','W  ']
                el.bc_con_el = [n-nSq, n+1, n+nSq, 0]
                el.bc_con_f = [3, 4, 1, 0]
            elif (position == 'sq_internal'):
                el.bc = ['E  ','E  ','E  ','E  ']
                el.bc_con_el = [n-nSq, n+1, n+nSq, n-1]
                el.bc_con_f = [3, 4, 1, 2]
            else:
                print('position of element not found!')
                sys.exit(1)
        else:                       # we are in the outer onion like region :-)
            if ('on_up' in position):   # we are in the upper onion part
                if (position == 'on_up_south_sq_west_y'):
                    el.bc = ['E  ','E  ','E  ','W  ']
                    el.bc_con_el = [n-nSq, n+1, n+nSq, 0]
                    el.bc_con_f = [3, 4, 1, 0]
                elif (position == 'on_up_south_sq_east_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = [n-nSq, n+1, n+nSq, n-1]
                    el.bc_con_f = [3, 4, 1, 2]
                elif (position == 'on_up_east_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = [n-nSq, n+1, n+nSq, n-1]
                    el.bc_con_f = [3, 3, 1, 2]
                elif (position == 'on_up_south_sq'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = [n-nSq, n+1, n+nSq, n-1]
                    el.bc_con_f = [3, 3, 1, 2]
                elif (position == 'on_up_west_y_north_wall'):
                    el.bc = ['E  ','E  ','W  ','W  ']
                    el.bc_con_el = [n-nSq, n+1, 0, 0]
                    el.bc_con_f = [3, 4, 0, 0]
                elif (position == 'on_up_west_y'):
                    el.bc = ['E  ','E  ','E  ','W  ']
                    el.bc_con_el = [n-nSq, n+1, n+nSq, 0]
                    el.bc_con_f = [3, 4, 1, 0]
                elif (position == 'on_up_north_wall_east_edge'):
                    el.bc = ['E  ','E  ','W  ','E  ']
                    el.bc_con_el = [n-nSq, n+1, 0, n-1]
                    el.bc_con_f = [3, 3, 0, 2]
                elif (position == 'on_up_north_wall'):
                    el.bc = ['E  ','E  ','W  ','E  ']
                    el.bc_con_el = [n-nSq, n+1, 0, n-1]
                    el.bc_con_f = [3, 4, 0, 2]
                elif (position == 'on_up_intern'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = [n-nSq, n+1, n+nSq, n-1]
                    el.bc_con_f = [3, 4, 1, 2]
            elif ('on_low' in position):    # we are in the lower onion part
                if (position == 'on_low_west_sq_south_x'):
                    el.bc = ['W  ','E  ','E  ','E  ']
                    el.bc_con_el = [0, n+nSq*2, n-1, (k+1)*nSq]
                    el.bc_con_f = [0, 4, 1, 2]
                elif (position == 'on_low_west_sq_north_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = [n+1, n+nSq*2, n-1, (k+1)*nSq]
                    el.bc_con_f = [3, 4, 2, 2]
                elif (position == 'on_low_west_sq'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = [n+1, n+nSq*2, n-1, (k+1)*nSq]
                    el.bc_con_f = [3, 4, 1, 2]
                elif (position == 'on_low_south_x_east_wall'):
                    el.bc = ['W  ','W  ','E  ','E  ']
                    el.bc_con_el = [0, 0, n-1, n-nSq*2]
                    el.bc_con_f = [0, 0, 1, 2]
                elif (position == 'on_low_south_x'):
                    el.bc = ['W  ','E  ','E  ','E  ']
                    el.bc_con_el = [0, n+nSq*2, n-1, n-nSq*2]
                    el.bc_con_f = [0, 4, 1, 2]
                elif (position == 'on_low_east_wall_north_edge'):
                    el.bc = ['E  ','W  ','E  ','E  ']
                    el.bc_con_el = [n+1, 0, n-1, n-nSq*2]
                    el.bc_con_f = [3, 0, 2, 2]
                elif (position == 'on_low_east_wall'):
                    el.bc = ['E  ','W  ','E  ','E  ']
                    el.bc_con_el = [n+1, 0, n-1, n-nSq*2]
                    el.bc_con_f = [3, 0, 1, 2]
                elif (position == 'on_low_north_edge'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = [n+1, n+nSq*2, n-1, n-nSq*2]
                    el.bc_con_f = [3, 2, 2, 2]
                elif (position == 'on_low_intern'):
                    el.bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = [n+1, n+nSq*2, n-1, n-nSq*2]
                    el.bc_con_f = [3, 4, 1, 2]
            else:
                print('position assignment was not correct!')
                os.exit(3)

    
def check_position(element, nR, nSq):
    """ Check position of the given element within the square region. """

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
                elif (k==nSq-1):     # north edge
                    el.pos = 'on_low_north_edge'
                    return
                else:
                    el.pos = 'on_low_east_wall'
                    return
                        ## interal lower part
            else:
                el.pos = 'on_low_intern'
                return
        else:
            print('Error in Position onion region.')
            os.exit(2)


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
    f.write('  -1.00000     P054: fixed flow rate dir: |p54|=1,2,3=x,y,z\n')
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
