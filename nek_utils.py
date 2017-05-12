## A collection of important functions
#----------------------------------------------------------------------
# The idea is to collect all tasks in separate (small) functions to
# clearly arange the code. 

## Import modules
#----------------------------------------------------------------------
import sys
import pdb
import numpy as np
import math as m
import my_math
import elementclass

## Function definitions:
#----------------------------------------------------------------------

# Distribution of points:
#----------------------------------------------------------------------
# Necessary for setting the vertex positions
def sq_dist(nSq, dr, dr_sq_ratio):
    """ Set distribution of elements within the "square" region in 
    radial direction in a certain way.
    Use geometric progression.
    """

    fact_sq = dr_sq_ratio**(1/(nSq-1))
    fact_sq_sum = 0
    for i in range(0,nSq):
        fact_sq_sum = fact_sq_sum + fact_sq**i
    dr_sq_max = nSq*dr/fact_sq_sum
    dr_sq_min = dr_sq_max*dr_sq_ratio
    dr_sq = np.zeros(nSq)
    for i in range(0,nSq):
        dr_sq[i] = my_math.geom_prog(nSq, dr_sq_max, dr_sq_min, i)
    return dr_sq


def on_dist(nR, nSq, dr, dr_sq, distri_on):
    """ Set distribution of elements in radial direction along axis 
    within onion region in a certain way.
    Use geometric progression for a small increase and then cosine for 
    sharp decrease of element size.
    """
    
    n_on_1 = int(m.floor(distri_on*(nR-nSq)))     # increasing region
    dr_sq_min = min(dr_sq)
    dr_on_interface = dr_sq_min

    if (n_on_1<2):  
        # If it is too small then only shrink 
        # Note that condition of total length being (nR-nSq)*dr
        # still needs to be fulfilled. Besides, we keep the first
        # element's size at a certain percentage of the adjacent square region's 
        # elements.
        # Hence, the only free parameter is the smallest element
        # close to the wall.
        n_on_2 = (nR-nSq)
        first_shrink = 1.0
        first_el = dr_sq_min*first_shrink

#        def on_shrink(x):   
#            """ This functions purpose is to ensure that the total 
#            length of the onion region in radial direction is correct.
#            It is needed for Newton-Raphson later on.
#            """
#
#            sum_sin = 0
#            for i in range(0, n_on_2):
#                sum_sin = sum_sin + m.sin( m.pi/2 * i/(n_on_2-1) )
#            ret = -dr*(nR-nSq) + n_on_2*first_el - (first_el -x )*\
#                    sum_sin
#            return ret
#
#        def d_on_shrink(x):   
#            """ This is just the first derivative of on_shrink(x).
#            """
#            sum_sin = 0
#            for i in range(0, n_on_2):
#                sum_sin = sum_sin + m.sin( m.pi/2 * i/(n_on_2-1) )
#            ret = sum_sin
#
#            return ret
#        # initial guess is somewhat arbitrarily set to first_el/10
#        dr_on_wall = my_math.newton_raphson(first_el/10, on_shrink, d_on_shrink)
#        pdb.set_trace()
        sum_sin = 0
        for i in range(0,n_on_2):
            sum_sin = sum_sin + m.sin(m.pi/2 *i/(n_on_2-1))
        dr_on_wall = (dr*(nR-nSq) - n_on_2*first_el)/(sum_sin) + first_el
        dr_on = np.zeros(n_on_2)
        for i in range(0,n_on_2):
            dr_on[i] = my_math.sin_dist(n_on_2, first_el, dr_on_wall,i)
    else:
        # Use increasing region first and afterwards decreasing
        n_on_2 = (nR-nSq) - n_on_1              #decreasing region
        dr_sq_min = min(dr_sq)
        dr_on_interface = dr_sq_min
        
        # The critical part here is to use geometric progression first and 
        # cosine distribution afterwards. However, the overall length needs 
        # to be the same as with the nominal values (nR-nSq)*dr.
        # In order to fulfill this condition, the following two functions 
        # are needed.
        def x_transition(x):
            """ This function is defined by the requirement that the total length
            of both onion regions needs to be (nR-nSq)*dr and we end up at r = R.
            Note that the first element of the second region is not the same as 
            the last of the first region, i.e. both regions "share" dr_transition.
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
        
        # Find the size of the element in between the increasing and decreasing regions
        # by the requirement that the total size of the onion region is still
        # (nR-nSq)*dr by using newton raphson algorithm for finding roots.
        dr_on_transition = my_math.newton_raphson(dr,x_transition, x_transition_prime)
        dr_on_1 = np.zeros(n_on_1)  # size distribution in increasing region
        for i in range(0,n_on_1):
            dr_on_1[i] = my_math.geom_prog(n_on_1,dr_on_interface, dr_on_transition, i)
        dr_on_2 = np.zeros(n_on_2)  # size distribution in decreasing region
        for i in range(0,n_on_2):
            dr_on_2[i] = dr_on_transition*m.cos((i+1)/(n_on_2+1)*m.pi/2)
        dr_on = np.concatenate([dr_on_1, dr_on_2])

    # make sure that last element is exactly at R and not only close because
    # of rounding errors
    desired_radius = np.sum(dr_sq)+dr*(nR-nSq)
    actual_radius = np.sum(dr_sq)+np.sum(dr_on)
    if (abs(desired_radius - actual_radius)>0.1):
        print('WARNING: desired radius is {0:10.5f}, but actual radius is {1:10.5f}'\
                .format(desired_radius, actual_radius))
        sys.exit(1)
    else:
        # adjust last element's size
        dr_last = desired_radius - (np.sum(dr_sq[:])+np.sum(dr_on[:-1]))
    dr_on[-1] = dr_last

    return dr_on


# MAJOR FUNCTION for setting the vertex positions:
#----------------------------------------------------------------------
def set_vertices(elements, nR, nSq, dr, dz, dr_sq_ratio, dr_sq_int_ratio, \
        stretch_sq, distri_on, a_interf, tog_r_out_const, tog_a_on_dist):
    """ Set vertex location for each element. 

    The vertices are set in a special way. 
    The inner section, called "square" section is put together by 
    ellipses and straight lines as a modification of just a simple regular
    mesh consisting of squares. 
    The outer part, called "onion" region, is built up by ellipses and
    straight lines, too.
    In the upper onion region, the semi-major axis a is decreasing each 
    layer outwards so that a=1 in the outermost layer which gives a circle. 
    The constant on the right hand side is kept constant and semi-minor 
    and semi-major axis are varied.
    """
    
    # Variable definitions
    ntheta = nSq*2  # number of elements in one onion layer
    r_const = 1         # constant in ellipses
    a_on = np.zeros(2)    # semi-major axis in onion region
    a_row = np.zeros(2)   # semi-maj. ax. in sq. along rows 
    a_col = np.zeros(2)   # semi-maj. ax. in sq. along colums
    b_on = np.zeros(2)    # semi-minor axis in onion region
    b_row = np.zeros(2)   # semi-min. ax. in sq. along row
    b_col = np.zeros(2)   # semi-min. ax. in sq. along col
    slope_on = np.zeros(2)     # slope of straight lines in onion region
    slope_row = np.zeros(2)    # slope on east and west faces
    slope_col = np.zeros(2)    # slope on south and north faces
    y_interc = np.zeros(2)  # y interception of straight lines 
    # in onion region
    
    # Intersection between straight lines and ellipses
    # at the interface between sq. and onion region.
    # Row: east and west faces
    # Col: south and north faces
    x_inters_row = np.zeros(2)     
    x_inters_col = np.zeros(2)      
    y_inters_col = np.zeros(2)
    y_inters_row = np.zeros(2)

    # pts along wall of the pipe
    pt_wall_x = np.zeros(2)
    pt_wall_y = np.zeros(2)

    a_test_on = np.zeros(2)
 
    # Set distribution of elements in radial direction along axis in a certain
    # way for square region
    dr_sq_nominal = dr*stretch_sq  # stretch el in square region
    dr_sq = sq_dist(nSq, dr_sq_nominal, dr_sq_ratio)
    # Set distribution of elements in radial direction along axis in a certain 
    # way for onion region
    dr_on_nominal = (dr*nR - dr_sq_nominal*nSq)/(nR-nSq)
    dr_on = on_dist(nR, nSq, dr_on_nominal, dr_sq, distri_on)

    for el in elements:
        # Check in which cross section we are
        nel_cross_section = (nSq**2+(nR-nSq)*nSq*2)*4
        cross_sec = int(el.number/nel_cross_section)
        # Reduce current element number to be within the first cross section
        n = el.number - cross_sec*nel_cross_section

        # Define streamwise location 
        z0 = dz*cross_sec
        z1 = z0
        z2 = z0
        z3 = z0
        z4 = dz*(cross_sec+1)
        z5 = z4
        z6 = z4
        z7 = z4



        if (n<= nSq**2):   
            # This is the inner, "square" section
            #--------------------------------------------------
            # OPEN square
            #--------------------------------------------------
            i = (n-1)%nSq     # column number
            j = int((n-1)/nSq)    # row number

            # Determine the semi-major and semi-minor axis for interface between "square" 
            # and onion section
            #----------------------------------------------------------------------
            a_interface = a_interf
            b_interface = np.sum(dr_sq[:])

            # Idea: define ellipses in the inner region such that they coincide with the 
            # elements in the onion region.
            # OPEN def semi-major and semi-minor axis
            #--------------------------------------------------
            slope_row[0] = m.tan(m.pi/2*(j/ntheta))  # slope of the straight line on the bottom side
            slope_row[1] = m.tan(m.pi/2*((j+1)/ntheta)) # slope of the straight line on the top side
            slope_col[0] = m.tan(m.pi/2*((ntheta-i)/ntheta))  # slope of the straight line on the left side
            slope_col[1] = m.tan(m.pi/2*((ntheta-i-1)/ntheta)) # slope of the straight line on the right side

            b_row[0] = np.sum(dr_sq[:j])    # small semi-minor axis  
            b_row[1] = np.sum(dr_sq[:j+1])   # large semi-minor axis
            a_col[0] = np.sum(dr_sq[:i])    # small semi-major axis  
            a_col[1] = np.sum(dr_sq[:i+1])   # large semi-major axis

            # Set distribution of points along intersection
            # Shrink elements when approaching northeast corner of square region
            dr_sq_inters_ratio = dr_sq_int_ratio
            x_inters_corner = my_math.intersec_ellip_ellip(a_interface, b_interface, r_const**2,\
                    b_interface, a_interface, r_const**2)
            dr_sq_inters_nominal = x_inters_corner/nSq
            # distribution of elements along northern interface
            x_inters_col_dist = sq_dist(nSq, dr_sq_inters_nominal, dr_sq_inters_ratio)
            x_inters_col[0] = np.sum(x_inters_col_dist[:i])
            x_inters_col[1] = np.sum(x_inters_col_dist[:i+1])
            y_inters_col[0] = my_math.ellipse(a_interface, b_interface, r_const**2, x_inters_col[0])
            y_inters_col[1] = my_math.ellipse(a_interface, b_interface, r_const**2, x_inters_col[1])

            # distribution of elements along eastern row
            y_inters_row_dist = x_inters_col_dist
            y_inters_row[0] = np.sum(y_inters_row_dist[:j])
            y_inters_row[1] = np.sum(y_inters_row_dist[:j+1])
            x_inters_row[0] = my_math.ellipse(a_interface, b_interface, r_const**2, y_inters_row[0])
            x_inters_row[1] = my_math.ellipse(a_interface, b_interface, r_const**2, y_inters_row[1])

            # Find semi-major axis by the points at the intersection, r_const and semi-minor axis, 
            # which is defined by the vertical position 
            if (j==0):
                a_row[0] = 0   # this is reset later
                a_row[1] = x_inters_row[1]*b_row[1]/( (r_const**2*b_row[1]**2 - y_inters_row[1]**2)**0.5 )
            else: 
                a_row[0] = x_inters_row[0]*b_row[0]/( (r_const**2*b_row[0]**2 - y_inters_row[0]**2)**0.5 )
                a_row[1] = x_inters_row[1]*b_row[1]/( (r_const**2*b_row[1]**2 - y_inters_row[1]**2)**0.5 )

            # Find semi-minor axis by the points at the intersection, r_const and semi-major axis, 
            # which is defined by the horizontal position 
            # note that x and y need to be switched here
            if (i==0): 
                b_col[0] = 0   # this is reset later 
                b_col[1] = y_inters_col[1]*a_col[1]/( (r_const**2*a_col[1]**2 - x_inters_col[1]**2)**0.5 )
            else:
                b_col[0] = y_inters_col[0]*a_col[0]/( (r_const**2*a_col[0]**2 - x_inters_col[0]**2)**0.5 )
                b_col[1] = y_inters_col[1]*a_col[1]/( (r_const**2*a_col[1]**2 - x_inters_col[1]**2)**0.5 )

            # CLOSE def semi-major and semi-minor axis
            #--------------------------------------------------



            # Set vertex position and curvature
            if (j==0): # first row
                if (i==0):  # first col
                    x0 = np.sum(dr_sq[:i])
                    x1 = np.sum(dr_sq[:i+1])
                    x2 = my_math.intersec_ellip_ellip(a_row[1],b_row[1],r_const**2,\
                            a_col[1],b_col[1],r_const**2)
                    x3 = np.sum(dr_sq[:i])

                    y0 = np.sum(dr_sq[:j])
                    y1 = np.sum(dr_sq[:j])
                    y2 = my_math.ellipse(a_row[1],b_row[1],r_const**2,x2)
                    y3 = np.sum(dr_sq[:j+1])
                    
                    
                    # use the midpoint between two vertices to calculate the curvature 
                    c0 = 0
                    # note that for c1 and c2 we use the midpoint of y-coordinates and
                    # switch semi-major and semi-minor axis so that symmetry for
                    # the curvature is preserved, i.e. same curvature at interface 
                    # sq-onion in upper part and sq-onion lower part.
                    c1 = my_math.get_rad_ell(b_col[1],a_col[1],r_const**2, (y1+y2)/2)
                    c2 = my_math.get_rad_ell(a_row[1],b_row[1],r_const**2, (x2+x3)/2)
                    c3 = 0

                    # Account for 3D
                    x4 = x0
                    x5 = x1
                    x6 = x2
                    x7 = x3
                    y4 = y0
                    y5 = y1
                    y6 = y2
                    y7 = y3
                    c4 = c0
                    c5 = c1
                    c6 = c2
                    c7 = c3

                    el.x = np.array([x0, x1, x2, x3, x4, x5, x6, x7])
                    el.y = np.array([y0, y1, y2, y3, y4, y5, y6, y7])
                    el.z = np.array([z0, z1, z2, z3, z4, z5, z6, z7])
                    # Set concave sides to negative and convex positive
                    el.c = np.array([c0, c1, c2, c3, c4, c5, c6, c7])
                else:
                    x0 = np.sum(dr_sq[:i])
                    x1 = np.sum(dr_sq[:i+1])
                    x2 = my_math.intersec_ellip_ellip(a_row[1],b_row[1],r_const**2,\
                            a_col[1],b_col[1],r_const**2)
                    x3 = my_math.intersec_ellip_ellip(a_row[1],b_row[1],r_const**2,\
                            a_col[0],b_col[0],r_const**2)
                    
                    y0 = np.sum(dr_sq[:j])
                    y1 = np.sum(dr_sq[:j])
                    y2 = my_math.ellipse(a_row[1],b_row[1],r_const**2,x2)
                    y3 = my_math.ellipse(a_row[1],b_row[1],r_const**2,x3)

                    c0 = 0
                    c1 = my_math.get_rad_ell(b_col[1],a_col[1],r_const**2, (y1+y2)/2)
                    c2 = my_math.get_rad_ell(a_row[1],b_row[1],r_const**2, (x2+x3)/2)
                    c3 = my_math.get_rad_ell(b_col[0],a_col[0],r_const**2, (y0+y3)/2)

                    # Account for 3D
                    x4 = x0
                    x5 = x1
                    x6 = x2
                    x7 = x3
                    y4 = y0
                    y5 = y1
                    y6 = y2
                    y7 = y3
                    c4 = c0
                    c5 = c1
                    c6 = c2
                    c7 = c3

                    el.x = np.array([x0, x1, x2, x3, x4, x5, x6, x7])
                    el.y = np.array([y0, y1, y2, y3, y4, y5, y6, y7])
                    el.z = np.array([z0, z1, z2, z3, z4, z5, z6, z7])
                    # Set concave sides to negative and convex positive
                    el.c = np.array([c0, c1, c2, -c3, c4, c5, c6, -c7])

            elif (j>0 and i==0): # first col
                x0 = np.sum(dr_sq[:i])
                x1 = my_math.intersec_ellip_ellip(a_row[0],b_row[0],r_const**2,\
                        a_col[1],b_col[1],r_const**2)
                x2 = my_math.intersec_ellip_ellip(a_row[1],b_row[1],r_const**2,\
                        a_col[1],b_col[1],r_const**2)
                x3 = np.sum(dr_sq[:i])

                y0 = np.sum(dr_sq[:j])    
                y1 = my_math.ellipse(a_row[0],b_row[0],r_const**2,x1)
                y2 = my_math.ellipse(a_row[1],b_row[1],r_const**2,x2)
                y3 = np.sum(dr_sq[:j+1])

                c0 = my_math.get_rad_ell(a_row[0],b_row[0],r_const**2, (x0+x1)/2)
                c1 = my_math.get_rad_ell(b_col[1],a_col[1],r_const**2, (y1+y2)/2)
                c2 = my_math.get_rad_ell(a_row[1],b_row[1],r_const**2, (x2+x3)/2)
                c3 = 0

                # Account for 3D
                x4 = x0
                x5 = x1
                x6 = x2
                x7 = x3
                y4 = y0
                y5 = y1
                y6 = y2
                y7 = y3
                c4 = c0
                c5 = c1
                c6 = c2
                c7 = c3

                el.x = np.array([x0, x1, x2, x3, x4, x5, x6, x7])
                el.y = np.array([y0, y1, y2, y3, y4, y5, y6, y7])
                el.z = np.array([z0, z1, z2, z3, z4, z5, z6, z7])
                # Set concave sides to negative and convex positive
                el.c = np.array([-c0, c1, c2, c3, -c4, c5, c6, c7])
            elif (i> 0 and j>0):    # inside
                #find intersection between both ellipses
                x0 = my_math.intersec_ellip_ellip(a_row[0],b_row[0],r_const**2,\
                        a_col[0],b_col[0],r_const**2)
                x1 = my_math.intersec_ellip_ellip(a_row[0],b_row[0],r_const**2,\
                        a_col[1],b_col[1],r_const**2)
                x2 = my_math.intersec_ellip_ellip(a_row[1],b_row[1],r_const**2,\
                        a_col[1],b_col[1],r_const**2)
                x3 = my_math.intersec_ellip_ellip(a_row[1],b_row[1],r_const**2,\
                        a_col[0],b_col[0],r_const**2)

                y0 = my_math.ellipse(a_row[0],b_row[0],r_const**2,x0)
                y1 = my_math.ellipse(a_row[0],b_row[0],r_const**2,x1)
                y2 = my_math.ellipse(a_row[1],b_row[1],r_const**2,x2)
                y3 = my_math.ellipse(a_row[1],b_row[1],r_const**2,x3)

                c0 = my_math.get_rad_ell(a_row[0],b_row[0],r_const**2, (x0+x1)/2)
                c1 = my_math.get_rad_ell(b_col[1],a_col[1],r_const**2, (y1+y2)/2)
                c2 = my_math.get_rad_ell(a_row[1],b_row[1],r_const**2, (x2+x3)/2)
                c3 = my_math.get_rad_ell(b_col[0],a_col[0],r_const**2, (y0+y3)/2)

                # Account for 3D
                x4 = x0
                x5 = x1
                x6 = x2
                x7 = x3
                y4 = y0
                y5 = y1
                y6 = y2
                y7 = y3
                c4 = c0
                c5 = c1
                c6 = c2
                c7 = c3

                el.x = np.array([x0, x1, x2, x3, x4, x5, x6, x7])
                el.y = np.array([y0, y1, y2, y3, y4, y5, y6, y7])
                el.z = np.array([z0, z1, z2, z3, z4, z5, z6, z7])
                # Set concave sides to negative and convex positive
                el.c = np.array([-c0, c1, c2, -c3, -c4, c5, c6, -c7])
            else:
                sys.exit(1)
            #--------------------------------------------------
            # END square
            #--------------------------------------------------
        else:                       
            # This is the outer, "onion" section
            #--------------------------------------------------
            # OPEN onion
            #--------------------------------------------------
            i = ((n-1)-nSq**2)%(nSq*2)      # position in clockwise manner through each layer
            k = abs(i-((nSq*2)-1))                  # position in anticlockwise manner
            j = int(((n-1)-nSq**2)/(nSq*2)) # onion like layer number, inner one is first,
            # starting from j=0

            # Determine the semi-major and semi-minor axis for "onion" section
            #----------------------------------------------------------------------
            a_wall = 0.5    # semi-major axis at last layer (wall)
            
            # Semi minor-axis:
            b_on[0] = np.sum(dr_sq)+np.sum(dr_on[:j])
            b_on[1] = np.sum(dr_sq)+np.sum(dr_on[:j+1])

            # Difference between prescribed semi-major axis at interface square-onion
            # and semi-minor axis at the interface
            a_diff = a_interface - np.sum(dr_sq)
       
            # Toggle for outermost layer having constant radial size
#            tog_r_out_const = 0
            if (tog_r_out_const == 1):
                # Set semi-major axis decreasing from the interface value to the last value
                # of semi-minor axis before the wall and finally to a_wall.
                # This ensures that the outermost onion layer has a constant radial size.
                b_last = np.sum(dr_sq)+np.sum(dr_on[:-1])
                if (j < nR-nSq-1):   # we are in the inner layers
                    a_on[0] = b_on[0] + a_diff*my_math.sin_dist(nR-nSq, 1, 0, j)
                    a_on[1] = b_on[1] + a_diff*my_math.sin_dist(nR-nSq, 1, 0, j+1)
                else:                   # we are in the outermost layers
                    a_on[0] = b_last
                    a_on[1] = a_wall
            else:
#                tog_a_on_dist = 1
                if (tog_a_on_dist == 0):
                    a_on[0] = b_on[0] + a_diff*my_math.exp_dist(nR-nSq+1, 1, 0, j)
                    a_on[1] = b_on[1] + a_diff*my_math.exp_dist(nR-nSq+1, 1, 0, j+1)
                else:
                    a_on[0] = b_on[0] + a_diff*my_math.sin_dist(nR-nSq+1, 1, 0, j)
                    a_on[1] = b_on[1] + a_diff*my_math.sin_dist(nR-nSq+1, 1, 0, j+1)


            # Straight line defined by points on intersection square-onion and equidistantly
            # spaced points along circumference
            # OPEN
            #------------------------------
            pt_wall_x[0] = nR*dr*m.cos(m.pi/2*(k/ntheta))
            pt_wall_y[0] = nR*dr*m.sin(m.pi/2*(k/ntheta))
            pt_wall_x[1] = nR*dr*m.cos(m.pi/2*((k+1)/ntheta))
            pt_wall_y[1] = nR*dr*m.sin(m.pi/2*((k+1)/ntheta))


            # Get points along intersection for each element in onion region 
            # Note that i and k are different than in square region
            dr_sq_inters_ratio = dr_sq_int_ratio
            x_inters_corner = my_math.intersec_ellip_ellip(a_interface, b_interface, r_const**2,\
                    b_interface, a_interface, r_const**2)
            dr_sq_inters_nominal = x_inters_corner/nSq
            x_inters_col_dist = sq_dist(nSq, dr_sq_inters_nominal, dr_sq_inters_ratio)
            x_inters_col[0] = np.sum(x_inters_col_dist[:i])
            x_inters_col[1] = np.sum(x_inters_col_dist[:i+1])
            y_inters_col[0] = my_math.ellipse(a_interface, b_interface, r_const**2, x_inters_col[0])
            y_inters_col[1] = my_math.ellipse(a_interface, b_interface, r_const**2, x_inters_col[1])

            y_inters_row_dist = x_inters_col_dist
            y_inters_row[0] = np.sum(y_inters_row_dist[:k])
            y_inters_row[1] = np.sum(y_inters_row_dist[:k+1])
            x_inters_row[0] = my_math.ellipse(a_interface, b_interface, r_const**2, y_inters_row[0])
            x_inters_row[1] = my_math.ellipse(a_interface, b_interface, r_const**2, y_inters_row[1])

            if (i<nSq):    # upper onion part
                # slope of the straight line on the right side of the element
                slope_on[0], y_interc[0] = my_math.get_line_params(pt_wall_x[0], pt_wall_y[0],\
                        x_inters_col[1], y_inters_col[1])
                # slope of the straight line on the left side of the element
                slope_on[1], y_interc[1] = my_math.get_line_params(pt_wall_x[1], pt_wall_y[1],\
                        x_inters_col[0], y_inters_col[0])
            else:       # lower onion part
                # slope of the straight line on the bottom side of the element
                slope_on[0], y_interc[0] = my_math.get_line_params(pt_wall_x[0], pt_wall_y[0],\
                        x_inters_row[0], y_inters_row[0])
                # slope of the straight line on the top side of the element
                slope_on[1], y_interc[1] = my_math.get_line_params(pt_wall_x[1], pt_wall_y[1],\
                        x_inters_row[1], y_inters_row[1])

            #------------------------------
            # END 
            # Definition of straight lines


            # Set vertex position in onion region as well as curvature
            if (i <= (nSq-1)):  # upper part, including border /
                x0 = my_math.intersec_ellip_line(a_on[0],b_on[0],r_const**2,slope_on[1],y_interc[1])
                x1 = my_math.intersec_ellip_line(a_on[0],b_on[0],r_const**2,slope_on[0],y_interc[0])
                x2 = my_math.intersec_ellip_line(a_on[1],b_on[1],r_const**2,slope_on[0],y_interc[0])
                x3 = my_math.intersec_ellip_line(a_on[1],b_on[1],r_const**2,slope_on[1],y_interc[1])

                y0 = my_math.ellipse(a_on[0],b_on[0],r_const**2,x0)
                y1 = my_math.ellipse(a_on[0],b_on[0],r_const**2,x1)
                y2 = my_math.ellipse(a_on[1],b_on[1],r_const**2,x2)
                y3 = my_math.ellipse(a_on[1],b_on[1],r_const**2,x3)

                c0 = my_math.get_rad_ell(a_on[0],b_on[0],r_const**2, (x0+x1)/2)
                c1 = 0
                c2 = my_math.get_rad_ell(a_on[1],b_on[1],r_const**2, (x2+x3)/2)
                c3 = 0

                # Account for 3D
                x4 = x0
                x5 = x1
                x6 = x2
                x7 = x3
                y4 = y0
                y5 = y1
                y6 = y2
                y7 = y3
                c4 = c0
                c5 = c1
                c6 = c2
                c7 = c3

                el.x = np.array([x0, x1, x2, x3, x4, x5, x6, x7])
                el.y = np.array([y0, y1, y2, y3, y4, y5, y6, y7])
                el.z = np.array([z0, z1, z2, z3, z4, z5, z6, z7])
                # Set concave sides to negative and convex positive
                el.c = np.array([-c0, c1, c2, c3, -c4, c5, c6, c7])
            elif (i >= nSq):     # lower part, including border /
                # note that semi-major and semi-minor axis are switched
                x0 = my_math.intersec_ellip_line(b_on[0],a_on[0],r_const**2,slope_on[0],y_interc[0])
                x1 = my_math.intersec_ellip_line(b_on[1],a_on[1],r_const**2,slope_on[0],y_interc[0])
                x2 = my_math.intersec_ellip_line(b_on[1],a_on[1],r_const**2,slope_on[1],y_interc[1])
                x3 = my_math.intersec_ellip_line(b_on[0],a_on[0],r_const**2,slope_on[1],y_interc[1])

                y0 = my_math.line(slope_on[0],x0,y_interc[0])
                y1 = my_math.line(slope_on[0],x1,y_interc[0])
                y2 = my_math.line(slope_on[1],x2,y_interc[1])
                y3 = my_math.line(slope_on[1],x3,y_interc[1])
                
                c0 = 0
                # Note that we use y-coordinates midpoint again as above. 
                # We do not need to switch semi-major and semi-minor axis.
                c1 = my_math.get_rad_ell(a_on[1],b_on[1],r_const**2, (y1+y2)/2)
                c2 = 0
                c3 = my_math.get_rad_ell(a_on[0],b_on[0],r_const**2, (y0+y3)/2)

                # Account for 3D
                x4 = x0
                x5 = x1
                x6 = x2
                x7 = x3
                y4 = y0
                y5 = y1
                y6 = y2
                y7 = y3
                c4 = c0
                c5 = c1
                c6 = c2
                c7 = c3

                el.x = np.array([x0, x1, x2, x3, x4, x5, x6, x7])
                el.y = np.array([y0, y1, y2, y3, y4, y5, y6, y7])
                el.z = np.array([z0, z1, z2, z3, z4, z5, z6, z7])
                # Set concave sides to negative and convex positive
                el.c = np.array([c0, c1, c2, -c3, c4, c5, c6, -c7])


def compl_mesh(elements, nR, nSq):
    """ Complete the quarter mesh to a whole cross section """

    # Mirror the first quarter along y-axis to create mesh in 2nd quadrant
    # then for 3rd and finally 4th quadrant
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
        mirr_el.x = np.array([mirr_el.x[1], mirr_el.x[0], mirr_el.x[3], mirr_el.x[2],\
                mirr_el.x[5], mirr_el.x[4], mirr_el.x[7], mirr_el.x[6]])
        mirr_el.y = el.y
        mirr_el.y = np.array([mirr_el.y[1], mirr_el.y[0], mirr_el.y[3], mirr_el.y[2],\
                mirr_el.y[5], mirr_el.y[4], mirr_el.y[7], mirr_el.y[6]])
        mirr_el.z = el.z

        
        # Position
        mirr_el.pos = el.pos

        # Boundary condition
        # (This is reset later)
        mirr_el.fl_bc = ['E  ','E  ','E  ','E  ']

        mirr_el.bc_con_f = np.zeros(4)  # connected face: 1: south, 2:east, 3:north, 4:west
        mirr_el.bc_con_el = np.zeros(4)  # number of the connected element

        # Curvature
        mirr_el.c = np.array([el.c[0],el.c[3],el.c[2],el.c[1],\
                el.c[4],el.c[7],el.c[6],el.c[5]])

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
        mirr_el.x = np.array([mirr_el.x[2], mirr_el.x[3], mirr_el.x[0], mirr_el.x[1],\
                mirr_el.x[6], mirr_el.x[7], mirr_el.x[4], mirr_el.x[5]])
        mirr_el.y = el.y*(-1)
        mirr_el.y = np.array([mirr_el.y[2], mirr_el.y[3], mirr_el.y[0], mirr_el.y[1],\
                mirr_el.x[6], mirr_el.y[7], mirr_el.y[4], mirr_el.y[5]])
        mirr_el.z = el.z
        
        # Position
        mirr_el.pos = el.pos

        # Boundary condition
        # (This is reset later)
        mirr_el.fl_bc = ['E  ','E  ','E  ','E  ']

        mirr_el.bc_con_f = np.zeros(4)  # connected face: 1: south, 2:east, 3:north, 4:west
        mirr_el.bc_con_el = np.zeros(4)  # number of the connected element

         # Curvature
        mirr_el.c = np.array([el.c[2], el.c[3], el.c[0], el.c[1],\
                el.c[6], el.c[7], el.c[4], el.c[5]])
       
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
        mirr_el.x = np.array([mirr_el.x[3], mirr_el.x[2], mirr_el.x[1], mirr_el.x[0],\
                mirr_el.x[7], mirr_el.x[6], mirr_el.x[5], mirr_el.x[4]])
        mirr_el.y = el.y*(-1)
        mirr_el.y = np.array([mirr_el.y[3], mirr_el.y[2], mirr_el.y[1], mirr_el.y[0],\
                mirr_el.y[7], mirr_el.y[6], mirr_el.y[5], mirr_el.y[4]])
        mirr_el.z = el.z
        
        # Position
        mirr_el.pos = el.pos

        # Boundary condition
        # (This is reset later)
        mirr_el.fl_bc = ['E  ','E  ','E  ','E  ']

        mirr_el.bc_con_f = np.zeros(4)  # connected face: 1: south, 2:east, 3:north, 4:west
        mirr_el.bc_con_el = np.zeros(4)  # number of the connected element

         # Curvature
        mirr_el.c = np.array([el.c[2], el.c[1], el.c[0], el.c[3],\
                el.c[6], el.c[5], el.c[4], el.c[7]])
       
        # Add mirrored element to the list of elements
        el_list_4th.append(mirr_el)

    # Add all new elements to the element list
    elements.extend(el_list_2nd)
    elements.extend(el_list_3rd)
    elements.extend(el_list_4th)



# Set the boundary conditions for each quadrant separately by checking the element's
# position and writing the corresponding BCs into the element's attributes.
# (This might be a possible cause of bugs since each BC is written manually
# depending on its location. Even though I was very careful, this approach is prone to errors.)
def set_bc_q1(elements,nR,nSq):
    """ Set boundary conditions for each face. 
    
    This is for quadrant 1.
    """

    nel_quarter = nSq**2 + (nR-nSq)*2*nSq       # number of elements in one quarter

    for el in elements:
        # Check in which cross section we are
        nel_cross_section = (nSq**2+(nR-nSq)*nSq*2)*4
        cross_sec = int(el.number/nel_cross_section)
        # Reduce current element number to be within the first cross section
        n = el.number - cross_sec*nel_cross_section


        # only consider elements in the first quadrant
        if (n<=nel_quarter):
#        elements = elements[0:nel_quarter]

#    for el in elements:
#            n = el.number
            check_position(el,nR,nSq)
            position = el.pos
    
            if (n <= nSq**2):   # we are in the square section
                i = (n-1)%nSq     # column number
                j = int((n-1)/nSq)    # row number
    
                if (position == 'sq_low_left'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nel_quarter*3, n+1, n+nSq, n+nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_low_right'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nel_quarter*3, nSq**2+(2*nSq-j), n+nSq, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_up_left'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n+1, n+nSq, n+nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_up_right'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, nSq**2+(2*nSq-j), n+nSq, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_north_row'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n+1, n+nSq, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_east_col'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, nSq**2+(2*nSq-j), n+nSq, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_south_row'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nel_quarter*3, n+1, n+nSq, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_west_col'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n+1, n+nSq, n+nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_internal'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
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
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-nSq, n+1, n+2*nSq, n+nel_quarter])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_up_south_sq_east_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-nSq, n+1, n+2*nSq, n-1])
                        el.bc_con_f = np.array([3, 3, 1, 2])
                    elif (position == 'on_up_east_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-2*nSq, n+1, n+2*nSq, n-1])
                        el.bc_con_f = np.array([3, 3, 1, 2])
                    elif (position == 'on_up_south_sq'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-nSq, n+1, n+2*nSq, n-1])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_up_west_y_north_wall'):
                        el.fl_bc = ['E  ','E  ','W  ','E  ']
                        el.th_bc = ['E  ','E  ','f  ','E  ']
                        el.bc_con_el = np.array([n-2*nSq, n+1, 0, n+nel_quarter])
                        el.bc_con_f = np.array([3, 4, 0, 2])
                    elif (position == 'on_up_west_y'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-2*nSq, n+1, n+2*nSq, n+nel_quarter])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_up_north_wall_east_edge'):
                        el.fl_bc = ['E  ','E  ','W  ','E  ']
                        el.th_bc = ['E  ','E  ','f  ','E  ']
                        el.bc_con_el = np.array([n-2*nSq, n+1, 0, n-1])
                        el.bc_con_f = np.array([3, 3, 0, 2])
                    elif (position == 'on_up_north_wall'):
                        el.fl_bc = ['E  ','E  ','W  ','E  ']
                        el.th_bc = ['E  ','E  ','f  ','E  ']                   
                        el.bc_con_el = np.array([n-2*nSq, n+1, 0, n-1])
                        el.bc_con_f = np.array([3, 4, 0, 2])
                    elif (position == 'on_up_intern'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  '] 
                        el.bc_con_el = np.array([n-2*nSq, n+1, n+2*nSq, n-1])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                elif ('on_low' in position):    # we are in the lower onion part
                    if (position == 'on_low_west_sq_south_x'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+nel_quarter*3, n+nSq*2, n-1, (k+1)*nSq])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_low_west_sq_north_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+1, n+nSq*2, n-1, (k+1)*nSq])
                        el.bc_con_f = np.array([3, 4, 2, 2])
                    elif (position == 'on_low_west_sq'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+1, n+nSq*2, n-1, (k+1)*nSq])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_low_south_x_east_wall'):
                        el.fl_bc = ['E  ','W  ','E  ','E  ']
                        el.th_bc = ['E  ','f  ','E  ','E  ']
                        el.bc_con_el = np.array([n+nel_quarter*3, 0, n-1, n-nSq*2])
                        el.bc_con_f = np.array([3, 0, 1, 2])
                    elif (position == 'on_low_south_x'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+nel_quarter*3, n+nSq*2, n-1, n-nSq*2])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_low_east_wall_north_edge'):
                        el.fl_bc = ['E  ','W  ','E  ','E  ']
                        el.th_bc = ['E  ','f  ','E  ','E  ']
                        el.bc_con_el = np.array([n+1, 0, n-1, n-nSq*2])
                        el.bc_con_f = np.array([3, 0, 2, 2])
                    elif (position == 'on_low_east_wall'):
                        el.fl_bc = ['E  ','W  ','E  ','E  ']
                        el.th_bc = ['E  ','f  ','E  ','E  ']
                        el.bc_con_el = np.array([n+1, 0, n-1, n-nSq*2])
                        el.bc_con_f = np.array([3, 0, 1, 2])
                    elif (position == 'on_low_north_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+1, n+nSq*2, n-1, n-nSq*2])
                        el.bc_con_f = np.array([3, 4, 2, 2])
                    elif (position == 'on_low_intern'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
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


    for el in elements:
        # Check in which cross section we are
        nel_cross_section = (nSq**2+(nR-nSq)*nSq*2)*4
        cross_sec = int(el.number/nel_cross_section)
        # Reduce current element number to be within the first cross section
        n = el.number - cross_sec*nel_cross_section


    # only consider elements in the second quadrant
        if (n>nel_quarter and n<=2*nel_quarter):
#        elements = elements[0:nel_quarter]

    # only consider elements in the second quadrant
#    elements = elements[nel_quarter:(nel_quarter*2)]

#    for el in elements:
#            n = el.number
#            check_position(el,nR,nSq)
            position = el.pos
    
            if (n-nel_quarter <= nSq**2):   # we are in the square section
                i = (n-nel_quarter-1)%nSq     # column number
                j = int((n-nel_quarter-1)/nSq)    # row number
    
                if (position == 'sq_low_left'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nel_quarter, n-nel_quarter, n+nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_low_right'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nel_quarter, n-1, n+nSq, nSq**2+(2*nSq-j)+nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_up_left'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n-nel_quarter, n+nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_up_right'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n-1, n+nSq, nSq**2+(2*nSq-j)+nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_north_row'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n-1, n+nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_east_col'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n-1, n+nSq, nSq**2+(2*nSq-j)+nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_south_row'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nel_quarter, n-1, n+nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_west_col'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n-nel_quarter, n+nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_internal'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n-nSq, n-1, n+nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                else:
                    print('position of element not found!')
                    # Debug point
                    print(n, position)
                    print(el.number)
                    print(el.x)
                    print(el.y)
                    print(el.z)


                    sys.exit(1)
            else:                       # we are in the outer onion like region :-)
                i = ((n-nel_quarter-1)-nSq**2)%(nSq*2) # position in clockwise manner through each layer
                k = abs(i-((nSq*2)-1))                  # position in anticlockwise manner
                j = int(((n-nel_quarter-1)-nSq**2)/(nSq*2)) # onion like layer number, inner one is first
    
                if ('on_up' in position):   # we are in the upper onion part
                    if (position == 'on_up_south_sq_west_y'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-nSq, n-nel_quarter, n+2*nSq, n+1])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_up_south_sq_east_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-nSq, n-1, n+2*nSq, n+1])
                        el.bc_con_f = np.array([3, 4, 1, 3])
                    elif (position == 'on_up_east_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-2*nSq, n-1, n+2*nSq, n+1])
                        el.bc_con_f = np.array([3, 4, 1, 3])
                    elif (position == 'on_up_south_sq'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-nSq, n-1, n+2*nSq, n+1])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_up_west_y_north_wall'):
                        el.fl_bc = ['E  ','E  ','W  ','E  ']
                        el.th_bc = ['E  ','E  ','f  ','E  ']
                        el.bc_con_el = np.array([n-2*nSq, n-nel_quarter, 0, n+1])
                        el.bc_con_f = np.array([3, 4, 0, 2])
                    elif (position == 'on_up_west_y'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-2*nSq, n-nel_quarter, n+2*nSq, n+1])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_up_north_wall_east_edge'):
                        el.fl_bc = ['E  ','E  ','W  ','E  ']
                        el.th_bc = ['E  ','E  ','f  ','E  ']
                        el.bc_con_el = np.array([n-2*nSq, n-1, 0, n+1])
                        el.bc_con_f = np.array([3, 4, 0, 3])
                    elif (position == 'on_up_north_wall'):
                        el.fl_bc = ['E  ','E  ','W  ','E  ']
                        el.th_bc = ['E  ','E  ','f  ','E  ']
                        el.bc_con_el = np.array([n-2*nSq, n-1, 0, n+1])
                        el.bc_con_f = np.array([3, 4, 0, 2])
                    elif (position == 'on_up_intern'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-2*nSq, n-1, n+2*nSq, n+1])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                elif ('on_low' in position):    # we are in the lower onion part
                    if (position == 'on_low_west_sq_south_x'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+nel_quarter, (k+1)*nSq+nel_quarter, n-1, n+nSq*2])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_low_west_sq_north_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+1, (k+1)*nSq+nel_quarter, n-1, n+nSq*2])
                        el.bc_con_f = np.array([3, 4, 4, 2])
                    elif (position == 'on_low_west_sq'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+1, (k+1)*nSq+nel_quarter, n-1, n+nSq*2])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_low_south_x_east_wall'):
                        el.fl_bc = ['E  ','E  ','E  ','W  ']
                        el.th_bc = ['E  ','E  ','E  ','f  ']
                        el.bc_con_el = np.array([n+nel_quarter, n-nSq*2, n-1, 0])
                        el.bc_con_f = np.array([3, 4, 1, 0])
                    elif (position == 'on_low_south_x'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+nel_quarter, n-nSq*2, n-1, n+nSq*2])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_low_east_wall_north_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','W  ']
                        el.th_bc = ['E  ','E  ','E  ','f  ']
                        el.bc_con_el = np.array([n+1, n-nSq*2, n-1, 0])
                        el.bc_con_f = np.array([3, 4, 4, 0])
                    elif (position == 'on_low_east_wall'):
                        el.fl_bc = ['E  ','E  ','E  ','W  ']
                        el.th_bc = ['E  ','E  ','E  ','f  ']
                        el.bc_con_el = np.array([n+1, n-nSq*2, n-1, 0])
                        el.bc_con_f = np.array([3, 4, 1, 0])
                    elif (position == 'on_low_north_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+1, n-nSq*2, n-1, n+nSq*2])
                        el.bc_con_f = np.array([3, 4, 4, 2])
                    elif (position == 'on_low_intern'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
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


    for el in elements:
        # Check in which cross section we are
        nel_cross_section = (nSq**2+(nR-nSq)*nSq*2)*4
        cross_sec = int(el.number/nel_cross_section)
        # Reduce current element number to be within the first cross section
        n = el.number - cross_sec*nel_cross_section


        # only consider elements in the third quadrant
        if (n>2*nel_quarter and n<=3*nel_quarter):
    
    
    #    # only consider elements in the third quadrant
    #    elements = elements[nel_quarter*2:(nel_quarter*3)]
    
    #    for el in elements:
#            n = el.number
            position = el.pos
    
            if (n-nel_quarter*2 <= nSq**2):   # we are in the square section
                i = (n-nel_quarter*2-1)%nSq     # column number
                j = int((n-nel_quarter*2-1)/nSq)    # row number
    
                if (position == 'sq_low_left'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n+nel_quarter, n-nel_quarter, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_low_right'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n-1, n-nel_quarter, nSq**2+(2*nSq-j)+nel_quarter*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_up_left'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n+nel_quarter, n-nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_up_right'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n-1, n-nSq, nSq**2+(2*nSq-j)+nel_quarter*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_north_row'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n-1, n-nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_east_col'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n-1, n-nSq, nSq**2+(2*nSq-j)+nel_quarter*2])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_south_row'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n-1, n-nel_quarter, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_west_col'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n+nel_quarter, n-nSq, n+1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_internal'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
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
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+2*nSq, n+nel_quarter, n-nSq, n+1])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_up_south_sq_east_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+2*nSq, n-1, n-nSq, n+1])
                        el.bc_con_f = np.array([3, 4, 1, 1])
                    elif (position == 'on_up_east_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+2*nSq, n-1, n-2*nSq, n+1])
                        el.bc_con_f = np.array([3, 4, 1, 1])
                    elif (position == 'on_up_south_sq'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+2*nSq, n-1, n-nSq, n+1])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_up_west_y_north_wall'):
                        el.fl_bc = ['W  ','E  ','E  ','E  ']
                        el.th_bc = ['f  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([0, n+nel_quarter, n-2*nSq, n+1])
                        el.bc_con_f = np.array([0, 4, 1, 2])
                    elif (position == 'on_up_west_y'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+2*nSq, n+nel_quarter, n-2*nSq, n+1])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_up_north_wall_east_edge'):
                        el.fl_bc = ['W  ','E  ','E  ','E  ']
                        el.th_bc = ['f  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([0, n-1, n-2*nSq, n+1])
                        el.bc_con_f = np.array([0, 4, 1, 1])
                    elif (position == 'on_up_north_wall'):
                        el.fl_bc = ['W  ','E  ','E  ','E  ']
                        el.th_bc = ['f  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([0, n-1, n-2*nSq, n+1])
                        el.bc_con_f = np.array([0, 4, 1, 2])
                    elif (position == 'on_up_intern'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+2*nSq, n-1, n-2*nSq, n+1])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                elif ('on_low' in position):    # we are in the lower onion part
                    if (position == 'on_low_west_sq_south_x'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, (k+1)*nSq+nel_quarter*2, n-nel_quarter, n+nSq*2])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_low_west_sq_north_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, (k+1)*nSq+nel_quarter*2, n+1, n+nSq*2])
                        el.bc_con_f = np.array([4, 4, 1, 2])
                    elif (position == 'on_low_west_sq'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, (k+1)*nSq+nel_quarter*2, n+1, n+nSq*2])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_low_south_x_east_wall'):
                        el.fl_bc = ['E  ','E  ','E  ','W  ']
                        el.th_bc = ['E  ','E  ','E  ','f  ']                   
                        el.bc_con_el = np.array([n-1, n-nSq*2, n-nel_quarter, 0])
                        el.bc_con_f = np.array([3, 4, 1, 0])
                    elif (position == 'on_low_south_x'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, n-nSq*2, n-nel_quarter, n+nSq*2])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_low_east_wall_north_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','W  ']
                        el.th_bc = ['E  ','E  ','E  ','f  ']
                        el.bc_con_el = np.array([n-1, n-nSq*2, n+1, 0])
                        el.bc_con_f = np.array([4, 4, 1, 0])
                    elif (position == 'on_low_east_wall'):
                        el.fl_bc = ['E  ','E  ','E  ','W  ']
                        el.th_bc = ['E  ','E  ','E  ','f  ']
                        el.bc_con_el = np.array([n-1, n-nSq*2, n+1, 0])
                        el.bc_con_f = np.array([3, 4, 1, 0])
                    elif (position == 'on_low_north_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, n-nSq*2, n+1, n+nSq*2])
                        el.bc_con_f = np.array([4, 4, 1, 2])
                    elif (position == 'on_low_intern'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
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


    for el in elements:
        # Check in which cross section we are
        nel_cross_section = (nSq**2+(nR-nSq)*nSq*2)*4
        cross_sec = int(el.number/nel_cross_section)
        # Reduce current element number to be within the first cross section
        n = el.number - cross_sec*nel_cross_section


        # only consider elements in the fourth quadrant
        if (n>3*nel_quarter and n<=4*nel_quarter):


    # only consider elements in the fourth quadrant
#    elements = elements[nel_quarter*3:(nel_quarter*4)]

#    for el in elements:
#            n = el.number
            position = el.pos
    
            if (n-nel_quarter*3 <= nSq**2):   # we are in the square section
                i = (n-nel_quarter*3-1)%nSq     # column number
                j = int((n-nel_quarter*3-1)/nSq)    # row number
    
                if (position == 'sq_low_left'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n+1, n-nel_quarter*3, n-nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_low_right'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, nSq**2+(2*nSq-j)+nel_quarter*3, n-nel_quarter*3, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_up_left'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n+1, n-nSq, n-nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_up_right'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, nSq**2+(2*nSq-j)+nel_quarter*3, n-nSq, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_north_row'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n+1, n-nSq, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_east_col'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, nSq**2+(2*nSq-j)+nel_quarter*3, n-nSq, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_south_row'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n+1, n-nel_quarter*3, n-1])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_west_col'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
                    el.bc_con_el = np.array([n+nSq, n+1, n-nSq, n-nel_quarter])
                    el.bc_con_f = np.array([3, 4, 1, 2])
                elif (position == 'sq_internal'):
                    el.fl_bc = ['E  ','E  ','E  ','E  ']
                    el.th_bc = ['E  ','E  ','E  ','E  ']
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
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+2*nSq, n+1, n-nSq, n-nel_quarter])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_up_south_sq_east_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+2*nSq, n+1, n-nSq, n-1])
                        el.bc_con_f = np.array([3, 1, 1, 4])
                    elif (position == 'on_up_east_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+2*nSq, n+1, n-2*nSq, n-1])
                        el.bc_con_f = np.array([3, 1, 1, 4])
                    elif (position == 'on_up_south_sq'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+2*nSq, n+1, n-nSq, n-1])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_up_west_y_north_wall'):
                        el.fl_bc = ['W  ','E  ','E  ','E  ']
                        el.th_bc = ['f  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([0, n+1, n-2*nSq, n-nel_quarter])
                        el.bc_con_f = np.array([0, 4, 1, 2])
                    elif (position == 'on_up_west_y'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+2*nSq, n+1, n-2*nSq, n-nel_quarter])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_up_north_wall_east_edge'):
                        el.fl_bc = ['W  ','E  ','E  ','E  ']
                        el.th_bc = ['f  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([0, n+1, n-2*nSq, n-1])
                        el.bc_con_f = np.array([0, 1, 1, 4])
                    elif (position == 'on_up_north_wall'):
                        el.fl_bc = ['W  ','E  ','E  ','E  ']
                        el.th_bc = ['f  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([0, n+1, n-2*nSq, n-1])
                        el.bc_con_f = np.array([0, 4, 1, 2])
                    elif (position == 'on_up_intern'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n+2*nSq, n+1, n-2*nSq, n-1])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                elif ('on_low' in position):    # we are in the lower onion part
                    if (position == 'on_low_west_sq_south_x'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, n+nSq*2, n-nel_quarter*3, (k+1)*nSq+nel_quarter*3])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_low_west_sq_north_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, n+nSq*2, n+1, (k+1)*nSq+nel_quarter*3])
                        el.bc_con_f = np.array([2, 4, 1, 4])
                    elif (position == 'on_low_west_sq'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, n+nSq*2, n+1, (k+1)*nSq+nel_quarter*3])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_low_south_x_east_wall'):
                        el.fl_bc = ['E  ','W  ','E  ','E  ']
                        el.th_bc = ['E  ','f  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, 0, n-nel_quarter*3, n-nSq*2])
                        el.bc_con_f = np.array([3, 0, 1, 4])
                    elif (position == 'on_low_south_x'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, n+nSq*2, n-nel_quarter*3, n-nSq*2])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                    elif (position == 'on_low_east_wall_north_edge'):
                        el.fl_bc = ['E  ','W  ','E  ','E  ']
                        el.th_bc = ['E  ','f  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, 0, n+1, n-nSq*2])
                        el.bc_con_f = np.array([2, 0, 1, 1])
                    elif (position == 'on_low_east_wall'):
                        el.fl_bc = ['E  ','W  ','E  ','E  ']
                        el.th_bc = ['E  ','f  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, 0, n+1, n-nSq*2])
                        el.bc_con_f = np.array([3, 0, 1, 4])
                    elif (position == 'on_low_north_edge'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, n+nSq*2, n+1, n-nSq*2])
                        el.bc_con_f = np.array([2, 4, 1, 2])
                    elif (position == 'on_low_intern'):
                        el.fl_bc = ['E  ','E  ','E  ','E  ']
                        el.th_bc = ['E  ','E  ','E  ','E  ']
                        el.bc_con_el = np.array([n-1, n+nSq*2, n+1, n-nSq*2])
                        el.bc_con_f = np.array([3, 4, 1, 2])
                else:
                    print('position assignment was not correct!')
                    sys.exit(3)

    
def check_position(element, nR, nSq):
    """ Check position of the given element within the quarter region
    and save it as an element attribute. 
    
    This is needed for choosing the right boundary conditions. 
    """

    el = element
    # Check in which cross section we are
    nel_cross_section = (nSq**2+(nR-nSq)*nSq*2)*4
    cross_sec = int(el.number/nel_cross_section)
    # Reduce current element number to be within the first cross section
    n = el.number - cross_sec*nel_cross_section


#    n = el.number
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
        z = el.z
        n = el.number

        mesh.append('{0:>19s} {1:10d} {2:6s}{3:1s}{4:12s}'.format\
            ('ELEMENT',n,'[    1','a',']  GROUP  0\n'))
        mesh.append('{0: 10.6E}{1: 14.6E}{2: 14.6E}{3: 14.6E}   {4:s}'.format\
            (x[0], x[1], x[2], x[3], '\n'))   # x coordinates
        mesh.append('{0: 10.6E}{1: 14.6E}{2: 14.6E}{3: 14.6E}   {4:s}'.format\
            (y[0], y[1], y[2], y[3], '\n'))  # y coordinates
        mesh.append('{0: 10.6E}{1: 14.6E}{2: 14.6E}{3: 14.6E}   {4:s}'.format\
            (z[0], z[1], z[2], z[3], '\n'))  # z coordinates
        mesh.append('{0: 10.6E}{1: 14.6E}{2: 14.6E}{3: 14.6E}   {4:s}'.format\
            (x[4], x[5], x[6], x[7], '\n'))   # x coordinates
        mesh.append('{0: 10.6E}{1: 14.6E}{2: 14.6E}{3: 14.6E}   {4:s}'.format\
            (y[4], y[5], y[6], y[7], '\n'))  # y coordinates
        mesh.append('{0: 10.6E}{1: 14.6E}{2: 14.6E}{3: 14.6E}   {4:s}'.format\
            (z[4], z[5], z[6], z[7], '\n'))  # z coordinates



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

def write_curv(elements):
    """ Write curvature information to rea file. 
    
    Note that only curved sides are allowed to be printed here.
    """

    # Count all the curved edges
    num_curv = 0
    for el in elements:
        for f in range(0,8):
            if (abs(el.c[f]) > 1e-15):
                num_curv = num_curv+1
    curv = []
    n_tot = len(elements)
    curv.append('{0:10d} Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n'.format(num_curv))

    # Check number of curved sides for correct formatting of curved side data 
    # (see Nek5000 user doc. p. 20)
    if (n_tot < 1e3):
        format_str = '{iedge:3d}{current_el:3d}{curve1:14.6f}\
{curve2:14.6f}{curve3:14.6f}{curve4:14.6f}{curve5:14.6f} {ccurve:s}\
{newline:s}'
    elif (n_tot < 1e6):
        format_str = '{iedge:2d}{current_el:6d}{curve1:14.6f}\
{curve2:14.6f}{curve3:14.6f}{curve4:14.6f}{curve5:14.6f} {ccurve:s}\
{newline:s}'
    else:
        format_str = '{iedge:2d}{current_el:12}{curve1:14.6f}\
{curve2:14.6f}{curve3:14.6f}{curve4:14.6f}{curve5:14.6f} {ccurve:s}\
{newline:s}'

    for el in elements:
        for f in range(8):
            if (abs(el.c[f]) > 1e-15):
                curv.append(format_str.format(iedge=f+1, current_el=el.number,\
                curve1=el.c[f],curve2=0.0,curve3=0.0,curve4=0.0,curve5=0.0,\
                ccurve='C',newline='\n'))

    f = open('base.rea','r')
    contents = f.readlines()
    f.close()
 
    # find row for curv data
    line_curv = 0
    while (not 'CURVED SIDE DATA' in contents[line_curv]):
        line_curv = line_curv + 1
 
    # and write it to rea file
    contents[line_curv+1:line_curv+1] = curv
    contents = "".join(contents)
    f = open('base.rea','w')
    f.write(contents)
    f.close()


def write_fl_bc(elements, nR, nSq):
    """ Write fluid boundary conditions to rea file. """

    bc = []
    dig_n_tot = len(str(elements[-1].number))   # size of element number
    for el in elements:
        for f in range(4):
            bc.append(' {boundary:3s}{current_el:{digits_n_tot}d} {face:2d}   \
{con_el:<07.1f}{con_f:14.5f}{zero1:14.5f}{zero2:14.5f}{zero3:14.5f}    {newline:s}'\
            .format(boundary=el.fl_bc[f], current_el=el.number, digits_n_tot=dig_n_tot, face=(f+1),\
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


def write_th_bc(elements, nR, nSq):
    """ Write fluid boundary conditions to rea file. """

    bc = []
    dig_n_tot = len(str(elements[-1].number))   # size of element number
    for el in elements:
        for f in range(4):
            bc.append(' {boundary:3s}{current_el:{digits_n_tot}d} {face:2d}   \
{con_el:<07.1f}{con_f:14.5f}{zero1:14.5f}{zero2:14.5f}{zero3:14.5f}    {newline:s}'\
            .format(boundary=el.th_bc[f], current_el=el.number, digits_n_tot=dig_n_tot, face=(f+1),\
            con_el=el.bc_con_el[f], con_f=el.bc_con_f[f],\
            zero1=0.0,zero2=0.0,zero3=0.0,newline='\n'))

    f = open('base.rea','r')
    contents = f.readlines()
    f.close()
 
    # find row for bc data
    line_bc = 0
    while (not 'THERMAL BOUNDARY CONDITIONS' in contents[line_bc]):
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
    f.write('   1.00000     P065: #iofiles (eg, 0 or 64); <0 --> sep. dirs\n')
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
#    f.write('       0 Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n')
    f.write('  ***** BOUNDARY CONDITIONS *****\n')
    f.write('  ***** FLUID BOUNDARY CONDITIONS *****\n')
#    f.write('   ***** NO THERMAL BOUNDARY CONDITIONS *****\n')
    f.write('   ***** THERMAL BOUNDARY CONDITIONS *****\n')
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


def dump_input_vars(R, nR, nSq, N, Re_t, stretch_sq, dr_sq_ratio,\
        dr_sq_int_ratio, distri_on, a_interf,\
        tog_r_out_const, tog_a_on_dist):
    """ Print all the input variables so the output can be 
    saved and used to correctly recreate the mesh.
    """

    print('INPUT VARIABLES:')
    print('----------------')
    print('R                = {0:10.5f}'.format(R))
    print('nR               = {0:10.5f}'.format(nR))
    print('nSq              = {0:10.5f}'.format(nSq))
    print('N                = {0:10.5f}'.format(N))
    print('Re_t             = {0:10.5f}'.format(Re_t))
    print('stretch_sq       = {0:10.5f}'.format(stretch_sq))
    print('dr_sq_ratio      = {0:10.5f}'.format(dr_sq_ratio))
    print('dr_sq_int_ratio  = {0:10.5f}'.format(dr_sq_int_ratio))
    print('distri_on        = {0:10.5f}'.format(distri_on))
    print('a_interf         = {0:10.5f}'.format(a_interf))
    print('tog_r_out_const  = {0:10.5f}'.format(tog_r_out_const))
    print('tog_a_on_dist    = {0:10.5f}'.format(tog_a_on_dist))


def check_mesh_quality(elements, nR, nSq, R, N , Re_t):
    """ Find minimum and maximum radial and circumferential 
    element lengths and element angles (distortion from 90°).
    Find resolution of the generated mesh considering 
    Gauss-Lobatto-Legendre distribution of grid points 
    within each element.
    """
    
    # only check first quadrant
    nel_quarter = nSq**2 + (nR-nSq)*2*nSq       # number of elements in one quarter
    nPhi = 8*nSq

    elements = elements[0:nel_quarter]
    l_r_max = 0
    l_r_min = 1e5
    l_p_max = 0
    l_p_min = 1e5
    alph_max = 0
    alph_min = 1e5
    x_gll = np.zeros(N+1)   # Distribution of GLL points in reference element [-1,1]
    d_x_gll = np.zeros(N+1) # Length between two adjacent GLL points
    el_wall_ind = nSq**2+nSq*2*(nR-nSq-1)   # index of element at wall
    w_ind = el_wall_ind     # copy
    cum_el = 0              # radial length of cumulative elements before 10th pt from the wall
    dz_rec = 0              # recommended size of streamwise elements


    # Get size of elements themselves
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

    # Get size of the actual grid by considering GLL distribution of grid points
    # GLL distribution on reference element x in [-1, 1]
    x_gll = my_math.get_gll(N)
    # Distance between two points on ref. element
    d_x_gll = x_gll[0:-2] - x_gll[1:-1]
    
    r_plus_min = l_r_min*min(d_x_gll)*0.5*Re_t
    r_plus_max = l_r_max*max(d_x_gll)*0.5*Re_t

    el_wall_ind = w_ind

    # Resolution in radial direction calculated at y=0
    r_1_plus = ((elements[w_ind].y[3] - elements[w_ind].y[0]))\
            *min(d_x_gll)*0.5*Re_t

    # First 10 point away from the wall is in which element?
    away_from_wall = int(m.ceil(10 / (N+1)))
    while away_from_wall > 1:
        # Cumulative elements' radial length at the wall closer than 10th pt
        cum_el = cum_el + elements[w_ind].y[3] - elements[w_ind].y[0]            

        # Update 
        away_from_wall = away_from_wall -1
        w_ind = w_ind - 2*nSq

    # Remaining pts up to 10th point x-th element
    rem_pts = 10 % (N+1)
    r_10_plus = (cum_el + \
            (elements[w_ind].y[3] - elements[w_ind].y[0])* \
            np.sum(d_x_gll[:rem_pts])*0.5)*Re_t

    # Resolution in circumferential direction
    
    # First, we have to find the angle theta spanned by one element
    theta_el = m.pi/2/(2*nSq)
    # Only a portion of that is spanned between two adjacent grid points
    theta_max_gp = theta_el*(max(d_x_gll)*0.5)
    theta_min_gp = theta_el*(min(d_x_gll)*0.5)
    r_theta_max = R*theta_max_gp*Re_t
    r_theta_min = R*theta_min_gp*Re_t



    dz_rec = 10/(max(d_x_gll)*0.5*Re_t)
            

    # Write a little output to stdout
    print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
    print('Some information about the mesh size:')
    print('-------------------------------------')
    print('Total number of elements in plane: {0:d}'.format(len(elements)*4))
    print('Delta R max      = {0:10.5f} at {1:d}'.format(l_r_max, el_r_max))
    print('Delta R min      = {0:10.5f} at {1:d}'.format(l_r_min, el_r_min))   
    print('Delta phi max    = {0:10.5f} at {1:d}'.format(l_p_max, el_p_max))
    print('Delta phi min    = {0:10.5f} at {1:d}'.format(l_p_min, el_p_min))   
    print('R*phi max        = {0:10.5f}'.format(2*m.pi/nPhi*R))
    print('alpha max        = {0:10.5f}° at {1:d}'.format(alph_max*(180/m.pi), el_alph_max))
    print('alpha min        = {0:10.5f}° at {1:d}'.format(alph_min*(180/m.pi), el_alph_min))
    print('Note that curvature is not considered here!')
    print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
    print('RESOLUTION')
    print('----------')
    print('r+ min           = {0:10.5f}  (< 1  )'.format(r_plus_min))
    print('r+ max           = {0:10.5f}  (< 5  )'.format(r_plus_max))
    print('r1+              = {0:10.5f}  (< 1  )'.format(r_1_plus))
    print('r10+             = {0:10.5f}  (<10  )'.format(r_10_plus))
    print('R theta plus max = {0:10.5f}° (<  5°)'.format(r_theta_max))
    print('R theta plus min = {0:10.5f}° (<1.5°)'.format(r_theta_min))
    print('For z+ < 10, element length in streamwise < {0:10.5f}'.format(dz_rec))
    print('Radial resolution is evaluated at vertical axis.')
    print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
