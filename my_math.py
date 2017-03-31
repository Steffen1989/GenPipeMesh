# Some mathematics
import math as m
import numpy as np
import pdb


def geom_prog(N, a_start, a_end, j):
    """ Geometric progression a(j) = a_start * r**j 

    with a(0) = a_start 
    and a(N) = a_start*r**(N-1) = a_end

    N       : number of steps
    a_start : starting value at j=0
    a_end   : finale value at j=N-1
    j       : step
    """

    r = (a_end/a_start)**(1/(N-1))
    ret = a_start * r**j

    return ret

def sin_dist(N, start, end, i):
    """ Sine distribution with a clustering at the end

    N       : number of steps
    start   : starting value at i=0
    end     : final value at i=N-1

    """
    if (end > start):   # swap
        start, end = end, start
    
    r = start - (start-end)*m.sin(m.pi/2*i/(N-1))

    return r

def exp_dist(N, start, end, i):
    """ Exponential decay distribution 

    N       : number of steps
    start   : starting value at i=0
    end     : final value at i=N-1

    """

    if (end > start):   # swap
        start, end = end, start
    
    if (i == 0):    # catch division by zero
        r = start
    else:
        r = (end - start)/(1-m.exp(-(N-1)))*(1-m.exp(-i)) + start

    return r
    


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

def get_line_params(pt1_x, pt1_y, pt2_x, pt2_y):
    """ Find the slope and y-intercept for a line 
    going through pt1 and pt2

    y = m*x + c
    """

    m = (pt2_y - pt1_y)/(pt2_x - pt1_x)
    c = pt1_y - m*pt1_x
    return (m, c)


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

def newton_raphson(x0, func, func_prime):
    """ Newton-Raphson algorithm for determining zeros.

    x0          : initial guess
    func        : function
    func_prime  : first derivative
    """
    delta = 1e3
    eps = 1e-16
    x_old = x0

    while delta > eps:
        x_new = x_old - func(x_old)/func_prime(x_old)
        x_old = x_new
        delta = abs(x_new-x_old)

    return x_old


def vec_angle(vec1, vec2):
    """ Return the angle between two vectors """

    if (vec1[1] == -0.0):   # for elements in first row vec1[1] == -0.0
        vec1[1] = 0.0

    a1 = m.atan2(vec1[1],vec1[0])
    a2 = m.atan2(vec2[1],vec2[0])

    return m.fabs(a1-a2)

def get_rad_ell(a, b, c, x):
    """ Calculate the curvature at x with the formula
    kappa = y'' / (1+y'**2)**(3/2)
    and return its inverse the radius

    a : semi-major axis
    b : semi-minor axis
    c : "radius" or constant on the rhs
    x : dependent variable
    """

    y_pr = -b*x/a**2 * (c-x**2/a**2)**(-1/2)
    y_prpr = -b*x**2/a**4 * (c-x**2/a**2)**(-3/2) \
            - b/a**2 * (c-x**2/a**2)**(-1/2)
    kappa = abs(y_prpr)/(1+y_pr**2)**(3/2)
    return 1/kappa

def get_gll(N):
    """ Get the Gauss-Lobatto-Legendre distribution of grid points.

    The quadrature nodes (respectively the grid points) are given
    by (see Canuto et al. 2006: Spectral methods - Fundamentals 
    in single domains, ch. 2.3 , p. 75-76)
    x[0] = -1
    x[N] = 1
    x[j] = zeros of L_N' for j = 1, ..., N
    where L_N' is the first derivative of the n-th Legendre polynomial.
    """

    x = np.zeros(N+1)  # location of quadrature nodes
    V = np.zeros((N+1, N+1))  # Legendre Vandermonde matrix
    # This matrix evaluates a polynomial at a set of points so that 
    # Vx = y as in y = P(x)
    
    # Initial guess taken from Gauss-Lobatta-Chebychev quadrature nodes
    for i in range(0, N+1):
        x[i] = m.cos(i*m.pi/N)

    x_old = np.zeros(N+1)   # variable to store value of previous iteration
    eps = np.finfo(float).eps

    while max(abs(x - x_old)) > eps:

        x_old = x   # update old value to the previous value

        # Fill first two columns of Vandermonde matrix with starting
        # values of the recursion formula.
        V[:, 0] = 1
        V[:, 1] = x

        # Calculate the remaining columns of the Vandermonde matrix
        # with recursion formula
        for k in range(1, N):
            V[:,k+1] = (2*k+1)/(k+1)*x*V[:,k] - k/(k+1)*V[:,k-1]
            
        # Update x values with Newton-Raphson algorithm to find the 
        # roots of (1-x**2)*L_N'(x) = 0
        # Derivation of the numerator and denominator are a bit tricky.
        x = x_old - (x*V[:, N] - V[:, N-1])/(V[:, N]*(1+N))

    return x

