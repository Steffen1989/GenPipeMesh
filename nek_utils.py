# A collection of important functions
#import re
import sys
import pdb

def set_vertices(elements,nSq,dx,dy):
    """ Set vertex location for each element. """
    for el in elements:
        i = (el.number-1)%nSq     # column number
        j = int((el.number-1)/nSq)    # row number
#    for j in range(nSq):    # loop through each column
#        for i in range(nSq):    # loop through each row
#            el = elements[i+j*nSq]
        el.x = [i*dx, (i+1)*dx, (i+1)*dx, i*dx]  # set vertex coordinates
        el.y = [j*dy, j*dy, (j+1)*dy, (j+1)*dy]  # set vertex coordinates

def set_bc(elements,nSq):
    """ Set boundary conditions for each face. """

    for el in elements:
        n = el.number
        position = check_position(n, nSq)
        if (position == 'lowleft'):
            el.bc = ['W  ','E  ','E  ','W  ']
            el.bc_con_el = [0, n+1, n+nSq, 0]
            el.bc_con_f = [0, 4, 1, 0]
        elif (position == 'lowright'):
            el.bc = ['W  ','W  ','E  ','E  ']
            el.bc_con_el = [0, n+1, n+nSq, n-1] 
            el.bc_con_f = [0, 4, 1, 2]
        elif (position == 'upleft'):
            el.bc = ['E  ','E  ','W  ','W  ']
            el.bc_con_el = [n-nSq, n+1, 0, 0]           
            el.bc_con_f = [3, 4, 0, 0]
        elif (position == 'upright'):
            el.bc = ['E  ','W  ','W  ','E  ']
            el.bc_con_el = [n-nSq, 0, 0, n-1]
            el.bc_con_f = [3, 0, 0, 2]
        elif (position == 'northrow'):
            el.bc = ['E  ','E  ','W  ','E  ']
            el.bc_con_el = [n-nSq, n+1, 0, n-1]
            el.bc_con_f = [3, 4, 0, 2]
        elif (position == 'eastcol'):
            el.bc = ['E  ','W  ','E  ','E  ']
            el.bc_con_el = [n-nSq, 0, n+nSq, n-1]
            el.bc_con_f = [3, 0, 1, 2]
        elif (position == 'southrow'):
            el.bc = ['W  ','E  ','E  ','E  ']
            el.bc_con_el = [0, n+1, n+nSq, n-1]
            el.bc_con_f = [0, 4, 1, 2]
        elif (position == 'westcol'):
            el.bc = ['E  ','E  ','E  ','W  ']
            el.bc_con_el = [n-nSq, n+1, n+nSq, 0]
            el.bc_con_f = [3, 4, 1, 0]
        elif (position == 'internal'):
            el.bc = ['E  ','E  ','E  ','E  ']
            el.bc_con_el = [n-nSq, n+1, n+nSq, n-1]
            el.bc_con_f = [3, 4, 1, 2]
        else:
            print('position of element not found!')
            sys.exit(1)
          
def check_position(number, nSq):
    """ Check position of the given element. """
    n = number
    if (n == 1 or n == nSq or n == nSq**2-nSq+1 or n == nSq**2):    # corners
        if (n == 1):  # we are on the first element on the lower left
            pos = 'lowleft'
        elif (n == nSq):    # we are on the lower right corner
            pos = 'lowright'
        elif (n == nSq**2-nSq+1):     # we are on the upper left corner
            pos = 'upleft'
        elif (n == nSq**2):   # last element on the upper right
            pos = 'upright'
        return pos
    elif (n > nSq**2-nSq or n%nSq == 0 or n < nSq or n%(nSq+1) == 0):  # edges
        if (n > nSq**2-nSq):   # northern row
            pos = 'northrow'
        elif ((n%nSq) == 0): # eastern column
            pos = 'eastcol'
        elif (n < nSq):  # southern row
            pos = 'southrow'
        elif ((n%(nSq+1)) == 0):  #  western column
            pos = 'westcol'
        return pos
    else:   # interior
        pos = 'internal'
        return pos



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

def write_bc(elements):
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
