try:
    from src.grumpy.polyhedron2D import Polyhedron2D, Figure
except ImportError:
    from coinor.grumpy.polyhedron2D import Polyhedron2D, Figure
import numpy as np
from cylp.cy import CyClpSimplex
from cylp.py.modeling import CyLPArray, CyLPModel

import MIP6 as LP

def disp_polyhedron(A = None, b = None, points = None, rays = None, c = None, obj_val = None,
                    opt = None, loc = None):
    if loc is None and opt is not None:
        loc = opt
    f = Figure()
    f.add_polyhedron(p, label = 'Polyhedron $P$', color = 'red')
    f.set_xlim(p.plot_min[0], p.plot_max[0]+1)
    f.set_ylim(p.plot_min[1], p.plot_max[1]+2)
    if c is not None and obj_val is not None:
        f.add_line(c, obj_val, p.plot_max + [0.8, 1.8], p.plot_min + [0.2, 0.2], 
                   linestyle = 'dashed', color = 'black', label = "Objective Function")
    if opt is not None:
        f.add_point(opt, 0.04, 'red')
        f.add_text(loc[0], loc[1], r'$x^* = %s$' % str(opt))
    f.show()
    
if LP.A is not None and LP.b is not None:
    p = Polyhedron2D(A = LP.A + [[-1, 0], [0, -1]], b = LP.b + [0, 0])
elif LP.points is not None and LP.rays is not None: 
    p = Polyhedron2D(points = LP.points, rays = LP.rays)
else:
    print 'Error: Must specify either A and b or points and rays'
    p = None

if p is not None:

    lp = CyClpSimplex()
        
    A = np.matrix(p.hrep.A)
    b = CyLPArray(p.hrep.b)
            
    print A
    print b

    disp_polyhedron(A = A, b = b)

    x = lp.addVariable('x', LP.numVars)
        
    lp += A * x <= b
    lp += x >= 0

    c = CyLPArray(LP.c)
    # We are maximizing, so negate objective
    lp.objective = -c * x
    lp.logLevel = 0
    lp.primal(startFinishOptions = 'x')
    np.set_printoptions(precision = 2, linewidth = 200)
    print 'Basic variables: ', lp.basicVariables
    print "Current tableaux and reduced costs:"
    print lp.reducedCosts
    print np.around(lp.tableau, decimals = 3)
    print 'Right-hand side of optimal tableaux:'
    #There is a bug in CyLP and this is wrong
    #print lp.rhs
    print np.dot(lp.basisInverse, lp.constraintsUpper)
    print 'Inverse of optimal basis:'
    print np.around(lp.basisInverse, 3)
    obj_val = -lp.objectiveValue
    psol = np.around(lp.primalVariableSolution['x'], 2)
    print 'Optimal Value:', obj_val
    print 'Primal solution:', psol
    print 'Dual solution:', lp.dualConstraintSolution['R_1']

    disp_polyhedron(A = A, b = b, c = c, obj_val = obj_val,
                    opt = psol.tolist(), loc = (psol[0]+0.1, psol[1]-0.1))
     