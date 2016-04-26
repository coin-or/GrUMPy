try:
    from src.grumpy.polyhedron2D import Polyhedron2D, Figure
except ImportError:
    from coinor.grumpy.polyhedron2D import Polyhedron2D, Figure

CYLP_INSTALLED = True
try:
    import numpy as np
    from cylp.cy import CyClpSimplex
    from cylp.py.modeling import CyLPArray, CyLPModel
except ImportError:
    CYLP_INSTALLED = False
    
import LP9 as LP

def disp_polyhedron(A = None, b = None, points = None, rays = None, c = None, obj_val = None,
                    opt = None, loc = None):
    if loc is None and opt is not None:
        loc = opt
    f = Figure()
    f.add_polyhedron(p, label = 'Polyhedron $P$', color = 'red')
    f.set_xlim([p.xlim[0], p.xlim[1]+1])
    f.set_ylim([p.ylim[0], p.ylim[1]+2])
    if c is not None and obj_val is not None:
        f.add_line(c, obj_val, p.xlim + [0.2, 0.8], p.ylim + [0.2, 1.8], 
                   linestyle = 'dashed', color = 'black', label = "Objective Function")
    if opt is not None:
        f.add_point(opt, 0.04, 'red')
        f.add_text(loc, r'$x^* = %s$' % str(opt))
    f.show()
    
try:
    p = Polyhedron2D(A = LP.A, b = LP.b)
except AttributeError:
    try:
        p = Polyhedron2D(points = LP.points, rays = LP.rays)
    except AttributeError:
        print 'Error: Must specify either A and b or points and rays'
        p = None

if p is not None:

    if CYLP_INSTALLED:
        lp = CyClpSimplex()
        
        A = np.matrix(p.hrep.A)
        b = CyLPArray(p.hrep.b)
        
        print A
        print b
    
        if LP.numVars == 2:
            disp_polyhedron(A = A, b = b)
    
        x = lp.addVariable('x', LP.numVars)
            
        if LP.sense[0] == '>=':
            lp += A * x >= b
        else:
            lp += A * x <= b
        #lp += x >= 0
    
        c = CyLPArray(LP.c)
        # We are maximizing, so negate objective
        if LP.sense[1] == 'Min':
            lp.objective = c * x
        else:
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
        if LP.sense[0] == '<=':
            print np.dot(lp.basisInverse, lp.constraintsUpper)
        else:
            print np.dot(lp.basisInverse, lp.constraintsLower)
        print 'Inverse of optimal basis:'
        print np.around(lp.basisInverse, 3)
        if LP.sense[1] == 'Min':
            obj_val = lp.objectiveValue
        else:
            obj_val = -lp.objectiveValue
            
        psol = np.around(lp.primalVariableSolution['x'], 2)
        print 'Optimal Value:', obj_val
        print 'Primal solution:', psol
        print 'Dual solution:', lp.dualConstraintSolution['R_1']
    
        if LP.numVars == 2:
            disp_polyhedron(A = A, b = b, c = c, obj_val = obj_val,
                            opt = psol.tolist(), loc = (psol[0]+0.1, psol[1]-0.1))
    else:
        print LP.A
        print LP.b
    
        disp_polyhedron(A = p.hrep.A, b = p.hrep.b)
