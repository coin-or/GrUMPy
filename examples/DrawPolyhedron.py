try:
    from src.grumpy.polyhedron2D import Polyhedron2D, Figure
except ImportError:
    from coinor.grumpy.polyhedron2D import Polyhedron2D, Figure

import MIP5 as LP

def disp_polyhedron(A = None, b = None, points = None, rays = None, c = None, obj_val = None,
                    opt = None, loc = None):
    if loc is None and opt is not None:
        loc = opt
    f = Figure()
    f.add_polyhedron(p, label = 'Polyhedron $P$', color = 'red')
    f.set_xlim(p.plot_min[0], p.plot_max[0])
    f.set_ylim(p.plot_min[1], p.plot_max[1])
    if c is not None and obj_val is not None:
        f.add_line(c, obj_val, p.plot_max - [0.2, 0.2], p.plot_min + [0.2, 0.2], 
                   linestyle = 'dashed', color = 'black', label = "Objective Function")
    f.add_point(opt, 0.04, 'red')
    f.add_text(loc[0], loc[1], r'$x^* = %s$' % str(opt))
    f.show()
    

if LP.A is not None and LP.b is not None:
    p = Polyhedron2D(A = LP.A, b = LP.b)
elif LP.points is not None and LP.rays is not None: 
    p = Polyhedron2D(points = LP.points, rays = LP.rays)
else:
    print 'Error: Must specify either A and b or points and rays'
    p = None

if p is not None:
    disp_polyhedron(A = p.hrep.A, b = p.hrep.b, 
                    c = LP.c, obj_val = LP.obj_val,
                    opt = LP.opt, loc = (LP.opt[0]+0.1, LP.opt[1]-0.1))

    if LP.cuts is not None:
        Acuts = p.A[:]
        bcuts = p.b[:]
        for i in range(len(LP.cuts)):
            Acuts.append(LP.cuts[i])
            bcuts.append(LP.rhs[i])
            disp_polyhedron(A = Acuts, b = bcuts)        