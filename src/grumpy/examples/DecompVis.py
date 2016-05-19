try:
    from src.grumpy.polyhedron2D import Polyhedron2D, Figure
except ImportError:
    from coinor.grumpy.polyhedron2D import Polyhedron2D, Figure

numVars = 2
numCons = 4
#points = [[0, 0], [3, 4], [8, 6], [6, 1]]
points = None
rays = []
A1 = [
     [1, 1],
     [4, -10],
     [-2, -2],
     [-6, -2],
     [-1, 4]
     ]
A2 = [
     [-7, 1],
     [0, -1],
     [1, -1],
     [4, 1],
     [0, 1],
     [-1, 5]
     ]
b1 = [
     8,
     -3,
     -9,
     -19,
     12
     ]
b2 = [
     -13,
     -1,
     3,
     27,
     5,
     20
     ]
sense = ('Min', '<=')
integerIndices = [0, 1]
c = [1, 0]
cuts = None
rhs = None
obj_val = 2

p1 = Polyhedron2D(A = A1, b = b1)
p2 = Polyhedron2D(A = A2, b = b2)
sR = p2.make_integer_hull()

f = Figure()
f.add_polyhedron(p1, label = '$\mathcal{P}^{\,\prime}$', color = 'blue', show_int_points = True)
f.add_polyhedron(p2, label = '$\mathcal{P}^{\,\prime\prime}$', color = 'black', show_int_points = True)
f.add_polyhedron(sR, label = '$\operatorname{conv}(\mathcal{P}^{\,\prime\prime} \cap \mathbb{Z}^{\,2})$', color = 'red', linestyle = 'dashed')
f.set_xlim(p2.xlim + [0, 1])
f.set_ylim(p2.ylim + [0, 1])
opt = [2, 3.5]
f.add_point(opt, .05, 'red')
f.add_text([opt[0]-0.25, opt[1]], r'$x^*$')
if c is not None and obj_val is not None:
    f.add_line(c, obj_val, p2.xlim, p2.ylim + [0.5, - 1.5], 
               linestyle = 'dashed', color = 'green', label = "Objective Function")
f.show()
cuts = [[-3, 1]]
rhs = [-5]
opt = [[2.42, 2.25]]
obj_val = [2.42]

for i in range(len(cuts)):
    A1.append(cuts[i])
    b1.append(rhs[i])
    p1 = Polyhedron2D(A = A1, b = b1)
    f = Figure()
    f.add_polyhedron(p1, label = '$\mathcal{P}^{\,\prime}$', color = 'blue', show_int_points = True)
    f.add_polyhedron(p2, label = '$\mathcal{P}^{\,\prime\prime}$', color = 'black', show_int_points = True)
    f.add_polyhedron(sR, label = '$\operatorname{conv}(\mathcal{P}^{\,\prime\prime} \cap \mathbb{Z}^{\,2})$', color = 'red', linestyle = 'dashed')
    f.set_xlim(p2.xlim + [0, 1])
    f.set_ylim(p2.ylim + [0, 1])
    f.add_point(opt[i], .05, 'red')
    f.add_text([opt[i][0]+0.1, opt[i][1]], r'$x^D$')
    f.add_line(cuts[i], rhs[i], p2.xlim, p2.ylim + [0.5, - 1.5], 
               linestyle = 'dashed', color = 'orange', label = "Valid Inequality")
    if c is not None and obj_val is not None:
        f.add_line(c, obj_val[i], p2.xlim, p2.ylim + [0.5, - 1.5], 
                   linestyle = 'dashed', color = 'green', label = "Objective Function")
    f.show()

