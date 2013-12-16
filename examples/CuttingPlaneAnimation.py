from coinor.grumpy.polyhedron2D import Polyhedron2D, add_line
import matplotlib.pyplot as plt

from MIP2 import A, b, cuts, rhs, points, rays
try:
    from MIP2 import c, obj_val
except ImportError:
    c = None
fig = plt.figure()
ax = fig.add_subplot(111)
ax.grid()
if points is not None:
    p = Polyhedron2D(points = points, rays = rays)
else:
    p = Polyhedron2D(A = A, b = b)
p.draw(ax, color = 'blue', linestyle = 'solid')
ax.set_xlim(p.plot_min[0], p.plot_max[0])
ax.set_ylim(p.plot_min[1], p.plot_max[1])
pI = p.make_integer_hull()
pI.draw(ax, color = 'red', linestyle = 'dashed')
if c is not None:
    add_line(ax, c, obj_val, p.plot_max - [0.2, 0.2], p.plot_min + [0.2, 0.2], 
             linestyle = 'dashed')
plt.show()

if cuts is not None:
    Acuts = A[:]
    bcuts = b[:]
    for i in range(len(cuts)):
        Acuts.append(cuts[i])
        bcuts.append(rhs[i])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid()
        p = Polyhedron2D(A = A, b = b)
        p.draw(ax, color = 'blue', linestyle = 'solid')
        ax.set_xlim(p.plot_min[0], p.plot_max[0])
        ax.set_ylim(p.plot_min[1], p.plot_max[1])
        pI = p.make_integer_hull()
        pI.draw(ax, color = 'red', linestyle = 'dashed') 
        pC = Polyhedron2D(A = Acuts, b = bcuts)
        pC.draw(ax, color = 'black', linestyle = 'solid')
        plt.show()
