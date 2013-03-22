from Polyhedron2D import Polyhedron2D, add_line
import matplotlib.pyplot as plt

A = [[-0.2, -0.3], [-0.4, -0.2], [-0.3, -0.4], 
     [1, 0], [0, 1], [-1, 0], [0, -1]]
b = [-0.5, -1.5, -2.0, 9, 6, 0, 0]
fig = plt.figure()
ax = fig.add_subplot(111)
ax.grid()
p = Polyhedron2D(A = A, b = b)
p.draw(ax, color = 'blue', linestyle = 'solid')
ax.set_xlim(p.plot_min[0], p.plot_max[0])
ax.set_ylim(p.plot_min[1], p.plot_max[1])
pI = p.make_integer_hull()
pI.draw(ax, color = 'red', linestyle = 'dashed') 
add_line(ax, [20, 15], 100, p.plot_max - [0.2, 0.2], p.plot_min + [0.2, 0.2], 
         linestyle = 'dashed')
plt.show()

cuts = [[-4.8, -4.4], 
        [-7.73, -5.87],
        [-6, -4], 
        [-3, -3],
        ]
rhs = [-26, 
       -38, 
       -28, 
       -18, 
       ]
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

